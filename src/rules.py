"""Executable rules + the held-out scoring harness.

This is the part of the project that is actually yours. Everything upstream
is a pip install; this is where an LLM hypothesis becomes falsifiable.

A rule is a claim of the form:
    molecules matching CONDITION have higher/lower y than those that don't

Two condition kinds:
  descriptor : {"kind":"descriptor","name":"MolLogP","op":">","value":3.0}
  smarts     : {"kind":"smarts","pattern":"[CX3](=O)[OX2H1]"}   # COOH

Start with descriptor rules. They almost always run. Upgrade to SMARTS only
if you have time left -- getting an LLM to emit valid, correctly-matching
SMARTS is the single most likely thing to eat your afternoon.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from rdkit import Chem

from .counterfactuals import descriptors

# A small starter vocabulary. Hand these to the LLM so it names real groups
# instead of inventing SMARTS. Verified to parse.
GROUP_SMARTS = {
    "carboxylic_acid": "[CX3](=O)[OX2H1]",
    "ester": "[CX3](=O)[OX2H0][#6]",
    "hydroxyl": "[OX2H]",
    "primary_amine": "[NX3;H2;!$(NC=O)]",
    "amide": "[NX3][CX3](=[OX1])",
    "nitro": "[N+](=O)[O-]",
    "halogen": "[F,Cl,Br,I]",
    "aromatic_ring": "a1aaaaa1",
    "ether": "[OD2]([#6])[#6]",
    "ketone": "[#6][CX3](=O)[#6]",
    "sulfonamide": "[SX4](=O)(=O)[NX3]",
    "nitrile": "[NX1]#[CX2]",
}


@dataclass
class Rule:
    text: str                  # the LLM's natural-language claim
    condition: dict            # executable form
    direction: str             # "higher" | "lower"
    provenance: dict = field(default_factory=dict)

    def matches(self, smiles: str) -> bool:
        c = self.condition
        if c["kind"] == "descriptor":
            d = descriptors(smiles)
            v = d.get(c["name"])
            if v is None:
                return False
            return {">": v > c["value"], "<": v < c["value"],
                    ">=": v >= c["value"], "<=": v <= c["value"]}[c["op"]]
        if c["kind"] == "smarts":
            pat = Chem.MolFromSmarts(c["pattern"])
            m = Chem.MolFromSmiles(smiles)
            return bool(pat and m and m.HasSubstructMatch(pat))
        raise ValueError(f"unknown condition kind {c['kind']}")


def validate(rule: Rule) -> str | None:
    """Reject malformed rules before wasting an eval. Returns an error string."""
    c = rule.condition
    if rule.direction not in ("higher", "lower"):
        return f"bad direction {rule.direction!r}"
    if c.get("kind") == "smarts":
        if Chem.MolFromSmarts(c.get("pattern", "")) is None:
            return f"unparseable SMARTS {c.get('pattern')!r}"
    elif c.get("kind") == "descriptor":
        if c.get("name") not in descriptors("CCO"):
            return f"unknown descriptor {c.get('name')!r}"
        if c.get("op") not in (">", "<", ">=", "<="):
            return f"bad op {c.get('op')!r}"
    else:
        return f"unknown kind {c.get('kind')!r}"
    return None


@dataclass
class Control:
    """The descriptor that dominates the property (MolLogP for solubility).

    One descriptor often explains most of a property, so almost any rule
    that tracks it survives the raw test. The control asks whether a claim
    says anything beyond it. Both the line (slope, intercept; diagnostic
    only) and the bin edges are fitted on the TRAIN split only.
    """
    descriptor: str
    slope: float
    intercept: float
    edges: np.ndarray

    @classmethod
    def fit(cls, train, descriptor: str, n_bins: int = 6) -> "Control":
        x = np.array([descriptors(s)[descriptor] for s in train["smiles"]])
        slope, intercept = np.polyfit(x, train["y"].values, 1)
        edges = np.quantile(x, np.linspace(0, 1, n_bins + 1)[1:-1])
        return cls(descriptor, float(slope), float(intercept), edges)


# Verdict names before the control descriptor was configurable (ESOL runs
# 1-3 used MolLogP). Apply when reading old results.
LEGACY_VERDICTS = {"survives_beyond_logp": "survives_beyond_control",
                   "logp_only": "control_only",
                   "untestable_beyond_logp": "untestable_beyond_control"}


def _stratified(y, mask, bins, direction: str, min_per_side: int = 3,
                n_perm: int = 2000, seed: int = 0) -> dict:
    """Match-vs-rest effect within LogP bins, no line fitted.

    Effect is the weighted mean of within-bin differences, weights
    n_match*n_rest/(n_match+n_rest). p-value is a two-sided permutation test
    that shuffles the match labels within each bin. Bins without at least
    min_per_side molecules on each side are skipped.
    """
    used = [b for b in np.unique(bins)
            if min(mask[bins == b].sum(), (~mask[bins == b]).sum()) >= min_per_side]
    if not used:
        return dict(effect_strat=None, p_value_strat=None, holds_strat=False,
                    strat_bins_used=0, strat_coverage=0.0)

    def effect(m):
        num = den = 0.0
        for b in used:
            yb, mb = y[bins == b], m[bins == b]
            k, j = mb.sum(), (~mb).sum()
            w = k * j / (k + j)
            num += w * (yb[mb].mean() - yb[~mb].mean())
            den += w
        return num / den

    obs = effect(mask)
    rng = np.random.default_rng(seed)
    perm = mask.copy()
    hits = 0
    for _ in range(n_perm):
        for b in used:
            idx = np.flatnonzero(bins == b)
            perm[idx] = rng.permutation(mask[idx])
        hits += abs(effect(perm)) >= abs(obs)
    covered = np.isin(bins, used) & mask
    return dict(
        effect_strat=float(obs),
        p_value_strat=float((hits + 1) / (n_perm + 1)),
        holds_strat=bool((obs > 0) == (direction == "higher")),
        strat_bins_used=len(used),
        strat_coverage=float(covered.sum() / max(mask.sum(), 1)),
    )


def _compare(a, b, direction: str) -> dict:
    from scipy import stats

    effect = float(a.mean() - b.mean())
    pooled = np.sqrt((a.var(ddof=1) + b.var(ddof=1)) / 2) or 1e-9
    t = stats.ttest_ind(a, b, equal_var=False)
    return dict(
        effect=effect,
        cohens_d=float(effect / pooled),
        p_value=float(t.pvalue),
        mean_match=float(a.mean()),
        mean_other=float(b.mean()),
        holds=bool((effect > 0) == (direction == "higher")),
    )


def evaluate(rule: Rule, df, min_support: int = 10,
             control: Control | None = None) -> dict:
    """Score a rule on data the LLM never saw.

    effect   : mean(y | match) - mean(y | no match), in target units
    holds    : effect sign agrees with the claimed direction
    cohens_d : standardised effect, so you can rank across rules
    p_value  : Welch t-test. With many rules this needs correcting --
               use Benjamini-Hochberg in report(), not a bare 0.05 cutoff.

    With a Control, the same stats are also computed on the residual of the
    control's linear fit (_adj suffix, diagnostic only -- the line can be
    miscalibrated, e.g. ESOL at high LogP), and a line-free effect within
    bins of the control descriptor is stored with a _strat suffix. _strat
    is what the beyond-control verdict uses.
    """
    mask = np.array([rule.matches(s) for s in df["smiles"]])
    y = df["y"].values
    n_hit = int(mask.sum())
    res = {
        "text": rule.text,
        "condition": rule.condition,
        "direction": rule.direction,
        "support": n_hit,
        "support_frac": n_hit / len(df),
    }
    if n_hit < min_support or (len(df) - n_hit) < min_support:
        res.update(status="insufficient_support", holds=False)
        return res
    res.update(status="evaluated", **_compare(y[mask], y[~mask], rule.direction))
    if control is None:
        return res
    x = np.array([descriptors(s)[control.descriptor] for s in df["smiles"]])
    resid = y - (control.slope * x + control.intercept)
    adj = _compare(resid[mask], resid[~mask], rule.direction)
    res.update({f"{k}_adj": v for k, v in adj.items()})
    res.update(_stratified(y, mask, np.digitize(x, control.edges), rule.direction))
    return res


def _bh(rows: list[dict], p_key: str, alpha: float) -> None:
    """Benjamini-Hochberg step-up; sets bh_significant[_adj] on each row."""
    out_key = "bh_significant" + p_key.removeprefix("p_value")
    rows = sorted(rows, key=lambda r: r[p_key])
    m = len(rows)
    k = max((i for i, r in enumerate(rows, 1) if r[p_key] <= alpha * i / m), default=0)
    for i, r in enumerate(rows, 1):
        r[out_key] = i <= k


def report(results: list[dict], alpha: float = 0.05, noise_floor: float = 0.6,
           min_coverage: float = 0.5, label: str = "LLM hypotheses",
           control_label: str = "the control"):
    """Benjamini-Hochberg across all tested rules, plus a noise floor.

    noise_floor is in target units and comes from the dataset config -- for
    ESOL logS, experimental error is ~0.6, so an effect smaller than that is
    not a finding no matter what the
    p-value says. This single line is what separates your submission from a
    demo that launders noise into confident chemistry prose.

    Call once per family (LLM rules, textbook baseline) -- each call is its
    own BH family.
    """
    ev = [r for r in results if r.get("status") == "evaluated"]
    ev.sort(key=lambda r: r["p_value"])
    _bh(ev, "p_value", alpha)
    # Beyond-control is judged by the binned test, and only where it has
    # enough coverage to mean anything. The linear _adj fields stay in the
    # output as a diagnostic only.
    strat = [r for r in ev if r.get("p_value_strat") is not None
             and r["strat_coverage"] >= min_coverage]
    _bh(strat, "p_value_strat", alpha)

    def finding(r, sfx):
        return (r.get("bh_significant" + sfx)
                and abs(r["effect" + sfx]) >= noise_floor)

    for r in ev:
        r["survives"] = bool(r["holds"] and finding(r, ""))
        r["strat_testable"] = any(r is s for s in strat)
        r["survives_strat"] = bool(
            r["strat_testable"] and r["holds_strat"] and finding(r, "_strat"))
        # Wrong sign AND a real effect: the claim is contradicted, not just
        # unsupported. These rows are deliverable 3.
        r["contradicted"] = bool(
            (not r["holds"] and finding(r, "")) or
            (r["strat_testable"] and not r["holds_strat"] and finding(r, "_strat")))
        r["verdict"] = (
            "contradicted" if r["contradicted"] else
            "fails" if not r["survives"] else
            "survives_beyond_control" if r["survives_strat"] else
            "control_only" if r["strat_testable"] else
            "untestable_beyond_control")
    for r in results:
        if r.get("status") != "evaluated":
            r.update(bh_significant=False, survives=False, survives_strat=False,
                     strat_testable=False, contradicted=False,
                     verdict="insufficient_support")
    n = max(len(results), 1)
    n_sur = sum(r["survives"] for r in ev)
    n_st = sum(r["survives_strat"] for r in ev)
    n_con = sum(r["contradicted"] for r in ev)
    print(f"\n{n_sur}/{len(results)} {label} survived held-out testing ({n_sur / n:.0%})")
    print(f"{n_st}/{len(results)} survived beyond {control_label}, binned ({n_st / n:.0%}); "
          f"{len(ev) - len(strat)} untestable at coverage < {min_coverage:.0%}")
    print(f"{n_con}/{len(results)} contradicted (wrong sign, real effect)")
    # Return every rule, including insufficient-support ones, so the saved
    # artefact has the same denominator as the headline.
    return ev + [r for r in results if r.get("status") != "evaluated"]
