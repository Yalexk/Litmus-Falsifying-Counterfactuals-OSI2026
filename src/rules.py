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


def fit_logp_baseline(train) -> tuple[float, float]:
    """Fit y ~ MolLogP on the TRAIN split only. Returns (slope, intercept).

    LogP alone explains most of ESOL, so almost any polarity-flavoured rule
    survives the raw test. Scoring rules on the residual y - f(LogP) asks
    whether a claim says anything beyond "greasy molecules don't dissolve".
    """
    x = np.array([descriptors(s)["MolLogP"] for s in train["smiles"]])
    slope, intercept = np.polyfit(x, train["y"].values, 1)
    return float(slope), float(intercept)


def logp_bin_edges(train, n_bins: int = 6) -> np.ndarray:
    """Interior LogP quantile edges from the TRAIN split, for stratified tests."""
    x = np.array([descriptors(s)["MolLogP"] for s in train["smiles"]])
    return np.quantile(x, np.linspace(0, 1, n_bins + 1)[1:-1])


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


def evaluate(rule: Rule, df, min_support: int = 10, logp_fit=None,
             logp_edges=None) -> dict:
    """Score a rule on data the LLM never saw.

    effect   : mean(y | match) - mean(y | no match), in target units
    holds    : effect sign agrees with the claimed direction
    cohens_d : standardised effect, so you can rank across rules
    p_value  : Welch t-test. With many rules this needs correcting --
               use Benjamini-Hochberg in report(), not a bare 0.05 cutoff.

    With logp_fit=(slope, intercept) from fit_logp_baseline, the same stats
    are also computed on the LogP residual and stored with an _adj suffix.
    With logp_edges from logp_bin_edges, a line-free LogP-stratified effect
    is stored with a _strat suffix -- use it to confirm or kill an _adj
    sign reversal, since the linear fit is miscalibrated at high LogP.
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
    if logp_fit is None and logp_edges is None:
        return res
    logp = np.array([descriptors(s)["MolLogP"] for s in df["smiles"]])
    if logp_fit is not None:
        slope, intercept = logp_fit
        resid = y - (slope * logp + intercept)
        adj = _compare(resid[mask], resid[~mask], rule.direction)
        res.update({f"{k}_adj": v for k, v in adj.items()})
    if logp_edges is not None:
        bins = np.digitize(logp, logp_edges)
        res.update(_stratified(y, mask, bins, rule.direction))
    return res


def _bh(rows: list[dict], p_key: str, alpha: float) -> None:
    """Benjamini-Hochberg step-up; sets bh_significant[_adj] on each row."""
    out_key = "bh_significant" + p_key.removeprefix("p_value")
    rows = sorted(rows, key=lambda r: r[p_key])
    m = len(rows)
    k = max((i for i, r in enumerate(rows, 1) if r[p_key] <= alpha * i / m), default=0)
    for i, r in enumerate(rows, 1):
        r[out_key] = i <= k


def report(results: list[dict], alpha: float = 0.05, noise_floor: float = 0.6):
    """Benjamini-Hochberg across all tested rules, plus a noise floor.

    noise_floor is in target units -- for ESOL logS, experimental error is
    ~0.6, so an effect smaller than that is not a finding no matter what the
    p-value says. This single line is what separates your submission from a
    demo that launders noise into confident chemistry prose.
    """
    ev = [r for r in results if r.get("status") == "evaluated"]
    ev.sort(key=lambda r: r["p_value"])
    adjusted = bool(ev) and "p_value_adj" in ev[0]
    _bh(ev, "p_value", alpha)
    if adjusted:
        _bh(ev, "p_value_adj", alpha)
    strat = [r for r in ev if r.get("p_value_strat") is not None]
    _bh(strat, "p_value_strat", alpha)
    for r in ev:
        r["survives_strat"] = bool(
            r.get("holds_strat") and r.get("bh_significant_strat")
            and abs(r["effect_strat"]) >= noise_floor
        )
    for r in ev:
        r["survives"] = bool(
            r["holds"] and r["bh_significant"] and abs(r["effect"]) >= noise_floor
        )
        if adjusted:
            r["survives_adj"] = bool(
                r["holds_adj"] and r["bh_significant_adj"]
                and abs(r["effect_adj"]) >= noise_floor
            )
    for r in results:
        if r.get("status") != "evaluated":
            r["bh_significant"] = False
            r["survives"] = False
            r["survives_strat"] = False
            if adjusted:
                r["survives_adj"] = False
    n = max(len(results), 1)
    n_sur = sum(r["survives"] for r in ev)
    print(f"\n{n_sur}/{len(results)} LLM hypotheses survived held-out testing "
          f"({n_sur / n:.0%})")
    if adjusted:
        n_adj = sum(r["survives_adj"] for r in ev)
        print(f"{n_adj}/{len(results)} survived beyond LogP ({n_adj / n:.0%})")
    if strat:
        n_st = sum(r["survives_strat"] for r in ev)
        print(f"{n_st}/{len(results)} survived within LogP bins ({n_st / n:.0%})")
    # Return every rule, including insufficient-support ones, so the saved
    # artefact has the same denominator as the headline.
    return ev + [r for r in results if r.get("status") != "evaluated"]
