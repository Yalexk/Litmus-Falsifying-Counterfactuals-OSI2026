"""End-to-end run.

    python run_pipeline.py --n-bases 10 --skip-llm     # no API key needed
    python run_pipeline.py --n-bases 20                # full loop

Stages are checkpointed to out/ so you can rerun the LLM stage without
regenerating counterfactuals. On hackathon day you will do that a lot.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

from src import counterfactuals as cf
from src import data as D
from src.model import BlackBox
from src.rules import Rule, evaluate, fit_logp_baseline, logp_bin_edges, report

PROPERTY = "aqueous solubility, log mol/L (higher = more soluble)"


def main(a):
    Path("out").mkdir(exist_ok=True)

    df = D.load_esol(a.csv) if a.csv else D.load_esol()
    train, valid, test = D.scaffold_split(df)
    print(f"{len(df)} molecules -> train {len(train)} valid {len(valid)} test {len(test)}")

    box = BlackBox().fit(train)
    metrics = [box.report(train, "train"), box.report(valid, "valid"), box.report(test, "test")]
    Path("out/model_metrics.json").write_text(json.dumps(metrics, indent=2))

    cf_path = Path("out/counterfactuals.json")
    if cf_path.exists() and not a.refresh_cf:
        records = json.loads(cf_path.read_text())
        print(f"reusing {len(records)} cached counterfactual sets")
    else:
        bases = train["smiles"].head(a.n_bases).tolist()
        print(f"generating counterfactuals for {len(bases)} base molecules...")
        records = cf.run_batch(bases, box, str(cf_path),
                               threshold=a.threshold, num_samples=a.num_samples)

    if a.skip_llm:
        print("\n--skip-llm set; stopping before hypothesis generation.")
        return

    from src.hypotheses import propose, revise

    rules = propose(records, PROPERTY, n=a.n_rules)
    print(f"\n{len(rules)} well-formed hypotheses proposed")

    logp_fit = fit_logp_baseline(train)
    print(f"LogP baseline (train): y = {logp_fit[0]:.2f}*LogP + {logp_fit[1]:.2f}")
    logp_edges = logp_bin_edges(train)
    ev = dict(logp_fit=logp_fit, logp_edges=logp_edges)

    # Revision feedback comes from the valid split only. The test split is
    # scored exactly once per final rule and never shown to the LLM.
    results, final = [], []
    for r in rules:
        if a.revise:
            res_valid = evaluate(r, valid, **ev)
            if not (res_valid.get("holds") and res_valid.get("holds_adj")):
                r2 = revise(r, res_valid, PROPERTY)
                if r2:
                    res2 = evaluate(r2, test, **ev)
                    res2["revised_from"] = r.text
                    res2["parent_valid"] = res_valid
                    results.append(res2)
                    final.append(r2)
                    continue
        results.append(evaluate(r, test, **ev))
        final.append(r)

    scored = report(results, noise_floor=a.noise_floor)
    Path("out/rules.json").write_text(json.dumps(scored, indent=2))

    print("\nRules (raw / beyond-LogP / within-LogP-bins effect):")
    for r in scored:
        if r.get("status") != "evaluated":
            print(f"  [insufficient support, n={r['support']}] {r['text']}")
            continue
        tag = ("SURVIVES" if r["survives_adj"] else
               "raw only" if r["survives"] else "fails")
        st = r["effect_strat"]
        st = "n/a" if st is None else f"{st:+.2f}"
        print(f"  {tag:8} [{r['effect']:+.2f} / {r['effect_adj']:+.2f} / {st}, "
              f"n={r['support']}] {r['text']}")


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--csv", default=None, help="local dataset CSV (else fetch ESOL)")
    p.add_argument("--n-bases", type=int, default=10)
    p.add_argument("--n-rules", type=int, default=8)
    p.add_argument("--num-samples", type=int, default=1500)
    p.add_argument("--threshold", type=float, default=1.0,
                   help="min |delta| in target units to count as a counterfactual")
    p.add_argument("--noise-floor", type=float, default=0.6,
                   help="effects smaller than this are not findings")
    p.add_argument("--refresh-cf", action="store_true")
    p.add_argument("--skip-llm", action="store_true")
    p.add_argument("--revise", action="store_true")
    main(p.parse_args())
