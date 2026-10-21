"""End-to-end run for one dataset, described by a config in configs/.

    python run_pipeline.py --dataset configs/esol.json --n-bases 10 --skip-llm
    python run_pipeline.py --dataset configs/esol.json --n-bases 20 --run-name r1

Stages are checkpointed to out/<dataset>/ so you can rerun the LLM stage
without regenerating counterfactuals. Runs given --run-name are archived to
results/<dataset>/<run-name>/.
"""
from __future__ import annotations

import argparse
import json
import shutil
import subprocess
import warnings
from collections import Counter
from datetime import datetime
from pathlib import Path

from src import config as C
from src import counterfactuals as cf
from src import data as D
from src.baseline import baseline_rules
from src.model import BlackBox
from src.rules import Control, evaluate, report

# exmol calls the model ~1M times; this sklearn warning fires on every call
# and produced a 344 MB log in run2.
warnings.filterwarnings("ignore", message=".*sklearn.utils.parallel.delayed.*")


def main(a):
    cfg = C.load(a.dataset)
    # CLI values override the config only when given; the manifest records both.
    noise_floor = cfg["noise_floor"] if a.noise_floor is None else a.noise_floor
    threshold = cfg["cf_threshold"] if a.threshold is None else a.threshold
    out = Path("out") / cfg["name"]
    out.mkdir(parents=True, exist_ok=True)

    df = D.load_dataset(cfg)
    train, valid, test = D.scaffold_split(df)
    print(f"[{cfg['name']}] {len(df)} molecules -> "
          f"train {len(train)} valid {len(valid)} test {len(test)}")

    box = BlackBox().fit(train)
    metrics = [box.report(train, "train"), box.report(valid, "valid"), box.report(test, "test")]
    (out / "model_metrics.json").write_text(json.dumps(metrics, indent=2))

    cf_path = out / "counterfactuals.json"
    if cf_path.exists() and not a.refresh_cf:
        records = json.loads(cf_path.read_text())
        print(f"reusing {len(records)} cached counterfactual sets")
    else:
        bases = D.spread_bases(train, a.n_bases)
        print(f"generating counterfactuals for {len(bases)} base molecules...")
        records = cf.run_batch(bases, box, str(cf_path),
                               threshold=threshold, num_samples=a.num_samples)

    control = (Control.fit(train, cfg["control"]["descriptor"])
               if cfg["control"] else None)
    control_label = cfg["control"]["descriptor"] if cfg["control"] else "the control"
    rep = dict(noise_floor=noise_floor, min_coverage=a.min_coverage,
               control_label=control_label)

    # Textbook baseline: same harness, own BH family. Needs no LLM.
    base = report([evaluate(r, test, control=control) | {"ref": r.provenance["ref"]}
                   for r in baseline_rules(cfg, train)],
                  label="textbook rules", **rep)
    (out / "baseline.json").write_text(json.dumps(base, indent=2))
    show(base)

    if a.skip_llm:
        print("\n--skip-llm set; stopping before hypothesis generation.")
        archive(a, cfg, out, base, None, noise_floor, threshold)
        return

    from src.hypotheses import propose, revise

    rules = propose(records, cfg, n=a.n_rules)
    print(f"\n{len(rules)} well-formed hypotheses proposed")

    # Revision feedback comes from the valid split only. The test split is
    # scored exactly once per final rule and never shown to the LLM.
    results = []
    for r in rules:
        if a.revise:
            res_valid = evaluate(r, valid, control=control)
            beyond_control_fails = (res_valid.get("strat_coverage", 0) >= a.min_coverage
                                    and not res_valid.get("holds_strat"))
            # Too rare to test on valid is not evidence against the rule; revising
            # it just swapped specific claims for generic ones (ESOL final, 4/4).
            testable = res_valid.get("status") == "evaluated"
            if testable and (not res_valid.get("holds") or beyond_control_fails):
                try:
                    r2 = revise(r, res_valid, cfg)
                except Exception as e:  # one flaky CLI call shouldn't cost the run
                    print(f"  revise failed, keeping original: {e}")
                    res = evaluate(r, test, control=control)
                    res["revise_failed"] = str(e)[:300]
                    results.append(res)
                    continue
                if r2:
                    res2 = evaluate(r2, test, control=control)
                    res2["revised_from"] = r.text
                    res2["parent_valid"] = res_valid
                    # Scored for the record only, never shown to the LLM: did
                    # revision turn a right claim wrong, or a wrong one right?
                    res2["parent_test"] = evaluate(r, test, control=control)
                    results.append(res2)
                    continue
        results.append(evaluate(r, test, control=control))

    # Count each claim once: a rule with the same condition and direction as an
    # earlier one (revision can reproduce one verbatim) is set aside, not scored,
    # so it can't inflate the survival rate or the BH family.
    seen, unique, dupes = set(), [], []
    for res in results:
        key = (json.dumps(res["condition"], sort_keys=True), res["direction"])
        (dupes if key in seen else unique).append(res)
        seen.add(key)
    (out / "rules_duplicates.json").write_text(json.dumps(dupes, indent=2))
    if dupes:
        print(f"\n{len(dupes)} duplicate rule(s) set aside; scoring {len(unique)} unique claims")

    scored = report(unique, **rep)
    (out / "rules.json").write_text(json.dumps(scored, indent=2))
    show(scored)
    archive(a, cfg, out, base, scored, noise_floor, threshold)


def archive(a, cfg, out, base, scored, noise_floor, threshold):
    """Copy this run's artefacts from out/<dataset>/ into
    results/<dataset>/<run-name>/ with a manifest and a copy of the config.
    out/ is scratch and gets overwritten; results/ is committed."""
    if not a.run_name:
        return
    dest = Path("results") / cfg["name"] / a.run_name
    dest.mkdir(parents=True, exist_ok=True)
    files = ["model_metrics.json", "counterfactuals.json", "baseline.json"]
    if scored is not None:
        files += ["rules.json", "rules_duplicates.json"]
    for f in files:
        shutil.copy(out / f, dest / f)
    shutil.copy(cfg["_path"], dest / "config.json")

    def git(*cmd):
        return subprocess.run(["git", *cmd], capture_output=True, text=True).stdout.strip()

    def verdicts(rows):
        return None if rows is None else dict(Counter(r["verdict"] for r in rows))

    manifest = {
        "dataset": cfg["name"],
        "run_name": a.run_name,
        "timestamp": datetime.now().isoformat(timespec="seconds"),
        "git_commit": git("rev-parse", "--short", "HEAD"),
        "git_dirty": bool(git("status", "--porcelain", "--untracked-files=no")),
        # When was the config first committed? For a pre-registered dataset
        # this must predate the run.
        "config_first_commit": git("log", "--diff-filter=A", "--format=%h %cI",
                                   "--", cfg["_path"]),
        "preregistered": cfg["preregistered"],
        "args": vars(a),
        "noise_floor_used": noise_floor,
        "cf_threshold_used": threshold,
        "control_descriptor": cfg["control"]["descriptor"] if cfg["control"] else None,
        "baseline_verdicts": verdicts(base),
        "llm_verdicts": verdicts(scored),
    }
    if manifest["git_dirty"]:
        # The commit hash alone can't reproduce a dirty run; keep the diff.
        (dest / "uncommitted.diff").write_text(git("diff", "HEAD"))
    if scored is not None:
        from src import hypotheses as H
        manifest["prompt_version"] = H.PROMPT_VERSION
        manifest["llm_model"] = H.MODEL
        (dest / "system_prompt.txt").write_text(H.system_prompt(cfg))
    (dest / "manifest.json").write_text(json.dumps(manifest, indent=2))
    print(f"\narchived to {dest}/")


def show(scored):
    print("  verdict                      raw    binned  cov   n  rule")
    for r in scored:
        st = r.get("effect_strat")
        raw = f"{r['effect']:+.2f}" if "effect" in r else "  n/a"
        st = "  n/a" if st is None else f"{st:+.2f}"
        cov = f"{r.get('strat_coverage', 0):.0%}"
        print(f"  {r['verdict']:26} {raw:>6} {st:>6} {cov:>4} {r['support']:>3}  {r['text']}")


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--dataset", default="configs/esol.json",
                   help="dataset config (see src/config.py)")
    p.add_argument("--n-bases", type=int, default=10)
    p.add_argument("--n-rules", type=int, default=8)
    p.add_argument("--num-samples", type=int, default=1500)
    p.add_argument("--threshold", type=float, default=None,
                   help="override the config's cf_threshold (target units)")
    p.add_argument("--noise-floor", type=float, default=None,
                   help="override the config's noise_floor (target units)")
    p.add_argument("--min-coverage", type=float, default=0.5,
                   help="min fraction of matching molecules in usable control bins "
                        "for a rule to be testable beyond the control")
    p.add_argument("--run-name", default=None,
                   help="archive this run to results/<dataset>/<run-name>/")
    p.add_argument("--refresh-cf", action="store_true")
    p.add_argument("--skip-llm", action="store_true")
    p.add_argument("--revise", action="store_true")
    main(p.parse_args())
