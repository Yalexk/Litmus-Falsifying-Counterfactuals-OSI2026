"""Cross-dataset summary: one row per archived run.

    python -m src.summary                 # every run under results/
    python -m src.summary results/esol/final-v2-revise results/lipophilicity/r1

Writes results/summary.md. Counts come straight from each run's rules.json /
baseline.json, so the table can't drift from the artefacts.
"""
from __future__ import annotations

import json
import sys
from collections import Counter
from pathlib import Path

from .rules import LEGACY_VERDICTS


def _counts(rows: list[dict]) -> dict:
    v = Counter(LEGACY_VERDICTS.get(r["verdict"], r["verdict"]) for r in rows)
    n = len(rows)
    raw = sum(r.get("survives", False) for r in rows)
    return {"n": n, "raw": raw, "beyond": v["survives_beyond_control"],
            "untestable": v["untestable_beyond_control"],
            "contradicted": v["contradicted"]}


def _frac(k: int, n: int) -> str:
    return f"{k}/{n} ({k / n:.0%})" if n else "-"


def row(run: Path) -> str | None:
    if not (run / "rules.json").exists():
        return None
    man = json.loads((run / "manifest.json").read_text()) if (run / "manifest.json").exists() else {}
    test = next((m for m in json.loads((run / "model_metrics.json").read_text())
                 if m["split"] == "test"), {})
    rules = json.loads((run / "rules.json").read_text())
    llm = _counts(rules)
    base = _counts(json.loads((run / "baseline.json").read_text()))
    revised = sum("revised_from" in r for r in rules)
    return (f"| {run.parent.name} | {run.name} | {man.get('control_descriptor', 'MolLogP')} "
            f"| {test.get('n', '')} | {test.get('r2', float('nan')):.2f} "
            f"| {_frac(llm['raw'], llm['n'])} | {_frac(llm['beyond'], llm['n'])} "
            f"| {llm['untestable']} | {llm['contradicted']} | {revised} "
            f"| {_frac(base['raw'], base['n'])} | {_frac(base['beyond'], base['n'])} "
            f"| {'yes' if man.get('preregistered') else 'no'} |")


def main(paths: list[str]):
    runs = [Path(p) for p in paths] or sorted(
        p.parent for p in Path("results").glob("*/*/rules.json"))
    lines = ["# Cross-dataset summary", "",
             "LLM and textbook rules scored on each dataset's held-out scaffold test split. "
             "'beyond control' = survives within bins of the dataset's control descriptor "
             "at >=50% coverage.", "",
             "| dataset | run | control | test n | test R² | LLM raw | LLM beyond control "
             "| untestable | contradicted | revised | textbook raw | textbook beyond control "
             "| pre-registered |",
             "|---|---|---|---|---|---|---|---|---|---|---|---|---|"]
    lines += [r for r in map(row, runs) if r]
    Path("results/summary.md").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main(sys.argv[1:])
