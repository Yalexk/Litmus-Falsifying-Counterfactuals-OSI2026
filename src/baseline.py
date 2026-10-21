"""Known-chemistry baseline: textbook rules through the same harness.

This is the yardstick AND the bug detector. Textbook rules must survive the
raw test; if they don't, the harness is broken. How many survive beyond the
control descriptor is the comparison the LLM rules are judged against.

The rules live in each dataset's config ("baseline"), fixed a priori from
the literature with a "ref" each. A condition value of "train_median" is
replaced by that descriptor's median over the TRAIN split, so nothing is
tuned on test; the rule text may show it via "{value}" format fields.
"""
from __future__ import annotations

import numpy as np

from .counterfactuals import descriptors
from .rules import Rule, validate


def _median(train, name: str) -> float:
    return float(np.median([descriptors(s)[name] for s in train["smiles"]]))


def baseline_rules(cfg: dict, train) -> list[Rule]:
    rules = []
    for b in cfg["baseline"]:
        cond = dict(b["condition"])
        value = None
        if cond.get("value") == "train_median":
            value = cond["value"] = _median(train, cond["name"])
        text = b["text"].format(value=value) if value is not None else b["text"]
        r = Rule(text, cond, b["direction"], {"stage": "baseline", "ref": b["ref"]})
        err = validate(r)
        if err:
            raise ValueError(f"{cfg['name']} baseline rule {b['text']!r}: {err}")
        rules.append(r)
    return rules
