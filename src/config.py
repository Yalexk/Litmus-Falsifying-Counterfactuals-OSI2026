"""Dataset configs: everything dataset-specific lives in configs/<name>.json.

A config names the data (url or path, SMILES and target columns), describes
the property for the prompts, and fixes the analysis choices that must not
be tuned after seeing results: noise floor, counterfactual threshold,
control descriptor, textbook baseline.

Pre-registration: for a new dataset, commit its config BEFORE the first run
and set "preregistered": true. Each run's manifest records the commit, which
is what proves the order.
"""
from __future__ import annotations

import json
from pathlib import Path

REQUIRED = {
    "name": str, "preregistered": bool, "source": dict, "smiles_col": str,
    "target_col": str, "property": str, "units": str, "unit_word": str,
    "higher_means": str, "noise_floor": (int, float), "noise_floor_source": str,
    "cf_threshold": (int, float), "control": (dict, type(None)), "baseline": list,
}
CONTROL_REQUIRED = ("descriptor", "name", "test_name", "proxies", "rationale", "ref")


def load(path: str) -> dict:
    cfg = json.loads(Path(path).read_text())
    for k, t in REQUIRED.items():
        if k not in cfg:
            raise ValueError(f"{path}: missing {k!r}")
        if not isinstance(cfg[k], t):
            raise ValueError(f"{path}: {k!r} should be {t}, got {type(cfg[k]).__name__}")
    if not (cfg["source"].get("url") or cfg["source"].get("path")):
        raise ValueError(f"{path}: source needs a url or a path")
    if cfg["control"] is not None:
        from .counterfactuals import descriptors
        missing = [k for k in CONTROL_REQUIRED if k not in cfg["control"]]
        if missing:
            raise ValueError(f"{path}: control missing {missing}")
        if cfg["control"]["descriptor"] not in descriptors("CCO"):
            raise ValueError(f"{path}: unknown control descriptor "
                             f"{cfg['control']['descriptor']!r}")
    for i, b in enumerate(cfg["baseline"]):
        for k in ("text", "condition", "direction", "ref"):
            if k not in b:
                raise ValueError(f"{path}: baseline[{i}] missing {k!r}")
    cfg["_path"] = str(path)
    return cfg


def property_line(cfg: dict) -> str:
    """One-line property description handed to the LLM."""
    return f"{cfg['property']}, {cfg['units']} (higher = {cfg['higher_means']})"
