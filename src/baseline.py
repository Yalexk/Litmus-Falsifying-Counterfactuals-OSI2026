"""Known-chemistry baseline: textbook solubility rules through the same harness.

This is the yardstick AND the bug detector. Textbook rules must survive the
raw test; if they don't, the harness is broken. How many survive beyond LogP
is the comparison the LLM rules are judged against.

Rules and directions are fixed a priori from the literature, not from any
result in this repo:
  Delaney 2004 (ESOL), J. Chem. Inf. Comput. Sci. 44, 1000:
    logS = 0.16 - 0.63 cLogP - 0.0062 MW + 0.066 RB - 0.74 AP
    -> LogP down, MW down, aromatic proportion down.
  Lovering et al. 2009, J. Med. Chem. 52, 6752 ("Escape from Flatland"):
    higher Fsp3 -> higher solubility.
  H-bond donors interact with water -> higher solubility (general textbook).
Continuous thresholds are TRAIN-split medians, so nothing is tuned on test.
"""
from __future__ import annotations

import numpy as np

from .counterfactuals import descriptors
from .rules import Rule


def _median(train, name: str) -> float:
    return float(np.median([descriptors(s)[name] for s in train["smiles"]]))


def textbook_rules(train) -> list[Rule]:
    def d(name, op, value):
        return {"kind": "descriptor", "name": name, "op": op, "value": value}

    logp, mw, fsp3 = (_median(train, n) for n in ("MolLogP", "MolWt", "FractionCSP3"))
    src = {"stage": "baseline"}
    return [
        Rule(f"Molecules with LogP above the training median ({logp:.2f}) are less soluble.",
             d("MolLogP", ">", logp), "lower", {**src, "ref": "Delaney 2004"}),
        Rule(f"Molecules heavier than the training median ({mw:.0f} Da) are less soluble.",
             d("MolWt", ">", mw), "lower", {**src, "ref": "Delaney 2004"}),
        Rule("Molecules with two or more aromatic rings are less soluble.",
             d("NumAromaticRings", ">=", 2), "lower", {**src, "ref": "Delaney 2004"}),
        Rule(f"Molecules with sp3 fraction above the training median ({fsp3:.2f}) "
             "are more soluble.",
             d("FractionCSP3", ">", fsp3), "higher", {**src, "ref": "Lovering 2009"}),
        Rule("Molecules with at least one hydrogen-bond donor are more soluble.",
             d("NumHDonors", ">=", 1), "higher", {**src, "ref": "textbook"}),
    ]
