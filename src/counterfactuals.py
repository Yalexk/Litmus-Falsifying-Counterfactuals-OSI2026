"""exmol / MMACE wrapper: base molecule -> nearby molecules with a flipped
prediction, serialised to JSON for the hypothesis stage.

MMACE: Wellawatte, Seshadri & White, Chem. Sci. 2022, DOI 10.1039/D1SC05259D.
It perturbs SELFIES strings via STONED, so every candidate is syntactically
valid -- but valid is not the same as synthesisable or sensible. Expect some
counterfactuals to be chemical nonsense. Keeping a few of those in the demo
is honest and makes the validation stage look necessary.
"""
from __future__ import annotations

import json
from pathlib import Path

import exmol
from rdkit import Chem
from rdkit.Chem import Descriptors


def counterfactuals_for(
    base_smiles: str,
    model,
    threshold: float = 1.0,
    num_samples: int = 2000,
    max_cf: int = 5,
):
    """exmol 3.3's cf_explain treats a counterfactual as any sample whose
    yhat differs from the base's. That is a classification assumption, so for
    a regression target we hand it a DISCRETISED label (-1/0/+1) relative to
    the base prediction, and keep the real values ourselves.

    `threshold` is in target units. For logS, keep it at or above experimental
    noise (~0.6 log units) or you are explaining measurement error.
    """
    base_yhat = model(base_smiles)
    real: dict[str, float] = {}

    def labelled(smiles, _selfies=None):
        v = model(smiles)
        real[smiles] = v
        if v >= base_yhat + threshold:
            return 1
        if v <= base_yhat - threshold:
            return -1
        return 0

    samples = exmol.sample_space(
        base_smiles, labelled, batched=False, num_samples=num_samples, quiet=True
    )
    cfs = exmol.cf_explain(samples, nmols=max_cf)
    base = cfs[0]
    out = []
    for c in cfs[1:]:
        v = real.get(c.smiles, model(c.smiles))
        out.append(
            {
                "smiles": c.smiles,
                "yhat": float(v),
                "similarity": float(c.similarity),
                "delta": float(v - base_yhat),
                "descriptors": descriptors(c.smiles),
            }
        )
    return {
        "base_smiles": base.smiles,
        "base_yhat": float(base_yhat),
        "base_descriptors": descriptors(base.smiles),
        "counterfactuals": out,
    }


def descriptors(smiles: str) -> dict:
    """The interpretable handles the LLM is allowed to reason over.

    MolLogP      lipophilicity (octanol/water). High = greasy, low solubility.
    TPSA         topological polar surface area, A^2. >~90 kills BBB entry.
    NumHDonors   OH / NH counts. More H-bonding to water = more soluble.
    NumHAcceptors  N / O counts.
    MolWt, NumRotatableBonds, NumAromaticRings, FractionCSP3 (sp3 fraction;
    flat aromatic molecules stack in crystals and dissolve badly).
    """
    m = Chem.MolFromSmiles(smiles)
    if m is None:
        return {}
    return {
        "MolWt": Descriptors.MolWt(m),
        "MolLogP": Descriptors.MolLogP(m),
        "TPSA": Descriptors.TPSA(m),
        "NumHDonors": Descriptors.NumHDonors(m),
        "NumHAcceptors": Descriptors.NumHAcceptors(m),
        "NumRotatableBonds": Descriptors.NumRotatableBonds(m),
        "NumAromaticRings": Descriptors.NumAromaticRings(m),
        "FractionCSP3": Descriptors.FractionCSP3(m),
        "RingCount": Descriptors.RingCount(m),
    }


def run_batch(base_smiles_list, model, out_path="out/counterfactuals.json", **kw):
    results = []
    for i, s in enumerate(base_smiles_list, 1):
        try:
            r = counterfactuals_for(s, model, **kw)
            if r["counterfactuals"]:
                results.append(r)
            print(f"  [{i}/{len(base_smiles_list)}] {s} -> {len(r['counterfactuals'])} cfs")
        except Exception as e:  # STONED occasionally chokes; keep going
            print(f"  [{i}] {s} FAILED: {e}")
    Path(out_path).parent.mkdir(parents=True, exist_ok=True)
    Path(out_path).write_text(json.dumps(results, indent=2))
    return results
