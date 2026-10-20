"""The black box being explained.

Deliberately boring. Counterfactuals explain the MODEL, not nature, so the
only thing that matters here is that you can state its held-out performance
on slide one. If the model is bad, every explanation downstream is an
explanation of a bad model, and a judge will ask.
"""
from __future__ import annotations

import numpy as np
from rdkit import Chem, DataStructs
from rdkit.Chem import rdFingerprintGenerator
from sklearn.ensemble import RandomForestRegressor
from sklearn.metrics import mean_absolute_error, r2_score

_GEN = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)


def ecfp(smiles: str) -> np.ndarray:
    """ECFP4 / Morgan fingerprint. 2048 bits, each bit = 'some circular
    substructure of radius<=2 hashed here'. Bits collide; a bit is not a
    functional group. Do not let the LLM pretend otherwise."""
    m = Chem.MolFromSmiles(smiles)
    if m is None:
        return np.zeros(2048, dtype=np.int8)
    fp = _GEN.GetFingerprint(m)
    arr = np.zeros((2048,), dtype=np.int8)
    DataStructs.ConvertToNumpyArray(fp, arr)
    return arr


def featurise(smiles_list) -> np.ndarray:
    return np.stack([ecfp(s) for s in smiles_list])


class BlackBox:
    def __init__(self, **kw):
        self.model = RandomForestRegressor(
            n_estimators=kw.pop("n_estimators", 500),
            n_jobs=-1,
            random_state=kw.pop("random_state", 0),
            **kw,
        )

    def fit(self, df):
        self.model.fit(featurise(df["smiles"]), df["y"].values)
        return self

    def predict(self, smiles_list):
        return self.model.predict(featurise(smiles_list))

    def __call__(self, smiles):
        """Signature exmol.sample_space expects (one SMILES in, scalar out)."""
        return float(self.predict([smiles])[0])

    def report(self, df, name="test"):
        p = self.predict(df["smiles"])
        r2, mae = r2_score(df["y"], p), mean_absolute_error(df["y"], p)
        print(f"[{name}] n={len(df)}  R2={r2:.3f}  MAE={mae:.3f}")
        return {"split": name, "n": len(df), "r2": r2, "mae": mae}
