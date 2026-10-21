"""Dataset loading and splitting.

Scaffold splitting matters more than you think. A random split puts close
analogues of the same molecule in both train and test, which inflates your
metric and makes the counterfactuals look better than they are.
"""
from __future__ import annotations

import io
import urllib.request
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem.Scaffolds import MurckoScaffold

RDLogger.DisableLog("rdApp.*")

ESOL_URL = (
    "https://deepchemdata.s3-us-west-1.amazonaws.com/datasets/delaney-processed.csv"
)


def load_dataset(cfg: dict) -> pd.DataFrame:
    """Load any dataset described by a config (see src/config.py).

    Reads the CSV from source.path if given, else source.url (cached under
    data/<name>.csv so reruns don't refetch), with any source.read_csv
    options passed to pandas, keeps the SMILES and target
    columns as smiles / y, then canonicalises and dedupes via clean().
    """
    src = cfg["source"]
    opts = src.get("read_csv", {})  # passed to pandas, e.g. sep / comment / names
    if src.get("path"):
        df = pd.read_csv(src["path"], **opts)
    else:
        cache = Path("data") / f"{cfg['name']}.csv"
        if not cache.exists():
            cache.parent.mkdir(exist_ok=True)
            # Some hosts (Harvard Dataverse) refuse urllib's default agent.
            req = urllib.request.Request(src["url"], headers={"User-Agent": "osi-falsifying-cf/0.1"})
            with urllib.request.urlopen(req, timeout=120) as r:
                cache.write_bytes(r.read())
        df = pd.read_csv(cache, **opts)
    missing = {cfg["smiles_col"], cfg["target_col"]} - set(df.columns)
    if missing:
        raise ValueError(f"{cfg['name']}: columns {missing} not in {list(df.columns)[:20]}")
    df = df.rename(columns={cfg["smiles_col"]: "smiles", cfg["target_col"]: "y"})
    df["y"] = pd.to_numeric(df["y"], errors="coerce")
    return clean(df[["smiles", "y"]])


def load_esol(path: str | None = None) -> pd.DataFrame:
    """ESOL / Delaney aqueous solubility. ~1128 molecules, target = log mol/L.

    Columns returned: smiles, y
    Experimental noise on logS is roughly 0.5-0.7 log units. Any 'effect'
    your pipeline reports below that is noise. Say so in the demo.
    """
    if path:
        df = pd.read_csv(path)
    else:
        with urllib.request.urlopen(ESOL_URL) as r:
            df = pd.read_csv(io.StringIO(r.read().decode()))
    df = df.rename(
        columns={
            "smiles": "smiles",
            "measured log solubility in mols per litre": "y",
        }
    )
    return clean(df[["smiles", "y"]])


def load_bigsoldb(path: str) -> pd.DataFrame:
    """BigSolDB 2.0 (download the CSV from Zenodo yourself, it's not a URL grab).

    Expect columns along the lines of SMILES / Solvent / Solubility / T.
    Rename to: smiles, solvent_smiles, y, temperature.

    The extra axis here is the point: you can generate counterfactuals over
    the SOLVENT as well as the solute, which MMACE does not do.
    """
    df = pd.read_csv(path)
    cols = {c.lower().strip(): c for c in df.columns}
    rename = {}
    for want, cands in {
        "smiles": ["smiles", "solute smiles", "smiles_solute"],
        "solvent_smiles": ["solvent smiles", "smiles_solvent", "solvent"],
        "y": ["solubility", "logs", "solubility_g_g", "solubility_mol_mol"],
        "temperature": ["t,k", "temperature", "t"],
    }.items():
        for c in cands:
            if c in cols:
                rename[cols[c]] = want
                break
    df = df.rename(columns=rename)
    missing = {"smiles", "y"} - set(df.columns)
    if missing:
        raise ValueError(
            f"Could not find {missing}. Columns present: {list(df.columns)[:20]}"
        )
    return clean(df)


def clean(df: pd.DataFrame) -> pd.DataFrame:
    """Canonicalise SMILES, drop anything RDKit refuses to parse, dedupe."""
    keep, canon = [], []
    for s in df["smiles"]:
        m = Chem.MolFromSmiles(s) if isinstance(s, str) else None
        keep.append(m is not None)
        canon.append(Chem.MolToSmiles(m) if m else None)
    df = df.loc[keep].copy()
    df["smiles"] = [c for c in canon if c is not None]
    df = df.dropna(subset=["smiles", "y"])
    df = df.groupby("smiles", as_index=False).first()
    return df.reset_index(drop=True)


def murcko_scaffold(smiles: str) -> str:
    try:
        return MurckoScaffold.MurckoScaffoldSmiles(smiles=smiles, includeChirality=False)
    except Exception:
        return ""


def scaffold_split(df: pd.DataFrame, frac_train=0.7, frac_valid=0.15, seed=0):
    """Deterministic scaffold split: largest scaffold groups go to train.

    Returns (train, valid, test). The test set is what your LLM-generated
    rules get scored on, and the LLM must never see it.
    """
    groups: dict[str, list[int]] = defaultdict(list)
    for i, s in enumerate(df["smiles"]):
        groups[murcko_scaffold(s)].append(i)
    sets = sorted(groups.values(), key=lambda g: (-len(g), g[0]))
    n = len(df)
    n_tr, n_va = int(frac_train * n), int(frac_valid * n)
    tr, va, te = [], [], []
    for g in sets:
        if len(tr) + len(g) <= n_tr:
            tr += g
        elif len(va) + len(g) <= n_va:
            va += g
        else:
            te += g
    rng = np.random.default_rng(seed)
    for part in (tr, va, te):
        rng.shuffle(part)
    return (
        df.iloc[tr].reset_index(drop=True),
        df.iloc[va].reset_index(drop=True),
        df.iloc[te].reset_index(drop=True),
    )


def spread_bases(df: pd.DataFrame, n: int, seed: int = 0) -> list[str]:
    """n base SMILES, one per distinct scaffold, scaffolds chosen at random.

    A plain random sample is dominated by the commonest scaffolds (benzene,
    acyclic), so the LLM would see the same chemistry over and over.
    """
    groups: dict[str, list[int]] = defaultdict(list)
    for i, s in enumerate(df["smiles"]):
        groups[murcko_scaffold(s)].append(i)
    rng = np.random.default_rng(seed)
    keys = sorted(groups)
    rng.shuffle(keys)
    return [df["smiles"].iloc[int(rng.choice(groups[k]))] for k in keys[:n]]


if __name__ == "__main__":
    d = load_esol()
    tr, va, te = scaffold_split(d)
    print(f"{len(d)} molecules -> train {len(tr)} valid {len(va)} test {len(te)}")
    print(d.head())
