# Falsifiable counterfactuals

Counterfactual explanations of a molecular property model are turned into
natural-language hypotheses by an LLM, converted into executable predicates,
and then **tested on held-out data the LLM never saw**. The headline output
is a number: what fraction of confident LLM structure–property claims actually
survive falsification.

Generating counterfactuals is a `pip install` ([exmol / MMACE](https://doi.org/10.1039/D1SC05259D)).
Generating prose explanations from them has been done ([XpertAI](https://doi.org/10.1038/s42004-024-01393-y)).
**Measuring whether those explanations are true has not.** That is this repo.

## Pipeline

```
dataset (SMILES, y)
  └─ scaffold split                     src/data.py
      └─ black-box model (ECFP + RF)    src/model.py
          └─ counterfactuals (exmol)    src/counterfactuals.py
              └─ LLM → hypotheses       src/hypotheses.py
                  └─ rule → predicate → held-out eval   src/rules.py
                      └─ BH correction + noise floor → scored table
```

## Setup

```bash
git clone <your repo> && cd cf-chem
./setup_git.sh                       # sets your git identity for this repo only
python -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt
export ANTHROPIC_API_KEY=...
```

## Run

```bash
python run_pipeline.py --n-bases 10 --skip-llm    # everything except the LLM stage
python run_pipeline.py --n-bases 20 --revise      # full loop
```

Outputs land in `out/`: `model_metrics.json`, `counterfactuals.json`,
`rules.json`. The counterfactual stage is cached — rerun the LLM stage freely
without regenerating it.

## Read this before writing code

[`CHEMISTRY.md`](CHEMISTRY.md) — the domain knowledge required to run this
honestly. Fifteen minutes, and it covers the three things that will otherwise
sink the submission: experimental noise floors, scaffold splits, and the fact
that these explain the *model*, not chemistry.

## Swapping the dataset

`src/data.py` ships an ESOL loader and a [BigSolDB 2.0](https://doi.org/10.5281/zenodo.15094979)
loader. BigSolDB gives a second counterfactual axis — perturb the *solvent*,
not just the solute, which MMACE does not do. Download the CSV from Zenodo to
`data/` and pass `--csv`. If it fights you in the first hour, fall back to
ESOL and move on.

## Known soft spots

- **SMARTS generation is the likely failure point.** Descriptor-based rules
  run reliably; LLM-authored SMARTS often parse but match the wrong thing.
  Everything is validated before evaluation, and rejected rules are logged
  rather than silently dropped — the rejection rate is itself a result.
- `--revise` allows one retry per failed rule. Resist turning it into a loop;
  unbounded revision just overfits the test set and destroys the headline number.
- No applicability-domain filter yet. Adding nearest-training-neighbour
  Tanimoto per counterfactual is the highest-value hour of work left.
