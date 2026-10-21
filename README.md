<h1 align="center">Falsifying Counterfactuals</h1>

<p align="center"><em>Do an LLM's explanations of a molecular property model survive contact with held-out data?</em></p>

When a model predicts a chemical property, counterfactual explanations show
which small structural changes would flip its prediction. An LLM can read
those and write confident rules: *"hydroxyl groups raise solubility"*,
*"carboxylic acids lower logD at pH 7.4"*. This project turns every such rule
into an executable test and checks it against measurements the LLM never saw,
then asks the published literature for a second, independent opinion.

What the pipeline does:

- **Explains a black-box property model** with counterfactuals (exmol), using
  only training molecules.
- **Has an isolated LLM propose rules**, each a sentence plus an executable
  condition (an RDKit descriptor threshold or a SMARTS pattern) and a
  direction.
- **Falsifies every rule on a held-out scaffold split**: raw effect, effect
  *beyond* the descriptor that already dominates the property, false-discovery
  control, and the dataset's experimental noise floor.
- **Benchmarks against textbook chemistry**: published rules go through the
  identical harness.
- **Cross-checks each claim against the literature**: Valency search, PubMed
  abstracts, and a blind judge whose quotes are verified word for word.
- **Runs on any molecular regression dataset** described by a pre-registered
  config file. Four ship with the repo.

<p align="center">
  <img src="docs/workflow.png" alt="Pipeline: explain the model, falsify on held-out data, cross-check against the literature" width="900">
</p>

---

## Goals

1. **One honest number per dataset**: what fraction of an LLM's confident
   structure–property claims survive falsification, and how many say anything
   beyond the obvious.
2. **The claims that fail, quoted verbatim**, next to the effect size that
   contradicts them.
3. **A harness, not a one-off study**: the same code and the same rules of
   evidence on every dataset, with analysis choices fixed before results exist.

### Where this sits

Hypothesis evaluation has been benchmarked in general domains: HypoBench
(arXiv 2504.11524) reports that the best model/method recovers 38.8% of
ground-truth hypotheses on synthetic tasks. Chemistry explanations have been
generated from XAI output and validated by expert judgement and
self-consistency ([XpertAI](https://doi.org/10.1038/s42004-024-01393-y)). This
project combines the two: counterfactual-derived chemistry hypotheses,
falsified against held-out measurements and cross-checked against literature.

---

## Results

**In short:** most of the LLM's claims pass a naive held-out test; far fewer
say anything beyond the descriptor that already dominates the property.
Asking it to revise weak rules made things worse: it reversed correct claims
into confident wrong ones with an invented mechanism. And the literature
check and the data check catch different failures.

### Headline numbers

One final run per dataset (prompt v2, 20 base molecules, one round of
revision). Each rule is tested on the held-out scaffold test split.

| Dataset | Test molecules | Model R² | LLM: holds raw | LLM: holds beyond control | LLM: contradicted | Textbook: holds beyond control |
|---|---|---|---|---|---|---|
| ESOL | 169 | 0.11 | 5/8 | 3/8 | 0 | 0/5 |
| Lipophilicity (logD 7.4) | 630 | 0.33 | 5/7 | **4/7** | 0 | 2/4 |
| FreeSolv | 148 | 0.65 | 2/8 | 1/8 | 1 | 2/4 |
| AqSolDB | 1,497 | 0.08 | 2/7 | 1/7 | **3** | 1/5 |

"Beyond control" means the effect survives within bins of the control
descriptor (calculated LogP; TPSA for FreeSolv), after Benjamini–Hochberg and
above the experimental noise floor.

- **How much a rule says beyond the obvious depends on the property.** logD
  has chemistry that calculated LogP misses (ionisation: carboxylic acids
  −1.77 log units at matched cLogP), and the LLM finds it. For solubility,
  most rules restate lipophilicity. Textbook rules follow the same pattern,
  so this reflects the chemistry rather than a quirk of the test.
- **LLM output varies between runs.** On ESOL, three development runs scored
  1/8, 0/8 and 3/8 beyond control. Treat each row as one sample, not a
  stable rate.
- **The LLM ignores what the data can test.** On FreeSolv, 5 of 8 rules
  matched too few test molecules to evaluate (thiones, N-halogens …).

### Exhibit 1 · Revision reverses correct claims

Shown validation feedback, the LLM rewrote correct rules as their opposites,
justified by crystal-lattice energy:

| Original (correct) | Revised | Effect on test |
|---|---|---|
| Primary amines → **more** soluble | "…**less** soluble, because NH₂ groups strengthen hydrogen bonding in the crystal lattice, raising lattice energy" | +1.05 log units, p = 4×10⁻¹², n = 175 (AqSolDB) |
| >3 H-bond donors → **more** soluble | "…**lower** solubility … high donor counts favour strong intermolecular H-bonding in the crystal lattice" | +1.08, p = 1×10⁻⁸, n = 181 (AqSolDB) |
| H-bond donors → more soluble | ">1 donor → lower solubility" | +1.61, p = 3×10⁻⁸ (ESOL development run 2) |

The originals are still scored on test for the record: both AqSolDB
originals were right. The literature agrees with the data: *"The
hydrogen-bond propensity of a compound always leads to an increase in
solubility"* ([Abraham & Le, 1999](https://pubmed.ncbi.nlm.nih.gov/10479348/)).
One caveat: the revise prompt explicitly allows proposing the opposite
direction. The lattice-energy mechanism is the model's own.

### Exhibit 2 · A rule that survives beyond the control

*"Molecules with a polar surface area below 10 Å² are less soluble than those
with higher TPSA at matched lipophilicity"* (ESOL): −3.70 raw, **−2.17 within
LogP bins**, n = 26. It appeared independently in two runs. The literature
supports TPSA as a solubility descriptor beyond logP in general, but nothing
at this threshold. It could be real chemistry, or Crippen LogP
under-predicting the hydrophobicity of hydrocarbons and halocarbons. That
question needs a chemist.

### Exhibit 3 · Literature and data disagree

Over the final runs the blind judge agreed with the data on 20 claims,
disagreed on 3, and had no usable evidence for 25; 166 of 169 quotes were
verified verbatim. The disagreements show why both checks are needed. The
judge rules on a claim's stated *mechanism*, the data on its *effect*:

- *"Zero H-bond donors → higher logD"* holds beyond cLogP in the data, but
  the literature rejects the stated reason (donors are less solvated than
  acceptors).
- *"Hydrazones/imines → less soluble, because planar groups stabilise
  crystal packing"* is contradicted by data (+0.74, p = 0.01) yet "supported"
  by literature, because the generic mechanism is well documented.

Per-dataset tables: `results/<dataset>/literature/adjudication.md`.

Full table: [`results/summary.md`](results/summary.md) (regenerate with
`python -m src.summary`). Every run keeps its counterfactuals, all rules with
their statistics, the textbook baseline, the exact system prompt, a copy of
the config and a manifest recording the git commit, under
`results/<dataset>/<run>/`.

---

## Quick start

**Needs:** Python 3.11+ and either an `ANTHROPIC_API_KEY` or the
[Claude Code](https://docs.anthropic.com/claude-code) CLI logged in. No other
keys.

```bash
git clone <this repo> && cd OSI-Falsifying-Counterfactuals
python -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt

# data, model, counterfactuals and textbook baseline; no LLM needed
python run_pipeline.py --dataset configs/esol.json --n-bases 20 --skip-llm

# the full loop, archived to results/esol/my-run/
python run_pipeline.py --dataset configs/esol.json --n-bases 20 --revise --run-name my-run

# one row per run, across datasets
python -m src.summary
```

Datasets download on first use and are cached in `data/`. Counterfactuals
are cached in `out/<dataset>/`, so rerunning the LLM stage takes minutes.

---

## How it works

### 1 · Explain the model

| Step | Module | What happens |
|---|---|---|
| Config | `src/config.py` | Loads and validates `configs/<name>.json` |
| Split | `src/data.py` | Canonicalise and dedupe SMILES; Murcko-scaffold split 70/15/15 |
| Model | `src/model.py` | Random forest on 2048-bit ECFP4, trained on train only |
| Counterfactuals | `src/counterfactuals.py` | exmol/STONED look-alikes for 20 training molecules, one per scaffold |
| Rules | `src/hypotheses.py` | The LLM proposes up to 8 rules; each is validated before scoring |

### 2 · Falsify on held-out data

`src/rules.py` scores each rule on the test split, which no LLM ever sees:

1. **Raw test**: matching vs non-matching molecules, Welch t-test, in the
   claimed direction.
2. **Beyond the control**: the same comparison within bins of the config's
   control descriptor (e.g. calculated LogP), with a within-bin permutation
   test. A rule that only restates the control shows nothing here.
3. **Rules of evidence**: Benjamini–Hochberg across the rule family, and the
   effect must exceed the experimental noise floor.

Each rule ends with one verdict:

| Verdict | Meaning |
|---|---|
| `survives_beyond_control` | Holds raw and at matched control descriptor |
| `control_only` | Holds raw, explained by the control |
| `untestable_beyond_control` | Holds raw, too collinear with the control to check |
| `fails` | No real effect in the claimed direction |
| `contradicted` | A real effect in the **opposite** direction |
| `insufficient_support` | Too few matching molecules to test |

With `--revise`, a rule with evidence against it on the validation split gets
**one** revision, from validation feedback only. The original is still scored
on test for the record, never shown to the LLM. Duplicate rules are set aside
so each claim counts once.

### 3 · Cross-check against the literature

`src/literature.py` runs as a separate stage:

```bash
python -m src.literature configs/<name>.json export results/<name>/<run>
#   retrieval: in an interactive Claude Code session with the Valency
#   connector, one semantic search per claim (query = claim text verbatim,
#   top 5), recorded in results/<name>/literature/hits.json
python -m src.literature configs/<name>.json fetch   # PubMed abstracts by ID
python -m src.literature configs/<name>.json judge   # blind, zero-tool judge
python -m src.literature configs/<name>.json table   # adjudication.md
```

The judge sees only the claim and the abstracts, never the data verdict, and
every stance must rest on a quote that appears verbatim in an abstract.
Disagreements between data and literature are the most interesting rows.

### What keeps it honest

- **Information barrier.** Counterfactuals come from training molecules; the
  test split reaches no LLM call.
- **Isolated LLM calls.** Without an API key the pipeline runs `claude -p` in
  an empty directory with `--tools "" --strict-mcp-config` and no MCP servers.
- **Pre-registration.** Noise floor, control descriptor and textbook rules are
  committed in the config before a dataset's first run; manifests record the
  config's first commit.
- **Scaffold split.** Test molecules have scaffolds the model never trained on.

---

## Datasets

| Config | Property | Molecules | Control | Noise floor | Pre-registered |
|---|---|---|---|---|---|
| `esol` | aqueous solubility (log mol/L) | 1,117 | MolLogP | 0.6 | no (development set) |
| `lipophilicity` | logD at pH 7.4 | 4,200 | MolLogP | 0.3 | yes |
| `freesolv` | hydration free energy (kcal/mol) | 624 | TPSA | 0.6 | yes |
| `aqsoldb` | aqueous solubility (log mol/L) | 9,980 | MolLogP | 0.6 | yes |

Every noise floor, control and textbook rule cites its source inside the
config.

### Add your own

Any CSV with a SMILES column and one continuous target:

1. Copy a file from `configs/` and set the source, columns and property
   description, then the analysis choices: `noise_floor` (with
   `noise_floor_source`), `cf_threshold` (the shipped configs use 5/3 × noise
   floor), `control` (or `null`), and `baseline` rules with references.
2. Set `"preregistered": true` and **commit the config before the first run**.
3. `python run_pipeline.py --dataset configs/<name>.json --skip-llm` checks the
   config, split and baseline.
4. Run the full loop, then the literature stage.

---

## Project structure

```
OSI-Falsifying-Counterfactuals/
├── configs/                 # one JSON per dataset: source, property, pre-registered choices
├── run_pipeline.py          # end-to-end run for one dataset
├── src/
│   ├── config.py            # config loading and validation
│   ├── data.py              # loading, cleaning, scaffold split, base selection
│   ├── model.py             # black-box predictor (RF on ECFP4)
│   ├── counterfactuals.py   # exmol wrapper, RDKit descriptors
│   ├── hypotheses.py        # isolated LLM: propose and revise
│   ├── rules.py             # the harness: evaluation, control bins, BH, verdicts
│   ├── baseline.py          # textbook rules from the config
│   ├── literature.py        # Valency -> PubMed -> blind judge -> adjudication
│   └── summary.py           # cross-dataset table
├── results/                 # committed artefacts per dataset and run
├── docs/                    # pipeline figures
└── CHEMISTRY.md             # the chemistry background needed to read the results
```

---

## Troubleshooting

- **`claude CLI failed` or a timeout during the LLM stage.** Calls get 20
  minutes and one automatic retry; a failed revision keeps the original rule.
  Rerun the same command: counterfactuals are cached.
- **A dataset download is refused.** Some hosts block default clients; the
  loader sends a user agent. Or download the file yourself and set
  `source.path` in the config.
- **A rule is `insufficient_support`.** It matched fewer than 10 test
  molecules (or nearly all of them). That is a finding about the claim's
  scope, not a crash.
- **Changing a pre-registered config.** Don't edit it after its first run;
  copy it under a new name, commit, and rerun.

---

## Limitations

- Continuous targets only.
- The predictor is deliberately simple and generalises modestly to unseen
  scaffolds. The rules explain *this model*, which is why every rule is
  tested against measurements rather than predictions.
- Public datasets may sit in any LLM's pretraining data; isolation stops the
  model reading the test split during a run, not remembering chemistry.
- AqSolDB likely overlaps ESOL; FreeSolv's validation split is small (40).

---

## Built with

- **Chemistry:** RDKit, exmol (STONED/SELFIES)
- **Modelling and statistics:** scikit-learn, NumPy, SciPy, pandas
- **LLM:** Claude, via the Anthropic SDK or the Claude Code CLI
- **Literature:** Valency search, NCBI E-utilities (PubMed)

## Acknowledgements

- Wellawatte, Seshadri & White for MMACE/exmol counterfactuals.
- The dataset authors: Delaney (ESOL), AstraZeneca via ChEMBL and MoleculeNet
  (Lipophilicity), Mobley & Guthrie (FreeSolv), Sorkun, Khetan & Er (AqSolDB).
- The literature behind each config's noise floor and textbook rules (cited
  inside the configs).
