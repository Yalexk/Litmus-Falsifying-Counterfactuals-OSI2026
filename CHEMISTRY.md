# Chemistry you need, and nothing more

Enough to not say something false on stage. This is not chemistry education,
it is the minimum model of the domain required to run this pipeline honestly.

## 1. How a molecule becomes data

A molecule is a graph: atoms are nodes, bonds are edges. Three text encodings
matter here.

**SMILES** — the standard. `CCO` is ethanol, `c1ccccc1` is benzene (lowercase
means aromatic), `CC(=O)O` is acetic acid. Branches in parentheses, ring
closures as digits. The same molecule has many valid SMILES strings, so you
*must* canonicalise before deduplicating or joining datasets. `data.clean()`
does this. Skipping it is a classic silent duplicate-leakage bug.

**SELFIES** — every syntactically valid string maps to a valid molecule.
This is why MMACE/STONED perturb SELFIES rather than SMILES: random edits to
SMILES mostly produce unparseable junk. Important caveat: *valid* does not
mean *stable*, *synthesisable*, or *sensible*. Some of your counterfactuals
will be molecules no chemist would ever make. Show one of these in your demo
rather than hiding it — it is evidence that the validation stage is needed.

**SMARTS** — a query language for substructures, like regex for molecules.
`[CX3](=O)[OX2H1]` matches a carboxylic acid. This is how a verbal hypothesis
becomes an executable predicate. See `GROUP_SMARTS` in `src/rules.py`.

## 2. Featurisation

**ECFP / Morgan fingerprints** (what the model uses): for each atom, hash the
circular neighbourhood out to radius *r*, set the corresponding bit in a
fixed-length vector. ECFP4 = radius 2, usually 2048 bits.

The thing to internalise: **bits collide and a bit is not a functional group.**
Several unrelated substructures land on the same bit. So you cannot read a
feature importance off this model and call it chemistry. That gap is precisely
why counterfactuals exist — and why your rule layer reasons over named
descriptors and SMARTS rather than over bits.

**Descriptors** (what the rules use) — interpretable scalars:

| Descriptor | Meaning | Rough intuition |
|---|---|---|
| `MolWt` | molecular weight (Da) | bigger is usually less soluble |
| `MolLogP` | octanol/water partition, computed | high = greasy/lipophilic = poorly water-soluble |
| `TPSA` | topological polar surface area (Å²) | polarity proxy; >~90 Å² ≈ won't cross blood-brain barrier |
| `NumHDonors` | OH / NH groups | H-bonds to water → more soluble |
| `NumHAcceptors` | N / O atoms | same direction |
| `NumRotatableBonds` | flexibility | >10 correlates with poor oral absorption |
| `NumAromaticRings` | flat ring count | flat molecules stack in crystals → dissolve badly |
| `FractionCSP3` | fraction of sp3 (non-flat) carbons | higher = more 3D = often more soluble |

**Lipinski's rule of five** (MW<500, logP<5, HBD≤5, HBA≤10) is the folk
heuristic these descriptors come from. Know it exists; it is a rough guide for
oral drugs, not a law, and quoting it as one will get you corrected.

## 3. Functional groups worth recognising on sight

| Group | SMILES fragment | Effect on water solubility |
|---|---|---|
| hydroxyl | `-OH` | ↑ strongly (H-bond donor *and* acceptor) |
| carboxylic acid | `-C(=O)OH` | ↑↑ (ionises to `-COO⁻` at pH 7) |
| amine | `-NH2` | ↑ (protonates to `-NH3⁺` at pH 7) |
| amide | `-C(=O)N-` | ↑ moderately |
| ester | `-C(=O)O-C` | ↓ vs the acid it came from — the ionisable H is gone |
| ether | `-O-` | weak acceptor, small effect |
| nitro | `-NO2` | polar but often ↓ solubility; also a tox-alert group |
| halogen | `-F -Cl -Br -I` | ↓ (adds lipophilicity and mass) |
| aromatic ring | `c1ccccc1` | ↓ (flat, greasy) |
| alkyl chain | `-CH2CH2CH3` | ↓ per carbon added |

**Ionisation is the one you'll miss.** A carboxylic acid is ~1000× more
soluble as its anion at physiological pH than as the neutral form. Descriptors
like `MolLogP` are computed for the *neutral* species, so a rule phrased on
logP alone silently ignores charge state. If the LLM proposes something about
acids, this is the first place to look for why it fails.

## 4. The properties

**Aqueous solubility (ESOL/Delaney)** — target is `logS`, log₁₀ of solubility
in mol/L, typically −8 to 0. More negative = less soluble. **Experimental
error is roughly 0.5–0.7 log units.** This number sets `--noise-floor`. Any
"effect" you report below it is not a finding.

**Lipophilicity (logD)** — partition between octanol and water at pH 7.4.
Measured, unlike computed `MolLogP`.

**BBB permeability (BBBP)** — binary, does the compound cross into the brain.
Strongly driven by TPSA and HBD count.

**Potency (BACE, IC50)** — concentration giving 50% inhibition. Usually
modelled as `pIC50 = −log₁₀(IC50 in M)`, so **higher pIC50 = more potent**.
Sign errors here are the single most common way to embarrass yourself.

## 5. Two concepts that will decide whether your evaluation is honest

**Scaffold splits.** The Bemis–Murcko scaffold is the ring system of a
molecule with side chains stripped. Chemical datasets are full of analogue
series — thirty variations on one core. A random split puts near-copies in
both train and test, and your R² is then measuring memorisation. Split by
scaffold (`data.scaffold_split`) and expect your numbers to drop noticeably.
That drop is the truth. Report it.

**Activity cliffs.** Pairs of very similar molecules with very different
property values. This is exactly what a counterfactual is, which cuts both
ways: it means counterfactuals surface real, interesting chemistry, and it
means a model with poor coverage of cliffs will produce counterfactuals that
are artefacts of its own smoothness rather than facts about molecules.

**Applicability domain.** Your model is only meaningful near its training
distribution. STONED happily wanders outside it. A counterfactual sitting in
a region where the model has seen nothing is a hallucination with a confidence
score attached. Worth logging the Tanimoto similarity of each counterfactual
to its nearest training molecule — a cheap, defensible check nobody else at
the hackathon will have bothered with.

## 6. The sentence to repeat whenever asked

> These are explanations of **the model**, not statements about chemistry —
> which is exactly why every one of them is tested on held-out data before we
> report it.

## Further reading, in priority order

1. [dmol.pub](https://dmol.pub) — White's free *Deep Learning for Molecules
   and Materials*. Chapters on representations and on XAI are the directly
   relevant ones.
2. MMACE paper: Wellawatte, Seshadri & White, *Chem. Sci.* 2022,
   [10.1039/D1SC05259D](https://doi.org/10.1039/D1SC05259D).
3. [exmol docs](https://ur-whitelab.github.io/exmol/).
4. XpertAI (XAI + LLM for structure–property relationships):
   [Commun. Chem. 2024](https://doi.org/10.1038/s42004-024-01393-y) — read so
   you can say precisely how your work differs from it, because someone will ask.
