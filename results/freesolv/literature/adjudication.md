# Data vs literature

Data: verdict on the held-out test split. Literature: blind judge over the top-5 Valency abstracts per claim, stances kept only with a verified verbatim quote. **DISAGREE** rows are the interesting ones.

| run | family | data verdict | data says | raw | binned | literature | agree? | claim | papers |
|---|---|---|---|---|---|---|---|---|---|
| final | textbook | untestable_beyond_control | supported | -5.22 | -1.82 | no_evidence | n/a | Molecules with polar surface area above the training median (17.1 A^2) have more favourable hydration free energy. |  |
| final | textbook | survives_beyond_control | supported | -3.98 | -1.31 | supported | agree | Molecules with at least one hydrogen-bond donor have more favourable hydration free energy. | [pubmed:29561152](https://pubmed.ncbi.nlm.nih.gov/29561152/) supports |
| final | textbook | survives_beyond_control | supported | -3.32 | -1.20 | no_evidence | n/a | Molecules containing a hydroxyl group have more favourable (more negative) hydration free energy. |  |
| final | textbook | untestable_beyond_control | supported | -3.35 | +1.02 | supported | agree | Molecules with two or more hydrogen-bond acceptors have more favourable hydration free energy. | [chemrxiv:10.26434/chemrxiv-2022-0mzxq](https://doi.org/10.26434/chemrxiv-2022-0mzxq) supports |
| final | llm | control_only | not supported beyond control | -2.02 | +0.40 | no_evidence | n/a | Molecules containing a primary amine not attached to a carbonyl have more favorable hydration free energy than molecules without a primary amine at matched TPSA. |  |
| final | llm | survives_beyond_control | supported | +2.51 | +1.12 | no_evidence | n/a | Molecules with MolLogP above 3 have less favorable hydration free energy than molecules below this threshold, independently of polar surface area. |  |
| final | llm | contradicted | contradicted | -1.33 | +0.15 | contradicted | agree | Molecules with more than one rotatable bond have less favorable hydration free energy than more rigid molecules of the same polar surface area. | [pubmed:18171044](https://pubmed.ncbi.nlm.nih.gov/18171044/) contradicts |
| final | llm | insufficient_support | unclear |  |  | not judged | n/a | Molecules containing an aromatic N-H group have more favorable hydration free energy than molecules without one, even at matched TPSA. |  |
| final | llm | insufficient_support | unclear |  |  | not judged | n/a | Molecules with FractionCSP3 above 0.8 are less favorably hydrated than less-saturated molecules of the same TPSA. |  |
| final | llm | insufficient_support | unclear |  |  | not judged | n/a | Molecules containing an N-halogen bond have less favorable hydration free energy than N-alkyl counterparts of identical polar surface area. |  |
| final | llm | insufficient_support | unclear |  |  | not judged | n/a | Molecules containing a thioketone or thione C=S group have less favorable hydration free energy than corresponding oxo-carbonyl compounds at matched TPSA. |  |
| final | llm | insufficient_support | unclear |  |  | not judged | n/a | Molecules containing a nitrile group have less favorable hydration free energy than other molecules of comparable TPSA. |  |

agree 3 · DISAGREE 0 · n/a 9
