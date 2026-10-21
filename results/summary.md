# Cross-dataset summary

LLM and textbook rules scored on each dataset's held-out scaffold test split. 'beyond control' = survives within bins of the dataset's control descriptor at >=50% coverage.

| dataset | run | control | test n | test R² | LLM raw | LLM beyond control | untestable | contradicted | revised | duplicates | textbook raw | textbook beyond control | pre-registered |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| esol | final-v2-revise | MolLogP | 169 | 0.11 | 6/8 (75%) | 0/8 (0%) | 0 | 0 | 4 | 1 | 5/5 (100%) | 0/5 (0%) | no |
| esol | run1-5bases | MolLogP | 169 | 0.11 | 8/8 (100%) | 1/8 (12%) | 4 | 0 | 0 | 0 | 5/5 (100%) | 0/5 (0%) | no |
| esol | run2-20bases | MolLogP | 169 | 0.11 | 6/8 (75%) | 1/8 (12%) | 1 | 2 | 2 | 0 | 5/5 (100%) | 0/5 (0%) | no |
| esol | run3-prompt-v2 | MolLogP | 169 | 0.11 | 5/8 (62%) | 1/8 (12%) | 0 | 0 | 0 | 0 | 5/5 (100%) | 0/5 (0%) | no |
| lipophilicity | r1-v2-revise | MolLogP | 630 | 0.33 | 5/8 (62%) | 5/8 (62%) | 0 | 0 | 3 | 1 | 3/4 (75%) | 2/4 (50%) | yes |
