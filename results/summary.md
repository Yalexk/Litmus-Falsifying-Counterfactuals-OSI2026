# Cross-dataset summary

LLM and textbook rules scored on each dataset's held-out scaffold test split. 'beyond control' = survives within bins of the dataset's control descriptor at >=50% coverage.

| dataset | run | control | test n | test R² | LLM raw | LLM beyond control | untestable | contradicted | revised | duplicates | textbook raw | textbook beyond control | pre-registered |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| aqsoldb | final | MolLogP | 1497 | 0.08 | 2/7 (29%) | 1/7 (14%) | 0 | 3 | 4 | 1 | 4/5 (80%) | 1/5 (20%) | yes |
| esol | final | MolLogP | 169 | 0.11 | 5/8 (62%) | 3/8 (38%) | 0 | 0 | 2 | 0 | 5/5 (100%) | 0/5 (0%) | no |
| esol | final-v2-revise | MolLogP | 169 | 0.11 | 6/8 (75%) | 0/8 (0%) | 0 | 0 | 4 | 1 | 5/5 (100%) | 0/5 (0%) | no |
| esol | run1-5bases | MolLogP | 169 | 0.11 | 8/8 (100%) | 1/8 (12%) | 4 | 0 | 0 | 0 | 5/5 (100%) | 0/5 (0%) | no |
| esol | run2-20bases | MolLogP | 169 | 0.11 | 6/8 (75%) | 1/8 (12%) | 1 | 2 | 2 | 0 | 5/5 (100%) | 0/5 (0%) | no |
| esol | run3-prompt-v2 | MolLogP | 169 | 0.11 | 5/8 (62%) | 1/8 (12%) | 0 | 0 | 0 | 0 | 5/5 (100%) | 0/5 (0%) | no |
| freesolv | final | TPSA | 148 | 0.65 | 2/8 (25%) | 1/8 (12%) | 0 | 1 | 0 | 0 | 4/4 (100%) | 2/4 (50%) | yes |
| lipophilicity | final | MolLogP | 630 | 0.33 | 5/7 (71%) | 4/7 (57%) | 0 | 0 | 0 | 1 | 3/4 (75%) | 2/4 (50%) | yes |
| lipophilicity | r1-v2-revise | MolLogP | 630 | 0.33 | 5/8 (62%) | 5/8 (62%) | 0 | 0 | 3 | 1 | 3/4 (75%) | 2/4 (50%) | yes |
