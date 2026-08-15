# Support-zero exclusion sensitivity amendment

**Specified:** 2026-08-15, before running this sensitivity analysis

**Status:** reviewer-requested post hoc sensitivity analysis
**Primary analysis remains unchanged.**

## Rationale

Every gene in the selected panel has at least one retained RNALocate record. Therefore, a gene-context pair with other-context support equal to zero must be observed in the candidate context. This creates a deterministic cohort-selection boundary. The sensitivity analysis tests whether the main comparison persists after removing that boundary.

## Analysis set

Exclude every gene-context pair with `other_context_support == 0` from training, validation, and test sets. Keep the original gene-level sequence-component split assignments, outcome definition, sequence features, model family, regularization grid, and random seed.

## Models and outcomes

Refit four prespecified models:

1. context only (C);
2. other-context support only (S);
3. context plus support (C+S);
4. context plus support plus the complete 70-feature low-level sequence vector (C+S+Q).

Select logistic-regression `C` from `{0.01, 0.1, 1, 10}` using validation average precision (AP), refit on training plus validation data, and evaluate once on the reduced held-out test set. Report AP, AUROC, Brier score, and 95% intervals from 1,000 bootstrap resamples of sequence components. Report paired AP differences for C+S minus C and C+S+Q minus C+S.

## Interpretation gate

The main structural interpretation is supported only if C+S remains clearly above C after exclusion and the sequence increment after C+S remains small. This analysis addresses the deterministic support-zero boundary but does not remove other selection, measurement, publication, or curation processes.
