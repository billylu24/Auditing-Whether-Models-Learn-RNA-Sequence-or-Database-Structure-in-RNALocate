# Frozen Analysis Plan: RNALocate Context-Shortcut Audit

**Frozen before the new high-school-paper analyses were run:** 2026-08-14  
**Status:** supplemental post hoc analysis plan; existing upstream results were already known  
**Audience:** general scientific and high-school research readers  
**Publisher-specific formatting:** not yet selected or audited

## Research question

To what extent can apparent performance in RNALocate-based mRNA localization
annotation recovery be explained by cellular context and database-coverage
structure, and how much incremental information is supplied by low-level RNA
sequence features?

The study concerns database records. It does not treat an absent record as a
verified biological negative and does not estimate a causal effect of context,
sequence, measurement, publication, or curation.

## Analysis A: record availability

### Population and outcome

- Population: the existing panel-complete grid of 18,753 eligible genes crossed
  with 16 retained contexts.
- Outcome: whether a gene-context pair contains at least one retained RNALocate
  target-location record.
- Split: the existing all-sequence MMseqs2 graph-component train, validation,
  and test assignment.
- Other-context support: the count of contexts with a retained record for the
  same gene, excluding the candidate context.
- Sequence representation: the existing 70-feature vector containing 64
  normalized 3-mer frequencies, five nucleotide/ambiguous-base fractions, and
  normalized log transcript length.

### Primary model matrix

Fit an intercept-only prevalence baseline and all seven non-empty combinations
of context (C), other-context support (S), and sequence (Q): C, S, Q, C+S,
C+Q, S+Q, and C+S+Q. Use the existing logistic-regression tuning grid
`C = {0.01, 0.1, 1, 10}` selected by validation average precision (AP), refit
on train plus validation, and evaluate once on held-out components.

Primary metric: AP. Secondary metrics: AUROC and Brier score. Use 1,000 paired
component-bootstrap replicates for model intervals and the following AP
contrasts:

1. C+Q minus C;
2. C+S minus C;
3. C+S+Q minus C+S;
4. C+S+Q minus C+Q.

### Coverage-gradient analysis

Report the observed test-set fraction and Wilson interval for each support
value from 0 to 15, together with its denominator. Display context-specific
fractions where estimable. Support zero is a selection boundary: because every
eligible gene has at least one retained record, zero records in all other
contexts implies a record in the candidate context. It must be displayed and
explained rather than interpreted as a causal dose-response point.

### Shuffled-support negative control

Within each context and split, permute support among genes. Do not move values
between train, validation, and test. For each of 200 fixed-seed replicates,
repeat hyperparameter selection, refitting, and held-out evaluation for the
context-plus-support model. This tests whether predictive gain depends on the
same-gene linkage after preserving context-specific support and outcome
margins. It does not identify a curation or measurement mechanism.

### Context heterogeneity

For each of 16 contexts, report sample size, prevalence, AP, AP minus
prevalence, normalized AP lift `(AP - prevalence) / (1 - prevalence)`, AUROC,
and Brier score for the prevalence, context-only, context-plus-support, and
full models. Undefined metrics remain missing rather than being set to zero.

### Sequence-view sensitivity

Using the same split and model family, compare composition plus length,
3-mer-only, and the complete 70-feature representation in sequence-only and
context-plus-support-plus-sequence models. This bounds claims to tested
low-level sequence representations; it is not a test of all possible sequence
models.

### Calibration

For context-only, context-plus-support, and the full model, construct ten
equal-frequency test-set probability bins. Report predicted and observed
fractions, bin sizes, and Wilson intervals. Bins are descriptive and do not
certify calibration.

## Analysis B: consequence for annotation recovery

Reuse the existing held-out results for prevalence, low-level sequence,
context, and context-plus-sequence annotation-recovery models. Reuse the
existing within-context comparison among context prevalence, sequence, and
context-plus-sequence. These outcomes are recorded localization labels among
already observed rows, not record availability on the complete grid.

## Figure plan

1. Database coverage by context and other-context support.
2. Complete factorial model comparison, incremental effects, and shuffled-support control.
3. Context heterogeneity, calibration, and low-level sequence sensitivity.
4. Consequences for localization annotation recovery, including within-context evidence.

All figures will be exported as editable SVG, opaque-background PNG, and PDF,
with source tables, captions, alt text, hashes, dimensions, and software
versions. No target-journal compliance is claimed until a journal and
submission phase are selected.

## Interpretation gates

- If shuffled support retains the real-support performance, do not attribute
  the gain to same-gene linkage.
- If performance is concentrated in a small number of contexts, restrict the
  claim to those contexts rather than the selected panel as a whole.
- If an alternative low-level sequence view materially increases AP beyond
  context plus support, revise the claim that tested sequence increments are
  small.
- Never describe coverage predictability as biological localization accuracy,
  causation, or direct evidence of curator bias.
