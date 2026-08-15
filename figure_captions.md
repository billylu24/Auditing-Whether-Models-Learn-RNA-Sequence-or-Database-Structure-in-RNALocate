# Figure captions and accessible descriptions

These are provisional general-publication figures. No target-journal or submission-phase compliance is claimed.

## Figure 1 | RNALocate database coverage is structured by context and same-gene cross-context support

**Caption.** (A) Fraction of candidate gene pairs with at least one retained RNALocate record in each of 16 contexts. The dashed line is the unweighted mean across contexts. (B) Test-set record fraction by the number of records for the same gene in the other 15 contexts. Shading shows Wilson 95% confidence intervals for binomial proportions. Support zero is marked as a cohort-selection boundary: every eligible gene has at least one retained record, so a gene with zero records in all other contexts must have a record in the candidate context. (C) Context-specific record fractions across support values; gray cells are combinations absent from the test set. These panels describe database coverage co-occurrence, not biological localization or a causal dose response.

**Alt text.** Three panels show substantial variation in RNALocate record coverage. Context fractions range from roughly one third to nearly all candidate genes. Except for the structurally constrained support-zero point, record availability generally rises as the same gene is recorded in more other contexts, with visible differences among contexts.

## Figure 2 | Context and same-gene support explain most predictable record availability

**Caption.** (A) Held-out average precision for the prevalence baseline and all seven non-empty combinations of context, target-excluding other-context support, and low-level sequence features. Error bars are 95% intervals from 1,000 bootstrap resamples of all-sequence graph components. The displayed nonzero axis is used for point estimates and intervals, not lengths. (B) Paired AP increments from the same bootstrap replicates. (C) Average precision from 200 controls in which support was shuffled among genes separately within every context and split, compared with the real same-gene context-plus-support model. Shuffling preserves context-specific support distributions but breaks the gene-support linkage. It does not identify measurement, publication, or curation mechanisms.

**Alt text.** The complete model has AP about 0.958. Context plus support is nearly identical, whereas context plus sequence is about 0.882. Adding sequence after context plus support changes AP by less than 0.001. All 200 shuffled-support models cluster near 0.811, far below the real context-plus-support result near 0.957.

## Figure 3 | Coverage prediction is heterogeneous but robust to tested low-level sequence views

**Caption.** (A) Prevalence and context-plus-support AP for each context. Connecting segments visualize AP lift without treating it as a causal effect. (B) Per-context AP difference after adding the complete low-level sequence vector to context plus support; green and vermillion encode positive and negative values with direction also represented by position relative to zero. (C) Equal-frequency reliability bins for context, context plus support, and the full model; vertical intervals are Wilson 95% confidence intervals for observed fractions. (D) Sensitivity to composition-plus-length, normalized 3-mer-only, and complete 70-feature sequence views in sequence-only and full models. The tested low-level views do not represent all possible RNA sequence models.

**Alt text.** Context-plus-support AP exceeds prevalence in every context, although performance is much weaker in IMR-90 than elsewhere. Adding sequence produces small and context-dependent changes. Context plus support and the full model track the calibration diagonal more closely than context alone. Stronger low-level sequence views improve sequence-only AP but barely alter the full model.

## Figure 4 | Database structure has a clear consequence for localization annotation recovery

**Caption.** (A) Macro average precision for recovering four recorded localization labels among already observed gene-context rows. Error bars are existing 95% graph-component bootstrap intervals. (B) Per-label average precision; the high nucleus prevalence and very low endoplasmic-reticulum prevalence make endpoint-specific results essential. (C) Mean average precision across the five eligible non-nucleus context-label cells with at least ten positive and ten negative test rows. A context-only model is constant within a context, whereas low-level sequence can rank genes. (D) Context and low-level-sequence macro AP under the primary exhaustive sequence graph and three MMseqs2 clustering stress tests. These analyses concern database annotations and do not validate localization for unrecorded pairs.

**Alt text.** Context strongly outperforms low-level sequence in aggregate annotation recovery, while combining both is best among the four simple models. Results differ sharply by localization label. Within fixed contexts, sequence improves mean AP over the constant context baseline. Context dominance is stable across four tested sequence-component partitions.
