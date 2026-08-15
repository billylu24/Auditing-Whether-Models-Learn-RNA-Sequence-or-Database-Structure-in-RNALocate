# Auditing Whether Models Learn RNA Sequence or Database Structure in RNALocate

## Abstract

RNA localization databases provide large training sets, but their coverage is uneven. Models may learn which genes and contexts are well documented rather than only sequence signals. We tested this possibility in RNALocate using 18,753 human genes and 16 cellular contexts, creating 300,048 gene-context pairs. We predicted whether each pair had a retained record, treating unrecorded pairs as unknown rather than biological negatives. Logistic-regression models used context, low-level RNA sequence features, and other-context support—the number of other contexts with a record for the same gene. Genes were separated by sequence-similarity components before training and testing. Context alone achieved held-out average precision (AP) of 0.8036, sequence alone 0.8234, and context plus support 0.9572. Adding sequence to context plus support increased AP by only 0.00087 (95% component-bootstrap confidence interval 0.00072–0.00105). In 200 negative controls, shuffling support among genes within each context and split reduced mean AP to 0.8108. The pattern appeared across the context panel and remained stable across three low-level sequence representations, although performance varied. In a separate annotation-recovery analysis, context was much stronger than 3-mer features overall, while sequence retained modest value within fixed contexts. Database structure can therefore create powerful predictive shortcuts. RNA localization benchmarks should report their observation process, split unit, and model inputs before interpreting high scores as transferable sequence biology.

**Keywords:** RNA localization; database coverage; shortcut learning; annotation bias; machine learning; RNALocate

## Introduction

RNA molecules are transported to particular parts of a cell, where their location can influence translation, stability, and regulation{1}.

Public resources such as RNALocate combine localization evidence from many studies, cell types, cell lines, RNA classes, and experimental methods{2}. This breadth is useful, but the database is not a single experiment in which every gene was measured in every context. A record can appear only after a gene and context are measured, reported, and curated. Consequently, the absence of a record may mean “not observed in this database” rather than “the RNA is biologically absent from this location.”

Computational models attempt to predict subcellular RNA localization from sequence. Approaches have included recurrent neural networks{3}, attention-based models{4}, transfer learning{5}, and pretrained RNA language models{6}. Their scores are often interpreted as evidence that sequence contains strong and transferable localization information.

A second concern is shortcut learning. A predictive model uses the easiest reliable pattern in its training data, even when that pattern is different from the scientific mechanism researchers intended to study{7}. Biological machine-learning studies can also obtain overly optimistic or misleading results when related sequences cross data partitions or when the available metadata reveal the target indirectly{8}. GraphPart{9} and DataSAIL{10} illustrate sequence-aware partitioning approaches developed to reduce these risks.

This study asks a narrow question: how much of RNALocate record availability can be predicted from database and context structure, and how much extra information is supplied by simple RNA sequence features? We compare every combination of three information sources: cellular context, other-context support for the same gene, and low-level sequence. We then test the central association using a shuffled-support negative control, context-specific analyses, calibration, and sequence-representation sensitivity. Finally, we examine whether the same structural signal affects recovery of recorded localization labels. The goal is not to build the most accurate RNA model. It is to determine what a high score means in this particular database benchmark.

## Methods

### Study design and claim boundary

We performed a computational audit of previously collected, publicly accessible RNALocate records. The primary task predicted record availability on a selected gene-context grid. It did not predict the true biological localization of an unrecorded RNA. All conclusions therefore concern structure within the selected database panel, not causal effects of sequence, context, publication, or curation.

### RNALocate panel and outcomes

RNALocate v3.0 archives were obtained from the official database{2}. We retained experimental *Homo sapiens* mRNA records for four exact locations: Nucleus, Cytoplasm, Endoplasmic reticulum, and Mitochondrion. Contexts with at least 1,000 unique observed genes were retained without using their four-location label distributions. The resulting panel contained 18,753 eligible genes and 16 contexts.

Crossing every eligible gene with every retained context produced 300,048 candidate pairs. A pair was marked observed if it contained at least one retained target-location record. There were 208,181 observed pairs and 91,867 unrecorded pairs. An unrecorded pair was not treated as an experimentally verified negative. A separate annotation-recovery task used only observed rows and predicted which of the four location labels had been recorded.

For each gene-context pair, other-context support was the number of the other 15 contexts containing a record for the same gene. The candidate context was excluded, so the feature did not directly read the outcome being predicted. Because every eligible gene had at least one retained record, support zero is a selection boundary: if a gene has no records in the other 15 contexts, its one required record must be in the candidate context.

### Sequence representation and partitioning

Eligible RefSeq coding accessions began with `NM_` or `XM_`. One representative transcript per gene was chosen by preferring `NM_`, then the longest sequence, with accession as a deterministic tie-break. The complete low-level sequence representation contained 70 features: 64 normalized 3-mer frequencies, five nucleotide or ambiguous-base fractions, and normalized log transcript length. This representative transcript is a reproducible gene-level proxy and may not be the isoform measured in the source experiment.

To limit train-test similarity, all sequences were compared once using MMseqs2 nucleotide search{11}. Edges required at least 80% nucleotide identity and 80% coverage in both directions. Connected sequence components, rather than individual rows, were assigned approximately 70% to training, 15% to validation, and 15% to testing. The primary graph contained 18,437 components. No retained threshold edge crossed partitions. The held-out test set contained 45,152 gene-context pairs from 2,764 components.

### Models and evaluation

The primary model family was logistic regression. We fit an intercept-only prevalence baseline and all seven non-empty combinations of context (C), other-context support (S), and sequence (Q): C, S, Q, C+S, C+Q, S+Q, and C+S+Q. Context was one-hot encoded and continuous features were standardized. The regularization parameter was chosen from 0.01, 0.1, 1, and 10 by validation AP. Each selected model was refit on training plus validation data and evaluated once on held-out components.

Average precision was the primary metric because the observed and unrecorded classes were imbalanced; AP should be interpreted relative to the positive prevalence{12}. The labels are also related to the positive-unlabeled setting, in which known positives and unlabelled examples do not form an ordinary positive-versus-negative dataset{13}. We also report area under the receiver-operating-characteristic curve (AUROC) and Brier score. Ninety-five percent percentile confidence intervals and paired AP differences used 1,000 bootstrap resamples of sequence components.

For the within-context annotation-recovery analysis, an eligible context-label cell was defined as a non-nucleus cell with at least 10 positive and at least 10 negative held-out rows. Five cells met this rule. The threshold was applied before averaging AP across cells.

### Negative control and robustness analyses

For the shuffled-support control, support values were randomly permuted among genes separately within every context and split. This preserved each context's support distribution, prevalence, and split membership while breaking the link between a gene and its actual cross-context coverage. Model selection and refitting were repeated for 200 fixed-seed permutations. The empirical one-sided probability was calculated as `(1 + number of shuffled AP values at least as large as the real AP) / 201`.

We calculated record fraction and Wilson 95% intervals for support values 0-15. For each context, we reported prevalence, AP, AP minus prevalence, normalized AP lift, AUROC, and Brier score. Calibration was summarized in ten equal-frequency test bins. Sequence sensitivity compared composition plus length (6 features), 3-mers only (64 features), and the complete representation (70 features). These are deliberately simple sequence views and do not test every possible sequence model.

Because support zero is a deterministic cohort-selection boundary, we specified a reviewer-requested post hoc sensitivity analysis before running it. Every pair with support zero was excluded from training, validation, and testing. Context, support, context-plus-support, and the full context-plus-support-plus-sequence model were then retuned and refit using the original split, model family, regularization grid, and 1,000 component-bootstrap replicates.

### Ethics and reproducibility

The study used previously published molecular database records and involved no direct interaction with human participants, identifiable clinical information, or new animal experiments. Analysis scripts, tests, the frozen analysis plan and sensitivity amendment, figure source tables, software versions, and file hashes are included in the project package. Redistribution of source-derived RNALocate records must follow the database terms.

## Results

### Database coverage varied by cellular context

The selected grid contained 300,048 pairs, of which 69.38% had at least one retained record. Coverage differed substantially among the 16 contexts (Figure 1A). In the test set, prevalence ranged from 0.3629 in IMR-90 to 0.9745 in HEK293T. Except for the structurally constrained support-zero group, record availability generally increased when the same gene had records in more other contexts (Figure 1B–C). This is a coverage co-occurrence pattern and not a biological dose response.

![Figure 1](figures/figure_1_database_coverage.png)

**Figure 1 | RNALocate database coverage is structured by context and same-gene cross-context support.** (A) Fraction of candidate gene pairs with at least one retained RNALocate record in each of 16 contexts. The dashed line is the unweighted mean across contexts. (B) Test-set record fraction by the number of records for the same gene in the other 15 contexts. Shading shows Wilson 95% confidence intervals. Support zero is marked as a cohort-selection boundary. (C) Context-specific record fractions across support values; gray cells are combinations absent from the test set. These panels describe database coverage co-occurrence, not biological localization or a causal dose response.

### Context and support explained most predictable record availability

The complete factorial comparison showed that context alone reached AP 0.8036 and the low-level sequence representation reached 0.8234 (Table 1; Figure 2A). Other-context support alone was stronger at 0.8996. Context plus sequence reached 0.8823, whereas context plus support reached 0.9572. The full model reached 0.9580.

**Table 1 | Held-out performance for the complete context-support-sequence model matrix.** Confidence intervals are 1,000 component-bootstrap intervals for AP. Lower Brier scores are better.

| Model | AP (95% CI) | AUROC | Brier score |
|---|---:|---:|---:|
| Prevalence baseline | 0.6876 (0.6754–0.6991) | 0.5000 | 0.2148 |
| Context | 0.8036 (0.7946–0.8127) | 0.6575 | 0.1980 |
| Other-context support | 0.8996 (0.8946–0.9040) | 0.8517 | 0.1326 |
| Low-level sequence | 0.8234 (0.8109–0.8344) | 0.7153 | 0.1854 |
| Context + sequence | 0.8823 (0.8736–0.8898) | 0.7839 | 0.1669 |
| Context + support | 0.9572 (0.9541–0.9597) | 0.9154 | 0.1061 |
| Support + sequence | 0.9071 (0.9015–0.9118) | 0.8541 | 0.1323 |
| Context + support + sequence | 0.9580 (0.9550–0.9606) | 0.9158 | 0.1058 |

Paired component bootstraps clarified the incremental information (Figure 2B). Adding sequence after context increased AP by 0.0787 (95% CI 0.0719-0.0863), whereas adding support after context increased AP by 0.1535 (0.1458-0.1616). Once context and support were included, adding sequence changed AP by only 0.00087 (0.00072-0.00105). Conversely, support still added 0.0757 (0.0695-0.0823) after context and sequence.

### Excluding the support-zero boundary did not change the conclusion

The support-zero exclusion removed 992 of 300,048 pairs across all splits, including 149 test pairs. The reduced test set contained 45,003 pairs with prevalence 0.6866. Context AP was 0.8014 (95% CI 0.7923-0.8115), support AP was 0.9006 (0.8959-0.9051), context-plus-support AP was 0.9573 (0.9544-0.9601), and full-model AP was 0.9582 (0.9552-0.9610). Context plus support remained 0.1559 AP above context (0.1475-0.1641), whereas sequence added 0.00086 after context plus support (0.00069-0.00104). Thus, the structural result did not depend on the deterministic support-zero pairs. This post hoc sensitivity addresses that boundary but not other forms of database selection.

![Figure 2](figures/figure_2_factorial_and_permutation.png)

**Figure 2 | Context and same-gene support explain most predictable record availability.** (A) Held-out AP for the prevalence baseline and all seven non-empty combinations of context, target-excluding other-context support, and low-level sequence features. Error bars are 95% intervals from 1,000 bootstrap resamples of sequence graph components. (B) Paired AP increments from the same bootstrap replicates. (C) AP from 200 controls in which support was shuffled among genes within every context and split, compared with the real context-plus-support model.

### Shuffled support removed the large gain

Across 200 shuffled-support controls, context plus shuffled support produced mean AP 0.8108; its 2.5th–97.5th percentile range was 0.8103–0.8113. The real context-plus-support model reached 0.9572, and none of the 200 shuffled results was as high, giving an empirical one-sided probability of 0.00498 (Figure 2C). Thus, the gain depended on linking each support value to the correct gene. The control does not reveal whether the pattern arose from measurement choices, publication, curation, biology, or a mixture of these processes.

### The pattern was broad but context-dependent

Context-plus-support AP exceeded prevalence in all 16 contexts (Figure 3A). Performance was especially weak in IMR-90, where AP was 0.4768 compared with prevalence 0.3629, but it was above 0.875 in every other context. Adding sequence to context plus support produced small, context-dependent changes (Figure 3B). Calibration plots showed that context-plus-support and full-model probabilities followed the observed fractions more closely than context alone (Figure 3C), consistent with their lower Brier scores.

The sequence-sensitivity analysis changed sequence-only AP from 0.7870 for composition plus length to 0.8139 for 3-mers and 0.8234 for the complete 70-feature view. In contrast, full-model AP ranged only from 0.9577 to 0.9580 (Figure 3D). The main conclusion was therefore not specific to one of these three low-level sequence encodings.

![Figure 3](figures/figure_3_robustness_and_calibration.png)

**Figure 3 | Coverage prediction is heterogeneous but robust to tested low-level sequence views.** (A) Prevalence and context-plus-support AP for each context. (B) Per-context AP difference after adding the complete low-level sequence vector to context plus support. (C) Equal-frequency reliability bins for context, context plus support, and the full model; vertical intervals are Wilson 95% confidence intervals. (D) Sensitivity to composition-plus-length, normalized 3-mer-only, and complete 70-feature sequence views.

### Database structure also affected annotation recovery

The separate annotation-recovery task predicted four recorded localization labels among already observed gene-context rows. Label-frequency AP was 0.2930, 3-mer sequence AP was 0.3054, context AP was 0.7198, and context plus 3-mers reached 0.7567 (Figure 4A). Results varied sharply by localization label (Figure 4B), so the overall average should not be treated as uniform performance.

Within a fixed context, a context-only model cannot rank genes. Across the five non-nucleus context-label cells with at least 10 positive and 10 negative held-out rows, mean within-context AP was 0.4523 for the context baseline, 0.5281 for 3-mers, and 0.5352 for context plus 3-mers (Figure 4C). This shows that the tested sequence features retained modest ranking information even though context dominated the global average. Context dominance was stable across the primary sequence graph and three MMseqs2 clustering stress tests (Figure 4D).

![Figure 4](figures/figure_4_annotation_recovery_consequence.png)

**Figure 4 | Database structure has a clear consequence for localization annotation recovery.** (A) Macro AP for recovering four recorded localization labels among already observed gene-context rows. (B) Per-label AP. (C) Mean AP across the five non-nucleus context-label cells with at least 10 positive and 10 negative held-out rows. (D) Context and low-level-sequence macro AP under the primary sequence graph and three clustering stress tests. These analyses concern database annotations and do not validate localization for unrecorded pairs.

## Discussion

### Principal finding

The main result is simple: RNALocate record availability is highly structured. Context and target-excluding same-gene support together reached AP 0.9572, and a low-level sequence vector added less than 0.001 AP after those variables were present. Shuffling support while preserving its context-specific distribution removed most of the gain. A model can therefore obtain a very high score by learning regularities in where database records exist, even without solving the intended biological localization problem.

### What the result does and does not mean

The study does not show that sequence is irrelevant. Sequence alone performed above prevalence, adding sequence after context improved AP, and sequence ranked genes within fixed contexts in the annotation task. The narrower conclusion is that database structure explained most of the predictable record-availability signal in this panel before sequence was considered.

The study also does not prove “curator bias.” Other-context support could reflect many linked processes: some genes are studied more often; some contexts receive deeper experiments; abundant or technically accessible transcripts may be easier to measure; biologically broad localization may generate more records; and database selection rules can preserve those differences. The present observational analysis cannot separate these explanations.

### Why this matters for benchmark design

High AP is meaningful only when the outcome matches the scientific question. In the availability task, AP measures how accurately a model ranks observed versus unrecorded database pairs. It does not measure whether an RNA truly localizes to a compartment. In the annotation task, performance is conditional on a row already having at least one retained record. These distinctions should be stated whenever a database benchmark is used.

Researchers can reduce misleading interpretations by reporting four items. First, define what an absent label means. Second, split by biological or sequence-related units rather than individual rows when related examples could cross partitions. Third, compare sequence models with simple metadata and database-structure baselines. Fourth, measure the incremental value of sequence after those baselines, not only the score of a combined model. A shuffled-feature control can further test whether a gain depends on the intended entity-level link.

### Limitations

The panel was selected after retaining genes and contexts with sufficient records, so its prevalence is not the prevalence of all human mRNA localization. Support zero is structurally constrained by gene eligibility; excluding it in a post hoc sensitivity analysis did not materially change the results. The selected representative transcript may not match the localized isoform. The primary split controls detected sequence links at one identity and coverage threshold but cannot rule out every distant relationship. Logistic regression and three low-level sequence views do not represent all possible RNA models; a richer representation could discover additional biological signal. Context-specific estimates came from one selected database panel, and the poorest context showed that the shortcut is not equally strong everywhere. Bootstrap intervals condition on the fixed dataset, preprocessing, model family, and split. Finally, no experiment here identifies the causal origin of the observed coverage structure.

### Future work and conclusion

A stronger biological follow-up would use a prospectively designed experiment in which the same genes are measured across multiple cellular contexts with one consistent protocol and verified positive and negative outcomes. Sequence models could then be evaluated on entirely new sequence groups and new contexts. Until such data are available, database benchmarks should separate record prediction from biological localization prediction.

In conclusion, context and same-gene database coverage can act as powerful shortcuts in RNALocate. The complete factorial comparison, shuffled-support control, context analysis, calibration, and sequence sensitivity all support the same interpretation: most of the tested record-availability performance came from database structure, while low-level sequence supplied a smaller residual contribution. Clear outcome definitions and baseline audits are essential before high localization scores are treated as evidence of transferable RNA biology.

### Data and code availability

Code, analysis scripts, aggregate results, paper-facing source tables, and figures are available at https://github.com/billylu24/Auditing-Whether-Models-Learn-RNA-Sequence-or-Database-Structure-in-RNALocate. RNALocate v3.0 is available from its official database{2}.

## References

1. AR Buxbaum, G Haimovich, RH Singer. In the right place at the right time: visualizing and understanding mRNA localization. Nature Reviews Molecular Cell Biology. Vol. 16, pg. 95–109, 2015, DOI: 10.1038/nrm3918.
2. L Wu, L Wang, S Hu, G Tang, J Chen, Y Yi, H Xie, J Lin, M Wang, D Wang, B Yang, Y Huang. RNALocate v3.0: advancing the repository of RNA subcellular localization with dynamic analysis and prediction. Nucleic Acids Research. Vol. 53, pg. D284–D292, 2025, DOI: 10.1093/nar/gkae872.
3. Z Yan, E Lécuyer, M Blanchette. Prediction of mRNA subcellular localization using deep recurrent neural networks. Bioinformatics. Vol. 35, pg. i333–i342, 2019, DOI: 10.1093/bioinformatics/btz337.
4. D Wang, Z Zhang, Y Jiang, Z Mao, D Wang, H Lin, D Xu. DM3Loc: multi-label mRNA subcellular localization prediction and analysis based on multi-head self-attention mechanism. Nucleic Acids Research. Vol. 49, pg. e46, 2021, DOI: 10.1093/nar/gkab016.
5. J Wang, M Horlacher, L Cheng, O Winther. DeepLocRNA: an interpretable deep learning model for predicting RNA subcellular localisation with domain-specific transfer-learning. Bioinformatics. Vol. 40, pg. btae065, 2024, DOI: 10.1093/bioinformatics/btae065.
6. M Zeng, X Zhang, Y Li, C Lu, R Yin, F Guo, M Li. RNALoc-LM: RNA subcellular localization prediction using pre-trained RNA language model. Bioinformatics. Vol. 41, pg. btaf127, 2025, DOI: 10.1093/bioinformatics/btaf127.
7. R Geirhos, JH Jacobsen, C Michaelis, R Zemel, W Brendel, M Bethge, FA Wichmann. Shortcut learning in deep neural networks. Nature Machine Intelligence. Vol. 2, pg. 665–673, 2020, DOI: 10.1038/s42256-020-00257-z.
8. J Bernett, DB Blumenthal, DG Grimm, F Haselbeck, R Joeres, OV Kalinina, M List. Guiding questions to avoid data leakage in biological machine learning applications. Nature Methods. Vol. 21, pg. 1444–1453, 2024, DOI: 10.1038/s41592-024-02362-y.
9. F Teufel, MH Gíslason, JJA Almagro Armenteros, AR Johansen, O Winther, H Nielsen. GraphPart: homology partitioning for biological sequence analysis. NAR Genomics and Bioinformatics. Vol. 5, pg. lqad088, 2023, DOI: 10.1093/nargab/lqad088.
10. R Joeres, DB Blumenthal, OV Kalinina. Data splitting to avoid information leakage with DataSAIL. Nature Communications. Vol. 16, pg. 3337, 2025, DOI: 10.1038/s41467-025-58606-8.
11. M Steinegger, J Söding. MMseqs2 enables sensitive protein sequence searching for the analysis of massive data sets. Nature Biotechnology. Vol. 35, pg. 1026–1028, 2017, DOI: 10.1038/nbt.3988.
12. T Saito, M Rehmsmeier. The precision-recall plot is more informative than the ROC plot when evaluating binary classifiers on imbalanced datasets. PLOS ONE. Vol. 10, pg. e0118432, 2015, DOI: 10.1371/journal.pone.0118432.
13. F Li, S Dong, A Leier, M Han, X Guo, J Xu, X Wang, S Pan, C Jia, Y Zhang, GI Webb, LJM Coin, C Li, J Song. Positive-unlabeled learning in bioinformatics and computational biology: a brief review. Briefings in Bioinformatics. Vol. 23, pg. bbab461, 2022, DOI: 10.1093/bib/bbab461.
