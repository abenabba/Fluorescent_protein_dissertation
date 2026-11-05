# Fluorescent Protein Dissertation Analysis

This repository contains scripts and intermediate analysis files used for the computational analyses described in Chapter 4 of my doctoral dissertation:

**Engineering Brighter Fluorescent Proteins: High-Diversity Libraries and Machine-Learning-Guided Design**

The Chapter 4 workflow combines multibin fluorescence-activated cell sorting, barcode sequencing, PacBio barcode-to-protein mappings, FlowSeq-based fluorescence inference, protein language model representations, AlphaFold2-derived structural representations, chromophore-local structural analysis, and machine-learning-guided candidate prioritization.

## Repository status

This repository is intended primarily as an archival and reproducibility resource for the analyses described in the dissertation.

The original analysis was developed iteratively across local workstations and the University of Oregon high-performance computing environment. Some scripts used during intermediate stages of the analysis have not yet been recovered or added to this repository. Where the original script is not currently available, the corresponding analysis step and its inputs/outputs are documented below.

Some scripts also contain file paths corresponding to the original analysis directories. These paths will need to be changed before the scripts are run on another system.

Large raw sequencing files, AlphaFold2 structure files, and some intermediate arrays are not included in the repository.

---

# Analysis workflow

The overall analysis proceeded approximately as follows:

```text
Multibin FACS
    ↓
Illumina barcode sequencing
    ↓
PacBio barcode-to-sequence mapping
    ↓
Protein × fluorescence-bin count matrix
    ↓
FlowSeq fractional-median inference
    ↓
MEFL fluorescence calibration
    ↓
Quality filtering / barcode-support analysis
    ↓
ESM2 sequence representations
    ↓
Whole-sequence brightness modeling
    ↓
AlphaFold2 structural representations
    ↓
Chromophore-local feature extraction
    ↓
Structural clustering and within-cluster analysis
    ↓
Residual-brightness selection
    ↓
Two-stage ProtGPT2 fine-tuning
    ↓
Sequence and structure-based candidate prioritization
```

---

# 1. Barcode processing and barcode-to-sequence mapping

### `true_bc_make.mak`

Makefile used during barcode-processing and mapping steps.

This file contains commands used to coordinate portions of the barcode-processing workflow. It should be run or inspected before downstream R analysis when reconstructing the workflow from raw barcode data.

The upstream sequencing workflow included:

1. paired-end Illumina read merging,
2. extraction of the 20-bp variant barcode using conserved flanking sequences,
3. barcode error correction/clustering,
4. comparison of barcode clusters against PacBio-derived barcode-to-protein mappings.

The exact raw-data preprocessing scripts used for all of these steps have not yet been completely recovered.

### `pacbiomapping_totalcount.R`

Used to summarize PacBio barcode-to-sequence mappings and associated barcode observations at the protein-sequence level.

This step links experimentally observed barcodes to the corresponding full-length amino-acid sequence and generates sequence-level mapping/count information used in downstream analyses.

**Use before:** barcode-support calculation and FlowSeq dataset construction.

---

# 2. Barcode-cluster support

### `barcode_cluster_support_bright_bc3_bc4.R`

Calculates barcode-cluster support for protein variants in the brighter-library sequencing data.

Barcode-cluster support was used as a measure of confidence in variant identity. Multiple distinct barcode clusters associated with the same protein sequence provide stronger evidence for that sequence than repeated sequencing of a single barcode.

The resulting sequence-level support measurement was later merged with the inferred-brightness dataset.

### Output

`bright_bc3_bc4_barcode_cluster_support.csv`

Contains protein sequences and their associated barcode-cluster support values, including the variable used in downstream analyses as:

```text
n_barcode_clusters_total
```

Support thresholds used in Chapter 4 included:

- no support threshold beyond sequence filtering,
- barcode support ≥ 3,
- barcode support ≥ 5.

---

# 3. FlowSeq brightness inference

### `inferred_brightness_calculation_two_libs.R`

Primary R workflow used to infer fluorescence distributions from multibin sequencing data.

This script performs the core FlowSeq calculations, including:

1. normalization of the medium- and brighter-blue presort populations by sequencing depth,
2. construction of the 50:50 medium/brighter presort baseline,
3. normalization of sequencing counts within each sorted fluorescence bin,
4. correction for the number of FACS events collected in each bin,
5. calculation of enrichment relative to the presort population,
6. conversion of enrichment values into a normalized per-variant distribution across fluorescence bins,
7. calculation of the fractional median fluorescence-bin position for each variant.

The number of cells collected in each fluorescence population is reported in Appendix Table C1 of the dissertation.

### Intermediate/source data

`Samples_Brighter_plus_Medium_combined_*`

Combined medium- and brighter-library count information used during construction of the FlowSeq dataset.

The exact filename currently present in the repository should be used in place of the wildcard above.

---

# 4. Conversion from inferred bin position to calibrated fluorescence

### `mefl_convert_excludezero.py`

Used after fractional-median FlowSeq inference to convert inferred fluorescence-bin positions into calibrated fluorescence measurements.

The analysis uses fluorescence calibration values derived from rainbow calibration beads and interpolates between neighboring bin medians to assign each variant a quantitative fluorescence estimate in molecules of equivalent fluorescein (MEFL).

Variants associated with the fluorescence floor or otherwise unresolved zero-valued measurements were handled separately from resolved non-floor measurements.

### Main output

`brighter_flowseq_inferred_brightness_updated.csv`

This is one of the primary downstream Chapter 4 datasets.

It contains the inferred brightness measurements together with sequence and FlowSeq-related metadata and was subsequently used for sequence filtering, representation generation, and machine-learning analyses.

---

# 5. ESM2 sequence representations

### `embed_full_csv_esm2.py`

Generates protein-language-model representations using:

```text
ESM2_t30_150M_UR50D
```

Protein sequences shorter than 150 amino acids were removed before embedding.

Final-layer residue representations were mean-pooled across sequence positions to generate one fixed-length **640-dimensional ESM2 representation per protein**.

The length-filtered dataset contained **15,228 protein sequences**.

The embeddings were saved as a numerical matrix with corresponding sequence metadata and were subsequently exported into CSV-compatible form for analyses performed in R.

### `pca_umap.R`

Performs exploratory dimensionality-reduction analyses of the ESM2 representations.

Analyses include:

- principal-component analysis (PCA),
- Uniform Manifold Approximation and Projection (UMAP),
- visualization by inferred brightness,
- barcode support,
- protein length,
- and other sequence/phenotype metadata.

The UMAP analysis used 30 nearest neighbors, minimum distance = 0.2, and Euclidean distance.

This analysis was used primarily to determine whether fluorescence phenotypes occupied visibly distinct regions of ESM2 representation space.

### `ridge_bright_ge3.py`

Performs ridge-regression benchmarking using protein representations and inferred fluorescence.

This script was used for analysis of the higher-confidence subset containing variants supported by at least three barcode clusters.

The target variable was log10-transformed inferred MEFL brightness.

Model performance was evaluated on a held-out test set using metrics including:

- R²,
- RMSE,
- Pearson correlation.

The ESM2 analyses showed little held-out predictive relationship between global mean-pooled sequence representations and fluorescence brightness.

---

# 6. Barcode-support and sequence-feature model benchmarks

The Chapter 4 benchmark analysis compared several sources of predictive information:

- barcode-cluster support,
- simple sequence-composition features,
- ESM2 representations,
- AlphaFold2 representations,
- combinations of support and sequence-derived variables.

Simple sequence features included:

- protein length,
- aromatic residue fraction,
- charged residue fraction,
- polar residue fraction,
- Gly/Pro fraction,
- hydrophobic residue fraction,
- candidate chromophore-motif indicators.

Support-only, sequence-based, residual, ridge-regression, and XGBoost models were compared.

### `all_model_performance_summaries_original*`

Summary table containing performance metrics from the model-comparison analyses.

This is an output/intermediate analysis file rather than a primary processing script.

---

# 7. AlphaFold2-derived whole-protein representations

Protein structures were predicted using ColabFold/AlphaFold2 with the `alphafold2_ptm` model in single-sequence mode.

Internal AlphaFold2 `single_repr` representations were retained. Each residue has a 256-dimensional representation, and these residue-level values were mean-pooled across the sequence to generate one 256-dimensional vector per protein.

The final pooled dataset contained **4,954 protein sequences**.

### `bright_partial_single_repr_pooled_4954emb.csv`

Pooled AlphaFold2 `single_repr` representation dataset containing the 4,954 variants used in the structure-representation analysis.

### `af2_brighter_embedding_master_analysis.R`

Primary downstream analysis of the AlphaFold2-derived whole-protein representations.

This analysis includes dimensionality reduction and brightness-model benchmarking using the pooled `single_repr` features.

For the regression analysis, the mean-pooled representations were reduced to 20 principal components before model fitting.

The scripts originally used to run the complete ColabFold/AlphaFold2 prediction workflow and extract/pool all `single_repr` arrays are **not currently included in this repository**.

The procedure used was:

1. select variants meeting the sequence-length and barcode-support criteria,
2. generate AlphaFold2 predictions using ColabFold,
3. retain internal `single_repr` features,
4. mean-pool residue-level representations,
5. combine representations with protein-sequence and fluorescence metadata.

---

# 8. Chromophore-motif identification

### `chromophore_detect_small.R`

Used to identify candidate chromophore-forming sequence motifs within fluorescent protein variants.

This step provides the sequence anchor used for subsequent structure-localized analysis.

Candidate motifs were evaluated within the protein sequence and used to identify the approximate location of the putative chromophore in the corresponding AlphaFold2-predicted structure.

Chromophore identification is used **before local structural feature extraction**.

---

# 9. Chromophore-local structural feature extraction

### `bright_af2_local_chromophore_features_par*`

Parallelized structural-feature extraction workflow for the AlphaFold2-predicted proteins.

Replace the wildcard with the exact filename currently present in the repository.

The analysis characterizes the predicted structural environment surrounding the putative chromophore rather than representing the entire protein.

Features were calculated from a local neighborhood around the chromophore motif, including measurements related to:

- local aromatic-residue content,
- polar-residue content,
- charged-residue content,
- local residue/packing density,
- distances between the chromophore and neighboring residues,
- local structural span,
- chromophore-motif span,
- local AlphaFold2 confidence,
- chromophore-motif pLDDT,
- other local geometric and confidence measurements.

For several structural measurements, neighboring residues were defined using a **10 Å chromophore-centered neighborhood**.

The expanded final chromophore-local analysis contained **1,132 non-floor variants with complete structural measurements**.

---

# 10. Chromophore-local clustering

The 14 chromophore-local structural features were standardized and analyzed using principal-component analysis followed by k-means clustering.

Four chromophore-environment clusters were obtained:

```text
Cluster 1: 395 variants
Cluster 2: 126 variants
Cluster 3: 254 variants
Cluster 4: 357 variants
```

Brightness was not used to define the clusters.

The original standalone script used for the final PCA/k-means clustering step has not yet been definitively identified in the recovered repository files. The resulting cluster assignments are contained in the downstream chromophore-analysis datasets.

---

# 11. Residual-brightness and within-cluster feature analysis

Residual brightness was calculated to reduce the influence of technical/support-associated variation.

A linear model was fit using:

```text
log10(MEFL) ~ barcode-cluster support + protein length
```

Residual brightness was defined as:

```text
observed log10(MEFL) - model-predicted log10(MEFL)
```

Positive residuals therefore represent proteins that were brighter than expected from barcode support and protein length alone.

### `within_cluster_feature_effects_residual_all14*`

Used to compare chromophore-local structural features between higher- and lower-residual-brightness variants within each structural cluster.

The analysis:

1. ranks variants within each structural cluster by residual brightness,
2. defines upper and lower residual-brightness quartiles,
3. compares all 14 chromophore-local features,
4. calculates standardized effect sizes,
5. performs statistical comparisons,
6. applies Benjamini-Hochberg multiple-testing correction.

This analysis identified a relatively strong association between local aromatic-residue number and residual brightness in one structural context, although the association was not significant after correction across all 56 cluster-feature comparisons.

### `figure_4_within_cluster_effects.R`

Generates the figure summarizing within-cluster chromophore-feature effects.

This script is used **after** the residual-brightness and within-cluster feature comparisons have been calculated.

---

# 12. Additional source and intermediate data files

### `ML_FACS_allstats_good3_v2_noCol.xlsx`

Legacy/source analysis spreadsheet retained for provenance from the ML/FACS analysis workflow.

This file contains intermediate statistics used during development and curation of the fluorescent-protein analysis.

It is not a standalone executable analysis.

### `bright_bc3_bc4_barcode_cluster_support.csv`

Output of the barcode-cluster support analysis.

### `bright_partial_single_repr_pooled_4954emb.csv`

Pooled AlphaFold2 whole-protein `single_repr` representations.

### `brighter_flowseq_inferred_brightness_updated.csv`

Primary inferred-brightness dataset used in downstream Chapter 4 analyses.

### `all_model_performance_summaries_original*`

Model-benchmark performance summary.

### `ml blue candidates (1).csv`

Candidate-tracking table used during prioritization of blue fluorescent protein sequences generated after two-stage ProtGPT2 fine-tuning.

---

# 13. Two-stage ProtGPT2 fine-tuning

The original training and generation scripts used for the final two-stage ProtGPT2 experiment have **not yet been added to this repository**.

The analysis proceeded as follows.

## Stage 1

ProtGPT2 was fine-tuned on:

```text
364 unique protein sequences
```

These variants were:

- non-floor variants,
- length filtered,
- drawn from the structurally characterized fluorescent protein dataset.

The goal of the first stage was to expose the model to sequence patterns represented within the experimentally characterized shuffled fluorescent protein library.

## Stage 2

The Stage-1 model was further fine-tuned using:

```text
95 variants
```

These sequences were selected from the **highest quartile of residual brightness across the chromophore-environment clusters**.

Residual brightness represented fluorescence not explained by barcode-cluster support and protein length.

The goal was therefore to bias the generative model toward proteins that were brighter than expected from these covariates rather than directly train a regression model to predict fluorescence.

## Sequence generation and filtering

The two-stage model generated:

```text
2,000 sequences
```

Candidate filtering then proceeded approximately as follows:

```text
2,000 generated sequences
        ↓
304 unique sequences between 170 and 306 aa
        ↓
150 sequences containing chromophore-like motifs
        ↓
additional filtering by motif plausibility,
sequence complexity, and sequence similarity
        ↓
9 candidates advanced to AlphaFold2 prediction
        ↓
3 candidates selected for synthesis
```

The three final candidates were:

| Candidate | Generation ID | Length |
|---|---|---:|
| Candidate 1 | `gen_7` | 260 aa |
| Candidate 2 | `gen_1983` | 234 aa |
| Candidate 3 | `gen_1760` | 218 aa |

The final candidates were selected using both sequence-level criteria and AlphaFold2-predicted structural plausibility.

### Candidate output

`ml blue candidates (1).csv`

Contains candidate-level information retained during the final prioritization workflow.

---

# Scripts not yet recovered

The following analysis steps were performed for Chapter 4 but the complete original scripts are not currently included in this repository.

## Raw Illumina barcode preprocessing

The original workflow:

1. merged paired-end reads using BBMerge,
2. located the expected barcode-flanking motifs,
3. extracted the intervening 20-bp barcode,
4. clustered near-identical barcodes using Starcode with Levenshtein distance 1,
5. compared barcode clusters with PacBio-derived barcode-to-sequence reference mappings.

Some of these operations may also be represented in `true_bc_make.mak`, but the complete original preprocessing workflow has not yet been reconstructed as a single standalone script.

## Flow-cytometry bead calibration

Rainbow calibration beads were used to convert fluorescence detector measurements into MEFL units.

The calibration workflow generated calibrated fluorescence values for each sorting bin, which were subsequently used by the FlowSeq/MEFL interpolation workflow.

The original standalone bead-calibration script is not currently included.

## AlphaFold2 prediction and `single_repr` extraction

The original workflow generated ColabFold/AlphaFold2 predictions and retained internal `single_repr` features before mean pooling.

The final pooled representation matrix is retained as:

```text
bright_partial_single_repr_pooled_4954emb.csv
```

## Final chromophore-feature PCA/k-means clustering

Four clusters were generated from 14 standardized chromophore-local features.

The final cluster assignments are retained in downstream data, but the standalone script responsible for the final clustering run has not yet been definitively identified.

## ProtGPT2 training and candidate-generation scripts

The original scripts used for:

- Stage-1 ProtGPT2 fine-tuning,
- Stage-2 residual-brightness fine-tuning,
- sequence generation,
- generated-sequence filtering,
- candidate ranking,

have not yet been added.

The exact selection logic and dataset sizes are documented above so that the analysis can be reconstructed even before the original scripts are recovered.

---

# Suggested execution order

For reconstruction of the Chapter 4 analysis using the currently recovered files:

```text
1. true_bc_make.mak
       ↓
2. pacbiomapping_totalcount.R
       ↓
3. barcode_cluster_support_bright_bc3_bc4.R
       ↓
4. inferred_brightness_calculation_two_libs.R
       ↓
5. mefl_convert_excludezero.py
       ↓
6. embed_full_csv_esm2.py
       ↓
7. pca_umap.R
       ↓
8. ridge_bright_ge3.py
       ↓
9. AlphaFold2/ColabFold single_repr generation
   [original generation scripts not yet recovered]
       ↓
10. af2_brighter_embedding_master_analysis.R
       ↓
11. chromophore_detect_small.R
       ↓
12. bright_af2_local_chromophore_features_par*
       ↓
13. chromophore-feature PCA / k-means clustering
    [standalone final script not yet identified]
       ↓
14. within_cluster_feature_effects_residual_all14*
       ↓
15. figure_4_within_cluster_effects.R
       ↓
16. two-stage ProtGPT2 training and generation
    [training scripts not yet added]
       ↓
17. candidate prioritization / ml blue candidates (1).csv
```

Not every exploratory analysis was used in the final dissertation, and some intermediate scripts were superseded by later versions. The dissertation Methods describe the final analyses used for interpretation.

---

# Dissertation correspondence

The repository primarily supports the following Chapter 4 sections:

| Dissertation section | Analysis |
|---|---|
| 4.3–4.5 | Multibin sorting, barcode sequencing, barcode-to-sequence mapping |
| 4.6 | FlowSeq brightness inference and MEFL calibration |
| 4.7 | Length filtering, barcode support, cross-library confidence |
| 4.8.1 | ESM2 sequence representations |
| 4.8.2 | AlphaFold2-derived protein representations |
| 4.8.3 | PCA and UMAP |
| 4.8.4–4.9 | Supervised brightness prediction and benchmark models |
| 4.10.2–4.10.4 | ESM2 and barcode-support analyses |
| 4.10.5–4.10.7 | Chromophore-local structural analyses |
| 4.10.8 | Two-stage ProtGPT2 fine-tuning and candidate prioritization |

---

# Reproducibility notes

Several scripts contain absolute paths corresponding to the original analysis environment. Before rerunning the analyses:

1. change the input/output paths to the local repository location;
2. verify that required intermediate CSV/NumPy/structure files are available;
3. inspect package imports at the beginning of each script;
4. preserve the filtering thresholds described in the dissertation when attempting to reproduce the final reported analyses.

The repository reflects an exploratory dissertation research workflow rather than a fully packaged software pipeline. Intermediate files have therefore been retained where they provide useful provenance for the final analyses.

---

# Code availability

This repository accompanies Chapter 4 of:

**Anissa Benabbas. _Engineering Brighter Fluorescent Proteins: High-Diversity Libraries and Machine-Learning-Guided Design_. University of Oregon, 2026.**

Repository:

`https://github.com/abenabba/Fluorescent_protein_dissertation`

Raw sequencing data and other large experimental datasets are described separately in the dissertation's data-availability sections.
