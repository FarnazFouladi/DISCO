
# DISCO

<!-- README.md is generated from README.Rmd. Please edit README.Rmd. -->

DISCO (**DIS**ruption in **CO**rrelations) identifies features whose
correlation patterns differ between two groups, such as healthy controls
and people with a disease. It combines evidence from pairwise
differential-correlation tests into a score, p-value, and adjusted
p-value (q-value) for each feature. This helps identify features whose
relationships with other features change between groups, regardless of
changes in average abundances of features between the two groups.

DISCO supports one modality (for example, correlations among microbial
taxa) or two modalities (for example, correlations between microbes and
metabolites). It can be applied to any high dimensional data such as
microbiome, transcriptome, and metabolome. Data should be normalized
appropriate to the data type before applying DISCO.

DISCO uses Spearman correlations with Fisher transformations and
[Bonett–Wright standard errors
(2000)](https://doi.org/10.1007/BF02294183), then combines
differential-correlation p-values using the [Cauchy combination test
(Liu and Xie, 2020)](https://doi.org/10.1080/01621459.2018.1554485).

## Installation

``` r
install.packages("remotes") # If needed
remotes::install_github("FarnazFouladi/DISCO")
library(DISCO)
```

## Inputs and arguments

The main function, `disco()`, compares `data1` (group 1, such as
controls) with `data2` (group 2, such as cases). Each input is a list
containing **one or two numeric matrices or data frames**, with
**samples in rows and features in columns**.

- **One modality:** DISCO compares correlations among features and
  assigns a score to each feature based on its differential correlations
  with the others.
- **Two modalities:** DISCO compares cross-modality correlations and
  scores features in the **first modality**. For example, list microbes
  first and metabolites second to score microbes by changes in their
  metabolite associations. Use the same modality order in both groups.

Features are matched by column name separately within each modality.
DISCO keeps features shared between groups and reorders group 2 to match
group 1; you do not need to align columns manually.

For two modalities, row names must identify samples consistently. Within
each group, the modalities must share samples, but their sample sets and
row order may differ. DISCO aligns rows by sample name and represents
measurements from an unavailable modality as missing. Each feature pair
must have at least `cooccur` paired, nonmissing observations **in each
group** to be tested.

**Normalize inputs appropriately before running DISCO.** The function
does not normalize data or adjust for covariates. Supply normalized and,
when needed, covariate-adjusted measurements. In the examples below we
use a robust centered log-ratio (CLR) transformation for microbiome
data: each sample is centred over its observed positive features, and
abundance zeros are represented as `NA`. DISCO then uses complete
observations separately for each feature pair; it does not require a
sample to be observed for every feature. Other data types may require
different preprocessing. DISCO itself treats numeric zeros as observed
values; missing values are excluded separately for each pair.

| Argument | Default | Meaning |
|----|---:|----|
| `cooccur` | `20` | Minimum pairwise non-missing observations in **each** group; integer ≥10. |
| `empirically_derive_threshold` | `TRUE` | Select the correlation filtering threshold using label permutations. |
| `filter_threshold` | `1` | Fixed correlation p-value threshold when empirical selection is disabled. |
| `q_value` | `0.2` | Feature-level q-value cutoff used to count discoveries in empirical mode. |
| `global_p_cutoff` | `0.05` | Select a threshold by the real-to-null ratio only when the global permutation p-value is below this cutoff; otherwise use `1`. |
| `n_resamples` | `100` | Number of size-matched label permutations. |
| `random_seed` | `1L` | Reproducible permutation seed; the caller’s RNG state is restored. |
| `adj_method` | `"BY"` | Multiple-testing adjustment passed to `p.adjust()`. |
| `workers` | `1L` | Number of parallel workers for permutations; capped at `n_resamples`. |
| `output_dir`, `prefix` | Required | Results directory and filename prefix. |

### Correlation filtering and threshold selection

DISCO estimates a Spearman correlation for each eligible pair in each
group, tests whether the correlations differ, and combines the raw
differential-test p-values separately for each focal feature. In
high-dimensional data, screening out weak correlations can reduce the
number of uninformative pairs contributing to a score and increase the
power, although this could be at the expense of increasing false
positives.

`filter_threshold` screens pairs using their **single-group correlation
p-values**, not their differential-correlation p-values. A pair is
retained if its correlation p-value is below the threshold in **at least
one group**. For example, `filter_threshold = 0.05` excludes pairs whose
correlation p-values are ≥0.05 in both groups. A threshold of `1` retain
all correlations.

Since the choice of filter_threshold could be arbitrary, DISCO can
derive a data-driven threshold empirically. By default, DISCO
(`empirically_derive_threshold = TRUE`) compares the observed data with
size-matched random label permutations at different filtering thresholds
`c(0.001, 0.005, 0.01, 0.05, 0.1, 0.5, 1)`. At each threshold, it
calculates the percentage of tested features with DISCO `q < q_value` in
the observed data and in each permutation. Then, a global permutation
test compares the maximum observed discovery percentage across
thresholds with the corresponding maximum from each valid permutation:

```text
global p = (1 + number of null maxima >= real maximum) / (1 + valid permutations)
```

If `global p < global_p_cutoff`, DISCO selects the threshold that
maximizes `max(real percentage, 1) / max(mean null percentage, 1)`. The
floors are in percentage points and prevent ratios from being dominated
by near-zero values. Ties select the largest threshold. Otherwise, DISCO
uses threshold `1`.

For two modalities, permutations preserve each modality’s group sizes
and sample overlap. A shared sample receives the same shuffled group
label in both modalities.

To use a fixed threshold, set `empirically_derive_threshold = FALSE` and
supply `filter_threshold` (default `1`). This skips permutations and
diagnostic plots. For example, with already preprocessed matrices
`controls` and `cases`:

``` r
result <- disco(
  data1 = list(controls),
  data2 = list(cases),
  empirically_derive_threshold = FALSE,
  filter_threshold = 1,
  output_dir = "disco_results",
  prefix = "controls_vs_cases"
)
head(result$disco)
```

## Returned results

| Component | Contents |
|----|----|
| `cor_diff` | Feature pairs, correlations `r1`/`r2`, counts `n1`/`n2`, and differential-correlation `p`/`q`, sorted by increasing `q`. |
| `disco` | Feature-level scores (`t`), p-values, q-values, and comparison counts (`num_corr`), sorted by increasing `q`. |
| `selected_threshold` | Selected empirical threshold or supplied fixed threshold. |
| `global_test` | Global p-value, cutoff, signal decision, permutation counts, and selected threshold. |
| `threshold_summary` | Real/null discovery percentages, null SD, floored percentages, and ratios across thresholds. |
| `discovery_plot`, `ratio_plot` | Real versus null discovery curves and their ratio. |
| `feature_null_summary`, `feature_null_plot` | Per-feature null discovery frequencies versus real p-values, across thresholds. |
| `output_files` | Paths to saved CSV tables and PNG figures. |

In `feature_null_summary`, a column such as `q_less_than_0.2` indicates
whether the real-data feature q-value is below `0.2`. Its name follows
the supplied `q_value` by the user. Null discovery percentages use all
requested permutations as the denominator, counting absent or untested
features as non-discoveries.

`feature_null_plot` shows **−log10(real-data p-value)** on the x-axis
and the percentage of permutations in which the feature has
`q < q_value` on the y-axis. Orange points identify real-data
discoveries. Low null discovery frequencies provide a useful diagnostic,
but are not formal feature-level FDR estimates.

In fixed mode, `cor_diff`, `disco`, `selected_threshold`, and
`output_files` are returned; the diagnostic components are `NULL`. Pair
and feature tables are saved as `<prefix>_cor_diff.csv` and
`<prefix>_disco.csv`.

## Example: One modality — taxon–taxon correlations in type 2 diabetes

This example uses **QinJ_2012**, comparing type 2 diabetes (**T2D**)
with strictly healthy controls in `curatedMetagenomicData`. Both the
one-modality and two-modality examples use this same comparison; studies
are not pooled.

The workflow follows these steps:

1.  Select one independent sample per participant, excluding conflicting
    disease labels and samples missing **gender, BMI, or age category**.
    Retain the earliest covariate-complete sample, breaking
    collection-day ties by sample ID.
2.  Retain features present in **≥20% of each group**, testing
    prevalence separately in controls and T2D.
3.  Apply **robust CLR** within each sample, centring log positive
    abundances over the retained observed features. Abundance zeros and
    existing missing values remain `NA`; no zero replacement or matrix
    completion is performed.
4.  Fit **`~ gender + BMI + age_category` separately within controls and
    T2D**, using each feature’s observed values. Gender and age category
    are categorical; BMI is continuous. Preserve missing residuals.
    Models are fitted once in the original groups, before permutations.
5.  Run DISCO on the residuals, using complete observations separately
    for each pair and requiring at least **20 observations in each
    group**. A sample missing one feature can contribute to other pairs.

The saved run uses **1000 permutations**, 12 workers, seed 1, BY
adjustment, feature q \< 0.2, and global p \< 0.05. Age categories are
adult, schoolage, and senior. Each group has only one schoolage
participant; a feature with an unestimable adjustment model is left
missing and its model status is recorded.

### Retrieve and preprocess the data

Install the data dependencies if needed:

``` r
install.packages("BiocManager")
BiocManager::install(c("curatedMetagenomicData", "SummarizedExperiment"))
```

The following blocks are imported directly from the [analysis
script](analysis/run_qinj_readme_examples.R) and [preprocessing
helpers](analysis/qinj_readme_preprocessing.R).

``` r
settings <- list(
  workers = 12L,
  n_resamples = 1000L, random_seed = 1L, cooccur = 20L, adj_method = "BY",
  q_value = 0.2, global_p_cutoff = 0.05, prevalence_cutoff = 0.2,
  transformation = "robust CLR; zeros remain NA",
  prevalence_scope = "at least 20% in each of healthy and T2D",
  covariate_adjustment = "separate healthy and T2D models: ~ gender + BMI + age_category"
)
results_dir <- file.path(analysis_dir, "qinj_readme_results")
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
source_md5 <- readme_source_hashes(repo_dir)
```

The selection helper handles independent participants and complete
covariates:

``` r
# Select one independent, gender/BMI/age-category-complete sample per QinJ participant.
select_qinj_metadata <- function(metadata) {
  metadata <- as.data.frame(metadata)
  required <- c("study_name", "sample_id", "subject_id", "disease", "gender", "BMI", "age_category")
  if (!all(required %in% names(metadata))) stop("Required QinJ metadata columns are missing.")
  metadata <- metadata[!is.na(metadata$study_name) & metadata$study_name == "QinJ_2012" &
    metadata$disease %in% c("healthy", "T2D"), , drop = FALSE]
  if (!nrow(metadata) || anyNA(metadata$sample_id) || any(!nzchar(metadata$sample_id)) ||
      anyDuplicated(metadata$sample_id)) stop("Sample IDs must be non-missing and unique.")
  metadata$subject_id <- trimws(as.character(metadata$subject_id))
  if (anyNA(metadata$subject_id) || any(!nzchar(metadata$subject_id))) stop("Subject IDs are required.")
  if (!is.numeric(metadata$BMI)) stop("BMI must be numeric.")
  conflicting <- names(which(vapply(split(metadata$disease, metadata$subject_id),
    function(x) length(unique(x)) > 1L, logical(1))))
  conflict <- metadata[metadata$subject_id %in% conflicting, , drop = FALSE]
  conflict$exclusion_reason <- rep("Subject has conflicting disease groups", nrow(conflict))
  metadata <- metadata[!metadata$subject_id %in% conflicting, , drop = FALSE]
  complete <- stats::complete.cases(metadata[, c("gender", "BMI", "age_category")]) &
    is.finite(metadata$BMI) & nzchar(trimws(as.character(metadata$gender))) &
    nzchar(trimws(as.character(metadata$age_category)))
  excluded_covariates <- metadata[!complete, , drop = FALSE]
  excluded_covariates$exclusion_reason <- rep("Missing or invalid gender, BMI, or age category", nrow(excluded_covariates))
  metadata <- metadata[complete, , drop = FALSE]
  day <- if ("days_from_first_collection" %in% names(metadata))
    metadata$days_from_first_collection else rep(NA_real_, nrow(metadata))
  if (!is.numeric(day)) stop("Collection day must be numeric.")
  metadata <- metadata[order(metadata$subject_id, day, metadata$sample_id, na.last = TRUE), , drop = FALSE]
  repeated <- duplicated(metadata$subject_id)
  excluded_repeats <- metadata[repeated, , drop = FALSE]
  excluded_repeats$exclusion_reason <- rep("Repeated subject: retained earliest covariate-complete sample",
    nrow(excluded_repeats))
  metadata <- metadata[!repeated, , drop = FALSE]
  rownames(metadata) <- NULL # returnSamples() requires ordinary row indices.
  list(metadata = metadata, excluded_covariates = excluded_covariates,
    excluded_subjects = rbind(conflict, excluded_repeats))
}
```

``` r
selection <- select_qinj_metadata(curatedMetagenomicData::sampleMetadata)
metadata <- selection$metadata
excluded_covariates <- selection$excluded_covariates
sample_metadata <- metadata
rownames(sample_metadata) <- sample_metadata$sample_id
message("Excluded ", nrow(excluded_covariates), " sample(s) missing gender, BMI or age category; ",
        nrow(selection$excluded_subjects), " conflicting or repeated sample(s).")
taxa_se <- curatedMetagenomicData::returnSamples(
  metadata, dataType = "relative_abundance", counts = FALSE, rownames = "short"
)
taxa <- t(as.matrix(SummarizedExperiment::assay(taxa_se)))
taxa_ids <- metadata$sample_id[metadata$sample_id %in% rownames(taxa)]
taxa <- taxa[taxa_ids, , drop = FALSE]
if (any(!is.finite(taxa)) || any(taxa < 0)) stop("Invalid taxon abundances.")
taxa <- taxa[rowSums(taxa) > 0, , drop = FALSE]
group <- stats::setNames(as.character(metadata$disease), metadata$sample_id)
taxa_group <- group[rownames(taxa)]
if (any(table(factor(taxa_group, levels = c("healthy", "T2D"))) < settings$cooccur)) {
  stop("Each group must have at least cooccur usable samples.")
}
```

Prevalence filtering and robust CLR are applied independently to each
modality:

``` r
# Samples are rows; features are columns. Apply independently per modality.
preprocess_abundance <- function(x, group, prevalence_cutoff = 0.10) {
  stopifnot(
    is.matrix(x), is.numeric(x), nrow(x) > 0L, ncol(x) > 0L,
    all(is.na(x) | is.finite(x)), all(x >= 0, na.rm = TRUE),
    all(rowSums(x, na.rm = TRUE) > 0),
    !is.null(rownames(x)), !anyDuplicated(rownames(x)),
    !is.null(colnames(x)), !anyDuplicated(colnames(x)),
    length(prevalence_cutoff) == 1L, is.finite(prevalence_cutoff),
    prevalence_cutoff > 0, prevalence_cutoff <= 1,
    identical(rownames(x), names(group)), !anyNA(group),
    setequal(as.character(group), c("healthy", "T2D"))
  )

  # Require the prevalence cutoff separately in BOTH original groups.
  x <- sweep(x, 1, rowSums(x, na.rm = TRUE), "/")
  prevalence <- cbind(
    healthy = colMeans(x[group == "healthy", , drop = FALSE] > 0, na.rm = TRUE),
    T2D = colMeans(x[group == "T2D", , drop = FALSE] > 0, na.rm = TRUE)
  )
  keep <- rowSums(is.finite(prevalence) & prevalence >= prevalence_cutoff) == 2L
  if (sum(keep) < 2L) stop("Fewer than two features pass the prevalence filter.")
  retained <- x[, keep, drop = FALSE]

  # Robust CLR: centre each sample using only its observed positive features.
  # Zeros and existing missing values stay missing; do not impute them.
  log_x <- retained
  log_x[!is.na(log_x) & log_x == 0] <- NA_real_
  log_x <- log(log_x)
  centre <- rowMeans(log_x, na.rm = TRUE)
  centre[!is.finite(centre)] <- NA_real_
  clr <- sweep(log_x, 1, centre, "-")

  list(
    clr = clr,
    retained_relative_abundance = retained,
    transformation = "robust CLR; zeros remain NA",
    prevalence = data.frame(
      feature = rownames(prevalence), prevalence, retained = unname(keep),
      row.names = NULL
    )
  )
}
```

Each feature is adjusted separately within each group using its observed
values:

``` r
# Fit independent gender/BMI/age-category models within the original disease groups.
residualize_abundance <- function(clr, metadata, group) {
  stopifnot(identical(rownames(clr), rownames(metadata)),
            all(is.na(clr) | is.finite(clr)),
            is.numeric(metadata$BMI),
            all(is.finite(metadata$BMI)),
            identical(rownames(clr), names(group)), !anyNA(group),
            setequal(as.character(group), c("healthy", "T2D")))
  residuals <- matrix(NA_real_, nrow(clr), ncol(clr), dimnames = dimnames(clr))
  models <- list()
  for (group_name in c("healthy", "T2D")) {
    rows <- which(group == group_name)
    model_data <- metadata[rows, , drop = FALSE]
    model_data$gender <- factor(trimws(as.character(model_data$gender)))
    model_data$age_category <- factor(trimws(as.character(model_data$age_category)))
    if (anyNA(model_data$gender) || any(!nzchar(as.character(model_data$gender))) ||
        nlevels(model_data$gender) < 2L) {
      stop("At least two non-missing gender levels are required in ", group_name, ".")
    }
    if (anyNA(model_data$age_category) || any(!nzchar(as.character(model_data$age_category))) ||
        nlevels(model_data$age_category) < 2L) {
      stop("At least two non-missing age categories are required in ", group_name, ".")
    }
    design <- stats::model.matrix(~ gender + BMI + age_category, data = model_data)
    if (qr(design)$rank < ncol(design) || nrow(design) <= ncol(design)) {
      stop("Covariate model requires full rank and positive residual degrees of freedom in ",
           group_name, ".")
    }
    coefficients <- matrix(NA_real_, ncol(design), ncol(clr),
                           dimnames = list(colnames(design), colnames(clr)))
    n_observed <- residual_df <- stats::setNames(integer(ncol(clr)), colnames(clr))
    status <- stats::setNames(rep("ok", ncol(clr)), colnames(clr))
    for (i in seq_len(ncol(clr))) {
      observed <- !is.na(clr[rows, i])
      n_observed[i] <- sum(observed)
      feature_design <- design[observed, , drop = FALSE]
      if (sum(observed) <= ncol(design)) {
        status[i] <- "insufficient_observations"
        residual_df[i] <- NA_integer_
        next
      }
      if (qr(feature_design)$rank < ncol(design)) {
        status[i] <- "rank_deficient"
        residual_df[i] <- NA_integer_
        next
      }
      fit <- stats::lm.fit(feature_design, clr[rows[observed], i])
      coefficients[, i] <- fit$coefficients
      residuals[rows[observed], i] <- fit$residuals
      residual_df[i] <- fit$df.residual
    }
    models[[group_name]] <- list(
      coefficients = coefficients, design = design,
      sample_ids = rownames(model_data), residual_df = residual_df,
      n_observed = n_observed, status = status
    )
    if (any(status != "ok")) {
      warning(sum(status != "ok"), " feature model(s) could not be fitted in ",
              group_name, "; residuals left NA. See models[[group]]$status.", call. = FALSE)
    }
  }
  stopifnot(all(is.na(residuals) | is.finite(residuals)),
            all(is.na(residuals[is.na(clr)])))

  list(residuals = residuals, models = models)
}
```

Run taxa-only DISCO using the adjusted residuals:

``` r
taxa_prepared <- preprocess_abundance(taxa, taxa_group, settings$prevalence_cutoff)
taxa_adjusted <- residualize_abundance(
  taxa_prepared$clr, sample_metadata[rownames(taxa_prepared$clr), , drop = FALSE],
  taxa_group
)
single_dir <- file.path(results_dir, "taxa")
single_input <- list(
  data1 = list(taxa_adjusted$residuals[taxa_group == "healthy", , drop = FALSE]),
  data2 = list(taxa_adjusted$residuals[taxa_group == "T2D", , drop = FALSE]),
  settings = settings, source_md5 = source_md5
)
single <- readme_cached_result(single_dir, single_input, recompute)
if (is.null(single)) {
  single <- DISCO::disco(
    data1 = single_input$data1, data2 = single_input$data2,
    workers = settings$workers,
    n_resamples = settings$n_resamples, random_seed = settings$random_seed,
    q_value = settings$q_value, global_p_cutoff = settings$global_p_cutoff,
    cooccur = settings$cooccur, adj_method = settings$adj_method,
    output_dir = single_dir,
    prefix = "QinJ_2012_taxa"
  )
  saveRDS(single_input, file.path(single_dir, "disco_inputs.rds"))
  saveRDS(single, file.path(single_dir, "disco_analysis.rds"))
}
```

### DISCO results

The saved analysis excluded **19** samples with incomplete covariates
and **0** conflicting or repeated samples. Package versions: DISCO
0.0.0.9000 and curatedMetagenomicData 3.14.0.

| control | T2D | retained_taxa |
|--------:|----:|--------------:|
|     174 | 170 |           128 |

**7 taxa have DISCO q \< 0.2.** The table includes all significant taxa,
ranked by q and then raw p. The DISCO score combines
differential-correlation evidence; the p and q columns are feature-level
results, and “Pairs” counts contributing feature pairs.

| Taxon                        | DISCO score | DISCO p  | DISCO q  | Pairs |
|:-----------------------------|:------------|:---------|:---------|------:|
| Bilophila wadsworthia        | 84.75       | 3.76e-03 | 1.87e-01 |    12 |
| Bacteroides thetaiotaomicron | 74.01       | 4.30e-03 | 1.87e-01 |    14 |
| Bacteroides uniformis        | 72.58       | 4.39e-03 | 1.87e-01 |    17 |
| Alistipes putredinis         | 68.21       | 4.67e-03 | 1.87e-01 |    13 |
| Flavonifractor plautii       | 54.22       | 5.87e-03 | 1.87e-01 |    13 |
| Phocaeicola vulgatus         | 50.94       | 6.25e-03 | 1.87e-01 |     7 |
| Parabacteroides distasonis   | 45.60       | 6.98e-03 | 1.87e-01 |    16 |

These scores measure changes in taxon–taxon relationships after
group-specific covariate adjustment. They are not differential-abundance
or causal-effect tests.

### Diagnostics for threshold selection

The global permutation p-value is **0.01299**; the selected
correlation-filtering threshold is **0.01**. The global test passes, so
DISCO selects the threshold using the maximum real-to-null discovery
ratio, with the documented tie rule.

The discovery plot compares observed discovery percentages with the
permutation mean ± SD. The ratio plot shows the criterion used for
threshold selection.

<img src="man/figures/README-qinj-taxa-discovery_plot.png" alt="" width="100%" />

<img src="man/figures/README-qinj-taxa-ratio_plot.png" alt="" width="100%" />

The feature diagnostic shows −log10(raw DISCO p) against the percentage
of permutations with DISCO q \< 0.2, separately at each filtering
threshold. Orange points are real-data discoveries. Null discovery
percentages are diagnostics, not formal feature-level FDR estimates.

<img src="man/figures/README-qinj-taxa-feature_null_plot.png" alt="" width="100%" />

### Taxon–taxon network

The network highlights the **top five significant taxa**, ranked by
DISCO q, then p. Edges are selected by **raw differential-correlation p
\< 0.05** from the saved pair table after DISCO’s selected correlation
filter. This edge rule is separate from the feature-level DISCO q
cutoff. The panels share nodes, edges, and positions; red/blue indicate
positive/negative group correlations, and edge width shows \|Spearman
r\|. Edge presence does not imply within-group correlation significance.

<img src="man/figures/README-qinj-taxa-network.png" alt="" width="100%" />

### Example pairwise scatterplot

The highest-ranked retained pair is **Bacteroides thetaiotaomicron
versus Bacteroides uniformis**, ranked by differential-correlation q,
then p. The plot uses the exact adjusted inputs passed to DISCO and only
complete pairs. Its annotations report saved Spearman correlations and
paired sample counts; dashed linear fits are visual guides. The
displayed pair p/q are distinct from the feature-level DISCO results
above.

<img src="man/figures/README-qinj-taxa-scatterplot.png" alt="" width="100%" />

## Example: Two modalities — taxon–pathway correlations

The same T2D/control comparison now tests **taxon–pathway
correlations**. Taxa are the first modality, so DISCO assigns scores to
taxa according to changes in their relationships with pathways. Pathways
are unstratified community totals, excluding `UNMAPPED` and
`UNINTEGRATED`.

Restrict both modalities to common usable participants, then repeat the
≥20% per-group prevalence filter and robust CLR separately for taxa and
pathways. Adjust **both modalities** for gender, BMI, and age category
in separate control/T2D models, preserving missing values throughout.

``` r
pathway_se <- curatedMetagenomicData::returnSamples(metadata, dataType = "pathway_abundance")
pathways <- t(as.matrix(SummarizedExperiment::assay(pathway_se)))
if (any(!is.finite(pathways)) || any(pathways < 0)) stop("Invalid pathway abundances.")
pathway_ids <- colnames(pathways)
pathways <- pathways[, !grepl("|", pathway_ids, fixed = TRUE) &
                       !toupper(trimws(pathway_ids)) %in% c("UNMAPPED", "UNINTEGRATED"),
                     drop = FALSE]
common_ids <- rownames(taxa)[rownames(taxa) %in% rownames(pathways)]
common_ids <- common_ids[rowSums(pathways[common_ids, , drop = FALSE]) > 0]
if (!length(common_ids)) stop("No usable common samples.")
two_group <- group[common_ids]
if (any(table(factor(two_group, levels = c("healthy", "T2D"))) < settings$cooccur)) {
  stop("Each group requires at least cooccur samples.")
}
two_taxa <- preprocess_abundance(taxa[common_ids, , drop = FALSE], two_group,
                                 settings$prevalence_cutoff)
two_pathways <- preprocess_abundance(pathways[common_ids, , drop = FALSE], two_group,
                                     settings$prevalence_cutoff)
two_metadata <- sample_metadata[common_ids, , drop = FALSE]
two_taxa_adjusted <- residualize_abundance(two_taxa$clr, two_metadata, two_group)
two_pathways_adjusted <- residualize_abundance(two_pathways$clr, two_metadata, two_group)
two_dir <- file.path(results_dir, "taxa_pathways")
two_input <- list(
  data1 = list(two_taxa_adjusted$residuals[two_group == "healthy", , drop = FALSE],
               two_pathways_adjusted$residuals[two_group == "healthy", , drop = FALSE]),
  data2 = list(two_taxa_adjusted$residuals[two_group == "T2D", , drop = FALSE],
               two_pathways_adjusted$residuals[two_group == "T2D", , drop = FALSE]),
  settings = settings, source_md5 = source_md5
)
two <- readme_cached_result(two_dir, two_input, recompute)
if (is.null(two)) {
  message("Running two modalities: ", ncol(two_taxa$clr), " taxa, ",
          ncol(two_pathways$clr), " pathways; control n = ",
          sum(two_group == "healthy"), "; T2D n = ", sum(two_group == "T2D"), ".")
  two <- DISCO::disco(
    data1 = two_input$data1, data2 = two_input$data2,
    workers = settings$workers,
    n_resamples = settings$n_resamples, random_seed = settings$random_seed,
    q_value = settings$q_value, global_p_cutoff = settings$global_p_cutoff,
    cooccur = settings$cooccur, adj_method = settings$adj_method,
    output_dir = two_dir, prefix = "QinJ_2012_taxa_pathways"
  )
  saveRDS(two_input, file.path(two_dir, "disco_inputs.rds"))
  saveRDS(two, file.path(two_dir, "disco_analysis.rds"))
}
```

### DISCO results

| control | T2D | retained_taxa | retained_pathways |
|--------:|----:|--------------:|------------------:|
|     174 | 170 |           128 |               385 |

**14 taxa have DISCO q \< 0.2.** The table lists all significant taxa,
using the same columns and ordering as the one-modality analysis. Here
“Pairs” counts contributing taxon–pathway pairs.

| Taxon                           | DISCO score | DISCO p  | DISCO q  | Pairs |
|:--------------------------------|:------------|:---------|:---------|------:|
| Bilophila wadsworthia           | 852.73      | 3.73e-04 | 7.18e-02 |    59 |
| \[Ruminococcus\] gnavus         | 229.21      | 1.39e-03 | 8.95e-02 |   116 |
| Blautia wexlerae                | 221.60      | 1.44e-03 | 8.95e-02 |    62 |
| Agathobaculum butyriciproducens | 171.07      | 1.86e-03 | 8.95e-02 |    98 |
| Prevotella copri                | 89.78       | 3.55e-03 | 1.36e-01 |   133 |
| Eggerthella lenta               | 50.68       | 6.28e-03 | 1.63e-01 |    48 |
| Gemmiger formicilis             | 40.99       | 7.76e-03 | 1.63e-01 |    16 |
| Intestinibacter bartlettii      | 39.42       | 8.07e-03 | 1.63e-01 |     2 |
| Lachnospira pectinoschiza       | 38.81       | 8.20e-03 | 1.63e-01 |     8 |
| Ruminococcus bromii             | 34.48       | 9.23e-03 | 1.63e-01 |    35 |
| Bacteroides uniformis           | 34.12       | 9.33e-03 | 1.63e-01 |   121 |
| Bacteroides fragilis            | 29.02       | 1.10e-02 | 1.73e-01 |    50 |
| Parabacteroides distasonis      | 27.28       | 1.17e-02 | 1.73e-01 |   122 |
| Phocaeicola plebeius            | 25.15       | 1.27e-02 | 1.74e-01 |    49 |

### Diagnostics for threshold selection

The global permutation p-value is **0.000999**; the selected threshold
is **0.01**. The global test passes and DISCO selects the threshold
using the real-to-null discovery ratio. These diagnostics are calculated
separately from taxa-only DISCO.

<img src="man/figures/README-qinj-taxa_pathways-discovery_plot.png" alt="" width="100%" />

<img src="man/figures/README-qinj-taxa_pathways-ratio_plot.png" alt="" width="100%" />

<img src="man/figures/README-qinj-taxa_pathways-feature_null_plot.png" alt="" width="100%" />

### Taxon–pathway network

The bipartite network links the **top five significant DISCO taxa** to
pathways with **raw differential-correlation p \< 0.05**, using the
saved filtered pair table. Gold nodes are taxa and grey nodes are
pathways. Pathways are labelled by their IDs; full names are exported in
[network_nodes.csv](analysis/qinj_taxa_pathway_networks/network_nodes.csv).
Both panels use identical nodes, edges, and positions. Edge colors and
widths show the group-specific Spearman correlations, as in the
taxa-only network.

<img src="man/figures/README-qinj-taxa_pathways-network.png" alt="" width="100%" />

### Example pairwise scatterplot

The highest-ranked retained pair is **Bilophila wadsworthia versus
PWY-6612** (superpathway of tetrahydrofolate biosynthesis). The plotted
values are the separately adjusted robust CLR taxon and pathway
measurements, with complete pairs only. Saved group Spearman
correlations, sample counts, and pair-level differential p/q are shown;
dashed lines remain visual guides.

<img src="man/figures/README-qinj-taxa_pathways-scatterplot.png" alt="" width="100%" />

## Overall example results

Across the two analyses, **18 unique taxa** have DISCO q \< 0.2 in at
least one analysis: seven in taxa–taxa DISCO, 14 in taxon–pathway DISCO,
and three in both. The heatmap compares the **raw feature-level DISCO
p-values** for this union, including each taxon’s nonsignificant result
in the other column.

Colors show **−log10(raw DISCO p)** without row scaling. An asterisk
marks **DISCO q \< 0.2** in that cell. Rows are hierarchically clustered
by Euclidean distance on the displayed values using complete linkage;
the two analysis columns stay in a fixed order. Thus the stars summarize
adjusted significance, while the colors retain the raw p-value
magnitude.

<img src="man/figures/README-qinj-DISCO-heatmap.png" alt="" width="100%" />

### Rebuilding these examples

`README.Rmd` is the source for this README. Rendering reads the saved
QinJ summary and results under `analysis/qinj_readme_results/` and the
portable figures under `man/figures/`; it does not download data or
rerun permutations. Analysis blocks are imported from their source
scripts, and rendering checks that their hashes match the saved run.

``` r
rmarkdown::render("README.Rmd")
```

Regenerate the analysis with
`Rscript analysis/run_qinj_readme_examples.R` (or `source()` with its
full path). Set `DISCO_README_RECOMPUTE=true` to force recomputation;
matching cached analyses are otherwise reused. Recreate the additional
visualizations with:

``` sh
Rscript analysis/plot_qinj_taxa_networks.R
Rscript analysis/plot_qinj_top_differential_scatterplots.R
Rscript analysis/plot_qinj_taxa_pathway_networks.R
Rscript analysis/plot_qinj_top_taxa_pathway_scatterplots.R
Rscript analysis/plot_qinj_disco_heatmap.R
```

Refresh the portable copies used in this README after regenerating the
plots:

``` r
diagnostic_names <- paste0("README-qinj-",
  rep(c("taxa", "taxa_pathways"), each = 3), "-",
  rep(c("discovery_plot", "ratio_plot", "feature_null_plot"), 2), ".png")
assets <- stats::setNames(
  file.path("analysis/qinj_readme_results/figures", diagnostic_names), diagnostic_names)
assets <- c(assets,
  "README-qinj-taxa-network.png" = "analysis/qinj_taxa_networks/top5_taxa_combined_control_T2D_networks.png",
  "README-qinj-taxa-scatterplot.png" = "analysis/qinj_differential_scatterplots/top10_page01_preview.png",
  "README-qinj-taxa_pathways-network.png" = "analysis/qinj_taxa_pathway_networks/top5_taxa_pathways_control_T2D_networks.png",
  "README-qinj-taxa_pathways-scatterplot.png" = "analysis/qinj_taxa_pathway_scatterplots/top1_page01_preview.png",
  "README-qinj-DISCO-heatmap.png" = "analysis/qinj_disco_heatmap/significant_taxa_DISCO_heatmap.png")
stopifnot(all(file.exists(assets)))
file.copy(assets, file.path("man/figures", names(assets)), overwrite = TRUE)
rmarkdown::render("README.Rmd")
```

The scripts save PNG/PDF plots and the underlying CSV tables in their
respective `analysis/qinj_*` directories. The network scripts expose
`top_n`; the heatmap script exposes `cluster_rows` and `row_clusters`.
Numerical preprocessing and adjustment checks can be run with
`Rscript analysis/test_qinj_readme_helpers.R`.

## License

DISCO is licensed under the MIT License.

## Citation

If you use DISCO, please cite the paper that introduced the method:

Fouladi F, Chen Y, Bera S, Jarmusch AK, Van Tyne D, Palella FJ,
Margolick JB, Chew KW, Sun J, Martinson J, Rinaldo CR, Peddada SD. A
taxon-specific measurement of disruption in a multi-modal study of
microbiomes and metabolomes reveals system-wide dysbiosis preceding
HIV-1 infection. Nat Commun. 2025 Nov 20;16(1):10204.
[doi:10.1038/s41467-025-64822-z](https://doi.org/10.1038/s41467-025-64822-z).
PMID: 41266305; PMCID: PMC12635353.
