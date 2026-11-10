
# DISCO

<!-- README.md is generated from README.Rmd. Please edit README.Rmd. -->

DISCO (**DIS**ruption in **CO**rrelations) is an exploratory tool that
attempts to identify features whose correlation patterns differ between
two groups, such as healthy controls and people with a disease. It
combines evidence from pairwise differential-correlation tests into a
score, p-value, and adjusted p-value (q-value) for each feature. This
helps identify features whose relationships with other features differ
between groups, even when their average abundances are similar.

DISCO supports one modality (for example, correlations among microbial
taxa) or two modalities (for example, correlations between microbes and
metabolites). It can be applied to high-dimensional data, including
microbiome, transcriptomic, and metabolomic measurements. Normalize data
appropriately for the data type before applying DISCO.

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

DISCO matches features between groups by column name within each
modality. It keeps shared features and reorders group 2 to match group
1.

For two modalities, row names must identify samples consistently. Within
each group, the modalities must share samples, but their sample sets and
row order may differ. DISCO aligns rows by sample name and represents
measurements from an unavailable modality as missing. Each feature pair
must have at least `cooccur` paired, nonmissing observations **in each
group** to be tested.

**Normalize inputs appropriately before running DISCO.** The function
does not normalize data or adjust for covariates. Supply normalized and,
when needed, covariate-adjusted measurements. In the examples below, we
use a robust centred log-ratio (CLR) transformation for microbiome data:
each sample is centred over its observed positive features, and
abundance zeros are represented as `NA`. DISCO then uses complete
observations separately for each feature pair. Other data types may
require different preprocessing. DISCO itself treats numeric zeros as
observed values; missing values are excluded separately for each pair.
If zeros should be excluded from the correlation calculations, represent
them as `NA` before running DISCO.

| Argument | Default | Meaning |
|----|---:|----|
| `cooccur` | `20` | Minimum number of pairwise nonmissing observations in **each** group; integer ≥10. |
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
high-dimensional data, the number of feature pairs can be large.
Screening out weak correlations can reduce the number of uninformative
pairs contributing to a score and may increase power, although it can
also increase false positives. `filter_threshold` screens pairs using
their **single-group correlation p-values**, not their
differential-correlation p-values. A pair is retained if its correlation
p-value is below the threshold in **at least one group**. For example,
`filter_threshold = 0.05` excludes pairs whose correlation p-values are
≥0.05 in both groups. A threshold of `1` is the least restrictive
setting.

To choose a threshold empirically, DISCO
(`empirically_derive_threshold = TRUE`, the default) compares the
observed data with size-matched random label permutations at filtering
thresholds `c(0.001, 0.005, 0.01, 0.05, 0.1, 0.5, 1)`. At each
threshold, it calculates the percentage of tested features with DISCO
`q < q_value` in the observed data and in each permutation. A global
permutation test then compares the maximum observed discovery percentage
across thresholds with the corresponding maximum from each valid
permutation:

``` text
global p = (1 + number of null maxima >= real maximum) / (1 + valid permutations)
```

If `global p < global_p_cutoff`, DISCO selects the threshold that
maximizes `max(real percentage, 1) / max(mean null percentage, 1)`. The
floors are in percentage points and prevent ratios from being dominated
by near-zero values. Ties select the largest threshold. If
`global p >= global_p_cutoff`, DISCO uses threshold `1`.

For two modalities, permutations preserve each modality’s group sizes
and sample overlap. A shared sample receives the same shuffled group
label in both modalities.

To use a fixed threshold, set `empirically_derive_threshold = FALSE` and
supply `filter_threshold` (default `1`). This skips permutations and
diagnostic plots. Here is a self-contained toy example with 30 samples
per group and five features. Feature 1 and feature 2 are correlated in
controls but independent in cases. These simulated values represent
preprocessed data.

``` r
# A toy example
library(DISCO)

set.seed(1)
controls <- matrix(rnorm(30 * 5), nrow = 30, ncol = 5)
cases <- matrix(rnorm(30 * 5), nrow = 30, ncol = 5)
colnames(controls) <- colnames(cases) <- paste0("feature", 1:5)
controls[, 2] <- controls[, 1] + rnorm(30, sd = 0.2)

result <- disco(
  data1 = list(controls),
  data2 = list(cases),
  cooccur = 20,
  empirically_derive_threshold = FALSE, # empirical mode is turned off
  filter_threshold = 1,
  q_value = 0.2,       # Used only in empirical mode
  n_resamples = 100,  # Used only in empirical mode
  random_seed = 1L,   # Used only in empirical mode
  adj_method = "BY",
  output_dir = "disco_results",
  prefix = "controls_vs_cases",
  global_p_cutoff = 0.05, # Used only in empirical mode
  workers = 1L           # Used only in empirical mode
)
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
whether the real-data feature q-value is below `0.2`. The column name
reflects the supplied `q_value`. For each feature, its null discovery
frequency is the percentage of all requested permutations in which it
has DISCO q \< `q_value`. Permutations in which that feature is absent
or untested contribute zero to its discovery count.

`feature_null_plot` shows **−log10(real-data p-value)** on the x-axis
and the percentage of permutations in which the feature has
`q < q_value` on the y-axis. Orange points identify real-data
discoveries. Low null discovery frequencies provide a useful diagnostic,
but are not formal feature-level FDR estimates.

In fixed mode, `cor_diff`, `disco`, `selected_threshold`, and
`output_files` are returned; the diagnostic components are `NULL`.
Pair-level and feature-level tables are saved as `<prefix>_cor_diff.csv`
and `<prefix>_disco.csv`.

## Example of one data modality: taxon–taxon correlations in type 2 diabetes

This example uses **QinJ_2012** from [Qin et
al. (2012)](https://doi.org/10.1038/nature11450), comparing type 2
diabetes (**T2D**) with healthy controls in `curatedMetagenomicData`.
Both the one-modality and two-modality examples use this same
comparison.

Dataset reference: Qin J, Li Y, Cai Z, et al. **A metagenome-wide
association study of gut microbiota in type 2 diabetes.** *Nature*.
2012;490:55–60.
[doi:10.1038/nature11450](https://doi.org/10.1038/nature11450).

The workflow follows these steps:

1.  Select one independent sample per participant, excluding conflicting
    disease labels and samples missing **gender or BMI**.
2.  Retain features present in **≥50% of each group**, testing
    prevalence separately in controls and T2D.
3.  Apply **robust CLR** within each sample, centring log positive
    abundances over the retained observed features. Abundance zeros and
    existing missing values remain `NA`; no zero imputation is
    performed.
4.  Fit **`~ gender + BMI` separately within controls and T2D**, using
    each feature’s observed values. Gender is categorical and BMI is
    continuous. Preserve missing residuals.
5.  Run DISCO on the residuals, using complete observations separately
    for each pair and requiring at least **20 observations in each
    group**.

The saved run uses **1000 permutations**, 12 workers, seed 1, BY
adjustment, feature q \< 0.2, and global p \< 0.05.

### Retrieve and preprocess the data

Install the data dependencies if needed:

``` r
install.packages("BiocManager")
BiocManager::install(c("curatedMetagenomicData", "SummarizedExperiment"))
```

The following blocks are imported directly from the [analysis
script](analysis/run_qinj_readme_examples.R) and [preprocessing
helpers](analysis/qinj_readme_preprocessing.R).

Load the installed DISCO package (see Installation above):

``` r
library(DISCO)
```

``` r
settings <- list(
  workers = 12L,
  n_resamples = 1000L, random_seed = 1L, cooccur = 20L, adj_method = "BY",
  q_value = 0.2, global_p_cutoff = 0.05, prevalence_cutoff = 0.5,
  transformation = "robust CLR; zeros remain NA",
  prevalence_scope = "at least 50% in each of healthy and T2D",
  covariate_adjustment = "separate healthy and T2D models: ~ gender + BMI"
)
results_dir <- file.path(analysis_dir, "qinj_readme_results")
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
source_md5 <- readme_source_hashes(repo_dir)
```

The selection helper handles independent participants and complete
covariates:

``` r
# Select one independent, gender/BMI-complete sample per QinJ participant.
select_qinj_metadata <- function(metadata) {
  metadata <- as.data.frame(metadata)
  required <- c("study_name", "sample_id", "subject_id", "disease", "gender", "BMI")
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
  complete <- stats::complete.cases(metadata[, c("gender", "BMI")]) &
    is.finite(metadata$BMI) & nzchar(trimws(as.character(metadata$gender)))
  excluded_covariates <- metadata[!complete, , drop = FALSE]
  excluded_covariates$exclusion_reason <- rep("Missing or invalid gender or BMI", nrow(excluded_covariates))
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
message("Excluded ", nrow(excluded_covariates), " sample(s) missing gender or BMI; ",
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
# Fit independent gender/BMI models within the original disease groups.
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
    if (anyNA(model_data$gender) || any(!nzchar(as.character(model_data$gender))) ||
        nlevels(model_data$gender) < 2L) {
      stop("At least two non-missing gender levels are required in ", group_name, ".")
    }
    design <- stats::model.matrix(~ gender + BMI, data = model_data)
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

Package versions: DISCO 0.0.0.9000 and curatedMetagenomicData 3.14.0.
The saved analysis excluded **19** samples with incomplete covariates.

| control | T2D | retained_taxa |
|--------:|----:|--------------:|
|     174 | 170 |            62 |

**16 taxa have DISCO q \< 0.2.** The table lists these taxa. “Pairs”
counts the feature pairs contributing to each taxon’s DISCO score. These
taxa show evidence of altered correlation patterns between controls and
T2D cases.

| Taxon                            | DISCO p  | DISCO q  | Pairs |
|:---------------------------------|:---------|:---------|------:|
| Lachnospira eligens              | 1.08e-03 | 1.78e-01 |     8 |
| Bacteroides intestinalis         | 1.78e-03 | 1.78e-01 |     6 |
| Collinsella aerofaciens          | 2.57e-03 | 1.78e-01 |     7 |
| Prevotella copri                 | 4.24e-03 | 1.78e-01 |     6 |
| Anaerostipes hadrus              | 4.68e-03 | 1.78e-01 |    24 |
| Alistipes shahii                 | 5.85e-03 | 1.78e-01 |    12 |
| Bifidobacterium longum           | 6.63e-03 | 1.78e-01 |    11 |
| Bacteroides uniformis            | 7.03e-03 | 1.78e-01 |    28 |
| Bilophila wadsworthia            | 7.66e-03 | 1.78e-01 |    17 |
| Agathobaculum butyriciproducens  | 8.40e-03 | 1.78e-01 |    23 |
| Fusicatenibacter saccharivorans  | 8.61e-03 | 1.78e-01 |    15 |
| Butyricimonas virosa             | 8.66e-03 | 1.78e-01 |    14 |
| Bacteroides thetaiotaomicron     | 8.98e-03 | 1.78e-01 |    30 |
| Blautia wexlerae                 | 9.46e-03 | 1.78e-01 |    10 |
| Parasutterella excrementihominis | 9.50e-03 | 1.78e-01 |     3 |
| Streptococcus salivarius         | 1.13e-02 | 1.98e-01 |     7 |

### Diagnostics for threshold selection

The global permutation p-value is **0.001998**. The global test detects
excess discoveries, so DISCO selects the threshold using the maximum
real-to-null discovery ratio. The selected correlation-filtering
threshold is **0.005**.

The discovery plot compares observed discovery percentages with the
permutation mean ± SD. The ratio plot shows the criterion used for
threshold selection.

<img src="man/figures/README-qinj-taxa-discovery_plot.png" alt="" width="100%" />

<img src="man/figures/README-qinj-taxa-ratio_plot.png" alt="" width="100%" />

The following feature diagnostic shows −log10(raw DISCO p) against the
percentage of permutations with DISCO q \< 0.2, separately at each
filtering threshold. Orange points are real-data discoveries. Null
discovery percentages are diagnostics, not formal feature-level FDR
estimates. We can see that at the selected threshold **0.005**, the
significant DISCO taxa have a low null discovery percentage (less than
5%).

<img src="man/figures/README-qinj-taxa-feature_null_plot.png" alt="" width="100%" />

### Taxon–taxon network

For taxa with significant DISCO scores, the `cor_diff` output shows how
their correlations with other taxa differ between groups. The network
highlights the **top ten significant taxa**, ranked by DISCO q, then
p. Edges are selected by **raw differential-correlation p \< 0.05** from
the saved pair table after DISCO’s selected correlation filter. The
panels share nodes, edges, and positions across the two groups; red/blue
indicate positive/negative group correlations, and edge width shows
\|Spearman r\|. Edge presence does not imply within-group correlation
significance.

Correlation source:
[QinJ_2012_taxa_cor_diff.csv](analysis/qinj_readme_results/taxa/QinJ_2012_taxa_cor_diff.csv).
Its columns identify the pair (`t1`, `t2`), control/T2D Spearman
correlations (`r1`, `r2`), complete-pair sample counts (`n1`, `n2`), and
differential-correlation `p` and `q`. The exact subset drawn as edges is
saved in
[network_edges.csv](analysis/qinj_taxa_networks/network_edges.csv).

<img src="man/figures/README-qinj-taxa-network.png" alt="" width="100%" />

### Example pairwise scatterplot

Scatterplots show the group-specific relationship for an individual
feature pair. The highest-ranked retained pair is **Anaerostipes hadrus
versus Lachnospira eligens**, ranked by differential-correlation q, then
p. The plot uses the exact adjusted inputs passed to DISCO and only
complete pairs. Its annotations report the saved Spearman correlations
and paired sample counts from `cor_diff`. Dashed linear fits are visual
guides.

Pair statistics come from
[QinJ_2012_taxa_cor_diff.csv](analysis/qinj_readme_results/taxa/QinJ_2012_taxa_cor_diff.csv);
the selected pairs and their annotations are exported in
[top10_pairs.csv](analysis/qinj_differential_scatterplots/top10_pairs.csv).
The plotted sample values are exported in
[scatterplot_data.csv](analysis/qinj_differential_scatterplots/scatterplot_data.csv),
extracted from the adjusted matrices in
[disco_inputs.rds](analysis/qinj_readme_results/taxa/disco_inputs.rds).

<img src="man/figures/README-qinj-taxa-scatterplot.png" alt="" width="100%" />

## Example of two data modalities: taxon–pathway correlations

The same T2D/control comparison now tests **taxon–pathway
correlations**. Taxa are the first modality, so DISCO assigns scores to
taxa according to changes in their relationships with pathways. Pathways
are unstratified community totals, excluding `UNMAPPED` and
`UNINTEGRATED`.

Restrict both modalities to common usable participants, then repeat the
≥50% per-group prevalence filter and robust CLR separately for taxa and
pathways. Adjust **both modalities** for gender and BMI in separate
control and T2D models, preserving missing values throughout.

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
|     174 | 170 |            62 |               351 |

**18 taxa have DISCO q \< 0.2.** The table lists all significant taxa,
using the same columns and ordering as the one-modality analysis. Here,
“Pairs” counts contributing taxon–pathway pairs.

| Taxon                           | DISCO p  | DISCO q  | Pairs |
|:--------------------------------|:---------|:---------|------:|
| Bilophila wadsworthia           | 9.30e-04 | 1.21e-01 |    40 |
| Collinsella aerofaciens         | 1.47e-03 | 1.21e-01 |    43 |
| Agathobaculum butyriciproducens | 1.54e-03 | 1.21e-01 |    60 |
| \[Ruminococcus\] gnavus         | 1.65e-03 | 1.21e-01 |    92 |
| Butyricimonas virosa            | 2.15e-03 | 1.26e-01 |    27 |
| Blautia wexlerae                | 2.88e-03 | 1.40e-01 |    47 |
| Bacteroides caccae              | 5.29e-03 | 1.69e-01 |    31 |
| Ruminococcus bromii             | 5.35e-03 | 1.69e-01 |     2 |
| Bifidobacterium longum          | 5.49e-03 | 1.69e-01 |    62 |
| Bacteroides uniformis           | 6.40e-03 | 1.69e-01 |   115 |
| Lachnospira pectinoschiza       | 6.81e-03 | 1.69e-01 |     9 |
| Prevotella copri                | 6.93e-03 | 1.69e-01 |   116 |
| Gordonibacter pamelaeae         | 7.90e-03 | 1.78e-01 |    11 |
| Flavonifractor plautii          | 8.51e-03 | 1.78e-01 |    84 |
| Phocaeicola plebeius            | 1.03e-02 | 1.85e-01 |    41 |
| Oscillibacter sp. 57_20         | 1.04e-02 | 1.85e-01 |    12 |
| Alistipes shahii                | 1.08e-02 | 1.85e-01 |    17 |
| Eggerthella lenta               | 1.20e-02 | 1.95e-01 |    29 |

### Diagnostics for threshold selection

The global permutation p-value is **0.000999**; the selected threshold
is **0.005**. The global test passes and DISCO selects the threshold
using the real-to-null discovery ratio. These diagnostics are calculated
independently for the taxon–pathway analysis.

<img src="man/figures/README-qinj-taxa_pathways-discovery_plot.png" alt="" width="100%" />

<img src="man/figures/README-qinj-taxa_pathways-ratio_plot.png" alt="" width="100%" />

<img src="man/figures/README-qinj-taxa_pathways-feature_null_plot.png" alt="" width="100%" />

### Taxon–pathway network

The following bipartite network links the **top five significant DISCO
taxa** to pathways with **raw differential-correlation p \< 0.05**,
using the saved filtered pair table. Gold nodes are taxa and grey nodes
are pathways. Pathways are labelled by their IDs; full names are
exported in
[network_nodes.csv](analysis/qinj_taxa_pathway_networks/network_nodes.csv).
Both panels use identical nodes, edges, and positions. Edge colors and
widths show the group-specific Spearman correlations, as in the
taxa-only network.

Correlation source:
[QinJ_2012_taxa_pathways_cor_diff.csv](analysis/qinj_readme_results/taxa_pathways/QinJ_2012_taxa_pathways_cor_diff.csv).
Here `t1` is the taxon and `t2` is the pathway; `r1`/`r2`, `n1`/`n2`,
and `p`/`q` have the same meanings as in the taxa-only pair table. The
exact plotted edge subset is exported in
[network_edges.csv](analysis/qinj_taxa_pathway_networks/network_edges.csv).

<img src="man/figures/README-qinj-taxa_pathways-network.png" alt="" width="100%" />

### Example pairwise scatterplot

The highest-ranked retained pair is **\[Ruminococcus\] gnavus versus
PWY-241** (C4 photosynthetic carbon assimilation cycle, NADP-ME type).
The plotted values are the separately adjusted robust CLR taxon and
pathway measurements, with complete pairs only. Saved group-specific
Spearman correlations, sample counts, and pair-level
differential-correlation p-values and q-values are shown; dashed lines
remain visual guides.

Pair statistics come from
[QinJ_2012_taxa_pathways_cor_diff.csv](analysis/qinj_readme_results/taxa_pathways/QinJ_2012_taxa_pathways_cor_diff.csv).
The selected example and its annotations are saved in
[top1_pairs.csv](analysis/qinj_taxa_pathway_scatterplots/top1_pairs.csv),
and its sample values are exported in
[scatterplot_data.csv](analysis/qinj_taxa_pathway_scatterplots/scatterplot_data.csv).
Those values are extracted from the adjusted taxon and pathway matrices
in
[disco_inputs.rds](analysis/qinj_readme_results/taxa_pathways/disco_inputs.rds).

<img src="man/figures/README-qinj-taxa_pathways-scatterplot.png" alt="" width="100%" />

## Overall example results

Across the two analyses, **25 unique taxa** have DISCO q \< 0.2 in at
least one analysis: 16 in taxa–taxa DISCO, 18 in taxon–pathway DISCO,
and 9 in both. The heatmap compares the **raw feature-level DISCO
p-values** for this union, including each taxon’s nonsignificant result
in the other column.

Colors show **−log10(raw DISCO p)** without row scaling. An asterisk
marks **DISCO q \< 0.2** in that cell. Rows with both results are
hierarchically clustered by Euclidean distance on the displayed values
using complete linkage. Taxa without a testable result in one analysis
are shown in a separate unclustered block, with unavailable values grey;
no p-values are imputed. The two analysis columns stay in a fixed order.
Thus the stars summarize adjusted significance, while the colors retain
the raw p-value magnitude.

<img src="man/figures/README-qinj-DISCO-heatmap.png" alt="" width="100%" />

### Rebuilding these examples

Install DISCO as shown above, then clone this repository or download its
GitHub ZIP to obtain the example scripts (they are not included in the
installed package). Run the following commands from the repository root.

Install the example dependencies in R:

``` r
install.packages(c("BiocManager", "igraph", "ggrepel", "pheatmap", "rmarkdown", "knitr"))
BiocManager::install(c("curatedMetagenomicData", "SummarizedExperiment"))
```

Rerun the analysis and generate the plots:

``` sh
Rscript analysis/run_qinj_readme_examples.R
Rscript analysis/plot_qinj_taxa_networks.R
Rscript analysis/plot_qinj_top_differential_scatterplots.R
Rscript analysis/plot_qinj_taxa_pathway_networks.R
Rscript analysis/plot_qinj_top_taxa_pathway_scatterplots.R
Rscript analysis/plot_qinj_disco_heatmap.R
```

The analysis uses `library(DISCO)` and saves results in
`analysis/qinj_readme_results/`. Adjust `settings` in the analysis
script as needed. Matching cached results are reused; set
`DISCO_README_RECOMPUTE=true` to force recomputation. Each script also
supports `source()` with its full path.

To refresh the figures embedded in this README and regenerate
`README.md`, run:

``` sh
Rscript analysis/render_qinj_readme.R
```

This last step requires Pandoc and reads saved results without rerunning
DISCO.

## Citation

If you use DISCO, please cite the paper that introduced the method:

Fouladi F, Chen Y, Bera S, Jarmusch AK, Van Tyne D, Palella FJ,
Margolick JB, Chew KW, Sun J, Martinson J, Rinaldo CR, Peddada SD. A
taxon-specific measurement of disruption in a multi-modal study of
microbiomes and metabolomes reveals system-wide dysbiosis preceding
HIV-1 infection. Nat Commun. 2025 Nov 20;16(1):10204.
[doi:10.1038/s41467-025-64822-z](https://doi.org/10.1038/s41467-025-64822-z).
PMID: 41266305; PMCID: PMC12635353.
