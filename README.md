
# DISCO

<!-- README.md is generated from README.Rmd. Please edit README.Rmd. -->

DISCO (**DIS**ruption in **CO**rrelations) identifies features whose
correlation patterns differ between two groups. It evaluates
relationships within one data modality or between two modalities and
combines pairwise tests into a feature-level DISCO score.

DISCO can be applied to **any high-dimensional data represented by
quantitative features**, including microbial abundances, gene
expression, proteins, and metabolites. It uses Spearman correlations
with Bonett–Wright standard errors and Cauchy combination of
differential-correlation p-values.

## Installation

``` r
install.packages("remotes") # If needed
remotes::install_github("FarnazFouladi/DISCO")
library(DISCO)
```

To use a local checkout during development, run
`devtools::load_all(".")` from the repository root.

## Inputs and arguments

`data1` and `data2` represent the two groups. Each is a list containing
one or two numeric matrices or data frames, with **samples in rows and
features in columns**. Supply unique feature names and, especially for
two modalities, unique sample row names. Corresponding modalities are
restricted to common features and put in a consistent order.

**Normalize inputs appropriately for the data type before running
DISCO.** The function does not normalize or adjust for covariates. CLR
is used for the compositional abundance examples below; other data types
require their own appropriate preprocessing. Supply covariate-adjusted
data if your analysis requires it. Zeros are observed values; only `NA`
reduces the pairwise-complete sample count.

| Argument | Default | Meaning |
|----|---:|----|
| `cooccur` | `10` | Minimum pairwise non-missing observations in **each** group; integer ≥4. |
| `empirically_derive_threshold` | `TRUE` | Select the correlation filtering threshold using label permutations. |
| `filter_threshold` | `0.05` | Fixed correlation p-value threshold when empirical selection is disabled. |
| `q_value` | `0.2` | Feature-level DISCO q cutoff for real and null discovery counts and plot colors. |
| `global_p_cutoff` | `0.05` | Global permutation p-value must be below this cutoff to select by ratio; otherwise use threshold `1`. |
| `n_resamples` | `100` | Number of size-matched label permutations. |
| `random_seed` | `1L` | Reproducible permutation seed; the caller’s RNG state is restored. |
| `adj_method` | `"BY"` | Multiple-testing adjustment passed to `p.adjust()`. |
| `output_dir`, `prefix` | Required | Results directory and filename prefix. |

Empirical mode evaluates `c(0.001, 0.005, 0.01, 0.05, 0.1, 0.5, 1)`. A
pair passes a threshold when its correlation p-value is below it in at
least one group. DISCO first compares the maximum real discovery
percentage across thresholds with the maximum from each permutation:

``` text
global p = (1 + number of null maxima >= real maximum) / (1 + valid permutations)
```

If the global test passes, DISCO selects the threshold maximizing
`max(real %, 1) / max(mean null %, 1)`; ties select the largest
threshold. Otherwise it uses threshold `1`. With 100 permutations, the
smallest possible global p-value is `1/101`. Thresholds without a real
discovery percentage are excluded from both maxima; permutations missing
any remaining threshold are excluded from the global test. Counts used
are reported in `global_test`. The function stops if no thresholds or
valid permutations remain.

For a fixed threshold, set `empirically_derive_threshold = FALSE` and
supply `filter_threshold`. This skips permutations and diagnostic plots.

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

In two-modality mode, DISCO tests **cross-modality pairs** and scores
features in the **first modality**. Samples may differ between
modalities; permutations preserve each modality’s group sizes and
overlap. A shared sample receives one shuffled label for both
modalities.

In `feature_null_summary`, `q_less_than_0.2` indicates a real-data DISCO
`q < 0.2`; its name follows `q_value`. Null discovery percentages use
all permutations as the denominator, counting untested features as
non-discoveries. These frequencies are not formal feature-level FDR
estimates. The plots display real p-values below `1e-8` at
`-log10(p) = 8`; original p-values are retained. The global test
accounts for searching across thresholds but does not adjust
feature-level q-values for threshold selection.

Fixed mode returns the first three components and `output_files`; the
diagnostic components are `NULL`. Pair and feature tables are saved as
`<prefix>_cor_diff.csv` and `<prefix>_disco.csv`.

## Yachida example: taxa

We use healthy controls and CRC samples from [Yachida et
al. (2019)](https://www.nature.com/articles/s41591-019-0458-7), accessed
through
[curatedMetagenomicData](https://waldronlab.io/curatedMetagenomicData/reference/returnSamples.html).
This example uses:

- Relative abundances (`counts = FALSE`).
- Features present in **≥10% of samples in each group**: a feature must
  pass separately in both healthy controls and CRC.
- Zero replacement with **half each feature’s minimum positive
  abundance**, estimated across that cohort after prevalence filtering.
- **CLR normalization** within each sample, separately for each
  modality.
- After CLR, each feature is residualized against **age, gender, and
  BMI** using **separate models within healthy controls and CRC**. Each
  group has its own intercept and covariate coefficients. Models are
  fitted once within the original groups; the resulting residuals feed
  the real comparison and all label permutations, without refitting the
  covariate models after shuffling.
- 100 label permutations, q cutoff 0.2, and global p-value cutoff 0.05.

Install the example dependencies if needed:

``` r
install.packages("BiocManager")
BiocManager::install(c("curatedMetagenomicData", "SummarizedExperiment"))
```

### Retrieve and preprocess the data

The blocks below are read directly from the [analysis
script](analysis/run_yachida_readme_examples.R) and its
[helpers](analysis/yachida_readme_preprocessing.R). The script
initializes the local package and repository paths (`analysis_dir`,
`repo_dir`), then runs these exact blocks. Cached runs are reused only
when inputs, settings, and source-code hashes match.

``` r
settings <- list(
  n_resamples = 100L, random_seed = 1L, cooccur = 10L, adj_method = "BY",
  q_value = 0.2, global_p_cutoff = 0.05, prevalence_cutoff = 0.1,
  prevalence_scope = "at least 10% in each of healthy and CRC",
  covariate_adjustment = "separate healthy and CRC models: ~ age + gender + BMI"
)
results_dir <- file.path(analysis_dir, "yachida_readme_groupwise_10pct_results")
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
source_md5 <- readme_source_hashes(repo_dir)
```

Samples missing age, gender, or BMI are excluded **before prevalence
filtering and normalization**. The excluded-sample table is saved with
the analysis.

``` r
metadata <- curatedMetagenomicData::sampleMetadata
metadata <- metadata[
  !is.na(metadata$study_name) & metadata$study_name == "YachidaS_2019" &
    metadata$disease %in% c("healthy", "CRC"), , drop = FALSE
]
if (!is.numeric(metadata$age) || !is.numeric(metadata$BMI)) {
  stop("Age and BMI must be numeric.")
}
complete_covariates <- stats::complete.cases(metadata[, c("age", "gender", "BMI")]) &
  is.finite(metadata$age) & is.finite(metadata$BMI) &
  nzchar(trimws(as.character(metadata$gender)))
excluded_covariates <- as.data.frame(metadata[!complete_covariates, , drop = FALSE])
metadata <- as.data.frame(metadata[complete_covariates, , drop = FALSE])
if (!nrow(metadata) || anyNA(metadata$sample_id) || anyDuplicated(metadata$sample_id)) {
  stop("Covariate-complete sample IDs must be non-missing and unique.")
}
rownames(metadata) <- NULL
sample_metadata <- metadata
rownames(sample_metadata) <- sample_metadata$sample_id
message("Excluded ", nrow(excluded_covariates),
        " sample(s) missing age, gender or BMI.")
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
if (any(table(factor(taxa_group, levels = c("healthy", "CRC"))) < settings$cooccur)) {
  stop("Each group must have at least cooccur usable samples.")
}
```

The same preprocessing helper is used for taxa and pathways:

``` r
# Samples are rows; features are columns. Apply independently per modality.
preprocess_abundance <- function(x, group, prevalence_cutoff = 0.10) {
  stopifnot(
    is.matrix(x), is.numeric(x), nrow(x) > 0L, ncol(x) > 0L,
    all(is.finite(x)), all(x >= 0), all(rowSums(x) > 0),
    !is.null(rownames(x)), !anyDuplicated(rownames(x)),
    !is.null(colnames(x)), !anyDuplicated(colnames(x)),
    length(prevalence_cutoff) == 1L, is.finite(prevalence_cutoff),
    prevalence_cutoff > 0, prevalence_cutoff <= 1,
    identical(rownames(x), names(group)), !anyNA(group),
    setequal(as.character(group), c("healthy", "CRC"))
  )

  # Require the prevalence cutoff separately in BOTH original groups.
  x <- sweep(x, 1, rowSums(x), "/")
  prevalence <- cbind(
    healthy = colMeans(x[group == "healthy", , drop = FALSE] > 0),
    CRC = colMeans(x[group == "CRC", , drop = FALSE] > 0)
  )
  keep <- apply(prevalence >= prevalence_cutoff, 1, all)
  if (sum(keep) < 2L) stop("Fewer than two features pass the prevalence filter.")
  retained <- x[, keep, drop = FALSE]

  # One half-minimum per feature, estimated across the combined cohort.
  half_min <- apply(retained, 2, function(z) min(z[z > 0]) / 2)
  imputed <- retained
  zero_positions <- which(imputed == 0, arr.ind = TRUE)
  if (nrow(zero_positions) > 0L) {
    imputed[zero_positions] <- half_min[zero_positions[, "col"]]
  }
  log_x <- log(imputed)
  clr <- sweep(log_x, 1, rowMeans(log_x), "-")

  list(
    clr = clr,
    retained_relative_abundance = retained,
    imputed_abundance = imputed,
    prevalence = data.frame(
      feature = rownames(prevalence), prevalence, retained = unname(keep),
      row.names = NULL
    ),
    half_minimum = data.frame(feature = names(half_min), replacement = unname(half_min))
  )
}
```

The covariate-adjustment helper fits a separate response model for each
feature within each group. Residuals keep the original sample and
feature names.

``` r
# Fit independent age/gender/BMI models within the original disease groups.
residualize_abundance <- function(clr, metadata, group) {
  stopifnot(identical(rownames(clr), rownames(metadata)), all(is.finite(clr)),
            is.numeric(metadata$age), is.numeric(metadata$BMI),
            all(is.finite(metadata$age)), all(is.finite(metadata$BMI)),
            identical(rownames(clr), names(group)), !anyNA(group),
            setequal(as.character(group), c("healthy", "CRC")))
  residuals <- matrix(NA_real_, nrow(clr), ncol(clr), dimnames = dimnames(clr))
  models <- list()
  for (group_name in c("healthy", "CRC")) {
    rows <- which(group == group_name)
    model_data <- metadata[rows, , drop = FALSE]
    model_data$gender <- factor(trimws(as.character(model_data$gender)))
    if (anyNA(model_data$gender) || any(!nzchar(as.character(model_data$gender))) ||
        nlevels(model_data$gender) < 2L) {
      stop("At least two non-missing gender levels are required in ", group_name, ".")
    }
    design <- stats::model.matrix(~ age + gender + BMI, data = model_data)
    fit <- stats::lm.fit(design, clr[rows, , drop = FALSE])
    if (fit$rank < ncol(design) || nrow(design) <= fit$rank) {
      stop("Covariate model requires full rank and positive residual degrees of freedom in ",
           group_name, ".")
    }
    residuals[rows, ] <- fit$residuals
    models[[group_name]] <- list(
      coefficients = fit$coefficients, design = design,
      sample_ids = rownames(model_data), residual_df = nrow(design) - fit$rank
    )
  }
  stopifnot(all(is.finite(residuals)))
  list(residuals = residuals, models = models)
}
```

``` r
taxa_prepared <- preprocess_abundance(taxa, taxa_group, settings$prevalence_cutoff)
taxa_adjusted <- residualize_abundance(
  taxa_prepared$clr, sample_metadata[rownames(taxa_prepared$clr), , drop = FALSE],
  taxa_group
)
single_dir <- file.path(results_dir, "taxa")
single_input <- list(
  data1 = list(taxa_adjusted$residuals[taxa_group == "healthy", , drop = FALSE]),
  data2 = list(taxa_adjusted$residuals[taxa_group == "CRC", , drop = FALSE]),
  settings = settings, source_md5 = source_md5
)
single <- readme_cached_result(single_dir, single_input, recompute)
if (is.null(single)) {
  single <- DISCO::disco(
    data1 = single_input$data1, data2 = single_input$data2,
    n_resamples = settings$n_resamples, random_seed = settings$random_seed,
    q_value = settings$q_value, global_p_cutoff = settings$global_p_cutoff,
    cooccur = settings$cooccur, adj_method = settings$adj_method,
    output_dir = single_dir,
    prefix = "YachidaS_2019_taxa"
  )
  saveRDS(single_input, file.path(single_dir, "disco_inputs.rds"))
  saveRDS(single, file.path(single_dir, "disco_analysis.rds"))
}
```

### Results

These outputs come from the completed analysis, saved before rendering
this README. Package version: 0.0.0.9000; curatedMetagenomicData
version: 3.12.0. Downloading a different resource version can change
results. The complete-case restriction excluded **1** sample(s) with
missing covariates.

| control | CRC | retained_taxa |
|--------:|----:|--------------:|
|     245 | 258 |           184 |

Global p-value: **0.009901**. Selected threshold: **0.01**. Taxa with
DISCO q \<0.2 at the selected threshold: **80**.

Top 10 taxa by DISCO q-value, breaking ties by p-value:

| Taxon                      | DISCO p  | DISCO q  | Comparisons |
|:---------------------------|:---------|:---------|------------:|
| Parvimonas micra           | 3.43e-07 | 3.31e-04 |         121 |
| Fusobacterium ulcerans     | 1.08e-06 | 3.31e-04 |         133 |
| Actinomyces sp. S6-Spd3    | 1.08e-06 | 3.31e-04 |         134 |
| Actinomyces sp. HMSC035G02 | 1.24e-06 | 3.31e-04 |         128 |
| Clostridium disporicum     | 5.48e-06 | 1.17e-03 |         110 |
| Lactobacillus rogosae      | 8.43e-06 | 1.36e-03 |          95 |
| Actinomyces sp. ICM47      | 8.93e-06 | 1.36e-03 |         119 |
| Acidaminococcus intestini  | 1.13e-05 | 1.50e-03 |          87 |
| \[Eubacterium\] sulci      | 1.26e-05 | 1.50e-03 |          98 |
| Clostridium sp. CAG:299    | 1.60e-05 | 1.71e-03 |         116 |

The first plot compares real discovery percentages with the null mean ±
SD. The second shows the ratio used to select the threshold.

<img src="man/figures/README-yachida-taxa-discovery_plot.png" width="100%" />

<img src="man/figures/README-yachida-taxa-ratio_plot.png" width="100%" />

Each point below is a taxon. Orange denotes real-data DISCO q \<0.2;
each panel reports its orange-taxon count. The y-axis is the percentage
of null permutations in which that taxon had DISCO q \<0.2.

<img src="man/figures/README-yachida-taxa-feature_null_plot.png" width="100%" />

## Yachida example: taxa and pathways

We now compare **taxon–pathway correlations** between CRC and healthy
controls. The first modality is taxa, so the DISCO scores describe taxa
whose relationships with pathways differ between groups. Pathways are
unstratified community totals; `UNMAPPED` and `UNINTEGRATED` rows are
removed. Both modalities are restricted to common usable samples, then
prevalence-filtered and normalized separately using the same ≥10%,
half-minimum, CLR workflow above. **Both taxa and pathways** are then
residualized against age, gender, and BMI separately within each group,
and DISCO compares correlations between those residuals.

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
if (any(table(factor(two_group, levels = c("healthy", "CRC"))) < settings$cooccur)) {
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
  data2 = list(two_taxa_adjusted$residuals[two_group == "CRC", , drop = FALSE],
               two_pathways_adjusted$residuals[two_group == "CRC", , drop = FALSE]),
  settings = settings, source_md5 = source_md5
)
two <- readme_cached_result(two_dir, two_input, recompute)
if (is.null(two)) {
  message("Running two modalities: ", ncol(two_taxa$clr), " taxa, ",
          ncol(two_pathways$clr), " pathways; control n = ",
          sum(two_group == "healthy"), "; CRC n = ", sum(two_group == "CRC"), ".")
  two <- DISCO::disco(
    data1 = two_input$data1, data2 = two_input$data2,
    n_resamples = settings$n_resamples, random_seed = settings$random_seed,
    q_value = settings$q_value, global_p_cutoff = settings$global_p_cutoff,
    cooccur = settings$cooccur, adj_method = settings$adj_method,
    output_dir = two_dir, prefix = "YachidaS_2019_taxa_pathways"
  )
  saveRDS(two_input, file.path(two_dir, "disco_inputs.rds"))
  saveRDS(two, file.path(two_dir, "disco_analysis.rds"))
}
```

### Results

| control | CRC | retained_taxa | retained_pathways |
|--------:|----:|--------------:|------------------:|
|     245 | 258 |           184 |               404 |

Global p-value: **1**. Selected threshold: **1**. Taxa with DISCO q
\<0.2 at the selected threshold: **0**.

The global test does not detect excess discoveries in this taxon–pathway
analysis, so DISCO falls back to threshold 1. The table below shows the
top 10 taxa by DISCO q-value, breaking ties by p-value. A ranking alone
does not establish significance; use the reported q-values.

| Taxon                          | DISCO p  | DISCO q | Comparisons |
|:-------------------------------|:---------|:--------|------------:|
| Parvimonas micra               | 1.05e-02 | 1e+00   |         404 |
| Bifidobacterium longum         | 1.32e-02 | 1e+00   |         404 |
| Actinomyces sp. HMSC035G02     | 1.86e-02 | 1e+00   |         404 |
| Odoribacter splanchnicus       | 3.45e-02 | 1e+00   |         404 |
| Erysipelatoclostridium ramosum | 3.56e-02 | 1e+00   |         404 |
| Mogibacterium diversum         | 4.18e-02 | 1e+00   |         404 |
| Enterocloster clostridioformis | 4.71e-02 | 1e+00   |         404 |
| Streptococcus sp. A12          | 4.86e-02 | 1e+00   |         404 |
| Adlercreutzia equolifaciens    | 5.21e-02 | 1e+00   |         404 |
| Enterocloster aldenensis       | 5.93e-02 | 1e+00   |         404 |

<img src="man/figures/README-yachida-taxa_pathways-discovery_plot.png" width="100%" />

<img src="man/figures/README-yachida-taxa_pathways-ratio_plot.png" width="100%" />

<img src="man/figures/README-yachida-taxa_pathways-feature_null_plot.png" width="100%" />

## Rebuilding these examples

`README.Rmd` is the source for this README. Rendering reads the portable
saved summaries in `README-results.rds` and figures under
`man/figures/`; it does not download data or rerun permutations.
Analysis code blocks are imported from the script with
`knitr::read_chunk()`, so they are never maintained separately in the
Rmd. Rendering stops if source-code hashes differ from the saved
results, preventing changed code from being shown with stale outputs.
From the repository root:

``` r
install.packages(c("rmarkdown", "knitr")) # If needed; Pandoc is also required
rmarkdown::render("README.Rmd")
```

To regenerate the analysis outputs and README assets with the current
local package, run:

``` sh
Rscript analysis/run_yachida_readme_examples.R
```

This requires the example dependencies and `pkgload`, may download
resources through ExperimentHub, and reuses existing results only when
inputs and settings and source hashes match. Set
`DISCO_README_RECOMPUTE=true` to recompute both examples. Results stay
under `analysis/yachida_readme_groupwise_10pct_results/`; previous
pooled-adjustment and unadjusted results are retained. The README bundle
contains only cohort summaries and top-result tables. Package checks can
be run with `devtools::test()` and `devtools::check()`. Preprocessing
and residualization helpers can be checked without downloading data
using `Rscript analysis/test_yachida_readme_helpers.R`.

## License

DISCO is licensed under the MIT License.
