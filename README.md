
# DISCO

<!-- README.md is generated from README.Rmd. Please edit README.Rmd. -->

DISCO (**DIS**ruption in **CO**rrelations) identifies features whose
correlation patterns differ between two groups, such as healthy controls
and people with a disease. It combines evidence from pairwise
differential-correlation tests into a score, p-value, and adjusted
p-value (q-value) for each feature. This helps identify features whose
relationships with other measurements change between groups, even when
their average abundance is not the focus of the analysis.

DISCO supports one modality (for example, correlations among microbial
taxa) or two modalities (for example, correlations between microbes and
metabolites). It can be applied to quantitative features such as
microbial abundances, gene expression, protein levels, and metabolites,
with preprocessing appropriate to the data type. Correlations describe
statistical associations, not necessarily direct biological
interactions.

The implementation uses Spearman correlations with Fisher
transformations and [Bonett–Wright standard errors
(2000)](https://doi.org/10.1007/BF02294183), then combines
differential-correlation p-values using the [Cauchy combination test
(Liu and Xie, 2020)](https://doi.org/10.1080/01621459.2018.1554485).
DISCO was introduced by [Fouladi et
al. (2025)](https://doi.org/10.1038/s41467-025-64822-z).

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
group 1; you do not need to align columns manually. Supply unique,
nonempty feature names.

For two modalities, row names must identify samples consistently. Within
each group, the modalities must share samples, but their sample sets and
row order may differ. DISCO aligns rows by sample name and represents
measurements from an unavailable modality as missing. Each feature pair
must have at least `cooccur` paired, nonmissing observations **in each
group** to be tested.

**Normalize inputs appropriately before running DISCO.** The function
does not normalize data or adjust for covariates. Supply normalized and,
when needed, covariate-adjusted measurements. The abundance examples
below use a centered log-ratio (CLR) transformation; other data types
may require different preprocessing. Zeros are treated as observed
values, while missing values are excluded separately for each pair.

| Argument | Default | Meaning |
|----|---:|----|
| `cooccur` | `10` | Minimum pairwise non-missing observations in **each** group; integer ≥10. |
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
number of uninformative pairs contributing to a score.

`filter_threshold` screens pairs using their **single-group correlation
p-values**, not their differential-correlation p-values. A pair is
retained if its correlation p-value is below the threshold in **at least
one group**. For example, `filter_threshold = 0.05` excludes pairs whose
correlation p-values are both ≥0.05. A threshold of `1` is the least
restrictive setting; minimum sample counts and valid correlation tests
are still required. The implementation uses a strict `<` comparison, so
a pair with p-values exactly `1` in both groups is also excluded.

By default, `empirically_derive_threshold = TRUE`. DISCO compares the
observed data with size-matched random label permutations at thresholds
`c(0.001, 0.005, 0.01, 0.05, 0.1, 0.5, 1)`. At each threshold, it
calculates the percentage of tested features with DISCO `q < q_value` in
the observed data and in each permutation.

A global permutation test compares the maximum observed discovery
percentage across thresholds with the corresponding maximum from each
valid permutation:

```text
global p = (1 + number of null maxima >= real maximum) / (1 + valid permutations)
```

If `global p < global_p_cutoff`, DISCO selects the threshold that
maximizes `max(real percentage, 1) / max(mean null percentage, 1)`. The
floors are in percentage points and prevent ratios from being dominated
by near-zero values. Ties select the largest threshold. Otherwise, DISCO
uses threshold `1`. This global test accounts for searching across
thresholds; it does not further adjust the individual feature q-values
for threshold selection.

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
the supplied `q_value`. Null discovery percentages use all requested
permutations as the denominator, counting absent or untested features as
non-discoveries.

`feature_null_plot` shows **−log10(real-data p-value)** on the x-axis
and the percentage of permutations in which the feature has
`q < q_value` on the y-axis. Orange points identify real-data
discoveries. Low null discovery frequencies provide a useful diagnostic,
but are not formal feature-level FDR estimates.

In fixed mode, `cor_diff`, `disco`, `selected_threshold`, and
`output_files` are returned; the diagnostic components are `NULL`. Pair
and feature tables are saved as `<prefix>_cor_diff.csv` and
`<prefix>_disco.csv`.

## Example: One modality — microbial composition

For this example, we use healthy controls and CRC metagenomes from
[Yachida et
al. (2019)](https://www.nature.com/articles/s41591-019-0458-7), accessed
through
[curatedMetagenomicData](https://waldronlab.io/curatedMetagenomicData/reference/returnSamples.html).

The workflow preprocesses the data as follows:

1.  Retrieve microbial relative abundances from `curatedMetagenomicData`
    and exclude samples missing age, gender, or BMI.
2.  Retain features present in **≥10% of samples in each group**,
    requiring the cutoff separately in controls and CRC.
3.  Replace zeros with half of each feature’s minimum positive relative
    abundance, estimated across the combined retained cohort.
4.  Apply the **CLR transformation** within each sample.
5.  Regress each transformed feature on **age, gender, and BMI**,
    fitting separate models in controls and CRC. Use the residuals as
    DISCO inputs. These models are fitted once in the original groups,
    before permutations.

The example uses 100 label permutations, a feature q-value cutoff of
0.2, and a global p-value cutoff of 0.05.

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

The preprocessing helper applies the prevalence filter, half-minimum
zero replacement, and CLR transformation. It is used independently for
each modality.

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

Run DISCO on the covariate-adjusted measurements:

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
version: 3.14.0. Downloading a different resource version can change
results. The complete-case restriction excluded **1** sample(s) with
missing covariates.

| control | CRC | retained_taxa |
|--------:|----:|--------------:|
|     245 | 258 |           184 |

Global p-value: **0.009901**. Selected threshold: **0.01**. Taxa with
DISCO q \< 0.2 at the selected threshold: **80**.

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

The two highest-ranked taxa are *Parvimonas micra* and *Fusobacterium
ulcerans*. *P. micra* has been associated with CRC and shown to promote
colorectal tumorigenesis in experimental models ([Zhao et al.,
2022](https://pmc.ncbi.nlm.nih.gov/articles/PMC9439953/)). For *F.
ulcerans*, the evidence is more specific: a lineage including *F.
ulcerans* and *F. varium* was associated with lymphovascular invasion in
CRC tumours, although that lineage was not enriched in CRC overall ([Bi
et al., 2022](https://www.nature.com/articles/s41467-022-30957-6)).
These studies provide biological context; DISCO identifies changes in
correlation patterns rather than testing abundance differences or
causality.

The discovery plot compares observed discovery percentages with the
permutation mean ± SD. Here, the global p-value is **0.009901**, and the
selected threshold is **0.01**. The ratio plot shows the criterion used
to select that threshold.

<img src="man/figures/README-yachida-taxa-discovery_plot.png" alt="" width="100%" />

<img src="man/figures/README-yachida-taxa-ratio_plot.png" alt="" width="100%" />

In the next figure, each point is a taxon and each panel is a filtering
threshold. Orange points have real-data DISCO q \< 0.2. The x-axis shows
−log10 of the real-data p-value; the y-axis shows the percentage of
label permutations in which the taxon has DISCO q \< 0.2. Compare panels
to see how filtering changes both the observed discoveries and their
null discovery frequencies.

<img src="man/figures/README-yachida-taxa-feature_null_plot.png" alt="" width="100%" />

## Example: Two modalities — microbial composition and pathways

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

Global p-value: **1**. Selected threshold: **1**. Taxa with DISCO q \<
0.2 at the selected threshold: **0**.

The global test does not detect an excess of discoveries in this
taxon–pathway analysis, so DISCO uses the least restrictive threshold,
1.

The table below lists the top 10 taxa by DISCO q-value, breaking ties by
p-value. At the selected threshold, **0** taxa have q \< 0.2.

In this cohort and with this preprocessing, the one-modality analysis
detects changes in taxon–taxon correlation patterns. The taxon–pathway
analysis does not provide comparable evidence of a global excess of
discoveries. This result does not establish that microbial functions or
individual taxon–pathway relationships are unchanged.

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

<img src="man/figures/README-yachida-taxa_pathways-discovery_plot.png" alt="" width="100%" />

<img src="man/figures/README-yachida-taxa_pathways-ratio_plot.png" alt="" width="100%" />

<img src="man/figures/README-yachida-taxa_pathways-feature_null_plot.png" alt="" width="100%" />

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

## Citation

If you use DISCO, please cite the paper that introduced the method:

Fouladi F, Chen Y, Bera S, Jarmusch AK, Van Tyne D, Palella FJ,
Margolick JB, Chew KW, Sun J, Martinson J, Rinaldo CR, Peddada SD. A
taxon-specific measurement of disruption in a multi-modal study of
microbiomes and metabolomes reveals system-wide dysbiosis preceding
HIV-1 infection. Nat Commun. 2025 Nov 20;16(1):10204.
[doi:10.1038/s41467-025-64822-z](https://doi.org/10.1038/s41467-025-64822-z).
PMID: 41266305; PMCID: PMC12635353.
