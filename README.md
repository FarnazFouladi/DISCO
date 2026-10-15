# DISCO

DISCO (**DIS**ruption in **CO**rrelations) quantifies feature-specific changes
in correlation patterns between two groups. It can compare relationships among
features from one data modality or relationships between features from two
different modalities.

DISCO uses Spearman correlations, applies strict pairwise-complete sample-size
filtering in both groups, tests correlation differences with Fisher
transformations, and combines pairwise p-values into feature-level scores using
a Cauchy transformation. It can select the single-group correlation filtering
threshold empirically or use a fixed threshold supplied by the user.

## Installation

Install the package from the repository root with either:

```sh
R CMD INSTALL .
```

or:

```r
devtools::install()
```

Then load it with:

```r
library(DISCO)
```

## Input format

`disco()` requires `data1` and `data2`, representing the two groups being
compared. Each argument must be a list containing either one or two numeric
matrices or data frames.

- Rows are samples and columns are features.
- `data1` and `data2` must contain the same number of modalities.
- Corresponding modalities are restricted to their common feature names and
  ordered according to `data1`.
- Feature names must be non-empty and unique within each matrix.
- In two-modality analyses, sample sets may differ between modalities. Sample
  names must be unique within each modality. Label permutation is stratified
  by modality availability so each modality's group size and the overlap
  between modalities remain fixed.
- Sample identifiers do not need to match between the two groups.
- Inputs should already be normalized and, if desired, adjusted for covariates.
  `disco()` does not normalize or residualize the input.
- Zeros are observed values. Only `NA` values reduce the pairwise-complete
  sample count.
- A feature pair is tested only when its pairwise-complete count is at least
  `cooccur` in both groups. `cooccur` must be an integer of at least 4.

In one-modality mode, DISCO evaluates within-modality feature pairs. In
two-modality mode, it evaluates cross-modality pairs and calculates
feature-level scores for features in the first modality.

## One-modality example

```r
library(DISCO)

set.seed(2026)
n <- 40

group1_signal <- rnorm(n)
group2_signal <- rnorm(n)

group1 <- data.frame(
  feature_a = group1_signal,
  feature_b = group1_signal + rnorm(n, sd = 0.2),
  feature_c = rnorm(n)
)

group2 <- data.frame(
  feature_a = group2_signal,
  feature_b = -group2_signal + rnorm(n, sd = 0.2),
  feature_c = rnorm(n)
)

result <- disco(
  data1 = list(group1),
  data2 = list(group2),
  n_resamples = 10, # Use the default 100 for an analysis
  output_dir = file.path(tempdir(), "disco-one-modality"),
  prefix = "example"
)

result$cor_diff
result$disco
result$selected_threshold
result$discovery_plot
result$ratio_plot
result$feature_null_plot
result$feature_null_summary
```

To use a fixed threshold and skip null permutations:

```r
result <- disco(
  data1 = list(group1),
  data2 = list(group2),
  empirically_derive_threshold = FALSE,
  filter_threshold = 0.05,
  output_dir = file.path(tempdir(), "disco-fixed-threshold"),
  prefix = "example"
)
```

## Two-modality example

Sample identifiers are used to match modalities within each group. The two
groups can have different identifiers.

```r
library(DISCO)

set.seed(2026)
n <- 40

group1_signal <- rnorm(n)
group2_signal <- rnorm(n)

group1_modality1 <- data.frame(
  taxon_a = group1_signal + rnorm(n, sd = 0.2),
  taxon_b = rnorm(n),
  row.names = paste0("group1_sample_", seq_len(n))
)
group1_modality2 <- data.frame(
  metabolite_a = group1_signal + rnorm(n, sd = 0.2),
  metabolite_b = rnorm(n),
  row.names = paste0("group1_sample_", seq_len(n))
)

group2_modality1 <- data.frame(
  taxon_a = group2_signal + rnorm(n, sd = 0.2),
  taxon_b = rnorm(n),
  row.names = paste0("group2_sample_", seq_len(n))
)
group2_modality2 <- data.frame(
  metabolite_a = -group2_signal + rnorm(n, sd = 0.2),
  metabolite_b = rnorm(n),
  row.names = paste0("group2_sample_", seq_len(n))
)

result <- disco(
  data1 = list(group1_modality1, group1_modality2),
  data2 = list(group2_modality1, group2_modality2),
  n_resamples = 10, # Use the default 100 for an analysis
  output_dir = file.path(tempdir(), "disco-two-modality"),
  prefix = "example"
)

result$cor_diff
result$disco
```

## Important parameters

| Argument | Default | Purpose |
|---|---:|---|
| `cooccur` | `10` | Minimum pairwise-complete count required in each group. |
| `empirically_derive_threshold` | `TRUE` | Select a threshold empirically; when `FALSE`, skip null permutations and plots. |
| `filter_threshold` | `0.05` | Fixed threshold used only when empirical selection is disabled. |
| `q_value` | `0.2` | Feature-level DISCO adjusted p-value cutoff used during empirical selection. |
| `global_p_cutoff` | `0.05` | Global permutation p-value must be strictly below this cutoff to select by ratio; otherwise use threshold `1`. |
| `n_resamples` | `100` | Number of size-matched random label permutations used during empirical selection. |
| `random_seed` | `1L` | Seed used for reproducible label permutations. |
| `adj_method` | `"BY"` | Adjustment method passed to `stats::p.adjust()`. |

Empirical mode evaluates the internal candidate set:

```r
c(0.001, 0.005, 0.01, 0.05, 0.1, 0.5, 1)
```

DISCO uses the Bonett–Wright standard error for Fisher-transformed Spearman
correlations. For numerical stability, the Cauchy calculation uses the nearest
p-value tail and fixes its minimum value at `1e-8`.

## Results

In empirical mode, `disco()` first compares the maximum real discovery
percentage across thresholds with the maximum from each label permutation.
The global p-value is `(1 + sum(null maxima >= real maximum)) / (B + 1)`,
where `B` is the number of valid permutations. If this p-value is below
`global_p_cutoff`, the threshold maximizing
`max(real percentage, 1) / max(mean null percentage, 1)` is selected.
If thresholds tie, the largest threshold is selected. If the global test does
not pass, the threshold is `1`.

Maxima use unfloored percentages and preserve each permutation across
thresholds. Thresholds with no finite real percentage are excluded from both
maxima, and permutations missing any remaining threshold are excluded from the
global test. The returned test summary reports the numbers used. The function
stops if there are no evaluable thresholds or complete permutations. With 100
permutations, the smallest possible global p-value is `1/101`. This global test
accounts for searching across thresholds; it does not adjust the feature-level
q-values for selecting a threshold. The function returns:

| Component | Contents |
|---|---|
| `cor_diff` | Retained feature pairs, group correlations (`r1`, `r2`), pairwise-complete counts (`n1`, `n2`), two-sample p-values (`p`), and adjusted p-values (`q`). |
| `disco` | Feature-level Cauchy statistics, p-values, adjusted p-values, and the number of pairwise comparisons included. |
| `selected_threshold` | Empirically selected correlation filtering threshold. |
| `threshold_summary` | Real and null discovery percentages, null standard deviations, both floored percentages, ratios, and the selected-threshold indicator. |
| `global_test` | Global p-value, cutoff, signal decision, maximum discovery summaries, valid permutation count, and selected threshold. |
| `discovery_plot` | Real and mean null discovery percentages across candidate thresholds. |
| `ratio_plot` | Real-to-null discovery ratio across candidate thresholds. |
| `feature_null_summary` | Real feature p-values/q-values, null testing/discovery counts, and percentage of all permutations declaring each feature significant, per threshold. |
| `feature_null_plot` | Panels of `-log10(real DISCO p)` versus feature null discovery percentage, colored by real `q < q_value`, with orange-feature counts. |
| `output_files` | Paths to the written CSV and PNG files. |

In fixed mode, `selected_threshold` is the supplied `filter_threshold`, and
`threshold_summary`, `global_test`, `discovery_plot`, `ratio_plot`,
`feature_null_summary`, and `feature_null_plot` are `NULL`.

The feature diagnostic reuses the same `n_resamples` label permutations. Its
denominator is all permutations, including those where a feature had no DISCO
result; these count as non-discoveries. The frequency is not a formal
feature-level FDR estimate. In two-modality mode it describes modality-1
features. Both colors and null discoveries use `q_value` (default `0.2`). Real
p-values below `1e-8` are plotted at 8 on the negative-log10 axis, while their
original values remain in `feature_null_summary`.

The following files are written to `output_dir`:

- `<prefix>_cor_diff.csv`
- `<prefix>_disco.csv`
- `<prefix>_threshold_summary.csv`
- `<prefix>_global_test.csv`
- `<prefix>_real_vs_null_thresholds.png`
- `<prefix>_real_to_null_ratio.png`
- `<prefix>_real_p_vs_null_discovery_by_threshold.csv`
- `<prefix>_real_p_vs_null_discovery_by_threshold.png`

The final six diagnostic files are created only in empirical mode.

## Development checks

```r
devtools::test()
devtools::check()
```

## License

DISCO is licensed under the MIT License.
