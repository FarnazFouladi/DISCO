# DISCO

DISCO (**DIS**ruption in **CO**rrelations) quantifies feature-specific changes
in correlation patterns between two groups. It can compare relationships among
features from one data modality or relationships between features from two
different modalities.

DISCO uses Spearman correlations, applies strict pairwise-complete sample-size
filtering in both groups, tests correlation differences with Fisher
transformations, and combines pairwise p-values into feature-level scores using
a Cauchy transformation.

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
- In two-modality analyses, sample sets may differ between modalities. DISCO
  reports the number of common samples and stops if there is no overlap.
- Sample identifiers do not need to match between the two groups.
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
  output_dir = file.path(tempdir(), "disco-one-modality"),
  prefix = "example"
)

result$cor_diff
result$cor_diff_sig
result$disco
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
  output_dir = file.path(tempdir(), "disco-two-modality"),
  prefix = "example"
)

result$cor_diff
result$cor_diff_sig
result$disco
```

## Important parameters

| Argument | Default | Purpose |
|---|---:|---|
| `cooccur` | `10` | Minimum pairwise-complete count required in each group. |
| `filt_threshold` | `0.001` | Retains a pair when its correlation meets this p-value threshold in at least one group. |
| `diff_threshold` | `0.3` | Minimum absolute difference between raw group correlations for a pair to be called significant. |
| `p_alpha` | `0.01` | Two-sample correlation-test p-value threshold. |
| `q_alpha` | `0.2` | Adjusted two-sample correlation-test p-value threshold. |
| `adj_method` | `"BY"` | Adjustment method passed to `stats::p.adjust()`. |

DISCO uses the Bonett–Wright standard error for Fisher-transformed Spearman
correlations. For numerical stability, the Cauchy calculation uses the nearest
p-value tail and fixes its minimum value at `1e-8`.

## Results

`disco()` returns a list containing three data frames:

| Component | Contents |
|---|---|
| `cor_diff` | Retained feature pairs, group correlations (`r1`, `r2`), pairwise-complete counts (`n1`, `n2`), two-sample p-values (`p`), and adjusted p-values (`q`). |
| `cor_diff_sig` | Pairs from `cor_diff` meeting `p_alpha`, `q_alpha`, and `diff_threshold`, with significant-correlation counts added. |
| `disco` | Feature-level Cauchy statistics, p-values, adjusted p-values, and the number of pairwise comparisons included. |

The same tables are written to `output_dir`:

- `<prefix>_cor_diff.csv`
- `<prefix>_cor_diff_sig.csv`
- `<prefix>_disco.csv`

## Development checks

```r
devtools::test()
devtools::check()
```

## License

DISCO is licensed under the MIT License.
