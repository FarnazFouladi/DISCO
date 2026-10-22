#' Calculate DISruption in COrrelations (DISCO)
#'
#' Quantifies feature-specific changes in Spearman correlation patterns between
#' two groups. Inputs may contain one data modality or two data modalities.
#' Corresponding modalities are aligned to their common features before the
#' correlations are compared. The correlation filtering threshold can either
#' be selected empirically from a fixed candidate set or supplied directly by
#' the user.
#'
#' Input values should already be normalized and, when desired, adjusted for
#' covariates. The function does not normalize or residualize the data.
#' Fisher-transformed Spearman correlations use the Bonett-Wright standard
#' error. Cauchy-transformed tail probabilities are bounded below at `1e-8`
#' for numerical stability.
#'
#' @param data1 A list containing one or two numeric matrices or data frames for
#'   the first group. Rows represent samples and columns represent features.
#' @param data2 A list with the same number of modalities as `data1`, containing
#'   data for the second group. Rows represent samples and columns represent
#'   features.
#' @param cooccur A single integer greater than or equal to 4. A feature pair is
#'   tested only when both groups have at least this many pairwise non-missing
#'   observations.
#' @param empirically_derive_threshold Logical. If `TRUE`, select the filtering
#'   threshold empirically by comparing observed discoveries with discoveries
#'   from size-matched random label permutations. If `FALSE`, use
#'   `filter_threshold` and do not generate null data or diagnostic plots.
#' @param filter_threshold Single-group correlation p-value threshold used when
#'   `empirically_derive_threshold = FALSE`. A pair is retained when its
#'   correlation is significant in at least one group at this threshold.
#' @param q_value Adjusted DISCO p-value cutoff used to declare a feature
#'   significant when estimating empirical real and null discovery
#'   percentages. Used only when `empirically_derive_threshold = TRUE`.
#' @param n_resamples Number of random, size-matched label permutations used to
#'   estimate the mean null discovery percentage. Used only when
#'   `empirically_derive_threshold = TRUE`.
#' @param random_seed Non-negative integer seed for reproducible label
#'   permutations. The caller's random-number state is restored on exit. Used
#'   only when `empirically_derive_threshold = TRUE`.
#' @param adj_method Multiple-testing adjustment method passed to
#'   [stats::p.adjust()].
#' @param output_dir Directory in which result files are written. Empirical mode
#'   writes CSV and PNG files; fixed mode writes only the result CSV files. The
#'   directory is created when it does not exist.
#' @param prefix Character prefix used for result filenames.
#' @param global_p_cutoff Significance cutoff for the global maximum-discovery
#'   permutation test in empirical mode. Defaults to `0.05`. Ratio selection
#'   proceeds only when the global p-value is strictly below this cutoff;
#'   otherwise the filtering threshold is `1`.
#' @param workers Positive integer number of parallel worker processes for null
#'   permutations. Defaults to `1` (sequential); used only in empirical mode.
#'   The number used is capped at `n_resamples`. Each worker holds its own
#'   correlation matrices, so start with `2` if memory is limited. Uses socket
#'   workers on all platforms. For a given seed, permutations and results do
#'   not depend on worker count. Workers are stopped on completion or error.
#'
#' @return A list containing `cor_diff`, `disco`, `selected_threshold`,
#'   `threshold_summary`, `global_test`, `discovery_plot`, `ratio_plot`,
#'   `feature_null_summary`, `feature_null_plot`, and `output_files`.
#'   The `cor_diff` table is sorted by ascending adjusted p-value (`q`).
#'   `feature_null_summary` contains real p-values and q-values and per-feature
#'   null testing and discovery counts at each candidate threshold.
#'   Its real significance column is named `q_less_than_<q_value>`, for example
#'   `q_less_than_0.2` when `q_value = 0.2`.
#'   `feature_null_plot` shows negative log10 real p-values versus the percentage
#'   of permutations declaring each feature significant, colored using
#'   `q_value`, with significant-feature counts in each panel.
#'   In empirical mode, `global_test` reports the global p-value,
#'   cutoff, signal decision, maximum discovery summaries, number of evaluable
#'   thresholds and valid permutations, and selected threshold. When the global
#'   test passes, selection maximizes `max(real percentage, 1) /
#'   max(mean null percentage, 1)`; ties select the largest threshold.
#'   Otherwise the threshold is `1`. In fixed mode, `threshold_summary`,
#'   `global_test`, `discovery_plot`, `ratio_plot`, `feature_null_summary`, and
#'   `feature_null_plot` are `NULL`.
#'
#' @details In empirical mode, samples from the two groups are pooled for each
#'   null resample, randomly permuted without replacement, and split into groups
#'   with the same sample sizes as `data1` and `data2`. For two-modality
#'   analyses, permutation is stratified by modality-availability pattern. This
#'   applies one coherent shuffled label per sample while preserving each
#'   modality's group sizes and the number of samples shared between modalities.
#'   The global statistic is the maximum discovery percentage across thresholds.
#'   Each permutation contributes its maximum across the same thresholds, with
#'   p-value `(1 + sum(null maxima >= real maximum)) / (1 + valid permutations)`.
#'   Thresholds without a finite real discovery percentage are excluded from
#'   both maxima; permutations missing any remaining threshold are excluded from
#'   the global test. An error is raised if no thresholds or permutations remain.
#'   The test uses unfloored percentages. It accounts for the search across
#'   thresholds but does not adjust individual feature-level DISCO q-values for
#'   threshold selection. With 100 permutations the smallest possible global
#'   p-value is `1/101`.
#'   Feature null discovery frequencies use all `n_resamples` permutations as
#'   their denominator; features absent or untested in a permutation count as
#'   non-discoveries. These frequencies are not formal feature-level FDR
#'   estimates. In two-modality mode the diagnostic covers modality-1 features.
#'   For plotting only, real p-values below `1e-8` are displayed at `1e-8`.
#' @export
disco <- function(data1,
                  data2,
                  cooccur = 10,
                  empirically_derive_threshold = TRUE,
                  filter_threshold = 0.05,
                  q_value = 0.2,
                  n_resamples = 100,
                  random_seed = 1L,
                  adj_method = "BY",
                  output_dir,
                  prefix,
                  global_p_cutoff = 0.05,
                  workers = 1L) {
  if (!is.logical(empirically_derive_threshold) ||
      length(empirically_derive_threshold) != 1L ||
      is.na(empirically_derive_threshold)) {
    stop("empirically_derive_threshold must be TRUE or FALSE", call. = FALSE)
  }

  if (empirically_derive_threshold) {
    return(run_empirical_disco(
      data1 = data1,
      data2 = data2,
      cooccur = cooccur,
      filter_thresholds = c(0.001, 0.005, 0.01, 0.05, 0.1, 0.5, 1),
      q_value = q_value,
      n_resamples = n_resamples,
      random_seed = random_seed,
      adj_method = adj_method,
      output_dir = output_dir,
      prefix = prefix,
      global_p_cutoff = global_p_cutoff,
      workers = workers
    ))
  }

  result <- calculate_fixed_threshold_disco(
    data1 = data1,
    data2 = data2,
    cooccur = cooccur,
    filt_threshold = filter_threshold,
    adj_method = adj_method,
    output_dir = output_dir,
    prefix = prefix
  )
  output_files <- c(
    cor_diff = file.path(output_dir, paste0(prefix, "_cor_diff.csv")),
    disco = file.path(output_dir, paste0(prefix, "_disco.csv"))
  )

  list(
    cor_diff = result$cor_diff,
    disco = result$disco,
    selected_threshold = filter_threshold,
    threshold_summary = NULL,
    global_test = NULL,
    discovery_plot = NULL,
    ratio_plot = NULL,
    feature_null_summary = NULL,
    feature_null_plot = NULL,
    output_files = output_files
  )
}
