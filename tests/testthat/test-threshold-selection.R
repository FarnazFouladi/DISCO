testthat::test_that("threshold selection validates its public settings", {
  testthat::expect_error(
    validate_threshold_selection(0, 0.2, 100, 1),
    "filter_thresholds"
  )
  testthat::expect_error(
    validate_threshold_selection(0.1, 0, 100, 1),
    "q_value"
  )
  testthat::expect_error(
    validate_threshold_selection(0.1, 0.2, 1.5, 1),
    "n_resamples"
  )
  testthat::expect_error(
    validate_threshold_selection(0.1, 0.2, 10, -1),
    "random_seed"
  )
  for (cutoff in list(0, -0.1, NA_real_, Inf, 1.1, c(0.01, 0.05))) {
    testthat::expect_error(
      validate_threshold_selection(0.1, 0.2, 100, 1, cutoff),
      "global_p_cutoff"
    )
  }
})

testthat::test_that("largest threshold is selected when maximum ratios tie", {
  summary <- data.frame(
    filter_threshold = c(0.001, 0.01, 0.1),
    real_to_null_ratio = c(1, 3, 3)
  )

  testthat::expect_identical(select_optimal_threshold(summary), 0.1)
})

testthat::test_that("global test preserves permutation maxima and includes ties", {
  summary <- data.frame(
    filter_threshold = c(0.001, 0.05, 1),
    real_percentage = c(20, 50, 30),
    real_to_null_ratio = c(2, 5, 3)
  )
  null <- rbind(c(10, 20, 30), c(60, 10, 0), c(0, 50, 10))
  global <- test_global_discoveries(summary, null, 0.05)
  testthat::expect_equal(global$global_p_value, 3 / 4)
  testthat::expect_equal(global$null_mean_max_percentage, mean(c(30, 60, 50)))
  testthat::expect_identical(global$n_valid_resamples, 3L)
  testthat::expect_identical(select_optimal_threshold(summary, global), 1)

  # Repeating three null maxima gives 99 complete permutations, all below 50.
  low_null <- matrix(rep(c(10, 20, 30), 99), ncol = 3, byrow = TRUE)
  significant <- test_global_discoveries(summary, low_null, 0.05)
  testthat::expect_equal(significant$global_p_value, 0.01)
  testthat::expect_true(significant$signal_detected)
  testthat::expect_identical(select_optimal_threshold(summary, significant), 0.05)
  at_cutoff <- test_global_discoveries(summary, low_null, 0.01)
  testthat::expect_false(at_cutoff$signal_detected)
  testthat::expect_identical(select_optimal_threshold(summary, at_cutoff), 1)
})

testthat::test_that("global test handles unevaluable thresholds consistently", {
  summary <- data.frame(real_percentage = c(50, NA_real_))
  null <- rbind(c(10, NA), c(60, 20), c(NA, 30))
  global <- test_global_discoveries(summary, null, 0.05)
  testthat::expect_equal(global$global_p_value, 2 / 3)
  testthat::expect_identical(global$n_thresholds, 1L)
  testthat::expect_identical(global$n_valid_resamples, 2L)
  testthat::expect_error(
    test_global_discoveries(data.frame(real_percentage = NA_real_),
                            matrix(10), 0.05),
    "No filtering threshold"
  )
  testthat::expect_error(
    test_global_discoveries(data.frame(real_percentage = 50),
                            matrix(NA_real_), 0.05),
    "No null permutation"
  )
})

testthat::test_that("both discovery percentages are floored only for the ratio", {
  make_result <- function(q) {
    list(disco = data.frame(p = rep(0.5, length(q)), q = q))
  }
  summary <- build_threshold_summary(
    real_results = list(make_result(c(0.5, 0.5)), make_result(c(0.1, 0.5))),
    null_percentages = rbind(c(0, 0.5), c(0, 0.5)),
    filter_thresholds = c(0.05, 1), q_value = 0.2
  )
  testthat::expect_equal(summary$real_percentage, c(0, 50))
  testthat::expect_equal(summary$null_mean_percentage, c(0, 0.5))
  testthat::expect_equal(summary$real_percentage_floor_1, c(1, 50))
  testthat::expect_equal(summary$real_to_null_ratio, c(1, 50))
})

testthat::test_that("fixed-threshold mode skips empirical calculations", {
  set.seed(99)
  n <- 15
  signal1 <- stats::rnorm(n)
  signal2 <- stats::rnorm(n)
  data1 <- data.frame(
    a = signal1,
    b = signal1 + stats::rnorm(n, sd = 0.2),
    c = stats::rnorm(n)
  )
  data2 <- data.frame(
    a = signal2,
    b = -signal2 + stats::rnorm(n, sd = 0.2),
    c = stats::rnorm(n)
  )
  output_dir <- tempfile("disco-fixed-threshold-")

  result <- disco(
    list(data1),
    list(data2),
    empirically_derive_threshold = FALSE,
    filter_threshold = 0.05,
    q_value = -1,
    n_resamples = 0,
    random_seed = -1,
    global_p_cutoff = -1,
    output_dir = output_dir,
    prefix = "fixed"
  )

  testthat::expect_identical(result$selected_threshold, 0.05)
  testthat::expect_null(result$threshold_summary)
  testthat::expect_null(result$global_test)
  testthat::expect_null(result$discovery_plot)
  testthat::expect_null(result$ratio_plot)
  testthat::expect_null(result$feature_null_summary)
  testthat::expect_null(result$feature_null_plot)
  testthat::expect_identical(
    sort(list.files(output_dir)),
    sort(c("fixed_cor_diff.csv", "fixed_disco.csv"))
  )
})

testthat::test_that("threshold mode and fixed threshold are validated", {
  testthat::expect_error(
    disco(
      list(),
      list(),
      empirically_derive_threshold = NA,
      output_dir = tempdir(),
      prefix = "invalid-mode"
    ),
    "empirically_derive_threshold must be TRUE or FALSE"
  )
  testthat::expect_error(
    disco(
      list(data.frame(a = seq_len(10), b = seq_len(10))),
      list(data.frame(a = seq_len(10), b = seq_len(10))),
      empirically_derive_threshold = FALSE,
      filter_threshold = 0,
      output_dir = tempdir(),
      prefix = "invalid-threshold"
    ),
    "filter_threshold must be a single number"
  )
})

testthat::test_that("disco empirically selects a threshold and returns plots", {
  set.seed(2026)
  n <- 20
  signal1 <- stats::rnorm(n)
  signal2 <- stats::rnorm(n)
  data1 <- data.frame(
    a = signal1,
    b = signal1 + stats::rnorm(n, sd = 0.2),
    c = stats::rnorm(n)
  )
  data2 <- data.frame(
    a = signal2,
    b = -signal2 + stats::rnorm(n, sd = 0.2),
    c = stats::rnorm(n)
  )
  # A shared strongly correlated pair keeps the discovery denominator defined
  # even at the strictest threshold in each small-sample permutation.
  data1$d <- data1$c + stats::rnorm(n, sd = 0.05)
  data2$d <- data2$c + stats::rnorm(n, sd = 0.05)
  seed_before <- .Random.seed
  output_dir <- tempfile("disco-threshold-selection-")

  result <- NULL
  testthat::expect_message(
    result <- disco(
      list(data1),
      list(data2),
      n_resamples = 2,
      random_seed = 17,
      output_dir = output_dir,
      prefix = "empirical"
    ),
    "Selected correlation filtering threshold"
  )

  testthat::expect_identical(.Random.seed, seed_before)
  testthat::expect_identical(
    names(result),
    c(
      "cor_diff", "disco", "selected_threshold", "threshold_summary", "global_test",
      "discovery_plot", "ratio_plot", "feature_null_summary",
      "feature_null_plot", "output_files"
    )
  )
  testthat::expect_true(result$selected_threshold %in% c(
    0.001, 0.005, 0.01, 0.05, 0.1, 0.5, 1
  ))
  testthat::expect_identical(nrow(result$threshold_summary), 7L)
  testthat::expect_identical(sum(result$threshold_summary$selected), 1L)
  # Two permutations cannot produce p < 0.05; this exercises the fallback.
  testthat::expect_identical(result$selected_threshold, 1)
  testthat::expect_false(result$global_test$signal_detected)
  fixed <- calculate_fixed_threshold_disco(
    list(data1), list(data2), cooccur = 10, filt_threshold = 1,
    adj_method = "BY", output_dir = tempfile("disco-fallback-reference-"),
    prefix = "reference"
  )
  testthat::expect_equal(result$disco, fixed$disco)
  testthat::expect_equal(result$cor_diff, fixed$cor_diff)
  testthat::expect_s3_class(result$discovery_plot, "ggplot")
  testthat::expect_s3_class(result$ratio_plot, "ggplot")
  testthat::expect_s3_class(result$feature_null_plot, "ggplot")
  testthat::expect_equal(
    result$feature_null_summary$null_discovery_percentage,
    100 * result$feature_null_summary$n_null_discoveries / 2
  )
  testthat::expect_true(all(result$feature_null_summary$n_null_discoveries <=
                              result$feature_null_summary$n_null_tested))
  testthat::expect_true(all(file.exists(result$output_files)))
  testthat::expect_identical(
    sort(list.files(output_dir)),
    sort(c(
      "empirical_cor_diff.csv",
      "empirical_disco.csv",
      "empirical_threshold_summary.csv",
      "empirical_global_test.csv",
      "empirical_real_vs_null_thresholds.png",
      "empirical_real_to_null_ratio.png",
      "empirical_real_p_vs_null_discovery_by_threshold.csv",
      "empirical_real_p_vs_null_discovery_by_threshold.png"
    ))
  )
  # A relaxed cutoff exercises the ratio path with these same permutations.
  signal_result <- suppressMessages(disco(
    list(data1), list(data2), n_resamples = 2, random_seed = 17,
    global_p_cutoff = 0.5, output_dir = tempfile("disco-global-signal-"),
    prefix = "signal"
  ))
  testthat::expect_true(signal_result$global_test$signal_detected)
  testthat::expect_identical(
    signal_result$selected_threshold,
    select_optimal_threshold(signal_result$threshold_summary)
  )
  saved_global <- utils::read.csv(signal_result$output_files[["global_test"]])
  testthat::expect_equal(saved_global, signal_result$global_test)
})

testthat::test_that("two-modality permutations preserve partial sample overlap", {
  sample_ids <- paste0("sample", seq_len(12))
  data1 <- list(
    matrix(
      stats::rnorm(24),
      nrow = 12,
      dimnames = list(sample_ids, c("a", "b"))
    ),
    matrix(
      stats::rnorm(24),
      nrow = 12,
      dimnames = list(c(sample_ids[-12], "other"), c("x", "y"))
    )
  )
  data2 <- list(
    matrix(
      stats::rnorm(24),
      nrow = 12,
      dimnames = list(paste0("control", seq_len(12)), c("a", "b"))
    ),
    matrix(
      stats::rnorm(24),
      nrow = 12,
      dimnames = list(paste0("control", seq_len(12)), c("x", "y"))
    )
  )

  prepared <- suppressMessages(prepare_disco_inputs(data1, data2))
  pooled <- pool_null_inputs(prepared)
  set.seed(11)
  null_groups <- permute_null_groups(pooled)
  for (modality_index in seq_along(data1)) {
    testthat::expect_identical(
      nrow(null_groups$data1[[modality_index]]),
      nrow(data1[[modality_index]])
    )
    testthat::expect_identical(
      nrow(null_groups$data2[[modality_index]]),
      nrow(data2[[modality_index]])
    )
  }
  testthat::expect_identical(
    length(intersect(
      rownames(null_groups$data1[[1]]),
      rownames(null_groups$data1[[2]])
    )),
    length(intersect(rownames(data1[[1]]), rownames(data1[[2]])))
  )

  result <- NULL
  testthat::expect_no_error(
    result <- suppressMessages(disco(
      data1,
      data2,
      n_resamples = 1,
      output_dir = tempfile("disco-mismatched-samples-"),
      prefix = "mismatch"
    ))
  )
  testthat::expect_identical(result$selected_threshold, 1)
})
