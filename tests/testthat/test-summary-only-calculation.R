make_summary_test_inputs <- function(two_modal = FALSE) {
  set.seed(271)
  make_matrix <- function(prefix) {
    x <- matrix(stats::rnorm(120), 24, 5,
                dimnames = list(paste0(prefix, seq_len(24)), letters[1:5]))
    x[1:2, 1] <- NA_real_
    x
  }
  data1 <- list(make_matrix("control"))
  data2 <- list(make_matrix("tumor"))
  if (two_modal) {
    data1[[2]] <- make_matrix("control")[3:24, , drop = FALSE]
    data2[[2]] <- make_matrix("tumor")[2:24, , drop = FALSE]
    colnames(data1[[2]]) <- colnames(data2[[2]]) <- paste0("path", 1:5)
  }
  list(data1 = data1, data2 = data2)
}

testthat::test_that("summary-only calculations preserve feature statistics", {
  for (two_modal in c(FALSE, TRUE)) {
    inputs <- make_summary_test_inputs(two_modal)
    prepared <- prepare_disco_inputs(inputs$data1, inputs$data2,
                                     report_overlap = FALSE)
    base <- calculate_disco_base(prepared, cooccur = 4)
    for (threshold in c(1e-100, 0.001, 0.05, 1)) {
      for (method in c("BY", "BH", "none")) {
        full <- calculate_disco_at_threshold(base, threshold, method)
        summary <- calculate_disco_at_threshold(
          base, threshold, method, include_cor_diff = FALSE
        )
        testthat::expect_null(summary$cor_diff)
        testthat::expect_identical(summary$disco, full$disco)
        testthat::expect_identical(
          disco_discovery_percentage(summary, 0.2),
          disco_discovery_percentage(full, 0.2)
        )
        pvals <- filter_correlation_tests(base$correlation_tests, threshold)$pvals
        testthat::expect_identical(
          build_cor_diff_table(base, pvals, method), full$cor_diff
        )
      }
    }
  }
})

testthat::test_that("empirical summaries match the full-table path with one output table", {
  inputs <- make_summary_test_inputs()
  calculate_original <- calculate_disco_at_threshold
  build_original <- build_cor_diff_table
  table_calls <- 0L
  testthat::local_mocked_bindings(
    build_cor_diff_table = function(...) {
      table_calls <<- table_calls + 1L
      build_original(...)
    }
  )
  run <- function() {
    run_empirical_disco(
      inputs$data1, inputs$data2, cooccur = 4,
      filter_thresholds = c(0.001, 0.05, 1), q_value = 0.2,
      n_resamples = 3, random_seed = 17, adj_method = "BY",
      output_dir = tempfile("summary-only-"), prefix = "test"
    )
  }
  optimized <- suppressMessages(run())
  testthat::expect_identical(table_calls, 1L)
  # Emulate the old path: materialize a pairwise table at every threshold.
  testthat::local_mocked_bindings(
    calculate_disco_at_threshold = function(base, filt_threshold, adj_method,
                                             include_cor_diff = TRUE) {
      calculate_original(base, filt_threshold, adj_method, include_cor_diff = TRUE)
    }
  )
  reference <- suppressMessages(run())
  for (name in c("cor_diff", "disco", "selected_threshold",
                 "threshold_summary", "global_test", "feature_null_summary")) {
    testthat::expect_identical(optimized[[name]], reference[[name]], info = name)
  }
})
