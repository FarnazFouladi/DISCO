testthat::test_that("null discoveries match feature names and the chosen q cutoff", {
  result <- list(disco = data.frame(
    species = c("c", "b", "a", "d"),
    p = c(0.01, 0.02, NA, 0.01), q = c(0.09, 0.1, 0.01, NA)
  ))
  flags <- null_feature_indicators(result, c("a", "b", "c", "d", "absent"), 0.1)
  testthat::expect_identical(flags$tested, c(FALSE, TRUE, TRUE, FALSE, FALSE))
  testthat::expect_identical(flags$discovered, c(FALSE, FALSE, TRUE, FALSE, FALSE))
  empty <- list(disco = result$disco[FALSE, ])
  testthat::expect_identical(
    null_feature_indicators(empty, c("a", "b"), 0.2)$discovered,
    c(FALSE, FALSE)
  )
})

testthat::test_that("feature summaries retain p-values and use all permutations", {
  results <- list(
    list(disco = data.frame(species = c("b", "a"), p = c(0, 0.05), q = c(0.05, 0.1))),
    list(disco = data.frame(species = c("a", "b"), p = c(0.1, 0.2), q = c(0.1, 0.3)))
  )
  discoveries <- matrix(c(0L, 20L, 10L, 0L), nrow = 2,
                        dimnames = list(c("a", "b"), NULL))
  tested <- matrix(c(60L, 50L, 100L, 90L), nrow = 2,
                   dimnames = list(c("a", "b"), NULL))
  summary <- build_feature_null_summary(results, discoveries, tested,
                                        c(0.05, 1), 0.1, 100)
  testthat::expect_equal(summary$null_discovery_percentage, c(20, 0, 10, 0))
  testthat::expect_equal(summary$n_null_tested, c(50, 60, 100, 90))
  testthat::expect_identical(summary$real_significant, c(TRUE, FALSE, FALSE, FALSE))
  testthat::expect_identical(summary$real_disco_p[[1L]], 0)
  testthat::expect_identical(summary$neg_log10_real_disco_p[[1L]], 8)
  testthat::expect_true(summary$p_value_floored_for_plot[[1L]])
  plot <- create_feature_null_plot(summary, c(0.05, 1), 0.1, 100)
  testthat::expect_s3_class(plot, "ggplot")
  testthat::expect_equal(plot$layers[[2L]]$data$annotation,
                         c("Real q < 0.1: n = 1", "Real q < 0.1: n = 0"))
  testthat::expect_equal(as.character(plot$data$real_q_category), c(
    "Real DISCO q < 0.1", rep("Real DISCO q >= 0.1", 3)
  ))
})

testthat::test_that("empty threshold panels report zero significant features", {
  empty <- list(disco = data.frame(species = character(), p = numeric(), q = numeric()))
  counts <- matrix(0L, nrow = 1, ncol = 1, dimnames = list("a", NULL))
  summary <- build_feature_null_summary(list(empty), counts, counts, 0.05, 0.2, 100)
  testthat::expect_equal(nrow(summary), 0L)
  plot <- create_feature_null_plot(summary, 0.05, 0.2, 100)
  testthat::expect_equal(plot$layers[[2L]]$data$annotation, "Real q < 0.2: n = 0")
  testthat::expect_no_error(ggplot2::ggplot_build(plot))
})
