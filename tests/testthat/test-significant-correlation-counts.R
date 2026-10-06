testthat::test_that("one-modality counts include both edge orientations", {
  significant_pairs <- data.frame(
    t1 = c("a", "b", "b"),
    t2 = c("b", "c", "d")
  )

  result <- add_significant_correlation_counts(
    significant_pairs,
    two_modal = FALSE
  )

  testthat::expect_identical(result$num_sig_corr_t1, c(1L, 3L, 3L))
  testthat::expect_identical(result$num_sig_corr_t2, c(3L, 1L, 1L))
})

testthat::test_that("two-modality counts remain side-specific", {
  significant_pairs <- data.frame(
    t1 = c("a", "a", "b"),
    t2 = c("x", "y", "x")
  )

  result <- add_significant_correlation_counts(
    significant_pairs,
    two_modal = TRUE
  )

  testthat::expect_identical(result$num_sig_corr_t1, c(2L, 2L, 1L))
  testthat::expect_identical(result$num_sig_corr_t2, c(2L, 1L, 2L))
})

testthat::test_that("empty significant results retain count columns", {
  significant_pairs <- data.frame(
    t1 = character(),
    t2 = character()
  )

  result <- add_significant_correlation_counts(
    significant_pairs,
    two_modal = FALSE
  )

  testthat::expect_identical(result$num_sig_corr_t1, integer())
  testthat::expect_identical(result$num_sig_corr_t2, integer())
  testthat::expect_equal(nrow(result), 0)
})
