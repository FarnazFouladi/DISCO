testthat::test_that("pair-level significance arguments are not public", {
  testthat::expect_false(any(
    c("diff_threshold", "p_alpha", "q_alpha") %in% names(formals(disco))
  ))
  testthat::expect_false("filt_threshold" %in% names(formals(disco)))
  testthat::expect_false("filter_thresholds" %in% names(formals(disco)))
  testthat::expect_identical(
    formals(disco)$empirically_derive_threshold,
    TRUE
  )
  testthat::expect_identical(formals(disco)$filter_threshold, 0.05)
  testthat::expect_identical(formals(disco)$q_value, 0.2)
})
