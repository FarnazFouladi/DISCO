testthat::test_that("pair-level significance arguments are not public", {
  testthat::expect_false(any(
    c("diff_threshold", "p_alpha", "q_alpha") %in% names(formals(disco))
  ))
})
