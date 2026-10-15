testthat::test_that("the Cauchy tolerance is not user-configurable", {
  testthat::expect_false("tol" %in% names(formals(disco)))
  testthat::expect_false("tol" %in% names(formals(safe_tan_cauchy)))
})

testthat::test_that("Cauchy transformation uses a fixed 1e-8 tail boundary", {
  boundary_term <- 1 / tan(pi * 1e-8)

  testthat::expect_equal(
    safe_tan_cauchy(0),
    boundary_term
  )
  testthat::expect_equal(
    safe_tan_cauchy(1),
    -boundary_term
  )
})

testthat::test_that("complementary extreme p-values cancel symmetrically", {
  testthat::expect_equal(
    safe_tan_cauchy(c(0, 1)),
    0
  )
})

testthat::test_that("stable and conventional Cauchy formulas agree away from boundaries", {
  p_values <- c(0.1, 0.2, 0.5, 0.8, 0.9)
  expected <- mean(tan((0.5 - p_values) * pi))

  testthat::expect_equal(safe_tan_cauchy(p_values), expected)
})

testthat::test_that("disco uses the fixed Cauchy boundary", {
  feature_a <- seq_len(20)
  data1 <- data.frame(a = feature_a, b = feature_a)
  data2 <- data.frame(a = feature_a, b = -feature_a)

  invisible(utils::capture.output(
    result <- disco(
      list(data1),
      list(data2),
      cooccur = 10,
      empirically_derive_threshold = FALSE,
      filter_threshold = 1,
      output_dir = tempfile("disco-fixed-cauchy-boundary-"),
      prefix = "default"
    )
  ))

  expected_p <- stats::pcauchy(
    1 / tan(pi * 1e-8),
    lower.tail = FALSE
  )
  testthat::expect_equal(result$disco$p, rep(expected_p, 2L))
})
