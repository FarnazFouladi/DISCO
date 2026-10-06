testthat::test_that("DISCO uses a 1e-16 default Cauchy tolerance", {
  testthat::expect_equal(formals(disco)$tol, 1e-16)
  testthat::expect_equal(formals(safe_tan_cauchy)$tol, 1e-16)
})

testthat::test_that("Cauchy transformation uses the supplied tolerance", {
  p_values <- c(1e-50, 0.2, 1 - 1e-50)
  tol <- 1e-16
  bounded <- pmax(pmin(p_values, 1 - tol), tol)
  expected <- 1 / length(bounded) * sum(tan((0.5 - bounded) * pi))

  testthat::expect_equal(
    safe_tan_cauchy(p_values, tol = tol),
    expected
  )

  larger_tol <- 1e-4
  bounded_larger <- pmax(pmin(p_values, 1 - larger_tol), larger_tol)
  expected_larger <- 1 / length(bounded_larger) *
    sum(tan((0.5 - bounded_larger) * pi))

  testthat::expect_equal(
    safe_tan_cauchy(p_values, tol = larger_tol),
    expected_larger
  )
})

testthat::test_that("invalid Cauchy tolerances are rejected", {
  invalid_tolerances <- list(
    0,
    -1e-16,
    0.5,
    1,
    NA_real_,
    Inf,
    c(1e-16, 1e-15),
    "1e-16"
  )

  for (tol in invalid_tolerances) {
    testthat::expect_error(
      validate_cauchy_tolerance(tol),
      "tol must be a single finite number"
    )
  }
})

testthat::test_that("disco propagates tol to Cauchy combination", {
  feature_a <- seq_len(20)
  data1 <- data.frame(a = feature_a, b = feature_a)
  data2 <- data.frame(a = feature_a, b = -feature_a)

  default_output <- tempfile("disco-tol-default-")
  larger_output <- tempfile("disco-tol-larger-")

  invisible(utils::capture.output(
    default_result <- disco(
      list(data1),
      list(data2),
      cooccur = 10,
      filt_threshold = 1,
      output_dir = default_output,
      prefix = "default"
    )
  ))
  invisible(utils::capture.output(
    larger_result <- disco(
      list(data1),
      list(data2),
      cooccur = 10,
      filt_threshold = 1,
      output_dir = larger_output,
      prefix = "larger",
      tol = 1e-4
    )
  ))

  testthat::expect_true(all(default_result$disco$p < larger_result$disco$p))
})
