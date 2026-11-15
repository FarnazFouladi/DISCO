testthat::test_that("Spearman standard error uses the Bonett-Wright formula", {
  R <- c(0, 0.5, 0.8)
  N <- c(20, 30, 40)

  testthat::expect_equal(
    spearman_fisher_se(R, N),
    sqrt((1 + R^2 / 2) / (N - 3))
  )
})

testthat::test_that("disco does not expose a Spearman-correction argument", {
  testthat::expect_false(
    "spearman_correction" %in% names(formals(disco))
  )
  testthat::expect_false(
    "correction" %in% names(formals(spearman_fisher_se))
  )
})

testthat::test_that("single-sample Fisher tests use Bonett-Wright", {
  R <- 0.6
  N <- 30
  z_bc <- atanh(R) - R / (2 * (N - 1))

  expected <- 2 * stats::pnorm(
    -abs(z_bc / sqrt((1 + R^2 / 2) / (N - 3)))
  )

  testthat::expect_equal(
    fisher_z_p(R, N),
    expected
  )
})

testthat::test_that("two-sample Fisher tests use group-specific Bonett-Wright SEs", {
  R1 <- 0.3
  N1 <- 25
  R2 <- 0.8
  N2 <- 35
  z1_bc <- atanh(R1) - R1 / (2 * (N1 - 1))
  z2_bc <- atanh(R2) - R2 / (2 * (N2 - 1))

  expected_se <- sqrt(
    (1 + R1^2 / 2) / (N1 - 3) +
      (1 + R2^2 / 2) / (N2 - 3)
  )

  testthat::expect_equal(
    fisher_z_p_2sample(R1, N1, R2, N2),
    2 * stats::pnorm(-abs((z1_bc - z2_bc) / expected_se))
  )
})

testthat::test_that("eligible comparisons use Bonett-Wright", {
  feature_names <- c("a", "b")
  R1 <- matrix(c(1, 0.3, 0.3, 1), nrow = 2,
               dimnames = list(feature_names, feature_names))
  R2 <- matrix(c(1, 0.8, 0.8, 1), nrow = 2,
               dimnames = list(feature_names, feature_names))
  C1 <- matrix(25, nrow = 2, ncol = 2,
               dimnames = list(feature_names, feature_names))
  C2 <- matrix(35, nrow = 2, ncol = 2,
               dimnames = list(feature_names, feature_names))

  comparison <- compare_eligible_correlations(
    R1, C1, R2, C2,
    cooccur = 10,
    filt_threshold = 1
  )

  testthat::expect_equal(
    comparison$pvals["a", "b"],
    fisher_z_p_2sample(0.3, 25, 0.8, 35)
  )
})

testthat::test_that("disco runs with the fixed Bonett-Wright correction", {
  set.seed(2026)
  group1_a <- stats::rnorm(30)
  group2_a <- stats::rnorm(30)
  data1 <- data.frame(
    a = group1_a,
    b = 0.7 * group1_a + stats::rnorm(30, sd = 0.6)
  )
  data2 <- data.frame(
    a = group2_a,
    b = -0.5 * group2_a + stats::rnorm(30, sd = 0.8)
  )

  invisible(utils::capture.output(
    result <- disco(
      list(data1),
      list(data2),
      cooccur = 10,
      empirically_derive_threshold = FALSE,
      filter_threshold = 1,
      output_dir = tempfile("disco-bonett-"),
      prefix = "bonett"
    )
  ))

  testthat::expect_true(nrow(result$cor_diff) > 0L)
  testthat::expect_true(all(is.finite(result$cor_diff$p)))
})
