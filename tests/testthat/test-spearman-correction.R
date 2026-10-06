testthat::test_that("Spearman standard-error corrections match their formulas", {
  R <- c(0, 0.5, 0.8)
  N <- c(20, 30, 40)

  testthat::expect_equal(
    spearman_fisher_se(R, N, correction = "bonett_wright"),
    sqrt((1 + R^2 / 2) / (N - 3))
  )
  testthat::expect_equal(
    spearman_fisher_se(R, N, correction = "fieller_1.06"),
    sqrt(1.06 / (N - 3))
  )
})

testthat::test_that("Bonett-Wright is the default correction", {
  R <- 0.7
  N <- 25

  testthat::expect_equal(
    spearman_fisher_se(R, N),
    spearman_fisher_se(R, N, correction = "bonett_wright")
  )
  testthat::expect_identical(
    eval(formals(disco)$spearman_correction),
    c("bonett_wright", "fieller_1.06")
  )
})

testthat::test_that("single-sample Fisher tests use the selected correction", {
  R <- 0.6
  N <- 30
  z_bc <- atanh(R) - R / (2 * (N - 1))

  expected_bonett <- 2 * stats::pnorm(
    -abs(z_bc / sqrt((1 + R^2 / 2) / (N - 3)))
  )
  expected_fieller <- 2 * stats::pnorm(
    -abs(z_bc / sqrt(1.06 / (N - 3)))
  )

  testthat::expect_equal(
    fisher_z_p(R, N, correction = "bonett_wright"),
    expected_bonett
  )
  testthat::expect_equal(
    fisher_z_p(R, N, correction = "fieller_1.06"),
    expected_fieller
  )
})

testthat::test_that("two-sample Fisher tests use group-specific corrections", {
  R1 <- 0.3
  N1 <- 25
  R2 <- 0.8
  N2 <- 35
  z1_bc <- atanh(R1) - R1 / (2 * (N1 - 1))
  z2_bc <- atanh(R2) - R2 / (2 * (N2 - 1))

  bonett_se <- sqrt(
    (1 + R1^2 / 2) / (N1 - 3) +
      (1 + R2^2 / 2) / (N2 - 3)
  )
  fieller_se <- sqrt(
    1.06 / (N1 - 3) + 1.06 / (N2 - 3)
  )

  bonett <- fisher_z_p_2sample(
    R1, N1, R2, N2, correction = "bonett_wright"
  )
  fieller <- fisher_z_p_2sample(
    R1, N1, R2, N2, correction = "fieller_1.06"
  )

  testthat::expect_equal(
    bonett,
    2 * stats::pnorm(-abs((z1_bc - z2_bc) / bonett_se))
  )
  testthat::expect_equal(
    fieller,
    2 * stats::pnorm(-abs((z1_bc - z2_bc) / fieller_se))
  )
})

testthat::test_that("invalid Spearman corrections are rejected", {
  testthat::expect_error(
    spearman_fisher_se(0.5, 20, correction = "pearson"),
    "should be one of"
  )
})

testthat::test_that("eligible comparisons propagate the selected correction", {
  feature_names <- c("a", "b")
  R1 <- matrix(c(1, 0.3, 0.3, 1), nrow = 2,
               dimnames = list(feature_names, feature_names))
  R2 <- matrix(c(1, 0.8, 0.8, 1), nrow = 2,
               dimnames = list(feature_names, feature_names))
  C1 <- matrix(25, nrow = 2, ncol = 2,
               dimnames = list(feature_names, feature_names))
  C2 <- matrix(35, nrow = 2, ncol = 2,
               dimnames = list(feature_names, feature_names))

  bonett <- compare_eligible_correlations(
    R1, C1, R2, C2,
    cooccur = 10,
    filt_threshold = 1,
    correction = "bonett_wright"
  )
  fieller <- compare_eligible_correlations(
    R1, C1, R2, C2,
    cooccur = 10,
    filt_threshold = 1,
    correction = "fieller_1.06"
  )

  testthat::expect_false(isTRUE(all.equal(
    bonett$pvals["a", "b"],
    fieller$pvals["a", "b"]
  )))
})

testthat::test_that("disco exposes both Spearman corrections", {
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
    bonett <- disco(
      list(data1),
      list(data2),
      cooccur = 10,
      filt_threshold = 1,
      output_dir = tempfile("disco-bonett-"),
      prefix = "bonett",
      spearman_correction = "bonett_wright"
    )
  ))
  invisible(utils::capture.output(
    fieller <- disco(
      list(data1),
      list(data2),
      cooccur = 10,
      filt_threshold = 1,
      output_dir = tempfile("disco-fieller-"),
      prefix = "fieller",
      spearman_correction = "fieller_1.06"
    )
  ))

  testthat::expect_gt(
    abs(
      bonett$cor_diff$p - fieller$cor_diff$p
    ),
    0
  )
})
