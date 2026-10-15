testthat::test_that("both groups must meet the co-occurrence threshold", {
  feature_names <- c("a", "b", "c")
  R1 <- matrix(
    c(
      1.0, 0.2, 0.4,
      0.2, 1.0, 0.3,
      0.4, 0.3, 1.0
    ),
    nrow = 3,
    dimnames = list(feature_names, feature_names)
  )
  R2 <- matrix(
    c(
      1.0, 0.6, 0.7,
      0.6, 1.0, 0.8,
      0.7, 0.8, 1.0
    ),
    nrow = 3,
    dimnames = list(feature_names, feature_names)
  )
  C1 <- matrix(
    c(
      20, 10, 9,
      10, 20, 15,
      9, 15, 20
    ),
    nrow = 3,
    dimnames = list(feature_names, feature_names)
  )
  C2 <- matrix(
    c(
      20, 12, 20,
      12, 20, 8,
      20, 8, 20
    ),
    nrow = 3,
    dimnames = list(feature_names, feature_names)
  )

  result <- compare_eligible_correlations(
    R1, C1, R2, C2,
    cooccur = 10,
    filt_threshold = 1
  )

  expected_eligible <- matrix(
    c(
      FALSE, TRUE, FALSE,
      TRUE, FALSE, FALSE,
      FALSE, FALSE, FALSE
    ),
    nrow = 3,
    dimnames = list(feature_names, feature_names)
  )

  testthat::expect_identical(result$eligible, expected_eligible)
  testthat::expect_false(is.na(result$pvals["a", "b"]))
  testthat::expect_false(is.na(result$pvals["b", "a"]))
  testthat::expect_true(is.na(result$pvals["a", "c"]))
  testthat::expect_true(is.na(result$pvals["b", "c"]))
  testthat::expect_true(all(is.na(diag(result$pvals))))
})

testthat::test_that("eligible pairs still undergo single-group screening", {
  feature_names <- c("a", "b")
  R1 <- matrix(c(1, 0, 0, 1), nrow = 2,
               dimnames = list(feature_names, feature_names))
  R2 <- R1
  C1 <- matrix(20, nrow = 2, ncol = 2,
               dimnames = list(feature_names, feature_names))
  C2 <- C1

  result <- compare_eligible_correlations(
    R1, C1, R2, C2,
    cooccur = 10,
    filt_threshold = 0.001
  )

  testthat::expect_true(result$eligible["a", "b"])
  testthat::expect_false(result$retained["a", "b"])
  testthat::expect_true(is.na(result$pvals["a", "b"]))
})

testthat::test_that("no eligible pairs return stable empty result schemas", {
  set.seed(123)
  feature_names <- c("a", "b", "c")
  data1 <- as.data.frame(matrix(
    stats::rnorm(18),
    nrow = 6,
    dimnames = list(NULL, feature_names)
  ))
  data2 <- as.data.frame(matrix(
    stats::rnorm(18),
    nrow = 6,
    dimnames = list(NULL, feature_names)
  ))

  prepared <- prepare_disco_inputs(list(data1), list(data2))
  base <- calculate_disco_base(prepared, cooccur = 7)
  result <- calculate_disco_at_threshold(
    base,
    filt_threshold = 1,
    adj_method = "BY"
  )

  testthat::expect_identical(
    names(result),
    c("cor_diff", "disco")
  )
  testthat::expect_equal(nrow(result$cor_diff), 0)
  testthat::expect_identical(
    names(result$disco),
    c("species", "t", "p", "num_corr", "q")
  )
  testthat::expect_equal(nrow(result$disco), 0)
})
