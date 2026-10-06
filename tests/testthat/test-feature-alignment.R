testthat::test_that("common features use data1 order", {
  data1 <- matrix(
    seq_len(12),
    nrow = 3,
    dimnames = list(NULL, c("feature_b", "feature_a", "data1_only", "feature_c"))
  )
  data2 <- matrix(
    seq_len(12),
    nrow = 3,
    dimnames = list(NULL, c("feature_c", "data2_only", "feature_b", "feature_a"))
  )

  testthat::expect_warning(
    aligned <- align_common_features(data1, data2, modality = 1),
    "dropped 1 feature\\(s\\) from data1 and 1 from data2; reordered data2"
  )

  expected_features <- c("feature_b", "feature_a", "feature_c")
  testthat::expect_identical(colnames(aligned$data1), expected_features)
  testthat::expect_identical(colnames(aligned$data2), expected_features)
  testthat::expect_equal(
    unname(aligned$data2),
    unname(data2[, expected_features, drop = FALSE])
  )
})

testthat::test_that("already aligned features are unchanged", {
  data1 <- matrix(seq_len(9), nrow = 3,
                  dimnames = list(NULL, c("a", "b", "c")))
  data2 <- matrix(seq_len(9) + 10, nrow = 3,
                  dimnames = list(NULL, c("a", "b", "c")))

  testthat::expect_no_warning(
    aligned <- align_common_features(data1, data2, modality = 1)
  )

  testthat::expect_identical(aligned$data1, data1)
  testthat::expect_identical(aligned$data2, data2)
})

testthat::test_that("invalid feature names produce informative errors", {
  named <- matrix(seq_len(6), nrow = 3,
                  dimnames = list(NULL, c("a", "b")))
  unnamed <- unname(named)
  duplicated <- named
  colnames(duplicated) <- c("a", "a")
  different <- named
  colnames(different) <- c("x", "y")

  testthat::expect_error(
    align_common_features(unnamed, named, modality = 1),
    "must have non-missing, non-empty feature names"
  )
  testthat::expect_error(
    align_common_features(duplicated, named, modality = 1),
    "Feature names in data1 modality 1 must be unique"
  )
  testthat::expect_error(
    align_common_features(named, different, modality = 1),
    "No common features were found for modality 1"
  )
})

testthat::test_that("each modality is aligned independently", {
  group1 <- list(
    matrix(seq_len(12), nrow = 4,
           dimnames = list(NULL, c("a", "b", "first_only"))),
    matrix(seq_len(8), nrow = 4,
           dimnames = list(NULL, c("x", "y")))
  )
  group2 <- list(
    matrix(seq_len(12), nrow = 4,
           dimnames = list(NULL, c("b", "a", "second_only"))),
    matrix(seq_len(8), nrow = 4,
           dimnames = list(NULL, c("y", "x")))
  )

  aligned <- lapply(seq_along(group1), function(i) {
    result <- NULL
    testthat::expect_warning(
      result <- align_common_features(group1[[i]], group2[[i]], modality = i)
    )
    result
  })

  testthat::expect_identical(colnames(aligned[[1]]$data1), c("a", "b"))
  testthat::expect_identical(colnames(aligned[[1]]$data2), c("a", "b"))
  testthat::expect_identical(colnames(aligned[[2]]$data1), c("x", "y"))
  testthat::expect_identical(colnames(aligned[[2]]$data2), c("x", "y"))
})
