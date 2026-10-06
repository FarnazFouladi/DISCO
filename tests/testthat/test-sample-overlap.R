testthat::test_that("unequal modality sample sets report their overlap", {
  modality1 <- matrix(
    seq_len(6),
    nrow = 3,
    dimnames = list(c("sample1", "sample2", "sample3"), c("a", "b"))
  )
  modality2 <- matrix(
    seq_len(6),
    nrow = 3,
    dimnames = list(c("sample2", "sample3", "sample4"), c("c", "d"))
  )

  testthat::expect_message(
    n_common <- report_sample_overlap(modality1, modality2, group = "data1"),
    "data1 modalities have 2 common sample\\(s\\)"
  )
  testthat::expect_equal(n_common, 2)
})

testthat::test_that("zero overlap between modalities is rejected", {
  modality1 <- matrix(
    seq_len(4),
    nrow = 2,
    dimnames = list(c("sample1", "sample2"), c("a", "b"))
  )
  modality2 <- matrix(
    seq_len(4),
    nrow = 2,
    dimnames = list(c("sample3", "sample4"), c("c", "d"))
  )

  testthat::expect_error(
    report_sample_overlap(modality1, modality2, group = "data1"),
    "data1 modalities have no common samples"
  )
})

testthat::test_that("disco_score checks modality overlap before merging", {
  data1 <- list(
    matrix(
      seq_len(8),
      nrow = 4,
      dimnames = list(paste0("a", 1:4), c("feature1", "feature2"))
    ),
    matrix(
      seq_len(8),
      nrow = 4,
      dimnames = list(paste0("b", 1:4), c("feature3", "feature4"))
    )
  )
  data2 <- list(
    matrix(
      seq_len(8),
      nrow = 4,
      dimnames = list(paste0("c", 1:4), c("feature1", "feature2"))
    ),
    matrix(
      seq_len(8),
      nrow = 4,
      dimnames = list(paste0("d", 1:4), c("feature3", "feature4"))
    )
  )

  testthat::expect_error(
    disco_score(
      data1,
      data2,
      output_dir = tempfile("disco-sample-overlap-"),
      prefix = "no-overlap"
    ),
    "data1 modalities have no common samples"
  )
})
