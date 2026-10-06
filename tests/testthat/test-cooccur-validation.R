testthat::test_that("cooccur accepts integer values of at least 4", {
  valid_values <- list(4, 4L, 10)

  for (cooccur in valid_values) {
    testthat::expect_equal(validate_cooccur(cooccur), cooccur)
  }
})

testthat::test_that("disco rejects invalid cooccur values", {
  invalid_values <- list(
    3,
    4.5,
    NA_real_,
    NaN,
    Inf,
    c(4, 5),
    "4",
    NULL
  )

  for (cooccur in invalid_values) {
    testthat::expect_error(
      disco(
        data1 = list(),
        data2 = list(),
        cooccur = cooccur,
        output_dir = tempfile("disco-cooccur-"),
        prefix = "invalid"
      ),
      "cooccur must be a single integer greater than or equal to 4"
    )
  }
})
