testthat::test_that("presentation changes preserve the previous numerical results", {
  reference <- readRDS(testthat::test_path("fixtures", "output-reference.rds"))
  result <- disco(reference$data1, reference$data2,
                  output_dir = tempfile("csv-output-"), prefix = "example")
  testthat::expect_null(result$global_test)
  testthat::expect_identical(result$selected_threshold, 1)
  testthat::expect_identical(names(result$disco),
    c("Feature", "Pairs", "DISCO_score", "DSICO_p", "DISCO_q"))
  testthat::expect_false(is.unsorted(-result$disco$DISCO_score))
  before <- reference$expected$disco
  index <- match(result$disco$Feature, before$species)
  testthat::expect_equal(result$disco$Pairs, before$num_corr[index])
  testthat::expect_equal(result$disco$DISCO_score, before$t[index], tolerance = 1e-14)
  testthat::expect_equal(result$disco$DSICO_p, before$p[index], tolerance = 1e-14)
  testthat::expect_equal(result$disco$DISCO_q, before$q[index], tolerance = 1e-14)
  expected_pairs <- reference$expected$cor_diff
  names(expected_pairs)[1:2] <- c("feature_1", "feature_2")
  testthat::expect_equal(result$cor_diff, expected_pairs, tolerance = 1e-14)
  testthat::expect_true(all(endsWith(result$output_files, ".csv")))
  testthat::expect_equal(utils::read.csv(result$output_files[["disco"]]),
                         result$disco)
})

testthat::test_that("CSV exports preserve quoted names, missing values and empty headers", {
  path <- tempfile(fileext = ".csv")
  data <- data.frame(Feature = c("a,b", 'c"d', "e"),
                     DISCO_score = c(5, NA, 1))
  write_result_csv(data, path)
  testthat::expect_equal(utils::read.csv(path), data)
  testthat::expect_identical(names(utils::read.csv(path)), names(data))
  write_result_csv(data[FALSE, ], path)
  testthat::expect_identical(names(utils::read.csv(path)), names(data))
  testthat::expect_equal(nrow(utils::read.csv(path)), 0)
})
