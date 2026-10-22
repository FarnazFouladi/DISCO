testthat::test_that("worker counts must be positive integers", {
  for (workers in list(0, -1, 1.5, NA_real_, Inf, TRUE, "2", c(1, 2))) {
    testthat::expect_error(validate_null_workers(workers), "workers")
  }
  testthat::expect_identical(validate_null_workers(2), 2L)
})

testthat::test_that("socket workers preserve sequential results and random state", {
  set.seed(729)
  make_input <- function(prefix) {
    x <- matrix(stats::rnorm(144), 24, 6,
                dimnames = list(paste0(prefix, 1:24), letters[1:6]))
    x[1:2, 1] <- NA_real_
    x
  }
  data1 <- list(make_input("control"))
  data2 <- list(make_input("tumor"))
  for (two_modal in c(FALSE, TRUE)) {
    if (two_modal) {
      data1[[2]] <- make_input("control")[3:24, , drop = FALSE]
      data2[[2]] <- make_input("tumor")[2:24, , drop = FALSE]
      colnames(data1[[2]]) <- colnames(data2[[2]]) <- paste0("path", 1:6)
    }
    run <- function(workers) {
      suppressMessages(run_empirical_disco(
        data1, data2, cooccur = 4, filter_thresholds = c(0.001, 0.05, 1),
        q_value = 0.2, n_resamples = 3, random_seed = 17,
        adj_method = "BY", output_dir = tempfile("parallel-disco-"),
        prefix = "test", workers = workers
      ))
    }
    seed_before <- .Random.seed
    sequential <- run(1)
    testthat::expect_identical(.Random.seed, seed_before)
    concurrent <- run(2)
    testthat::expect_identical(.Random.seed, seed_before)
    for (name in c("cor_diff", "disco", "selected_threshold", "global_test",
                   "threshold_summary", "feature_null_summary")) {
      testthat::expect_identical(sequential[[name]], concurrent[[name]], info = name)
    }
    # Pre-drawing assignments must also match the original serial label draws.
    prepared <- prepare_disco_inputs(data1, data2, report_overlap = FALSE)
    pooled <- pool_null_inputs(prepared)
    set.seed(17)
    original <- lapply(1:3, function(i) permute_null_groups(pooled))
    set.seed(17)
    assignments <- lapply(1:3, function(i) draw_null_assignment(pooled))
    reconstructed <- lapply(assignments, function(a) permute_null_groups(pooled, a))
    testthat::expect_identical(original, reconstructed)
  }
})

testthat::test_that("parallel errors propagate and release workers", {
  before <- nrow(showConnections())
  testthat::expect_error(
    run_null_permutations(
      list(list(data1 = "a", data2 = "b"), list(data1 = "a", data2 = "b")),
      pooled = list(pooled = list(matrix(1, 1, 1))), cooccur = 4,
      filter_thresholds = 1, adj_method = "BY", feature_names = "x",
      q_value = 0.2, workers = 2
    ),
    "dimnames"
  )
  testthat::expect_identical(nrow(showConnections()), before)
})
