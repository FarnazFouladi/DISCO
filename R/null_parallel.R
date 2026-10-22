# Keep parallelism opt-in: memory use grows with the number of workers.
validate_null_workers <- function(workers) {
  if (!is.numeric(workers) || length(workers) != 1L || is.na(workers) ||
      !is.finite(workers) || workers < 1 || workers != floor(workers) ||
      workers > .Machine$integer.max) {
    stop("workers must be a single positive integer", call. = FALSE)
  }
  as.integer(workers)
}

# One permutation returns small feature summaries, never correlation matrices.
calculate_null_summary <- function(assignment, pooled, cooccur,
                                    filter_thresholds, adj_method,
                                    feature_names, q_value) {
  groups <- permute_null_groups(pooled, assignment)
  prepared <- prepare_disco_inputs(groups$data1, groups$data2,
                                   report_overlap = FALSE)
  base <- calculate_disco_base(prepared, cooccur)
  discovered <- tested <- matrix(FALSE, length(feature_names),
                                  length(filter_thresholds))
  percentages <- numeric(length(filter_thresholds))
  for (i in seq_along(filter_thresholds)) {
    result <- calculate_disco_at_threshold(
      base, filter_thresholds[[i]], adj_method, include_cor_diff = FALSE
    )
    percentages[[i]] <- disco_discovery_percentage(result, q_value)$percentage
    indicators <- null_feature_indicators(result, feature_names, q_value)
    discovered[, i] <- indicators$discovered
    tested[, i] <- indicators$tested
  }
  list(percentages = percentages, discovered = discovered, tested = tested)
}

# Export only the calculation helpers and pooled inputs. Rebinding helpers to
# this small environment avoids capturing real-data matrices and ensures that
# devtools::load_all() uses the current code, not an older installed DISCO.
make_null_worker <- function(pooled, cooccur, filter_thresholds, adj_method,
                             feature_names, q_value) {
  worker_env <- new.env(parent = baseenv())
  helper_names <- c(
    "calculate_null_summary", "permute_null_groups", "prepare_disco_inputs",
    "align_common_features", "calculate_disco_base", "get_npairs", "get_corr",
    "calculate_correlation_tests", "fisher_z_p", "fisher_z_p_2sample",
    "spearman_fisher_se", "calculate_disco_at_threshold",
    "filter_correlation_tests", "safe_tan_cauchy",
    "disco_discovery_percentage", "null_feature_indicators"
  )
  for (name in helper_names) {
    helper <- get(name, envir = environment(make_null_worker))
    environment(helper) <- worker_env
    assign(name, helper, envir = worker_env)
  }
  worker_env$pnorm <- stats::pnorm
  worker_env$`%>%` <- dplyr::`%>%`
  worker_env$column_to_rownames <- tibble::column_to_rownames
  worker_env$pooled <- pooled
  worker_env$cooccur <- cooccur
  worker_env$filter_thresholds <- filter_thresholds
  worker_env$adj_method <- adj_method
  worker_env$feature_names <- feature_names
  worker_env$q_value <- q_value
  worker <- function(assignment) {
    calculate_null_summary(assignment, pooled, cooccur, filter_thresholds,
                           adj_method, feature_names, q_value)
  }
  environment(worker) <- worker_env
  worker
}

run_null_permutations <- function(assignments, pooled, cooccur,
                                  filter_thresholds, adj_method,
                                  feature_names, q_value, workers) {
  workers <- min(validate_null_workers(workers), length(assignments))
  if (workers == 1L) {
    return(lapply(assignments, calculate_null_summary, pooled = pooled,
                  cooccur = cooccur, filter_thresholds = filter_thresholds,
                  adj_method = adj_method, feature_names = feature_names,
                  q_value = q_value))
  }
  worker <- make_null_worker(pooled, cooccur, filter_thresholds, adj_method,
                             feature_names, q_value)
  cluster <- parallel::makePSOCKcluster(workers)
  on.exit(parallel::stopCluster(cluster), add = TRUE)
  parallel::clusterCall(cluster, function(paths) .libPaths(paths), .libPaths())
  parallel::parLapply(cluster, assignments, worker)
}
