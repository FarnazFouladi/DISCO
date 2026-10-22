# Validate empirical threshold-selection settings
validate_threshold_selection <- function(filter_thresholds, q_value,
                                         n_resamples, random_seed,
                                         global_p_cutoff = 0.05) {
  if (!is.numeric(filter_thresholds) || length(filter_thresholds) == 0L ||
      anyNA(filter_thresholds) || any(!is.finite(filter_thresholds)) ||
      any(filter_thresholds <= 0 | filter_thresholds > 1)) {
    stop(
      "filter_thresholds must contain finite numeric values greater than 0 and less than or equal to 1",
      call. = FALSE
    )
  }

  if (!is.numeric(q_value) || length(q_value) != 1L || is.na(q_value) ||
      !is.finite(q_value) || q_value <= 0 || q_value > 1) {
    stop("q_value must be a single number greater than 0 and less than or equal to 1",
         call. = FALSE)
  }

  if (!is.numeric(n_resamples) || length(n_resamples) != 1L ||
      is.na(n_resamples) || !is.finite(n_resamples) || n_resamples < 1 ||
      n_resamples != floor(n_resamples)) {
    stop("n_resamples must be a positive integer", call. = FALSE)
  }

  if (!is.numeric(random_seed) || length(random_seed) != 1L ||
      is.na(random_seed) || !is.finite(random_seed) || random_seed < 0 ||
      random_seed > .Machine$integer.max || random_seed != floor(random_seed)) {
    stop("random_seed must be a single non-negative integer", call. = FALSE)
  }

  if (!is.numeric(global_p_cutoff) || length(global_p_cutoff) != 1L ||
      is.na(global_p_cutoff) || !is.finite(global_p_cutoff) ||
      global_p_cutoff <= 0 || global_p_cutoff > 1) {
    stop("global_p_cutoff must be a single number greater than 0 and less than or equal to 1",
         call. = FALSE)
  }

  list(
    filter_thresholds = sort(unique(filter_thresholds)),
    q_value = q_value,
    n_resamples = as.integer(n_resamples),
    random_seed = as.integer(random_seed),
    global_p_cutoff = global_p_cutoff
  )
}

# Pool the observed groups while preserving sample correspondence by modality
pool_null_inputs <- function(prepared) {
  groups <- list(prepared$data1, prepared$data2)
  group_records <- vector("list", length(groups))

  for (group_index in seq_along(groups)) {
    group <- groups[[group_index]]
    modality_ids <- lapply(seq_along(group), function(modality_index) {
      ids <- rownames(group[[modality_index]])
      if (is.null(ids)) {
        ids <- paste0("sample_", seq_len(nrow(group[[modality_index]])))
      }
      if (anyNA(ids) || any(!nzchar(ids)) || anyDuplicated(ids)) {
        stop(
          "Each modality must have unique, non-missing sample row names for empirical threshold selection",
          call. = FALSE
        )
      }
      ids
    })

    all_ids <- unique(unlist(modality_ids, use.names = FALSE))
    internal_ids <- paste0(
      ".disco_group", group_index, "_sample", seq_along(all_ids)
    )
    id_map <- stats::setNames(internal_ids, all_ids)

    for (modality_index in seq_along(group)) {
      rownames(group[[modality_index]]) <- unname(
        id_map[modality_ids[[modality_index]]]
      )
    }
    groups[[group_index]] <- group

    availability <- vapply(modality_ids, function(ids) {
      all_ids %in% ids
    }, logical(length(all_ids)))
    if (is.null(dim(availability))) {
      availability <- matrix(availability, ncol = 1L)
    }
    pattern <- apply(availability, 1L, function(x) {
      paste(as.integer(x), collapse = "")
    })
    group_records[[group_index]] <- data.frame(
      sample_id = internal_ids,
      original_group = group_index,
      availability_pattern = pattern,
      stringsAsFactors = FALSE
    )
  }

  records <- do.call(rbind, group_records)
  rownames(records) <- NULL
  strata <- split(records, records$availability_pattern, drop = TRUE)
  pooled <- lapply(seq_along(groups[[1]]), function(modality_index) {
    rbind(groups[[1]][[modality_index]], groups[[2]][[modality_index]])
  })

  list(pooled = pooled, strata = strata)
}

# Permute labels within modality-availability strata
draw_null_assignment <- function(pooled) {
  assignments <- lapply(pooled$strata, function(stratum) {
    shuffled <- sample(stratum$sample_id, nrow(stratum), replace = FALSE)
    n_data1 <- sum(stratum$original_group == 1L)
    list(
      data1 = if (n_data1 > 0L) {
        shuffled[seq_len(n_data1)]
      } else {
        character()
      },
      data2 = if (n_data1 < length(shuffled)) {
        shuffled[n_data1 + seq_len(length(shuffled) - n_data1)]
      } else {
        character()
      }
    )
  })
  data1_ids <- unlist(lapply(assignments, `[[`, "data1"), use.names = FALSE)
  data2_ids <- unlist(lapply(assignments, `[[`, "data2"), use.names = FALSE)

  list(data1 = data1_ids, data2 = data2_ids)
}

permute_null_groups <- function(pooled, assignment = draw_null_assignment(pooled)) {

  subset_modalities <- function(ids) {
    lapply(pooled$pooled, function(x) {
      available_ids <- ids[ids %in% rownames(x)]
      x[available_ids, , drop = FALSE]
    })
  }

  list(
    data1 = subset_modalities(assignment$data1),
    data2 = subset_modalities(assignment$data2)
  )
}

# Percentage of tested modality-1 features declared significant
disco_discovery_percentage <- function(result, q_value) {
  scores <- result$disco
  tested <- !is.na(scores$p) & !is.na(scores$q)
  n_tested <- sum(tested)
  n_significant <- sum(tested & scores$q < q_value)

  data.frame(
    n_tested = n_tested,
    n_significant = n_significant,
    percentage = if (n_tested > 0L) {
      100 * n_significant / n_tested
    } else {
      NA_real_
    }
  )
}

mean_or_na <- function(x) {
  if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)
}

sd_or_na <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) < 2L) NA_real_ else stats::sd(x)
}

# Build the observed-versus-null summary used for threshold selection
build_threshold_summary <- function(real_results, null_percentages,
                                    filter_thresholds, q_value) {
  real_summaries <- lapply(real_results, disco_discovery_percentage,
                           q_value = q_value)
  real_summary <- do.call(rbind, real_summaries)

  null_mean <- apply(null_percentages, 2L, mean_or_na)
  null_sd <- apply(null_percentages, 2L, sd_or_na)
  null_floor <- pmax(null_mean, 1)
  real_floor <- pmax(real_summary$percentage, 1)
  ratio <- real_floor / null_floor

  data.frame(
    filter_threshold = filter_thresholds,
    real_n_tested = real_summary$n_tested,
    real_n_significant = real_summary$n_significant,
    real_percentage = real_summary$percentage,
    null_mean_percentage = null_mean,
    null_sd_percentage = null_sd,
    null_percentage_floor_1 = null_floor,
    real_percentage_floor_1 = real_floor,
    real_to_null_ratio = ratio,
    stringsAsFactors = FALSE
  )
}

# Compare the observed maximum with each permutation's maximum across the
# same evaluable thresholds, preserving the pairing between thresholds.
test_global_discoveries <- function(summary, null_percentages,
                                    global_p_cutoff) {
  evaluable <- is.finite(summary$real_percentage)
  if (!any(evaluable)) {
    stop("No filtering threshold produced a real DISCO discovery percentage",
         call. = FALSE)
  }
  null_values <- null_percentages[, evaluable, drop = FALSE]
  complete <- apply(null_values, 1L, function(x) all(is.finite(x)))
  if (!any(complete)) {
    stop("No null permutation produced discovery percentages at every evaluable threshold",
         call. = FALSE)
  }
  real_max <- max(summary$real_percentage[evaluable])
  null_max <- apply(null_values[complete, , drop = FALSE], 1L, max)
  global_p_value <- (1 + sum(null_max >= real_max)) / (length(null_max) + 1)

  data.frame(
    real_max_percentage = real_max,
    null_mean_max_percentage = mean(null_max),
    null_sd_max_percentage = sd_or_na(null_max),
    n_thresholds = sum(evaluable),
    n_valid_resamples = length(null_max),
    global_p_value = global_p_value,
    global_p_cutoff = global_p_cutoff,
    signal_detected = global_p_value < global_p_cutoff
  )
}

# Select the largest threshold among those attaining the maximum ratio
select_optimal_threshold <- function(summary, global_test = NULL) {
  if (!is.null(global_test) && !global_test$signal_detected) {
    return(1)
  }
  valid <- is.finite(summary$real_to_null_ratio)
  if (!any(valid)) {
    stop(
      "No filtering threshold produced both real and null DISCO discovery percentages",
      call. = FALSE
    )
  }

  maximum_ratio <- max(summary$real_to_null_ratio[valid])
  candidates <- summary$filter_threshold[
    valid & summary$real_to_null_ratio == maximum_ratio
  ]
  max(candidates)
}

# Create threshold-selection diagnostic plots
create_threshold_plots <- function(summary, q_value, n_resamples,
                                   selected_threshold, global_test) {
  selection_label <- paste0(
    "Global p = ", format.pval(global_test$global_p_value, digits = 3),
    "; cutoff = ", format(global_test$global_p_cutoff),
    "; ", if (global_test$signal_detected) "maximum-ratio selection" else
      "fallback to threshold 1"
  )
  thresholds <- summary$filter_threshold
  threshold_labels <- format(
    thresholds,
    scientific = FALSE,
    trim = TRUE,
    drop0trailing = TRUE
  )
  comparison_data <- rbind(
    data.frame(
      filter_threshold = thresholds,
      percentage = summary$real_percentage,
      lower = NA_real_,
      upper = NA_real_,
      result_type = "Real data"
    ),
    data.frame(
      filter_threshold = thresholds,
      percentage = summary$null_mean_percentage,
      lower = pmax(0, summary$null_mean_percentage - summary$null_sd_percentage),
      upper = pmin(100, summary$null_mean_percentage + summary$null_sd_percentage),
      result_type = "Null mean"
    )
  )
  null_data <- comparison_data[
    comparison_data$result_type == "Null mean",
    ,
    drop = FALSE
  ]

  discovery_plot <- ggplot2::ggplot(
    comparison_data,
    ggplot2::aes(
      x = filter_threshold,
      y = percentage,
      color = result_type,
      linetype = result_type,
      group = result_type
    )
  ) +
    ggplot2::geom_errorbar(
      data = null_data,
      ggplot2::aes(ymin = lower, ymax = upper),
      width = 0,
      alpha = 0.6,
      show.legend = FALSE
    )
  if (length(thresholds) > 1L) {
    discovery_plot <- discovery_plot +
      ggplot2::geom_line(linewidth = 1, na.rm = TRUE)
  }
  discovery_plot <- discovery_plot +
    ggplot2::geom_point(size = 2.8, na.rm = TRUE) +
    ggplot2::geom_vline(
      xintercept = selected_threshold,
      color = "#666666",
      linetype = "dotted"
    ) +
    ggplot2::scale_x_log10(breaks = thresholds, labels = threshold_labels) +
    ggplot2::scale_color_manual(values = c(
      "Real data" = "#0072B2",
      "Null mean" = "#D55E00"
    )) +
    ggplot2::labs(
      title = "Real and null DISCO discoveries across filtering thresholds",
      subtitle = paste0(
        "Null mean +/- SD across ", n_resamples,
        " size-matched label permutations\n", selection_label
      ),
      x = "Correlation filtering threshold (log10 scale)",
      y = paste0("Features with DISCO q < ", format(q_value), " (%)"),
      color = NULL,
      linetype = NULL
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      legend.position = "top"
    )

  ratio_plot <- ggplot2::ggplot(
    summary,
    ggplot2::aes(
      x = filter_threshold,
      y = real_to_null_ratio,
      group = 1L
    )
  ) +
    ggplot2::geom_hline(
      yintercept = 1,
      color = "#666666",
      linetype = "dashed"
    )
  if (length(thresholds) > 1L) {
    ratio_plot <- ratio_plot +
      ggplot2::geom_line(color = "#6A3D9A", linewidth = 1, na.rm = TRUE)
  }
  ratio_plot <- ratio_plot +
    ggplot2::geom_point(color = "#6A3D9A", size = 2.8, na.rm = TRUE) +
    ggplot2::geom_vline(
      xintercept = selected_threshold,
      color = "#E69F00",
      linetype = "dotted",
      linewidth = 0.9
    ) +
    ggplot2::scale_x_log10(breaks = thresholds, labels = threshold_labels) +
    ggplot2::labs(
      title = "Real-to-null DISCO discovery ratio",
      subtitle = paste0(
        "Selected threshold = ", format(selected_threshold),
        "; both percentages floored at 1\n", selection_label
      ),
      x = "Correlation filtering threshold (log10 scale)",
      y = "max(real percentage, 1) / max(mean null percentage, 1)"
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
    )

  list(discovery_plot = discovery_plot, ratio_plot = ratio_plot)
}

# Run empirical threshold selection and write the selected results
run_empirical_disco <- function(data1, data2, cooccur, filter_thresholds,
                                q_value, n_resamples, random_seed, adj_method,
                                output_dir, prefix, global_p_cutoff = 0.05,
                                workers = 1L) {
  workers <- validate_null_workers(workers)
  cooccur <- validate_cooccur(cooccur)
  settings <- validate_threshold_selection(
    filter_thresholds = filter_thresholds,
    q_value = q_value,
    n_resamples = n_resamples,
    random_seed = random_seed,
    global_p_cutoff = global_p_cutoff
  )
  filter_thresholds <- settings$filter_thresholds
  q_value <- settings$q_value
  n_resamples <- settings$n_resamples
  random_seed <- settings$random_seed

  prepared <- prepare_disco_inputs(data1, data2, report_overlap = TRUE)
  pooled <- pool_null_inputs(prepared)

  real_base <- calculate_disco_base(prepared, cooccur = cooccur)
  real_results <- lapply(filter_thresholds, function(threshold) {
    calculate_disco_at_threshold(
      real_base,
      filt_threshold = threshold,
      adj_method = adj_method,
      include_cor_diff = FALSE
    )
  })

  null_percentages <- matrix(
    NA_real_,
    nrow = n_resamples,
    ncol = length(filter_thresholds),
    dimnames = list(NULL, format(filter_thresholds, scientific = FALSE))
  )
  feature_names <- colnames(prepared$data1[[1L]])
  null_discoveries <- matrix(
    0L, nrow = length(feature_names), ncol = length(filter_thresholds),
    dimnames = list(feature_names, colnames(null_percentages))
  )
  null_tested <- null_discoveries

  seed_existed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (seed_existed) {
    previous_seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  }
  on.exit({
    if (seed_existed) {
      assign(".Random.seed", previous_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)
  set.seed(random_seed)

  message(
    "Estimating null discoveries using ", n_resamples,
    " size-matched label permutations with ", min(workers, n_resamples),
    " worker(s)."
  )
  # Draw only labels here, not copies of the expression matrices. This retains
  # the serial RNG sequence independently of scheduling and worker count.
  assignments <- lapply(seq_len(n_resamples), function(i) draw_null_assignment(pooled))
  null_summaries <- run_null_permutations(
    assignments, pooled, cooccur, filter_thresholds, adj_method,
    feature_names, q_value, workers
  )
  for (resample_index in seq_len(n_resamples)) {
    summary <- null_summaries[[resample_index]]
    null_percentages[resample_index, ] <- summary$percentages
    null_discoveries <- null_discoveries + summary$discovered
    null_tested <- null_tested + summary$tested
  }

  threshold_summary <- build_threshold_summary(
    real_results = real_results,
    null_percentages = null_percentages,
    filter_thresholds = filter_thresholds,
    q_value = q_value
  )
  global_test <- test_global_discoveries(
    threshold_summary, null_percentages, settings$global_p_cutoff
  )
  selected_threshold <- select_optimal_threshold(threshold_summary, global_test)
  global_test$selected_threshold <- selected_threshold
  threshold_summary$selected <-
    threshold_summary$filter_threshold == selected_threshold
  selected_index <- match(selected_threshold, filter_thresholds)
  selected_result <- real_results[[selected_index]]
  # Only the selected real-data threshold needs a pairwise output table.
  selected_pvals <- filter_correlation_tests(
    real_base$correlation_tests, filt_threshold = selected_threshold
  )$pvals
  selected_result$cor_diff <- build_cor_diff_table(
    real_base, selected_pvals, adj_method
  )

  plots <- create_threshold_plots(
    summary = threshold_summary,
    q_value = q_value,
    n_resamples = n_resamples,
    selected_threshold = selected_threshold,
    global_test = global_test
  )
  feature_null_summary <- build_feature_null_summary(
    real_results, null_discoveries, null_tested, filter_thresholds,
    q_value, n_resamples
  )
  feature_null_plot <- create_feature_null_plot(
    feature_null_summary, filter_thresholds, q_value, n_resamples
  )

  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  output_files <- c(
    cor_diff = file.path(output_dir, paste0(prefix, "_cor_diff.csv")),
    disco = file.path(output_dir, paste0(prefix, "_disco.csv")),
    threshold_summary = file.path(
      output_dir,
      paste0(prefix, "_threshold_summary.csv")
    ),
    global_test = file.path(output_dir, paste0(prefix, "_global_test.csv")),
    discovery_plot = file.path(
      output_dir,
      paste0(prefix, "_real_vs_null_thresholds.png")
    ),
    ratio_plot = file.path(
      output_dir,
      paste0(prefix, "_real_to_null_ratio.png")
    ),
    feature_null_summary = file.path(
      output_dir, paste0(prefix, "_real_p_vs_null_discovery_by_threshold.csv")
    ),
    feature_null_plot = file.path(
      output_dir, paste0(prefix, "_real_p_vs_null_discovery_by_threshold.png")
    )
  )

  utils::write.csv(
    selected_result$cor_diff,
    output_files[["cor_diff"]],
    quote = TRUE,
    row.names = FALSE
  )
  utils::write.csv(
    selected_result$disco,
    output_files[["disco"]],
    quote = TRUE,
    row.names = FALSE
  )
  utils::write.csv(
    threshold_summary,
    output_files[["threshold_summary"]],
    quote = TRUE,
    row.names = FALSE
  )
  utils::write.csv(global_test, output_files[["global_test"]], row.names = FALSE)
  utils::write.csv(feature_null_summary, output_files[["feature_null_summary"]],
                   row.names = FALSE)
  ggplot2::ggsave(
    output_files[["discovery_plot"]],
    plots$discovery_plot,
    width = 8.5,
    height = 6.5,
    units = "in",
    dpi = 300
  )
  ggplot2::ggsave(
    output_files[["ratio_plot"]],
    plots$ratio_plot,
    width = 8.5,
    height = 6.5,
    units = "in",
    dpi = 300
  )

  ggplot2::ggsave(
    output_files[["feature_null_plot"]], feature_null_plot,
    width = 13, height = 8.5, units = "in", dpi = 300
  )

  message("Global permutation p-value: ", format(global_test$global_p_value),
          if (global_test$signal_detected) "; signal detected." else
            "; using fallback threshold 1.")
  message("Selected correlation filtering threshold: ", selected_threshold)

  list(
    cor_diff = selected_result$cor_diff,
    disco = selected_result$disco,
    selected_threshold = selected_threshold,
    threshold_summary = threshold_summary,
    global_test = global_test,
    discovery_plot = plots$discovery_plot,
    ratio_plot = plots$ratio_plot,
    feature_null_summary = feature_null_summary,
    feature_null_plot = feature_null_plot,
    output_files = output_files
  )
}
