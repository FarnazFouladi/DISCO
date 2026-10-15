# Match by feature name; absent or untested features are not discoveries.
null_feature_indicators <- function(result, feature_names, q_value) {
  scores <- result$disco
  index <- match(feature_names, scores$species)
  tested <- is.finite(scores$p[index]) & is.finite(scores$q[index])
  list(tested = tested, discovered = tested & scores$q[index] < q_value)
}

# Null frequencies use all requested permutations as the denominator.
build_feature_null_summary <- function(real_results, null_discoveries,
                                       null_tested, filter_thresholds,
                                       q_value, n_resamples) {
  rows <- lapply(seq_along(filter_thresholds), function(i) {
    scores <- real_results[[i]]$disco
    scores <- scores[
      is.finite(scores$p) & scores$p >= 0 & scores$p <= 1 &
        is.finite(scores$q), , drop = FALSE
    ]
    index <- match(scores$species, rownames(null_discoveries))
    data.frame(
      filter_threshold = rep(filter_thresholds[[i]], nrow(scores)),
      species = scores$species,
      real_disco_p = scores$p,
      real_disco_q = scores$q,
      real_significant = scores$q < q_value,
      neg_log10_real_disco_p = -log10(pmax(scores$p, 1e-8)),
      p_value_floored_for_plot = scores$p < 1e-8,
      n_null_tested = null_tested[index, i],
      n_null_discoveries = null_discoveries[index, i],
      n_resamples = rep(n_resamples, nrow(scores)),
      null_discovery_percentage =
        100 * null_discoveries[index, i] / n_resamples,
      stringsAsFactors = FALSE
    )
  })
  summary <- do.call(rbind, rows)
  rownames(summary) <- NULL
  summary
}

# Plot the same q cutoff for real and null results, one panel per threshold.
create_feature_null_plot <- function(summary, filter_thresholds, q_value,
                                      n_resamples) {
  labels <- paste0("filter_threshold = ", format(
    filter_thresholds, scientific = FALSE, trim = TRUE, drop0trailing = TRUE
  ))
  summary$panel <- factor(
    labels[match(summary$filter_threshold, filter_thresholds)], levels = labels
  )
  significant_label <- paste0("Real DISCO q < ", format(q_value))
  other_label <- paste0("Real DISCO q >= ", format(q_value))
  summary$real_q_category <- factor(
    ifelse(summary$real_significant, significant_label, other_label),
    levels = c(significant_label, other_label)
  )
  annotations <- data.frame(
    panel = factor(labels, levels = labels),
    annotation = paste0("Real q < ", format(q_value), ": n = ", vapply(
      filter_thresholds,
      function(threshold) {
        sum(summary$filter_threshold == threshold & summary$real_significant)
      }, integer(1L)
    ))
  )

  ggplot2::ggplot(summary, ggplot2::aes(
    x = neg_log10_real_disco_p, y = null_discovery_percentage,
    color = real_q_category
  )) +
    ggplot2::geom_point(size = 1.8, alpha = 0.72) +
    ggplot2::geom_label(
      data = annotations,
      ggplot2::aes(x = Inf, y = Inf, label = annotation),
      inherit.aes = FALSE, hjust = 1.05, vjust = 1.2,
      color = "#D55E00", fill = "white", linewidth = 0.25, size = 3
    ) +
    ggplot2::facet_wrap(ggplot2::vars(panel), ncol = 4, drop = FALSE) +
    ggplot2::scale_x_continuous(limits = c(0, 8), breaks = seq(0, 8, 2)) +
    ggplot2::scale_y_continuous(
      limits = c(0, 100), breaks = seq(0, 100, 20),
      labels = function(x) paste0(x, "%")
    ) +
    ggplot2::scale_color_manual(
      values = stats::setNames(c("#D55E00", "#999999"),
                               c(significant_label, other_label)),
      limits = c(significant_label, other_label),
      drop = FALSE
    ) +
    ggplot2::labs(
      title = "Real DISCO p-values versus feature null discovery frequency",
      subtitle = paste0("Null percentages use all ", n_resamples,
                        " size-matched label permutations"),
      x = expression(-log[10]("Real-data DISCO p-value")),
      y = paste0("Null permutations with DISCO q < ", format(q_value), " (%)"),
      color = NULL,
      caption = paste0(
        "Each point is one modality-1 feature; orange counts are shown per panel.\n",
        "Untested features count as non-discoveries in a permutation. ",
        "This frequency is not a formal feature-level FDR estimate.\n",
        "Real p-values below 1e-8 are plotted at 8; original p-values are retained in the data."
      )
    ) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(face = "bold"),
      plot.title = ggplot2::element_text(face = "bold"),
      plot.title.position = "plot",
      plot.caption = ggplot2::element_text(hjust = 0),
      legend.position = "top"
    )
}
