script <- local({
  paths <- unlist(lapply(sys.frames(), function(x) x$ofile), use.names = FALSE)
  paths <- c(rev(paths), sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)))
  paths <- paths[basename(paths) == "plot_ibd_hallmark_significance.R" & file.exists(paths)]
  if (!length(paths)) stop("Source the saved script or invoke it with Rscript.")
  normalizePath(paths[1])
})
root <- dirname(dirname(script))
# Plot current pathway results as horizontal bars, faceted by comparison.
# Keep pathways with >=50% of tested genes having within-pathway BY q < 0.05
# in any comparison. All facets share the descending iCD pathway order.
# Optional options: ibd.plot_input_dir, ibd.pathway_min_percent (default 50),
# and ibd.pathway_label_max_chars (default 60).
if (!requireNamespace("ggplot2", quietly = TRUE)) stop("ggplot2 is required.")
source(file.path(dirname(script), "ibd_pathway_file_helpers.R"))
source(file.path(dirname(script), "ibd_pathway_plot_helpers.R"))
min_percent <- getOption("ibd.pathway_min_percent", 50)
label_max_chars <- getOption("ibd.pathway_label_max_chars", 60L)
args <- commandArgs(trailingOnly = TRUE)
if (length(args)) {
  input <- args[1]
} else if (!is.null(getOption("ibd.plot_input_dir", NULL))) {
  input <- getOption("ibd.plot_input_dir")
} else {
  input <- file.path(root, "results")
}
if (!is.character(input) || length(input) != 1L || is.na(input) || !dir.exists(input))
  stop("Provide an existing results directory via ibd.plot_input_dir or the first command-line argument.")
provenance <- readRDS(file.path(input, "provenance.rds"))
database_label <- provenance$settings$pathway_database
if (is.null(database_label)) database_label <- "Hallmark"
pathway_lookup <- ibd_read_pathway_lookup(input, provenance$pathways)
cutoff <- 0.05
scores <- read.csv(file.path(input, "all_pathway_gene_scores.csv"))
panels <- read.csv(file.path(input, "pathway_analysis_summary.csv"))
group_sizes <- unique(panels[c("comparison", "n_disease", "n_control")])
if (anyDuplicated(group_sizes$comparison))
  stop("Sample sizes differ between pathways within a comparison; cannot use a single column label.")
stopifnot(all(is.finite(group_sizes$n_disease)), all(is.finite(group_sizes$n_control)),
          all(group_sizes$n_disease > 0), all(group_sizes$n_control > 0))
comparison_labels <- stats::setNames(paste0("", group_sizes$comparison,
  "\nDisease n = ", group_sizes$n_disease, "\nNon-IBD n = ", group_sizes$n_control), group_sizes$comparison)
stopifnot(is.numeric(cutoff), length(cutoff) == 1, is.finite(cutoff), cutoff > 0, cutoff < 1,
          !anyDuplicated(panels[c("comparison", "pathway")]))
expected <- expand.grid(comparison = provenance$settings$comparisons, pathway = provenance$pathways)
key <- function(x) paste(x$comparison, x$pathway, sep = "::")
if (!setequal(key(panels), key(expected)) || any(panels$status != "completed"))
  stop("Expected all comparison-pathway comparisons to be completed successfully.")
# Aggregate once across all genes; preserve the existing within-pathway denominator.
tested <- scores$reason %in% "tested" & is.finite(scores$q_within_pathway)
counts <- rowsum(cbind(n_members = rep(1L, nrow(scores)), n_tested = as.integer(tested),
  n_significant = as.integer(tested & scores$q_within_pathway < cutoff)), key(scores), reorder = FALSE)
idx <- match(key(panels), rownames(counts))
stopifnot(!anyNA(idx), setequal(key(panels), rownames(counts)))
summary <- data.frame(comparison = panels$comparison, pathway = panels$pathway, q_cutoff = cutoff,
  n_pathway_members = counts[idx, "n_members"], n_tested = counts[idx, "n_tested"],
  n_untested = counts[idx, "n_members"] - counts[idx, "n_tested"],
  n_significant = counts[idx, "n_significant"], row.names = NULL)
summary$percent_significant <- with(summary, ifelse(n_tested > 0, 100 * n_significant / n_tested, NA_real_))
stopifnot(all(summary$n_tested == panels$n_tested_genes),
  all(is.na(summary$percent_significant) | summary$percent_significant >= 0 & summary$percent_significant <= 100))
prepared <- ibd_prepare_pathway_rows(summary, pathway_lookup, provenance$settings$comparisons,
  min_percent, label_max_chars)
row_order <- prepared$selected_pathways
output <- file.path(input, "figure")
dir.create(output, recursive = TRUE, showWarnings = FALSE)
utils::write.csv(pathway_lookup, file.path(output, "pathway_lookup.csv"), row.names = FALSE)
utils::write.csv(summary, file.path(output, "pathway_significance_percentages.csv"), row.names = FALSE)
utils::write.csv(prepared$rows, file.path(output, "pathway_selection.csv"), row.names = FALSE)
annotation <- prepared$rows[match(row_order, prepared$rows$pathway), , drop = FALSE]
annotation$row_from_top <- seq_along(row_order)
utils::write.csv(annotation, file.path(output, "pathway_order.csv"), row.names = FALSE)
comparisons <- provenance$settings$comparisons
summary <- summary[summary$pathway %in% row_order, , drop = FALSE]
summary$row_from_top <- match(summary$pathway, row_order)
summary <- summary[order(match(summary$comparison, comparisons), summary$row_from_top), ]
utils::write.csv(summary, file.path(output, "pathway_plot_data.csv"), row.names = FALSE)
summary$comparison <- factor(summary$comparison, levels = comparisons)
summary$pathway <- factor(summary$pathway, levels = rev(row_order))
summary$label <- ifelse(is.finite(summary$percent_significant),
  sprintf("%.1f%%", summary$percent_significant), "NA")
summary$label_x <- ifelse(is.finite(summary$percent_significant), summary$percent_significant, 0)
summary$label_hjust <- ifelse(summary$label_x >= 90, 1.1, -0.15)
y_labels <- stats::setNames(annotation$pathway_label, annotation$pathway)
plot <- ggplot2::ggplot(summary, ggplot2::aes(percent_significant, pathway, fill = comparison)) +
  ggplot2::geom_vline(xintercept = min_percent, linetype = "dashed", color = "grey65", linewidth = .4) +
  ggplot2::geom_col(data = summary[is.finite(summary$percent_significant), ], width = .75) +
  ggplot2::geom_text(ggplot2::aes(x = label_x, label = label, hjust = label_hjust),
    size = 2.8, color = "grey20") +
  ggplot2::facet_wrap(~comparison, nrow = 1,
    labeller = ggplot2::as_labeller(comparison_labels)) +
  ggplot2::scale_x_continuous(limits = c(0, 100), breaks = seq(0, 100, 25),
    expand = ggplot2::expansion(mult = c(0, .03))) +
  ggplot2::scale_y_discrete(labels = y_labels, drop = FALSE) +
  ggplot2::scale_fill_manual(values = c(iCD = "#2166AC", cCD = "#67A9CF",
    UC = "#B2182B", CD = "#4393C3")) +
  ggplot2::labs(title = paste("IBD", database_label, "pathway gene disruption"),
    subtitle = paste0("Within-pathway BY q < ", cutoff, " | At least ", min_percent,
      "% disrupted genes in any comparison | ", length(row_order), " pathways"),
    x = paste0("Genes with DISCO q < ", cutoff, " (%)"), y = NULL,
    caption = paste0("100 x disrupted genes / tested genes within each pathway and comparison; untested genes excluded.\n",
      "Pathways ordered by iCD percentage, highest to lowest; ties ordered alphabetically. Dashed line: ", min_percent, "%.")) +
  ggplot2::theme_minimal(base_size = 11) +
  ggplot2::theme(panel.grid.major.y = ggplot2::element_blank(),
    panel.grid.minor = ggplot2::element_blank(), axis.text.y = ggplot2::element_text(size = 9),
    strip.text = ggplot2::element_text(face = "bold"), legend.position = "none",
    plot.caption = ggplot2::element_text(hjust = 0), plot.title.position = "plot",
    panel.spacing = grid::unit(1, "lines"), plot.margin = ggplot2::margin(12, 15, 10, 12))
height <- max(6, .28 * length(row_order) + 2)
for (extension in c("png", "pdf")) {
  # Finish graphics locally before replacing a cloud-synced output file.
  temporary_plot <- tempfile(fileext = paste0(".", extension))
  ggplot2::ggsave(temporary_plot, plot, width = 16, height = height,
    units = "in", dpi = 150, bg = "white", limitsize = FALSE)
  destination <- file.path(output, paste0("pathway_significance_bars.", extension))
  if (!file.copy(temporary_plot, destination, overwrite = TRUE))
    stop("Could not replace plot: ", destination)
  unlink(temporary_plot)
}
saveRDS(list(input = input, min_percent = min_percent, selection_rule = "at least one comparison",
  label_max_chars = label_max_chars, count_order = comparisons, order_by = "iCD",
  selected_pathways = row_order, q_cutoff = cutoff, plot_type = "faceted_bars"),
  file.path(output, "plot_settings.rds"))
# Remove the replaced heatmap outputs so the current results contain one display.
unlink(file.path(output, c("pathway_significance_heatmap.png", "pathway_significance_heatmap.pdf",
  "heatmap_row_order.csv", "heatmap_pathway_selection.csv", "pathway_row_clustering.rds")))
message("Input: ", input, "\nSelected pathways: ", length(row_order), "\nSaved plots and counts: ", output)
