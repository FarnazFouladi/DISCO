script <- local({
  paths <- unlist(lapply(sys.frames(), function(x) x$ofile), use.names = FALSE)
  paths <- c(rev(paths), sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)))
  paths <- paths[basename(paths) == "plot_ibd_hallmark_significance.R" & file.exists(paths)]
  if (!length(paths)) stop("Source the saved script or invoke it with Rscript.")
  normalizePath(paths[1])
})
root <- dirname(dirname(script))
# Adapted from the corresponding BRCA plotting script; calculations and plot conventions preserved.
# Run from the repository root; optionally supply a results directory.
# In RStudio: options(ibd.plot_input_dir = "/path/to/completed/results") before source().
# Priority: command-line directory, ibd.plot_input_dir, current results folder.
# Percentages use tested genes with finite within-pathway q-values only.
# Row clustering defaults to TRUE; set options(ibd.pathway_cluster_rows = FALSE)
# before sourcing for alphabetical rows. Columns always retain comparison order.
# options(ibd.pathway_min_percent = 50) keeps pathways with >=50% significant
# genes in at least ONE comparison (default 50).
# options(ibd.pathway_label_max_chars = 60L) truncates display names.
# Full names remain in CSVs; right-hand counts follow the comparison-column order.
cluster_rows <- getOption("ibd.pathway_cluster_rows", TRUE)
stopifnot(is.logical(cluster_rows), length(cluster_rows) == 1L, !is.na(cluster_rows))
if (!requireNamespace("ggplot2", quietly = TRUE)) stop("ggplot2 is required.")
source(file.path(dirname(script), "ibd_pathway_file_helpers.R"))
source(file.path(dirname(script), "ibd_heatmap_helpers.R"))
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
cutoff <- provenance$settings$q_value
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
  all(summary$n_significant == panels$n_significant_genes),
  all(is.na(summary$percent_significant) | summary$percent_significant >= 0 & summary$percent_significant <= 100))
prepared <- ibd_prepare_heatmap_rows(summary, pathway_lookup, provenance$settings$comparisons,
  min_percent, label_max_chars)
selected_pathways <- prepared$rows$pathway[prepared$rows$selected]
output <- file.path(input, "figure")
dir.create(output, recursive = TRUE, showWarnings = FALSE)
unlink(file.path(output, "pathway_row_clustering.rds"))
utils::write.csv(pathway_lookup, file.path(output, "pathway_lookup.csv"), row.names = FALSE)
utils::write.csv(summary, file.path(output, "pathway_significance_percentages.csv"), row.names = FALSE)
utils::write.csv(prepared$rows, file.path(output, "heatmap_pathway_selection.csv"), row.names = FALSE)
summary <- summary[summary$pathway %in% selected_pathways, , drop = FALSE]
row_order <- selected_pathways
clustering <- NULL
if (cluster_rows && length(row_order) > 1L) {
  percentages <- prepared$percentages
  if (any(!is.finite(percentages)))
    stop("Row clustering requires finite percentages for each selected pathway and comparison; disable clustering for missing values.")
  clustering <- stats::hclust(stats::dist(percentages, method = "euclidean"), method = "complete")
  row_order <- rownames(percentages)[clustering$order]
  saveRDS(clustering, file.path(output, "pathway_row_clustering.rds"))
}
annotation <- prepared$rows[match(row_order, prepared$rows$pathway), , drop = FALSE]
annotation$row_from_top <- seq_along(row_order)
utils::write.csv(annotation, file.path(output, "heatmap_row_order.csv"), row.names = FALSE)
comparisons <- provenance$settings$comparisons
n_comparisons <- length(comparisons)
summary$comparison_column <- match(summary$comparison, comparisons)
summary$pathway <- factor(summary$pathway, levels = rev(row_order))
annotation$pathway <- factor(annotation$pathway, levels = rev(row_order))
summary$label <- ifelse(is.na(summary$percent_significant), "NA", sprintf("%.1f", summary$percent_significant))
y_labels <- stats::setNames(annotation$pathway_label, as.character(annotation$pathway))
plot <- ggplot2::ggplot(summary, ggplot2::aes(comparison_column, pathway, fill = percent_significant)) +
  ggplot2::geom_tile(color = "white", linewidth = .4, width = 1, height = 1) +
  ggplot2::geom_text(ggplot2::aes(label = label), size = 3) +
  ggplot2::geom_text(data = annotation, ggplot2::aes(x = n_comparisons + 1.6, y = pathway, label = tested_genes),
    inherit.aes = FALSE, size = 3) +
  ggplot2::scale_x_continuous(breaks = c(seq_len(n_comparisons), n_comparisons + 1.6),
    labels = c(unname(comparison_labels[comparisons]), paste0("Tested genes\n(", paste(comparisons, collapse = ", "), ")")),
    limits = c(.5, n_comparisons + 3.1), expand = c(0, 0)) +
  ggplot2::scale_y_discrete(labels = y_labels, expand = ggplot2::expansion(add = .5)) +
  ggplot2::scale_fill_gradient(low = "#F7FBFF", high = "#4292C6", limits = c(0, 100), na.value = "grey85",
    name = "Significant genes (%)") +
  ggplot2::labs(title = paste("IBD", database_label, "pathway DISCO significance"),
    subtitle = paste0("Within-pathway BY q < ", cutoff, " | At least ", min_percent,
      "% significant genes in any comparison | ", length(row_order), " pathways"),
    x = "Disease group (versus non-IBD control)", y = NULL,
    caption = paste0("100 x significant genes / tested genes within each pathway and comparison. Untested genes excluded.\n",
      "Right-hand counts follow the comparison columns; full names and display labels are saved in heatmap_row_order.csv.\n",
      if (!is.null(clustering)) "Rows: hierarchical clustering, Euclidean distance, complete linkage; no scaling."
      else "Rows: alphabetical order.")) +
  ggplot2::theme_minimal(base_size = 11) +
  ggplot2::theme(panel.grid = ggplot2::element_blank(), axis.text.y = ggplot2::element_text(size = 9),
    legend.position = "bottom", plot.caption = ggplot2::element_text(hjust = 0), plot.title.position = "plot",
    plot.margin = ggplot2::margin(12, 15, 10, 12))
height <- max(8, .22 * length(row_order) + 3)
for (extension in c("png", "pdf")) ggplot2::ggsave(
  file.path(output, paste0("pathway_significance_heatmap.", extension)), plot,
  width = 16, height = height, units = "in", dpi = 150, bg = "white", limitsize = FALSE)
saveRDS(list(input = input, min_percent = min_percent, selection_rule = "at least one comparison",
  label_max_chars = label_max_chars, count_order = comparisons, cluster_rows = cluster_rows,
  selected_pathways = row_order, q_cutoff = cutoff), file.path(output, "plot_settings.rds"))
message("Input: ", input, "\nSelected pathways: ", length(row_order), "\nSaved plots and counts: ", output)
