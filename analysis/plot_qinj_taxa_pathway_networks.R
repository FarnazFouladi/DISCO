# Bipartite taxon-pathway network containing the top N significant DISCO taxa
# and their significant differential-correlation partners, one panel per group.
# Run: Rscript analysis/plot_qinj_taxa_pathway_networks.R
# Optional arguments: input_directory output_directory q_cutoff top_n
# Also supports source("/path/to/DISCO/analysis/plot_qinj_taxa_pathway_networks.R").
# q_cutoff selects top DISCO taxa; the edge p cutoff is fixed at 0.05.
# Uses existing DISCO outputs only. Group 1 is control; group 2 is T2D.
# Both panels use the same undirected edges and deterministic node positions.
# Edges use raw differential-correlation p < 0.05, without a differential q filter.
# Selection uses the saved cor_diff table, which retains the DISCO run's
# selected single-group correlation filter; omitted pairs are not recomputed.
# Pathway labels use IDs; full names are retained in network_nodes.csv.
# Both panels include all qualifying edges to the selected taxa.

# Change this number to select how many significant DISCO taxa to highlight.
top_n <- 5L

script_file <- local({
  files <- as.character(unlist(lapply(sys.frames(), function(frame) frame$ofile),
                               use.names = FALSE))
  files <- files[basename(files) == "plot_qinj_taxa_pathway_networks.R"]
  if (length(files)) return(normalizePath(tail(files, 1L), mustWork = TRUE))
  args <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(args) == 1L)
    return(normalizePath(sub("^--file=", "", args[[1L]]), mustWork = TRUE))
  stop("Use Rscript or source() with the full path to this script.")
})
analysis_dir <- dirname(script_file)
sourced <- any(vapply(sys.frames(), function(frame) {
  !is.null(frame$ofile) && identical(basename(frame$ofile), "plot_qinj_taxa_pathway_networks.R")
}, logical(1)))
args <- if (sourced) character(0) else commandArgs(trailingOnly = TRUE)
if (length(args) > 4L) {
  stop("Usage: Rscript plot_qinj_taxa_pathway_networks.R [input_dir] [output_dir] [q_cutoff] [top_n]")
}
if (length(args) >= 4L) top_n <- suppressWarnings(as.numeric(args[[4L]]))
if (!is.numeric(top_n) || length(top_n) != 1L || !is.finite(top_n) ||
    top_n < 1 || top_n != floor(top_n) || top_n > .Machine$integer.max)
  stop("top_n must be a positive integer.")
top_n <- as.integer(top_n)
top_role <- paste("Top", top_n, "DISCO taxon")
input_dir <- if (length(args) >= 1L) args[[1L]] else file.path(
  analysis_dir, "qinj_readme_results", "taxa_pathways"
)
output_dir <- if (length(args) >= 2L) args[[2L]] else file.path(
  analysis_dir, "qinj_taxa_pathway_networks"
)
if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Install ggplot2 first.")
result <- readRDS(file.path(input_dir, "disco_analysis.rds"))
q_cutoff <- if (length(args) >= 3L) {
  suppressWarnings(as.numeric(args[[3L]]))
} else {
  readRDS(file.path(input_dir, "disco_inputs.rds"))$settings$q_value
}
if (!is.numeric(q_cutoff) || length(q_cutoff) != 1L || !is.finite(q_cutoff) ||
    q_cutoff <= 0 || q_cutoff >= 1) stop("q_cutoff must be between 0 and 1.")
if (!all(c("species", "p", "q") %in% names(result$disco)) ||
    !all(c("t1", "t2", "r1", "r2", "p", "q") %in% names(result$cor_diff))) {
  stop("Input must contain DISCO feature scores and a pairwise cor_diff table.")
}
# This script expects the saved QinJ taxon-pathway analysis: r1 is control,
# r2 is T2D, as defined by run_qinj_readme_examples.R.
scores <- result$disco
scores <- scores[is.finite(scores$p) & is.finite(scores$q) & scores$q < q_cutoff, ]
scores <- scores[order(scores$q, scores$p, scores$species), , drop = FALSE]
if (nrow(scores) < top_n) stop("Requested ", top_n, " taxa, but only ", nrow(scores),
  " have DISCO q < ", q_cutoff, ". Reduce top_n or change q_cutoff.")
top <- head(scores, top_n)
top$rank <- seq_len(nrow(top))
top$panel <- paste0(top$rank, ". ", top$species, "\nDISCO q = ",
                    formatC(top$q, format = "g", digits = 3))

pairs <- result$cor_diff
edge_p_cutoff <- 0.05
pairs <- pairs[is.finite(pairs$p) & pairs$p < edge_p_cutoff &
                 is.finite(pairs$r1) & is.finite(pairs$r2), , drop = FALSE]
inputs <- readRDS(file.path(input_dir, "disco_inputs.rds"))
if (length(inputs$data1) != 2L || length(inputs$data2) != 2L)
  stop("Expected a two-modality taxa/pathway analysis.")
if (!all(result$cor_diff$t1 %in% colnames(inputs$data1[[1L]])) ||
    !all(result$cor_diff$t2 %in% colnames(inputs$data1[[2L]])))
  stop("Expected t1 to identify taxa and t2 to identify pathways.")
edges <- pairs[pairs$t1 %in% top$species, , drop = FALSE]
if (!nrow(edges)) stop("No top-taxon pairs pass differential-correlation p < 0.05.")
# Each row is one taxon-pathway pair; use typed node IDs to avoid collisions.
key <- paste(edges$t1, edges$t2, sep = "\r")
edges <- edges[!duplicated(key), , drop = FALSE]
pathways <- unique(edges$t2)
# Arrange pathways by their first connected top taxon, then their full name.
first_taxon <- vapply(pathways, function(p) min(match(edges$t1[edges$t2 == p], top$species)), integer(1))
pathways <- pathways[order(first_taxon, pathways)]
spread <- function(n) if (n == 1L) 0 else seq(1, -1, length.out = n)
taxon_nodes <- data.frame(node_id = paste0("taxon::", top$species),
  feature = top$species, modality = "Taxon", role = top_role,
  x = -0.65, y = spread(nrow(top)) * 0.9,
  label = vapply(top$species, function(x) paste(strwrap(x, width = 24), collapse = "\n"), ""))
pathway_nodes <- data.frame(node_id = paste0("pathway::", pathways),
  feature = pathways, modality = "Pathway", role = "Pathway",
  x = 0.65, y = spread(length(pathways)), label = sub(":.*$", "", pathways))
nodes <- rbind(taxon_nodes, pathway_nodes)
edges$from <- paste0("taxon::", edges$t1)
edges$to <- paste0("pathway::", edges$t2)
start <- match(edges$from, nodes$node_id)
end <- match(edges$to, nodes$node_id)
edges$x <- nodes$x[start]; edges$y <- nodes$y[start]
edges$xend <- nodes$x[end]; edges$yend <- nodes$y[end]
plot_edges <- rbind(transform(edges, group = "Control", r = r1),
                    transform(edges, group = "T2D", r = r2))
plot_nodes <- rbind(transform(nodes, group = "Control"), transform(nodes, group = "T2D"))
plot_edges$group <- factor(plot_edges$group, levels = c("Control", "T2D"))
plot_nodes$group <- factor(plot_nodes$group, levels = c("Control", "T2D"))
plot_edges$abs_r <- abs(plot_edges$r)
plot <- ggplot2::ggplot() +
  ggplot2::geom_segment(data = plot_edges,
    ggplot2::aes(x = x, y = y, xend = xend, yend = yend, color = r, linewidth = abs_r), alpha = .55) +
  ggplot2::geom_point(data = plot_nodes,
    ggplot2::aes(x = x, y = y, fill = role, size = role),
    shape = 21, color = "#334155", stroke = .7) +
  ggplot2::geom_text(data = subset(plot_nodes, modality == "Taxon"),
    ggplot2::aes(x = x - .06, y = y, label = label),
    hjust = 1, fontface = "bold", size = 3.1, lineheight = .95) +
  ggplot2::geom_text(data = subset(plot_nodes, modality == "Pathway"),
    ggplot2::aes(x = x + .06, y = y, label = label), hjust = 0, size = 2.6) +
  ggplot2::facet_wrap(~ group, nrow = 1) +
  ggplot2::scale_fill_manual(name = NULL,
    values = c(stats::setNames("#FBBF24", top_role), "Pathway" = "#CBD5E1")) +
  ggplot2::scale_size_manual(values = c(stats::setNames(6, top_role), "Pathway" = 2.8), guide = "none") +
  ggplot2::scale_color_gradient2(name = "Spearman r", low = "#2563EB", mid = "#E2E8F0",
    high = "#DC2626", midpoint = 0, limits = c(-1, 1), breaks = c(-1, -.5, 0, .5, 1)) +
  ggplot2::scale_linewidth_continuous(name = "|Spearman r|", range = c(.2, 1.8),
    limits = c(0, 1), breaks = c(.25, .5, .75)) +
  ggplot2::coord_cartesian(xlim = c(-1.9, 1.9), ylim = c(-1.07, 1.07), expand = FALSE) +
  ggplot2::labs(title = paste("Taxon-pathway networks of the top", top_n, "DISCO taxa"),
    subtitle = paste0("Top taxa: DISCO q < ", q_cutoff,
      " | Edges: differential-correlation p < ", edge_p_cutoff,
      " | ", nrow(top), " taxa, ", length(pathways), " pathways, ", nrow(edges), " pairs per group"),
    caption = paste0(
      "Gold nodes: top taxa ranked by cross-modality DISCO q, then p. Grey nodes: pathways (IDs shown; full names in network_nodes.csv).\n",
      "Red: positive correlation; blue: negative correlation. Width: |r|. Identical nodes, edges and positions in both groups.\n",
      "Inputs: robust CLR (zeros left NA), adjusted separately within each group for gender and BMI.\n",
      "Raw edge p-values are unadjusted and do not imply within-group significance. Selected DISCO filter: ", result$selected_threshold, ".")) +
  ggplot2::theme_void(base_size = 12) +
  ggplot2::theme(strip.text = ggplot2::element_text(size = 15, face = "bold"),
    strip.background = ggplot2::element_rect(fill = "#F1F5F9", color = NA),
    panel.spacing = grid::unit(1, "lines"), legend.position = "bottom",
    plot.title = ggplot2::element_text(size = 20, face = "bold"),
    plot.subtitle = ggplot2::element_text(size = 10, margin = ggplot2::margin(b = 12)),
    plot.caption = ggplot2::element_text(hjust = 0, size = 9, margin = ggplot2::margin(t = 14)),
    plot.margin = ggplot2::margin(18, 18, 18, 18))
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
base_path <- file.path(output_dir, paste0("top", top_n, "_taxa_pathways_control_T2D_networks"))
# Scale height with the number of pathway labels; no qualifying edges are dropped.
plot_height <- max(11, 4 + .18 * length(pathways))
ggplot2::ggsave(paste0(base_path, ".png"), plot, width = 18, height = plot_height, dpi = 180, bg = "white")
ggplot2::ggsave(paste0(base_path, ".pdf"), plot, width = 18, height = plot_height, bg = "white")
utils::write.csv(nodes, file.path(output_dir, "network_nodes.csv"), row.names = FALSE)
utils::write.csv(edges, file.path(output_dir, "network_edges.csv"), row.names = FALSE)
utils::write.csv(top[, setdiff(names(top), "panel")], file.path(output_dir, "selected_taxa.csv"), row.names = FALSE)
message("Saved paired taxon-pathway networks in ", normalizePath(output_dir))
