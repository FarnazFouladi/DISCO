# Feature-level DISCO heatmap: union of taxa with q < q_cutoff in either mode.
# Rows: taxa; columns: taxa-taxa and taxon-pathway analyses.
# Colors: -log10(raw feature-level DISCO p), with no row scaling.
# Asterisks: feature-level DISCO q < q_cutoff (not raw p significance).
# Rscript analysis/plot_qinj_disco_heatmap.R [input_dir] [output_dir] [cluster_rows] [row_clusters]
# Also supports source("/path/to/DISCO/analysis/plot_qinj_disco_heatmap.R").

q_cutoff <- 0.2
cluster_rows <- TRUE
row_clusters <- NULL # Optional number of row groups, e.g. 3L; requires clustering.
row_clustering_method <- "complete"

script_file <- local({
  files <- as.character(unlist(lapply(sys.frames(), function(frame) frame$ofile), use.names = FALSE))
  files <- files[basename(files) == "plot_qinj_disco_heatmap.R"]
  if (length(files)) return(normalizePath(tail(files, 1L), mustWork = TRUE))
  args <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(args) == 1L)
    return(normalizePath(sub("^--file=", "", args[[1L]]), mustWork = TRUE))
  stop("Use Rscript or source() with the full path to this script.")
})
analysis_dir <- dirname(script_file)
sourced <- any(vapply(sys.frames(), function(frame) {
  !is.null(frame$ofile) && identical(basename(frame$ofile), "plot_qinj_disco_heatmap.R")
}, logical(1)))
args <- if (sourced) character(0) else commandArgs(trailingOnly = TRUE)
if (length(args) > 4L)
  stop("Usage: Rscript plot_qinj_disco_heatmap.R [input_dir] [output_dir] [cluster_rows] [row_clusters]")
input_dir <- if (length(args) >= 1L) args[[1L]] else file.path(analysis_dir, "qinj_readme_results")
output_dir <- if (length(args) >= 2L) args[[2L]] else file.path(analysis_dir, "qinj_disco_heatmap")
if (length(args) >= 3L) {
  flag <- tolower(args[[3L]])
  if (!flag %in% c("true", "false")) stop("cluster_rows must be TRUE or FALSE.")
  cluster_rows <- flag == "true"
}
if (length(args) >= 4L) row_clusters <- suppressWarnings(as.numeric(args[[4L]]))
if (!is.numeric(q_cutoff) || length(q_cutoff) != 1L || !is.finite(q_cutoff) ||
    q_cutoff <= 0 || q_cutoff >= 1) stop("q_cutoff must be between 0 and 1.")
if (!is.logical(cluster_rows) || length(cluster_rows) != 1L || is.na(cluster_rows))
  stop("cluster_rows must be TRUE or FALSE.")
if (!requireNamespace("pheatmap", quietly = TRUE)) stop("Install pheatmap first.")

modes <- c("taxa", "taxa_pathways")
column_labels <- c("Taxa-taxa", "Taxon-pathway")
scores <- lapply(modes, function(mode) {
  x <- readRDS(file.path(input_dir, mode, "disco_analysis.rds"))$disco
  if (!all(c("species", "p", "q") %in% names(x)) ||
      anyNA(x$species) || anyDuplicated(x$species)) stop("Invalid feature scores in ", mode)
  if (any(!is.na(x$p) & (!is.finite(x$p) | x$p < 0 | x$p > 1)) ||
      any(!is.na(x$q) & (!is.finite(x$q) | x$q < 0 | x$q > 1)))
    stop("Invalid p or q values in ", mode)
  if (any(is.finite(x$q) & x$q < q_cutoff & !is.finite(x$p)))
    stop("Significant feature has no valid raw p-value in ", mode)
  x
})
names(scores) <- column_labels
taxa <- sort(unique(unlist(lapply(scores, function(x) x$species[is.finite(x$q) & x$q < q_cutoff]))))
if (!length(taxa)) stop("No taxa have DISCO q < ", q_cutoff, " in either analysis.")
p_values <- q_values <- matrix(NA_real_, length(taxa), length(scores),
  dimnames = list(taxa, column_labels))
for (j in seq_along(scores)) {
  idx <- match(taxa, scores[[j]]$species)
  p_values[, j] <- scores[[j]]$p[idx]
  q_values[, j] <- scores[[j]]$q[idx]
}
# Guard against numerical p=0 only for the log transform; export original values.
heatmap_values <- -log10(pmax(p_values, .Machine$double.xmin))
stars <- matrix("", nrow(p_values), ncol(p_values), dimnames = dimnames(p_values))
stars[is.finite(q_values) & q_values < q_cutoff] <- "*"
if (!is.null(row_clusters) && (!cluster_rows || !is.numeric(row_clusters) ||
    length(row_clusters) != 1L || !is.finite(row_clusters) ||
    row_clusters < 1 || row_clusters != floor(row_clusters) || row_clusters > length(taxa)))
  stop("row_clusters must be an integer between 1 and the number of rows, with cluster_rows=TRUE.")
# Cluster only taxa with both DISCO results. Incomplete rows are displayed
# separately below the clustered block, with missing cells grey and no imputation.
complete_rows <- which(stats::complete.cases(heatmap_values))
main_rows <- if (cluster_rows && length(complete_rows)) complete_rows else seq_along(taxa)
partial_rows <- if (cluster_rows && length(complete_rows)) setdiff(seq_along(taxa), complete_rows) else integer(0)
if (!is.null(row_clusters) && row_clusters > length(main_rows))
  stop("row_clusters exceeds the number of rows available for clustering.")
row_tree <- FALSE
if (cluster_rows && length(main_rows) > 1L && length(complete_rows)) {
  row_tree <- stats::hclust(stats::dist(heatmap_values[main_rows, , drop = FALSE]),
    method = row_clustering_method)
}
upper <- max(heatmap_values, na.rm = TRUE)
if (!is.finite(upper) || upper <= 0) upper <- 1
make_heatmap <- function(rows, tree = FALSE, title, legend = TRUE, show_colnames = TRUE, groups = NULL) {
  pheatmap::pheatmap(heatmap_values[rows, , drop = FALSE],
    cluster_rows = tree, cluster_cols = FALSE,
    cutree_rows = if (is.null(groups)) NA_integer_ else groups,
    scale = "none", color = grDevices::colorRampPalette(c("#FFF7EC", "#FDBB84", "#D7301F", "#7F0000"))(100),
    breaks = seq(0, upper, length.out = 101), na_col = "#D1D5DB",
    display_numbers = stars[rows, , drop = FALSE], number_color = "black", fontsize_number = 16,
    fontsize_row = 10, fontsize_col = 11, angle_col = "0", border_color = "white",
    cellwidth = 105, cellheight = 22, main = title, legend = legend,
    show_colnames = show_colnames, labels_col = c("Taxa", "Pathway"), silent = TRUE)
}
heatmap <- make_heatmap(main_rows, tree = row_tree,
  title = paste0("QinJ_2012: significant DISCO taxa\nColor: -log10(raw DISCO p); *: DISCO q < ", q_cutoff),
  show_colnames = TRUE, groups = row_clusters)
partial_heatmap <- NULL
if (length(partial_rows)) {
  partial_heatmap <- make_heatmap(partial_rows,
    title = "Unclustered rows: grey cells indicate unavailable DISCO results", legend = FALSE)
  widths <- grid::unit.pmax(heatmap$gtable$widths, partial_heatmap$gtable$widths)
  heatmap$gtable$widths <- partial_heatmap$gtable$widths <- widths
}
row_order <- c(if (inherits(row_tree, "hclust")) main_rows[row_tree$order] else main_rows, partial_rows)
clusters <- stats::setNames(rep(NA_integer_, length(taxa)), taxa)
if (!is.null(row_clusters) && inherits(row_tree, "hclust"))
  clusters[taxa[main_rows]] <- stats::cutree(row_tree, k = row_clusters)
if (!dir.exists(output_dir) && !dir.create(output_dir, recursive = TRUE)) stop("Could not create output directory.")
base_path <- file.path(output_dir, "significant_taxa_DISCO_heatmap")
plot_height <- max(6, (if (length(partial_rows)) 3.5 else 2.5) + .31 * length(taxa))
save_heatmap <- function(type) {
  if (type == "pdf") grDevices::pdf(paste0(base_path, ".pdf"), width = 10, height = plot_height,
    useDingbats = FALSE) else grDevices::png(paste0(base_path, ".png"), width = 10,
    height = plot_height, units = "in", res = 180)
  on.exit(grDevices::dev.off(), add = TRUE)
  grid::grid.newpage()
  if (length(partial_rows)) {
    grid::pushViewport(grid::viewport(layout = grid::grid.layout(2, 1,
      heights = grid::unit(c(2 + .31 * length(main_rows), 1.5 + .31 * length(partial_rows)), "null"))))
    grid::pushViewport(grid::viewport(layout.pos.row = 1))
    grid::grid.draw(heatmap$gtable)
    grid::popViewport()
    grid::pushViewport(grid::viewport(layout.pos.row = 2))
    grid::grid.draw(partial_heatmap$gtable)
    grid::popViewport(2)
  } else grid::grid.draw(heatmap$gtable)
}
save_heatmap("png")
save_heatmap("pdf")
for (entry in list(list("raw_p_values", p_values), list("q_values", q_values),
    list("neg_log10_p_values", heatmap_values))) {
  table <- data.frame(taxon = taxa[row_order], entry[[2L]][row_order, , drop = FALSE], check.names = FALSE)
  utils::write.csv(table, file.path(output_dir, paste0(entry[[1L]], ".csv")), row.names = FALSE)
}
utils::write.csv(data.frame(taxon = taxa[row_order], row_order = seq_along(taxa),
  row_cluster = unname(clusters[taxa[row_order]])), file.path(output_dir, "row_order.csv"), row.names = FALSE)
saveRDS(list(raw_p = p_values, q = q_values, neg_log10_p = heatmap_values, stars = stars,
  row_tree = row_tree, row_order = taxa[row_order], row_clusters = clusters,
  clustered_taxa = if (inherits(row_tree, "hclust")) taxa[main_rows] else character(),
  unclustered_partial_taxa = taxa[partial_rows],
  settings = list(q_cutoff = q_cutoff, cluster_rows = cluster_rows,
    row_clusters = row_clusters, row_clustering_method = row_clustering_method,
    p_zero_floor = .Machine$double.xmin, column_clustering = FALSE,
    missing_value_handling = "no imputation; incomplete rows appended outside clustering")),
  file.path(output_dir, "heatmap_data.rds"))
message("Saved heatmap for ", length(taxa), " taxa in ", normalizePath(output_dir))
