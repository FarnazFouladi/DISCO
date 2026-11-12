# Top retained taxon-pathway differential correlations: one PDF page per pair,
# with separate Control and T2D panels using the exact matrices passed to DISCO.
# Run: Rscript analysis/diabetes_microbiome/code/plot_qinj_top_taxa_pathway_scatterplots.R
# Optional arguments: input_directory output_directory top_n
# Also supports source("/path/to/DISCO/analysis/diabetes_microbiome/code/plot_qinj_top_taxa_pathway_scatterplots.R").
# Ranks by differential q, then p, then feature names; no new tests are fitted.
# The displayed linear trend is a visual guide, not the Spearman test model.

# Change this setting to show a different number of top-ranked pairs.
top_n <- 1L

script_file <- local({
  files <- as.character(unlist(lapply(sys.frames(), function(frame) frame$ofile), use.names = FALSE))
  files <- files[basename(files) == "plot_qinj_top_taxa_pathway_scatterplots.R"]
  if (length(files)) return(normalizePath(tail(files, 1L), mustWork = TRUE))
  args <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(args) == 1L)
    return(normalizePath(sub("^--file=", "", args[[1L]]), mustWork = TRUE))
  stop("Use Rscript or source() with the full path to this script.")
})
analysis_dir <- dirname(script_file)
sourced <- any(vapply(sys.frames(), function(frame) {
  !is.null(frame$ofile) && identical(basename(frame$ofile), "plot_qinj_top_taxa_pathway_scatterplots.R")
}, logical(1)))
args <- if (sourced) character(0) else commandArgs(trailingOnly = TRUE)
if (length(args) > 3L)
  stop("Usage: Rscript plot_qinj_top_taxa_pathway_scatterplots.R [input_dir] [output_dir] [top_n]")
if (length(args) >= 3L) top_n <- suppressWarnings(as.numeric(args[[3L]]))
if (!is.numeric(top_n) || length(top_n) != 1L || !is.finite(top_n) ||
    top_n < 1 || top_n != floor(top_n) || top_n > .Machine$integer.max)
  stop("top_n must be a positive integer.")
top_n <- as.integer(top_n)
input_dir <- if (length(args) >= 1L) args[[1L]] else file.path(
  dirname(analysis_dir), "results", "taxa_pathways"
)
output_dir <- if (length(args) >= 2L) args[[2L]] else file.path(
  dirname(analysis_dir), "results", "taxa_pathway_scatterplots"
)
if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Install ggplot2 first.")
result <- readRDS(file.path(input_dir, "disco_analysis.rds"))
inputs <- readRDS(file.path(input_dir, "disco_inputs.rds"))
if (length(inputs$data1) != 2L || length(inputs$data2) != 2L) {
  stop("This script expects the QinJ taxon-pathway analysis (two modalities).")
}
matrices <- list(Control = inputs$data1, T2D = inputs$data2)
pairs <- result$cor_diff
required <- c("t1", "t2", "r1", "r2", "n1", "n2", "p", "q")
if (!all(required %in% names(pairs))) stop("Missing columns in cor_diff.")
pairs <- pairs[is.finite(pairs$p) & is.finite(pairs$q) &
                 is.finite(pairs$r1) & is.finite(pairs$r2), , drop = FALSE]
pairs <- pairs[order(pairs$q, pairs$p, pairs$t1, pairs$t2), , drop = FALSE]
if (nrow(pairs) < top_n) stop("Requested ", top_n, " pairs, but only ", nrow(pairs), " valid pairs are available.")
pairs <- head(pairs, top_n)
pairs$page <- seq_len(nrow(pairs))
format_p <- function(x) formatC(x, format = "e", digits = 2)
plot_data <- vector("list", nrow(pairs))
plots <- vector("list", nrow(pairs))
for (i in seq_len(nrow(pairs))) {
  pair <- pairs[i, ]
  group_data <- lapply(seq_along(matrices), function(j) {
    taxa <- as.matrix(matrices[[j]][[1L]])
    pathways <- as.matrix(matrices[[j]][[2L]])
    if (!pair$t1 %in% colnames(taxa) || !pair$t2 %in% colnames(pathways))
      stop("Taxon or pathway missing from saved inputs.")
    ids <- rownames(taxa)[rownames(taxa) %in% rownames(pathways)]
    values <- cbind(taxa[ids, pair$t1], pathways[ids, pair$t2])
    rownames(values) <- ids
    # Match DISCO's pairwise complete observations and validate saved outputs.
    if (any(is.infinite(values))) stop("Infinite measurement in saved inputs.")
    complete <- stats::complete.cases(values)
    values <- values[complete, , drop = FALSE]
    observed_r <- stats::cor(values[, 1], values[, 2], method = "spearman")
    expected_r <- if (j == 1L) pair$r1 else pair$r2
    expected_n <- if (j == 1L) pair$n1 else pair$n2
    if (nrow(values) != expected_n || !isTRUE(all.equal(observed_r, expected_r, tolerance = 1e-8))) {
      stop("Saved inputs and correlation outputs do not agree on page ", i, ".")
    }
    data.frame(
      page = i, t1 = pair$t1, t2 = pair$t2,
      sample_id = ids[complete], group = names(matrices)[j],
      x = values[, 1], y = values[, 2], stringsAsFactors = FALSE
    )
  })
  data <- do.call(rbind, group_data)
  data$group <- factor(data$group, levels = c("Control", "T2D"))
  panel_labels <- c(
    Control = sprintf("Control\nSpearman r = %.3f | n = %d", pair$r1, pair$n1),
    T2D = sprintf("T2D\nSpearman r = %.3f | n = %d", pair$r2, pair$n2)
  )
  plots[[i]] <- ggplot2::ggplot(data, ggplot2::aes(x = x, y = y, color = group)) +
    ggplot2::geom_hline(yintercept = 0, linewidth = 0.3, color = "#E2E8F0") +
    ggplot2::geom_vline(xintercept = 0, linewidth = 0.3, color = "#E2E8F0") +
    ggplot2::geom_point(size = 1.8, alpha = 0.6) +
    ggplot2::geom_smooth(method = "lm", formula = y ~ x, se = FALSE,
                          linewidth = 0.8, linetype = "dashed") +
    ggplot2::facet_wrap(~ group, nrow = 1, labeller = ggplot2::as_labeller(panel_labels)) +
    ggplot2::scale_color_manual(values = c(Control = "#2563EB", T2D = "#DC2626"), guide = "none") +
    ggplot2::labs(
      title = paste0("Differential correlation ", i, " of ", top_n),
      subtitle = paste0(paste(strwrap(paste(pair$t1, "vs", pair$t2), width = 95), collapse = "\n"), "\nDifferential p = ",
                         format_p(pair$p), " | adjusted q = ", format_p(pair$q)),
      x = paste0(pair$t1, "\nCovariate-adjusted robust CLR abundance"),
      y = paste0(paste(strwrap(pair$t2, width = 58), collapse = "\n"), "\nCovariate-adjusted robust CLR abundance"),
      caption = paste0(
        "Each point is one sample. Identical axis scales across groups within each page.\n",
        "Inputs are robust CLR abundances residualized against gender and BMI separately within each group.\n",
        "Only complete residual pairs are plotted; abundance zeros remain missing.\n",
        "Dashed lines: linear visual guides; annotations report saved Spearman correlations. ",
        "Pairs ranked by differential q, then p; selected DISCO filter = ", result$selected_threshold, "."
      )
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 18),
      plot.subtitle = ggplot2::element_text(size = 11, lineheight = 1.25),
      strip.text = ggplot2::element_text(face = "bold", size = 12, lineheight = 1.2),
      panel.grid.minor = ggplot2::element_blank(),
      panel.spacing = grid::unit(1.5, "lines"),
      axis.title = ggplot2::element_text(size = 11),
      plot.caption = ggplot2::element_text(hjust = 0, size = 9, margin = ggplot2::margin(t = 14)),
      plot.margin = ggplot2::margin(18, 18, 18, 18)
    )
  plot_data[[i]] <- data
}
# Validate all pages before opening the output device.
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
pdf_path <- file.path(output_dir, paste0("top", top_n, "_differential_correlations_control_T2D.pdf"))
save_pdf <- function() {
  # Complete the PDF locally before replacing a file in a cloud-synced folder.
  temporary_pdf <- tempfile(fileext = ".pdf")
  on.exit(unlink(temporary_pdf), add = TRUE)
  grDevices::pdf(temporary_pdf, width = 12, height = 8.5, onefile = TRUE, useDingbats = FALSE)
  tryCatch({
    for (plot in plots) print(plot)
  }, finally = grDevices::dev.off())
  if (!file.copy(temporary_pdf, pdf_path, overwrite = TRUE))
    stop("Could not replace PDF: ", pdf_path)
}
save_pdf()
ggplot2::ggsave(file.path(output_dir, paste0("top", top_n, "_page01_preview.png")), plots[[1L]],
                 width = 12, height = 8.5, dpi = 160, bg = "white")
utils::write.csv(pairs[, c("page", required)], file.path(output_dir, paste0("top", top_n, "_pairs.csv")), row.names = FALSE)
utils::write.csv(do.call(rbind, plot_data), file.path(output_dir, "scatterplot_data.csv"), row.names = FALSE)
message("Saved ", top_n, "-page PDF: ", normalizePath(pdf_path))
