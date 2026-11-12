# Refresh the current example figures and render README.md without rerunning DISCO.
# Rscript analysis/render_readme.R
# Also supports source("/path/to/DISCO/analysis/render_readme.R").
script_file <- local({
  files <- as.character(unlist(lapply(sys.frames(), function(frame) frame$ofile), use.names = FALSE))
  files <- files[basename(files) == "render_readme.R"]
  if (length(files)) return(normalizePath(tail(files, 1L), mustWork = TRUE))
  args <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(args) == 1L)
    return(normalizePath(sub("^--file=", "", args[[1L]]), mustWork = TRUE))
  stop("Use Rscript or source() with the full path to this script.")
})
repo_dir <- dirname(dirname(script_file))
if (!requireNamespace("rmarkdown", quietly = TRUE)) stop("Install rmarkdown first.")

# Refresh microbiome plots and their tables from the current saved analyses.
for (name in c("plot_qinj_taxa_networks.R", "plot_qinj_top_differential_scatterplots.R",
               "plot_qinj_taxa_pathway_networks.R", "plot_qinj_top_taxa_pathway_scatterplots.R",
               "plot_qinj_disco_heatmap.R")) {
  script <- file.path(repo_dir, "analysis/diabetes_microbiome/code", name)
  status <- system2(file.path(R.home("bin"), "Rscript"), c("--vanilla", shQuote(script)))
  if (status != 0L) stop("Figure generation failed: ", name, call. = FALSE)
}

# Rebuild the IBD figure from the current results, preserving saved display settings.
# This does not rerun DISCO.
plot_code <- paste(
  "args <- commandArgs(trailingOnly = TRUE);",
  "settings_file <- file.path(args[[1L]], 'figure', 'plot_settings.rds');",
  "if (file.exists(settings_file)) { settings <- readRDS(settings_file);",
  "options(ibd.pathway_min_percent = settings$min_percent,",
  "ibd.pathway_label_max_chars = settings$label_max_chars,",
  "ibd.pathway_cluster_rows = settings$cluster_rows) };",
  "source(args[[2L]])"
)
status <- system2(file.path(R.home("bin"), "Rscript"),
  args = c("--vanilla", "-e", shQuote(plot_code),
    shQuote(file.path(repo_dir, "analysis/ibd_transcriptomics/results")),
    shQuote(file.path(repo_dir, "analysis/ibd_transcriptomics/code/plot_ibd_hallmark_significance.R"))))
if (status != 0L) stop("IBD figure generation failed; README was not rendered.", call. = FALSE)

diagnostic_names <- paste0("README-qinj-",
  rep(c("taxa", "taxa_pathways"), each = 3), "-",
  rep(c("discovery_plot", "ratio_plot", "feature_null_plot"), 2), ".png")
assets <- stats::setNames(
  file.path(repo_dir, "analysis/diabetes_microbiome/results/figures", diagnostic_names), diagnostic_names)
assets <- c(assets,
  "README-qinj-taxa-network.png" = file.path(repo_dir, "analysis/diabetes_microbiome/results/taxa_networks/top10_taxa_combined_control_T2D_networks.png"),
  "README-qinj-taxa-scatterplot.png" = file.path(repo_dir, "analysis/diabetes_microbiome/results/differential_scatterplots/top10_page01_preview.png"),
  "README-qinj-taxa_pathways-network.png" = file.path(repo_dir, "analysis/diabetes_microbiome/results/taxa_pathway_networks/top5_taxa_pathways_control_T2D_networks.png"),
  "README-qinj-taxa_pathways-scatterplot.png" = file.path(repo_dir, "analysis/diabetes_microbiome/results/taxa_pathway_scatterplots/top1_page01_preview.png"),
  "README-qinj-DISCO-heatmap.png" = file.path(repo_dir, "analysis/diabetes_microbiome/results/disco_heatmap/significant_taxa_DISCO_heatmap.png"),
  "README-ibd-hallmark-significance.png" = file.path(repo_dir, "analysis/ibd_transcriptomics/results/figure/pathway_significance_heatmap.png"))
if (!all(file.exists(assets))) stop("Missing example figures; rerun the plotting scripts first.")
copied <- file.copy(assets, file.path(repo_dir, "man", "figures", names(assets)), overwrite = TRUE)
if (!all(copied)) stop("Could not refresh all README figures.")
# Render the saved results in a clean R process.
render_code <- paste(
  "args <- commandArgs(trailingOnly = TRUE);",
  "rmarkdown::render(args[[1L]], knit_root_dir = args[[2L]])"
)
status <- system2(file.path(R.home("bin"), "Rscript"),
  args = c("--vanilla", "-e", shQuote(render_code),
    shQuote(file.path(repo_dir, "README.Rmd")), shQuote(repo_dir)))
if (status != 0L) stop("README rendering failed; see the error above.", call. = FALSE)
