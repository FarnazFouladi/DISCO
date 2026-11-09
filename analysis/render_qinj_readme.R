# Refresh the QinJ README figures and render README.md without rerunning DISCO.
# Rscript analysis/render_qinj_readme.R
# Also supports source("/path/to/DISCO/analysis/render_qinj_readme.R").
script_file <- local({
  files <- as.character(unlist(lapply(sys.frames(), function(frame) frame$ofile), use.names = FALSE))
  files <- files[basename(files) == "render_qinj_readme.R"]
  if (length(files)) return(normalizePath(tail(files, 1L), mustWork = TRUE))
  args <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(args) == 1L)
    return(normalizePath(sub("^--file=", "", args[[1L]]), mustWork = TRUE))
  stop("Use Rscript or source() with the full path to this script.")
})
repo_dir <- dirname(dirname(script_file))
if (!requireNamespace("rmarkdown", quietly = TRUE)) stop("Install rmarkdown first.")

diagnostic_names <- paste0("README-qinj-",
  rep(c("taxa", "taxa_pathways"), each = 3), "-",
  rep(c("discovery_plot", "ratio_plot", "feature_null_plot"), 2), ".png")
assets <- stats::setNames(
  file.path(repo_dir, "analysis/qinj_readme_results/figures", diagnostic_names), diagnostic_names)
assets <- c(assets,
  "README-qinj-taxa-network.png" = file.path(repo_dir, "analysis/qinj_taxa_networks/top10_taxa_combined_control_T2D_networks.png"),
  "README-qinj-taxa-scatterplot.png" = file.path(repo_dir, "analysis/qinj_differential_scatterplots/top10_page01_preview.png"),
  "README-qinj-taxa_pathways-network.png" = file.path(repo_dir, "analysis/qinj_taxa_pathway_networks/top5_taxa_pathways_control_T2D_networks.png"),
  "README-qinj-taxa_pathways-scatterplot.png" = file.path(repo_dir, "analysis/qinj_taxa_pathway_scatterplots/top1_page01_preview.png"),
  "README-qinj-DISCO-heatmap.png" = file.path(repo_dir, "analysis/qinj_disco_heatmap/significant_taxa_DISCO_heatmap.png"))
if (!all(file.exists(assets))) stop("Missing example figures; rerun the plotting scripts first.")
copied <- file.copy(assets, file.path(repo_dir, "man", "figures", names(assets)), overwrite = TRUE)
if (!all(copied)) stop("Could not refresh all README figures.")
rmarkdown::render(file.path(repo_dir, "README.Rmd"), knit_root_dir = repo_dir)
