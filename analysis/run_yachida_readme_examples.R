# Generate the real-data results and portable figures used by README.Rmd.
# Run from any directory: Rscript /path/to/DISCO/analysis/run_yachida_readme_examples.R
# Reuses the validated taxa-only run when its CLR input and settings match.
# Set DISCO_README_RECOMPUTE=true to recompute both examples.

file_argument <- grep("^--file=", commandArgs(FALSE), value = TRUE)
if (length(file_argument) != 1L) stop("Run this script with Rscript.")
analysis_dir <- dirname(normalizePath(sub("^--file=", "", file_argument[[1L]])))
repo_dir <- dirname(analysis_dir)
source(file.path(analysis_dir, "yachida_readme_preprocessing.R"))
required <- c("pkgload", "curatedMetagenomicData", "SummarizedExperiment", "ggplot2")
missing <- required[!vapply(required, requireNamespace, logical(1L), quietly = TRUE)]
if (length(missing)) stop("Install: ", paste(missing, collapse = ", "))
pkgload::load_all(repo_dir, quiet = TRUE, export_all = FALSE)
recompute <- identical(tolower(Sys.getenv("DISCO_README_RECOMPUTE")), "true")

metadata <- curatedMetagenomicData::sampleMetadata
metadata <- metadata[
  !is.na(metadata$study_name) & metadata$study_name == "YachidaS_2019" &
    metadata$disease %in% c("healthy", "CRC"), , drop = FALSE
]
taxa_se <- curatedMetagenomicData::returnSamples(
  metadata, dataType = "relative_abundance", counts = FALSE, rownames = "short"
)
taxa <- t(as.matrix(SummarizedExperiment::assay(taxa_se)))
taxa_ids <- metadata$sample_id[metadata$sample_id %in% rownames(taxa)]
taxa <- taxa[taxa_ids, , drop = FALSE]
taxa <- taxa[rowSums(taxa) > 0, , drop = FALSE]
group <- stats::setNames(as.character(metadata$disease), metadata$sample_id)
taxa_group <- group[rownames(taxa)]
taxa_prepared <- preprocess_abundance(taxa)

single_dir <- file.path(analysis_dir, "yachida_local_disco_results",
                        "prevalence_20pct_half_min_clr")
single_file <- file.path(single_dir, "disco_analysis.rds")
single_input <- file.path(single_dir, "preprocessed_data.rds")
reuse_single <- !recompute && file.exists(single_file) && file.exists(single_input)
if (reuse_single) {
  old <- readRDS(single_input)
  reuse_single <- isTRUE(all.equal(old$taxa_clr, taxa_prepared$clr)) &&
    identical(old$group, taxa_group) && old$settings$n_resamples == 100 &&
    old$settings$q_value == 0.2 && old$settings$global_p_cutoff == 0.05 &&
    old$settings$random_seed == 1L && old$settings$cooccur == 10L &&
    old$settings$adj_method == "BY" && old$settings$covariate_adjustment == "none"
}
if (reuse_single) {
  message("Reusing the validated taxa-only Yachida results.")
  single <- readRDS(single_file)
} else {
  single <- DISCO::disco(
    data1 = list(taxa_prepared$clr[taxa_group == "healthy", , drop = FALSE]),
    data2 = list(taxa_prepared$clr[taxa_group == "CRC", , drop = FALSE]),
    n_resamples = 100, random_seed = 1L, q_value = 0.2, global_p_cutoff = 0.05,
    output_dir = file.path(analysis_dir, "yachida_readme_taxa_results"),
    prefix = "YachidaS_2019_taxa"
  )
}

pathway_se <- curatedMetagenomicData::returnSamples(metadata, dataType = "pathway_abundance")
pathways <- t(as.matrix(SummarizedExperiment::assay(pathway_se)))
pathway_ids <- colnames(pathways)
pathways <- pathways[, !grepl("|", pathway_ids, fixed = TRUE) &
                       !toupper(trimws(pathway_ids)) %in% c("UNMAPPED", "UNINTEGRATED"),
                     drop = FALSE]
common_ids <- rownames(taxa)[rownames(taxa) %in% rownames(pathways)]
common_ids <- common_ids[rowSums(pathways[common_ids, , drop = FALSE]) > 0]
if (!length(common_ids)) stop("No usable common samples.")
two_group <- group[common_ids]
if (any(table(factor(two_group, levels = c("healthy", "CRC"))) < 10L)) {
  stop("Each group requires at least 10 samples.")
}
two_taxa <- preprocess_abundance(taxa[common_ids, , drop = FALSE])
two_pathways <- preprocess_abundance(pathways[common_ids, , drop = FALSE])
two_dir <- file.path(analysis_dir, "yachida_readme_two_modality_results")
two_file <- file.path(two_dir, "disco_analysis.rds")
two_input <- file.path(two_dir, "preprocessed_data.rds")
settings <- list(n_resamples = 100L, random_seed = 1L, cooccur = 10L,
                 adj_method = "BY", q_value = 0.2, global_p_cutoff = 0.05,
                 prevalence_cutoff = 0.2, covariate_adjustment = "none")
reuse_two <- !recompute && file.exists(two_file) && file.exists(two_input)
if (reuse_two) {
  old <- readRDS(two_input)
  reuse_two <- identical(old$settings, settings) && identical(old$group, two_group) &&
    isTRUE(all.equal(old$taxa$clr, two_taxa$clr)) &&
    isTRUE(all.equal(old$pathways$clr, two_pathways$clr))
}
if (reuse_two) {
  two <- readRDS(two_file)
} else {
  dir.create(two_dir, recursive = TRUE, showWarnings = FALSE)
  saveRDS(list(taxa = two_taxa, pathways = two_pathways, group = two_group,
               settings = settings), two_input)
  message("Running two modalities: ", ncol(two_taxa$clr), " taxa, ",
          ncol(two_pathways$clr), " pathways; control n = ",
          sum(two_group == "healthy"), "; CRC n = ", sum(two_group == "CRC"), ".")
  two <- DISCO::disco(
    data1 = list(two_taxa$clr[two_group == "healthy", , drop = FALSE],
                 two_pathways$clr[two_group == "healthy", , drop = FALSE]),
    data2 = list(two_taxa$clr[two_group == "CRC", , drop = FALSE],
                 two_pathways$clr[two_group == "CRC", , drop = FALSE]),
    n_resamples = 100, random_seed = 1L, q_value = 0.2, global_p_cutoff = 0.05,
    output_dir = two_dir, prefix = "YachidaS_2019_taxa_pathways"
  )
  saveRDS(two, two_file)
}

summarize_example <- function(result, group, taxa, n_pathways = NA_integer_) {
  scores <- result$disco[order(result$disco$q, result$disco$p), , drop = FALSE]
  list(
    cohort = data.frame(control = sum(group == "healthy"), CRC = sum(group == "CRC"),
                         retained_taxa = ncol(taxa), retained_pathways = n_pathways),
    global_test = result$global_test,
    selected_threshold = result$selected_threshold,
    significant_taxa = sum(result$disco$q < 0.2, na.rm = TRUE),
    top_taxa = utils::head(scores[, c("species", "p", "q", "num_corr")], 10L)
  )
}
bundle <- list(
  settings = settings,
  taxa = summarize_example(single, taxa_group, taxa_prepared$clr),
  taxa_pathways = summarize_example(two, two_group, two_taxa$clr, ncol(two_pathways$clr)),
  versions = vapply(c("DISCO", "curatedMetagenomicData"), function(p) {
    as.character(utils::packageVersion(p))
  }, character(1L))
)
saveRDS(bundle, file.path(repo_dir, "README-results.rds"))
figure_dir <- file.path(repo_dir, "man", "figures")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
for (mode in c("taxa", "taxa_pathways")) {
  result <- if (mode == "taxa") single else two
  for (plot_name in c("discovery_plot", "ratio_plot", "feature_null_plot")) {
    is_feature <- plot_name == "feature_null_plot"
    ggplot2::ggsave(
      filename = file.path(figure_dir, paste0("README-yachida-", mode, "-", plot_name, ".png")),
      plot = result[[plot_name]], width = if (is_feature) 13 else 8.5,
      height = if (is_feature) 8.5 else 6.5, dpi = 150, bg = "white"
    )
  }
}
capture.output(utils::sessionInfo(), file = file.path(two_dir, "session_info.txt"))
message("Saved README-results.rds and six README figures. Render README.Rmd to update README.md.")
