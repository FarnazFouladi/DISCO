# Generate QinJ_2012 T2D/control results, following the Nielsen README workflow.
# Run from any directory: Rscript /path/to/DISCO/analysis/diabetes_microbiome/code/run_qinj_readme_examples.R
# Or in R: source("/path/to/DISCO/analysis/diabetes_microbiome/code/run_qinj_readme_examples.R")
# Both modalities use separate healthy/T2D gender/BMI models after robust CLR; zeros remain NA.
# Both analyses use the installed DISCO package, loaded with library(DISCO).
# Set DISCO_README_RECOMPUTE=true to recompute both examples.

script_file <- local({
  files <- as.character(unlist(lapply(sys.frames(), function(frame) frame$ofile),
                               use.names = FALSE))
  files <- files[basename(files) == "run_qinj_readme_examples.R"]
  if (length(files)) return(normalizePath(tail(files, 1L), mustWork = TRUE))
  args <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(args) == 1L)
    return(normalizePath(sub("^--file=", "", args[[1L]]), mustWork = TRUE))
  stop("Use Rscript or source() with the full path to this script.")
})
analysis_dir <- dirname(script_file)
repo_dir <- dirname(dirname(dirname(analysis_dir)))
source(file.path(analysis_dir, "qinj_readme_preprocessing.R"))
required <- c("DISCO", "curatedMetagenomicData", "SummarizedExperiment", "ggplot2")
missing <- required[!vapply(required, requireNamespace, logical(1L), quietly = TRUE)]
if (length(missing)) stop("Install: ", paste(missing, collapse = ", "))
## ---- readme-package ----
library(DISCO)

## ---- internal-readme-startup ----
expected <- c("workers", "empirically_derive_threshold", "global_p_cutoff")
if (!all(expected %in% names(formals(DISCO::disco))))
  stop("Install an updated DISCO package before running this example.")
recompute <- identical(tolower(Sys.getenv("DISCO_README_RECOMPUTE")), "true")

## ---- readme-settings ----
settings <- list(
  workers = 12L,
  n_resamples = 1000L, random_seed = 1L, cooccur = 20L, adj_method = "BY",
  q_value = 0.2, global_p_cutoff = 0.05, prevalence_cutoff = 0.5,
  transformation = "robust CLR; zeros remain NA",
  prevalence_scope = "at least 50% in each of healthy and T2D",
  covariate_adjustment = "separate healthy and T2D models: ~ gender + BMI"
)
results_dir <- file.path(dirname(analysis_dir), "results")
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
source_md5 <- readme_source_hashes(repo_dir)

## ---- retrieve-taxa ----
selection <- select_qinj_metadata(curatedMetagenomicData::sampleMetadata)
metadata <- selection$metadata
excluded_covariates <- selection$excluded_covariates
sample_metadata <- metadata
rownames(sample_metadata) <- sample_metadata$sample_id
message("Excluded ", nrow(excluded_covariates), " sample(s) missing gender or BMI; ",
        nrow(selection$excluded_subjects), " conflicting or repeated sample(s).")
taxa_se <- curatedMetagenomicData::returnSamples(
  metadata, dataType = "relative_abundance", counts = FALSE, rownames = "short"
)
taxa <- t(as.matrix(SummarizedExperiment::assay(taxa_se)))
taxa_ids <- metadata$sample_id[metadata$sample_id %in% rownames(taxa)]
taxa <- taxa[taxa_ids, , drop = FALSE]
if (any(!is.finite(taxa)) || any(taxa < 0)) stop("Invalid taxon abundances.")
taxa <- taxa[rowSums(taxa) > 0, , drop = FALSE]
group <- stats::setNames(as.character(metadata$disease), metadata$sample_id)
taxa_group <- group[rownames(taxa)]
if (any(table(factor(taxa_group, levels = c("healthy", "T2D"))) < settings$cooccur)) {
  stop("Each group must have at least cooccur usable samples.")
}

## ---- run-taxa ----
taxa_prepared <- preprocess_abundance(taxa, taxa_group, settings$prevalence_cutoff)
taxa_adjusted <- residualize_abundance(
  taxa_prepared$clr, sample_metadata[rownames(taxa_prepared$clr), , drop = FALSE],
  taxa_group
)
single_dir <- file.path(results_dir, "taxa")
single_input <- list(
  data1 = list(taxa_adjusted$residuals[taxa_group == "healthy", , drop = FALSE]),
  data2 = list(taxa_adjusted$residuals[taxa_group == "T2D", , drop = FALSE]),
  settings = settings, source_md5 = source_md5
)
single <- readme_cached_result(single_dir, single_input, recompute)
if (is.null(single)) {
  single <- DISCO::disco(
    data1 = single_input$data1, data2 = single_input$data2,
    workers = settings$workers,
    n_resamples = settings$n_resamples, random_seed = settings$random_seed,
    q_value = settings$q_value, global_p_cutoff = settings$global_p_cutoff,
    cooccur = settings$cooccur, adj_method = settings$adj_method,
    output_dir = single_dir,
    prefix = "QinJ_2012_taxa"
  )
  saveRDS(single_input, file.path(single_dir, "disco_inputs.rds"))
  saveRDS(single, file.path(single_dir, "disco_analysis.rds"))
}

## ---- run-two-modalities ----
pathway_se <- curatedMetagenomicData::returnSamples(metadata, dataType = "pathway_abundance")
pathways <- t(as.matrix(SummarizedExperiment::assay(pathway_se)))
if (any(!is.finite(pathways)) || any(pathways < 0)) stop("Invalid pathway abundances.")
pathway_ids <- colnames(pathways)
pathways <- pathways[, !grepl("|", pathway_ids, fixed = TRUE) &
                       !toupper(trimws(pathway_ids)) %in% c("UNMAPPED", "UNINTEGRATED"),
                     drop = FALSE]
common_ids <- rownames(taxa)[rownames(taxa) %in% rownames(pathways)]
common_ids <- common_ids[rowSums(pathways[common_ids, , drop = FALSE]) > 0]
if (!length(common_ids)) stop("No usable common samples.")
two_group <- group[common_ids]
if (any(table(factor(two_group, levels = c("healthy", "T2D"))) < settings$cooccur)) {
  stop("Each group requires at least cooccur samples.")
}
two_taxa <- preprocess_abundance(taxa[common_ids, , drop = FALSE], two_group,
                                 settings$prevalence_cutoff)
two_pathways <- preprocess_abundance(pathways[common_ids, , drop = FALSE], two_group,
                                     settings$prevalence_cutoff)
two_metadata <- sample_metadata[common_ids, , drop = FALSE]
two_taxa_adjusted <- residualize_abundance(two_taxa$clr, two_metadata, two_group)
two_pathways_adjusted <- residualize_abundance(two_pathways$clr, two_metadata, two_group)
two_dir <- file.path(results_dir, "taxa_pathways")
two_input <- list(
  data1 = list(two_taxa_adjusted$residuals[two_group == "healthy", , drop = FALSE],
               two_pathways_adjusted$residuals[two_group == "healthy", , drop = FALSE]),
  data2 = list(two_taxa_adjusted$residuals[two_group == "T2D", , drop = FALSE],
               two_pathways_adjusted$residuals[two_group == "T2D", , drop = FALSE]),
  settings = settings, source_md5 = source_md5
)
two <- readme_cached_result(two_dir, two_input, recompute)
if (is.null(two)) {
  message("Running two modalities: ", ncol(two_taxa$clr), " taxa, ",
          ncol(two_pathways$clr), " pathways; control n = ",
          sum(two_group == "healthy"), "; T2D n = ", sum(two_group == "T2D"), ".")
  two <- DISCO::disco(
    data1 = two_input$data1, data2 = two_input$data2,
    workers = settings$workers,
    n_resamples = settings$n_resamples, random_seed = settings$random_seed,
    q_value = settings$q_value, global_p_cutoff = settings$global_p_cutoff,
    cooccur = settings$cooccur, adj_method = settings$adj_method,
    output_dir = two_dir, prefix = "QinJ_2012_taxa_pathways"
  )
  saveRDS(two_input, file.path(two_dir, "disco_inputs.rds"))
  saveRDS(two, file.path(two_dir, "disco_analysis.rds"))
}

## ---- internal-readme-export ----
# Save this workflow separately; retain the previous pooled and unadjusted runs.
saveRDS(list(taxa = taxa_prepared, adjustment = taxa_adjusted, group = taxa_group,
             metadata = sample_metadata[rownames(taxa), , drop = FALSE], settings = settings),
        file.path(single_dir, "preprocessed_data.rds"))
saveRDS(list(taxa = two_taxa, pathways = two_pathways,
             taxa_adjustment = two_taxa_adjusted,
             pathway_adjustment = two_pathways_adjusted, group = two_group,
             metadata = two_metadata, settings = settings),
        file.path(two_dir, "preprocessed_data.rds"))
utils::write.csv(selection$excluded_subjects, file.path(results_dir, "excluded_subject_samples.csv"),
                 row.names = FALSE)
utils::write.csv(excluded_covariates, file.path(results_dir, "excluded_covariate_samples.csv"),
                 row.names = FALSE)
summarize_example <- function(result, group, taxa, n_pathways = NA_integer_) {
  scores <- result$disco[order(result$disco$q, result$disco$p), , drop = FALSE]
  list(
    cohort = data.frame(control = sum(group == "healthy"), T2D = sum(group == "T2D"),
                         retained_taxa = ncol(taxa), retained_pathways = n_pathways),
    global_test = result$global_test,
    selected_threshold = result$selected_threshold,
    significant_taxa = sum(result$disco$q < settings$q_value, na.rm = TRUE),
    top_taxa = utils::head(scores[, c("species", "p", "q", "num_corr")], 10L)
  )
}
bundle <- list(
  study = "QinJ_2012", comparison = "T2D versus healthy controls",
  DISCO_package_path = find.package("DISCO"),
  settings = settings,
  source_md5 = source_md5,
  excluded_covariate_samples = nrow(excluded_covariates),
  excluded_subject_samples = nrow(selection$excluded_subjects),
  taxa = summarize_example(single, taxa_group, taxa_prepared$clr),
  taxa_pathways = summarize_example(two, two_group, two_taxa$clr, ncol(two_pathways$clr)),
  versions = vapply(c("DISCO", "curatedMetagenomicData"), function(p) {
    as.character(utils::packageVersion(p))
  }, character(1L))
)
saveRDS(bundle, file.path(results_dir, "README-results.rds"))
figure_dir <- file.path(results_dir, "figures")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
for (mode in c("taxa", "taxa_pathways")) {
  result <- if (mode == "taxa") single else two
  for (plot_name in c("discovery_plot", "ratio_plot", "feature_null_plot")) {
    is_feature <- plot_name == "feature_null_plot"
    ggplot2::ggsave(
      filename = file.path(figure_dir, paste0("README-qinj-", mode, "-", plot_name, ".png")),
      plot = result[[plot_name]], width = if (is_feature) 13 else 8.5,
      height = if (is_feature) 8.5 else 6.5, dpi = 150, bg = "white"
    )
  }
}
capture.output(utils::sessionInfo(), file = file.path(two_dir, "session_info.txt"))
message("Saved QinJ_2012 summaries and six figures in ", results_dir, ".")
