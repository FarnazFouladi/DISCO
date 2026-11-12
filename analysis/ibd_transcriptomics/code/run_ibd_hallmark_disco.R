# Adapted for the DISCO README from DISCO_Analysis; scientific calculations unchanged.
# Rscript analysis/ibd_transcriptomics/code/run_ibd_hallmark_disco.R
# options(ibd.pathway_preprocess_only=TRUE) validates all inputs without DISCO.
# options(ibd.pathway_names="HALLMARK_TNFA_SIGNALING_VIA_NFKB") runs one pathway.
# options(ibd.disco_workers=6L) sets permutation workers (default 6).
# Default: saved Hallmark membership; no pathway-size restriction.
# Reactome: options(ibd.msigdb_collection="C2", ibd.msigdb_subcollection="CP:REACTOME",
#                   ibd.min_pathway_genes=50L)
# Hallmark: options(ibd.msigdb_collection="H", ibd.msigdb_subcollection=NULL,
#                   ibd.min_pathway_genes=0L, ibd.pathway_names=NULL)
# Minimum size counts unique database gene symbols BEFORE expression filtering.
# Output folders use MSigDB gs_id; pathway_lookup.csv retains the full names.
# Defaults: iCD/cCD/UC, each versus non-IBD controls; optional pooled CD.
# options(ibd.comparisons=c("CD", "iCD", "cCD", "UC")) adds pooled Crohn's disease.
# options(ibd.exploration_dir="...") selects the audited metadata and downloaded counts.
# Complete age-at-diagnosis and sex required. No batch IDs available; no batch matching.
# Reuses the established GENCODE v26 redundancy audit and Hallmark snapshot.
# Controls are reused across comparisons; pooled CD overlaps iCD/cCD.
# Default preprocessing: one complete-case four-group filter, TMM normalization,
# age/sex correction preserving all four groups, and valid variance in every group.
# options(ibd.preprocessing_mode="comparison_specific") restores the earlier workflow.
# Bundled metadata and reference inputs live in ../results/inputs relative to this code directory.
# No automatic downloads or automatic plotting; run plotting scripts after completion.
script <- local({
  paths <- as.character(unlist(lapply(sys.frames(), function(x) x$ofile), use.names = FALSE))
  paths <- c(rev(paths), sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)),
    "analysis/ibd_transcriptomics/code/run_ibd_hallmark_disco.R")
  paths <- paths[basename(paths) == "run_ibd_hallmark_disco.R" & file.exists(paths)]
  if (!length(paths)) stop("Source this script or use Rscript from the repository root.")
  normalizePath(paths[1])
})
root <- dirname(dirname(script))
for (p in c("DISCO", "recount3", "edgeR", "limma"))
  if (!requireNamespace(p, quietly = TRUE)) stop("Required package: ", p)
library(DISCO)
helpers <- file.path(dirname(script), c("ibd_hallmark_disco_helpers.R", "ibd_pathway_file_helpers.R",
                                      "ibd_result_files.R"))
for (helper in helpers) source(helper)
input <- getOption("ibd.exploration_dir", file.path(root, "results", "inputs"))
reference <- getOption("ibd.pathway_reference", file.path(root, "results", "inputs", "reference"))
comparisons <- ibd_validate_comparisons(getOption("ibd.comparisons", c("iCD", "cCD", "UC")))
preprocessing_mode <- getOption("ibd.preprocessing_mode", "joint_four_groups")
if (!is.character(preprocessing_mode) || length(preprocessing_mode) != 1L || is.na(preprocessing_mode) ||
    !preprocessing_mode %in% c("joint_four_groups", "comparison_specific"))
  stop("ibd.preprocessing_mode must be joint_four_groups or comparison_specific.")
count_input_files <- file.path(input, c("sample_metadata.csv", "raw/recount_metadata.rds",
  "raw/sra.gene_sums.SRP042228.G026.gz"))
redundancy_file <- file.path(reference, "redundancy_candidates_all_genes.csv")
files <- c(count_input_files, redundancy_file)
if (!all(file.exists(files))) stop("Missing input files: ", paste(files[!file.exists(files)], collapse = ", "))
collection <- getOption("ibd.msigdb_collection", "H")
subcollection <- getOption("ibd.msigdb_subcollection", NULL)
loaded <- ibd_load_pathway_membership(reference, collection, subcollection)
membership <- loaded$membership
selection <- ibd_select_pathways(membership, getOption("ibd.min_pathway_genes", 0L),
  getOption("ibd.pathway_names", NULL))
pathways <- selection$pathways
pathway_lookup <- ibd_pathway_lookup(membership, pathways)
files <- c(files, loaded$input_file)
hashes <- tools::md5sum(files)
database_label <- if (collection == "H") "Hallmark" else paste(c(collection, subcollection), collapse = ":")
prefix <- if (collection == "H" && is.null(subcollection)) "hallmark" else
  tolower(gsub("[^A-Za-z0-9]+", "_", paste(c(collection, subcollection), collapse = "_")))
settings <- list(cooccur = 10L, q_value = 0.05, empirical = FALSE, filter_threshold = 1,
  n_resamples = 100L, seed = 2026L, workers = getOption("ibd.disco_workers", 6L),
  preprocess_only = isTRUE(getOption("ibd.pathway_preprocess_only", FALSE)),
  comparisons = comparisons, batch = NULL, covariates = c("age_at_diagnosis", "sex"),
  preprocessing_mode = preprocessing_mode,
  preprocessing_groups = if (preprocessing_mode == "joint_four_groups") c("Non-IBD control", "iCD", "cCD", "UC") else NULL,
  variance_threshold = 1e-10,
  preserve = if (preprocessing_mode == "joint_four_groups") "four disjoint control/iCD/cCD/UC groups" else "non-IBD control/disease group",
  cohort = "SRP042228 / GSE57945 diagnostic ileal biopsies",
  control_group = "Non-IBD control", batch_policy = "No batch correction or matching; identifiers unavailable",
  collection = collection, subcollection = subcollection, pathway_database = database_label,
  min_pathway_genes = getOption("ibd.min_pathway_genes", 0L),
  pathway_size_basis = "unique native human gene symbols before expression filtering",
  pathway_file_layout = "msigdb_id_v1")
output <- getOption("ibd.output_dir", file.path(root, "results"))
if (settings$preprocess_only) output <- file.path(output, "preparation")
ibd_reset_results(output)
write_table <- function(x, file, directory = output) utils::write.csv(x, file.path(directory, paste0(file, ".csv")), row.names = FALSE)
saveRDS(membership, file.path(output, "pathway_membership.rds"))
write_table(selection$sizes, "pathway_size_selection")
write_table(pathway_lookup, "pathway_lookup")
saveRDS(list(settings = settings, pathways = pathways, input_files = files, input_md5 = hashes,
  source_md5 = tools::md5sum(c(script, helpers)),
  disco_version = as.character(utils::packageVersion("DISCO")),
  msigdb_version = unique(membership$db_version), msigdbr_version = loaded$package_version,
  membership_md5 = tools::md5sum(file.path(output, "pathway_membership.rds")),
  pathway_lookup_md5 = tools::md5sum(file.path(output, "pathway_lookup.csv"))), file.path(output, "provenance.rds"))
message("Output: ", output)
message("Database: ", database_label, "; pathways: ", length(pathways),
  "; minimum database genes: ", settings$min_pathway_genes)
warning("Exploratory inference: estimated age/sex correction is not fully accounted for by DISCO's ordinary Fisher tests or label permutations. Controls are reused across comparisons; CD includes iCD and cCD.", call. = FALSE)
exclusion_fun <- if (collection == "H" && is.null(subcollection)) ibd_hallmark_exclusions else ibd_pathway_exclusions
# Representatives use the whole selected collection, independently of the pathway-size cutoff.
exclusions <- exclusion_fun(read.csv(redundancy_file), unique(membership$ensembl_gene))
write_table(exclusions, "redundancy_exclusions")
loaded_counts <- ibd_load_count_inputs(input)
counts <- loaded_counts$counts
z <- loaded_counts$metadata
all_gene_ids <- rownames(counts)
stopifnot(all(exclusions$excluded_gene %in% all_gene_ids), all(exclusions$retained_gene %in% all_gene_ids))
counts <- counts[!rownames(counts) %in% exclusions$excluded_gene, , drop = FALSE]
rm(loaded_counts); gc()
joint <- NULL
if (settings$preprocessing_mode == "joint_four_groups") {
  joint <- ibd_prepare_joint_cohort(counts, z, settings$cooccur)
  joint_path <- file.path(output, "joint_preprocessing")
  dir.create(joint_path)
  saveRDS(joint, file.path(joint_path, "prepared_expression.rds"))
  write_table(joint$cohort_audit, "sample_selection_audit", joint_path)
  write_table(joint$metadata, "sample_metadata", joint_path)
  write_table(joint$filter, "gene_expression_filter", joint_path)
  write_table(joint$variances, "adjusted_gene_variance", joint_path)
  write_table(joint$libraries, "library_normalization", joint_path)
  write_table(data.frame(sample_id = joint$metadata$sample_id, joint$adjustment$design_matrix), "adjustment_design", joint_path)
  write_table(data.frame(gene_id = rownames(joint$expression)), "common_analysis_genes", joint_path)
  write_table(data.frame(n_input_genes = length(all_gene_ids), n_after_redundancy = nrow(counts),
    n_after_expression_filter = sum(joint$filter$retained), n_after_variance_filter = nrow(joint$expression),
    n_samples = nrow(joint$metadata), adjustment_residual_df = joint$adjustment$residual_df),
    "filtering_summary", joint_path)
}
summaries <- scores <- cohorts <- list()
for (comparison in settings$comparisons) {
  path <- file.path(output, paste0("comparison_", comparison, "_vs_control"))
  dir.create(path)
  prepared <- if (is.null(joint)) ibd_prepare_pathway_comparison(counts, z, comparison, settings$cooccur) else
    ibd_subset_joint_comparison(joint, comparison, settings$cooccur)
  saveRDS(prepared, file.path(path, "prepared_expression.rds"))
  write_table(prepared$cohort_audit, "sample_selection_audit", path)
  write_table(data.frame(sample_id = prepared$metadata$sample_id, prepared$adjustment$design_matrix),
    "adjustment_design", path)
  write_table(prepared$metadata, "sample_metadata", path)
  write_table(prepared$filter, "gene_expression_filter", path)
  write_table(prepared$variances, "adjusted_gene_variance", path)
  write_table(prepared$libraries, "library_normalization", path)
  mm <- prepared$metadata
  cohorts[[comparison]] <- data.frame(comparison,
    n_control = sum(mm$analysis_group == "Control"), n_disease = sum(mm$analysis_group == "Disease"),
    n_paired_patients = length(intersect(mm$patient_id[mm$analysis_group == "Control"], mm$patient_id[mm$analysis_group == "Disease"])),
    n_expressed_genes = sum(prepared$filter$retained), n_variable_genes = nrow(prepared$expression),
    preprocessing_mode = settings$preprocessing_mode,
    adjustment_fitted_n_samples = if (is.null(joint)) nrow(mm) else nrow(joint$metadata),
    adjustment_residual_df = prepared$adjustment$residual_df)
  write_table(do.call(rbind, cohorts), "comparison_sample_counts")
  for (pathway in pathways) {
    message("Comparison ", comparison, ": ", pathway)
    pathway_id <- pathway_lookup$pathway_id[match(pathway, pathway_lookup$pathway)]
    audit <- ibd_pathway_gene_audit(membership, pathway, all_gene_ids, exclusions, prepared)
    result <- ibd_run_pathway_panel(prepared, audit, comparison, pathway, settings,
                                      file.path(path, pathway_id))
    result$summary$pathway_id <- pathway_id
    result$scores$pathway_id <- pathway_id
    summaries[[length(summaries) + 1L]] <- result$summary
    scores[[length(scores) + 1L]] <- result$scores
    write_table(do.call(rbind, summaries), "pathway_analysis_summary")
    combined <- do.call(rbind, scores)
    combined$q_all_gene_pathway_comparison_tests <- stats::p.adjust(combined$p, "BY")
    write_table(combined, "all_pathway_gene_scores")
  }
}
stopifnot(identical(hashes, tools::md5sum(files)))
stopifnot(identical(readRDS(file.path(output, "pathway_membership.rds")), membership))
capture.output(sessionInfo(), file = file.path(output, "session_info.txt"))
saveRDS(list(status = "finished", completed_at = Sys.time(),
  n_analyses = length(summaries), n_failed = sum(do.call(rbind, summaries)$status == "failed")),
  file.path(output, "completion.rds"))
print(do.call(rbind, cohorts), row.names = FALSE)
if (any(do.call(rbind, summaries)$status == "failed")) warning("Some pathways failed; inspect pathway_analysis_summary.csv.")
message("Finished: ", output, if (settings$preprocess_only) " (preparation only; DISCO not run)" else "")
