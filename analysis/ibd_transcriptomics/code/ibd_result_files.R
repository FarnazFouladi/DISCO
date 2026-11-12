# Replace generated analysis outputs while preserving downloaded inputs.
ibd_reset_results <- function(output) {
  dir.create(output, recursive = TRUE, showWarnings = FALSE)
  generated <- c(
    paste0("comparison_", c("CD", "iCD", "cCD", "UC"), "_vs_control"),
    "joint_preprocessing", "figure", "pathway_membership.rds",
    "pathway_size_selection.csv", "pathway_lookup.csv", "provenance.rds",
    "redundancy_exclusions.csv", "comparison_sample_counts.csv",
    "pathway_analysis_summary.csv", "all_pathway_gene_scores.csv",
    "session_info.txt", "completion.rds"
  )
  files <- file.path(output, generated)
  unlink(files, recursive = TRUE)
  if (any(file.exists(files))) stop("Could not replace previous analysis outputs.")
  invisible(output)
}
