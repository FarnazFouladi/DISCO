# Refresh Excel exports from the saved examples without rerunning statistics.
# Run from the repository root, or source from analysis/render_readme.R.
source("analysis/format_saved_results.R")
source("R/output.R")
for (mode in c("taxa", "taxa_pathways")) {
  directory <- file.path("analysis/diabetes_microbiome/results", mode)
  result <- read_disco_result(file.path(directory, "disco_analysis.rds"))
  prefix <- paste0("QinJ_2012_", mode)
  tables <- c("cor_diff", "disco", "threshold_summary", "global_test", "feature_null_summary")
  suffixes <- c("cor_diff", "disco", "threshold_summary", "global_test",
                "real_p_vs_null_discovery_by_threshold")
  for (i in seq_along(tables)) {
    if (!is.null(result[[tables[i]]]))
      write_result_xlsx(result[[tables[i]]],
                         file.path(directory, paste0(prefix, "_", suffixes[i], ".xlsx")))
  }
}
