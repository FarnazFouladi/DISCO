# Run from the DISCO repository root after downloading the counts (see README).
# Adapted from DISCO_Analysis's example_ibd_workflow.R; only the Hallmark
# significance heatmap is generated here. This runs the full analysis.

## ---- ibd-workflow
library(DISCO)

options(
  ibd.exploration_dir = "analysis/ibd_transcriptomics/inputs",
  ibd.pathway_reference = "analysis/ibd_transcriptomics/inputs/reference",
  ibd.comparisons = c("iCD", "cCD", "UC"),
  ibd.preprocessing_mode = "joint_four_groups",
  ibd.msigdb_collection = "H",
  ibd.msigdb_subcollection = NULL,
  ibd.min_pathway_genes = 50L,
  ibd.pathway_names = NULL,
  ibd.pathway_preprocess_only = FALSE,
  ibd.disco_workers = 1L
)
source("analysis/ibd_transcriptomics/code/run_ibd_hallmark_disco.R")
ibd_results_dir <- output

## ---- ibd-plot
# To plot the bundled results without rerunning DISCO, use instead:
# ibd_results_dir <- "analysis/ibd_transcriptomics/results"
options(
  ibd.plot_input_dir = ibd_results_dir,
  ibd.pathway_min_percent = 50,
  ibd.pathway_label_max_chars = 60L,
  ibd.pathway_cluster_rows = TRUE
)
source("analysis/ibd_transcriptomics/code/plot_ibd_hallmark_significance.R")
