# Run from the DISCO repository root after downloading the counts (see README).
# Adapted from DISCO_Analysis's example_ibd_workflow.R; only the Hallmark
# significance bar plot is generated here. This runs the full analysis.
# Each run replaces the previous analysis in results/ and the bar plot in results/figure/.
# Script order: (1) fetch_srp042228_metadata.R, (2) explore_srp042228_metadata.R,
# (3) fetch_srp042228_counts.R, (4) run_ibd_hallmark_disco.R,
# (5) plot_ibd_hallmark_significance.R. Steps 1-2 can be skipped with bundled
# metadata. This wrapper configures and runs steps 4-5; helpers load automatically.

## ---- ibd-workflow
library(DISCO)

options(
  ibd.exploration_dir = "analysis/ibd_transcriptomics/results/inputs",
  ibd.pathway_reference = "analysis/ibd_transcriptomics/results/inputs/reference",
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
# To plot existing results without rerunning DISCO, use instead:
# ibd_results_dir <- "analysis/ibd_transcriptomics/results"
options(
  ibd.plot_input_dir = ibd_results_dir,
  ibd.pathway_min_percent = 50,
  ibd.pathway_label_max_chars = 60L
)
source("analysis/ibd_transcriptomics/code/plot_ibd_hallmark_significance.R")
