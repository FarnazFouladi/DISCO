# Maintainer command: rerun both examples, replace their outputs, and render README.md.
# Run from the repository root: Rscript analysis/rebuild_readme.R
if (!file.exists("DESCRIPTION") || !dir.exists("analysis"))
  stop("Run this script from the DISCO repository root.")

scripts <- c(
  "analysis/ibd_transcriptomics/code/fetch_srp042228_counts.R",
  "analysis/diabetes_microbiome/code/run_qinj_readme_examples.R",
  "analysis/ibd_transcriptomics/code/example_ibd_workflow.R",
  "analysis/render_readme.R"
)
for (script in scripts) {
  message("Running ", script)
  args <- c("--vanilla", shQuote(script))
  if (basename(script) == "fetch_srp042228_counts.R")
    args <- c(args, shQuote("analysis/ibd_transcriptomics/results/inputs"))
  status <- system2(file.path(R.home("bin"), "Rscript"), args)
  if (status != 0L) stop("Stopped because this script failed: ", script, call. = FALSE)
}
