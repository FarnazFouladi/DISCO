# Rscript analysis/ibd_transcriptomics/code/fetch_srp042228_metadata.R
# Metadata only: no expression matrix, FASTQ files, normalization or DISCO run.
# Optional first argument resumes downloads in an existing output directory.
for (pkg in c("recount3", "curl"))
  if (!requireNamespace(pkg, quietly = TRUE)) stop("Required package: ", pkg)
script <- local({
  sourced <- unlist(lapply(sys.frames(), function(x) x$ofile), use.names = FALSE)
  cli <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  paths <- c(rev(sourced), cli, "analysis/ibd_transcriptomics/code/fetch_srp042228_metadata.R")
  paths <- paths[basename(paths) == "fetch_srp042228_metadata.R" & file.exists(paths)]
  if (!length(paths)) stop("Run from the repository root or source this saved script.")
  normalizePath(paths[1])
})
root <- dirname(dirname(script))
args <- commandArgs(trailingOnly = TRUE)
output <- if (length(args)) args[1] else getOption("ibd.metadata_output", file.path(root, "results", "inputs"))
raw_dir <- file.path(output, "raw")
dir.create(raw_dir, recursive = TRUE, showWarnings = FALSE)
message("Output: ", output)
url <- "https://recount-opendata.s3.amazonaws.com/recount3/release"
urls <- recount3::locate_url("SRP042228", "data_sources/sra", type = "metadata",
  organism = "human", recount3_url = url)
if (!length(urls)) stop("No recount3 metadata URLs for SRP042228.")
files <- file.path(raw_dir, basename(urls))
names(files) <- names(urls)
for (i in seq_along(urls)) {
  if (!file.exists(files[i])) curl::curl_download(urls[i], files[i], quiet = FALSE)
}
metadata <- recount3::read_metadata(files)
saveRDS(metadata, file.path(raw_dir, "recount_metadata.rds"))
fields <- c("study_accession", "sample_accession", "secondary_sample_accession", "experiment_accession",
  "run_accession", "sample_alias", "sample_title", "experiment_alias", "experiment_title",
  "library_name", "library_strategy", "library_source", "library_selection", "library_layout",
  "instrument_model", "instrument_platform", "read_count", "base_count", "first_public", "last_updated")
other_urls <- c(
  "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE57nnn/GSE57945/soft/GSE57945_family.soft.gz",
  paste0("https://www.ebi.ac.uk/ena/portal/api/filereport?accession=SRP042228&result=read_run&fields=",
    paste(fields, collapse = ","), "&format=tsv"),
  "https://www.ebi.ac.uk/europepmc/webservices/rest/PMC4109533/fullTextXML",
  "https://www.ebi.ac.uk/europepmc/webservices/rest/PMC4362234/fullTextXML")
other_files <- file.path(raw_dir, c("GSE57945_family.soft.gz", "ena_runs.tsv", "PMC4109533.xml", "PMC4362234_corrigendum.xml"))
for (i in seq_along(other_urls)) {
  if (!file.exists(other_files[i])) {
    if (i <= 2L) curl::curl_download(other_urls[i], other_files[i], quiet = FALSE) else
      tryCatch(curl::curl_download(other_urls[i], other_files[i], quiet = FALSE),
        error = function(e) warning("Optional manuscript XML unavailable: ", conditionMessage(e)))
  }
}
manifest <- data.frame(file = basename(c(files, other_files)), url = c(urls, other_urls),
  available = file.exists(c(files, other_files)), md5 = unname(tools::md5sum(c(files, other_files))))
utils::write.csv(manifest, file.path(output, "source_manifest.csv"), row.names = FALSE)
saveRDS(list(study = "SRP042228", geo = "GSE57945", n_recount_rows = nrow(metadata),
  source_md5 = tools::md5sum(script), retrieved_at = Sys.time(), expression_downloaded = FALSE),
  file.path(output, "fetch_provenance.rds"))
capture.output(sessionInfo(), file = file.path(output, "fetch_session_info.txt"))
message("Finished metadata download: ", output)
