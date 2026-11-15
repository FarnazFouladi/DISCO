# Rscript analysis/ibd_transcriptomics/code/fetch_srp042228_metadata.R
# Download recount3, GEO and ENA metadata, then build fresh sample_metadata.csv.
# No bundled clinical metadata, expression matrix, FASTQ or article XML is used.
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
# These two tables supply sample/run identifiers and average mapped read lengths.
# Other recount3 metadata tables are not used by ibd_load_count_inputs().
file_names <- c("sra.sra.SRP042228.MD.gz", "sra.recount_qc.SRP042228.MD.gz")
base_url <- paste0("https://recount-opendata.s3.amazonaws.com/recount3/release/",
                   "human/data_sources/sra/metadata/28/SRP042228")
urls <- stats::setNames(paste0(base_url, "/", file_names), file_names)
files <- file.path(raw_dir, basename(urls))
names(files) <- names(urls)
for (i in seq_along(urls)) {
  if (!file.exists(files[i])) curl::curl_download(urls[i], files[i], quiet = FALSE)
}
metadata <- recount3::read_metadata(files)
saveRDS(metadata, file.path(raw_dir, "recount_metadata.rds"))
# GEO provides clinical annotations; ENA provides the sample/run cross-check.
fields <- c("study_accession", "sample_accession", "secondary_sample_accession", "experiment_accession",
  "run_accession", "sample_alias", "sample_title", "experiment_alias", "experiment_title",
  "library_name", "library_strategy", "library_source", "library_selection", "library_layout",
  "instrument_model", "instrument_platform", "read_count", "base_count", "first_public", "last_updated")
clinical_urls <- c(
  "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE57nnn/GSE57945/soft/GSE57945_family.soft.gz",
  paste0("https://www.ebi.ac.uk/ena/portal/api/filereport?accession=SRP042228&result=read_run&fields=",
         paste(fields, collapse = ","), "&format=tsv"))
clinical_files <- file.path(raw_dir, c("GSE57945_family.soft.gz", "ena_runs.tsv"))
for (i in seq_along(clinical_urls)) {
  if (!file.exists(clinical_files[i]))
    curl::curl_download(clinical_urls[i], clinical_files[i], quiet = FALSE)
}

# Always regenerate clinical metadata from these sources, even if a previous
# sample_metadata.csv exists. The audit checks GEO/ENA/recount3 correspondence.
audit_script <- file.path(dirname(script), "explore_srp042228_metadata.R")
audit_env <- new.env(parent = environment())
sys.source(audit_script, envir = audit_env)
audit_env$audit_srp042228(output)
sample_file <- file.path(output, "sample_metadata.csv")

all_files <- c(files, clinical_files)
manifest <- data.frame(file = basename(all_files), url = c(unname(urls), clinical_urls),
  available = file.exists(all_files), md5 = unname(tools::md5sum(all_files)))
utils::write.csv(manifest, file.path(output, "source_manifest.csv"), row.names = FALSE)
saveRDS(list(study = "SRP042228", geo = "GSE57945", n_recount_rows = nrow(metadata),
  source_md5 = tools::md5sum(c(script, audit_script)),
  sample_metadata_md5 = tools::md5sum(sample_file),
  clinical_metadata_source = "Regenerated from GEO and audited against ENA/recount3",
  retrieved_at = Sys.time(), expression_downloaded = FALSE),
  file.path(output, "fetch_provenance.rds"))
capture.output(sessionInfo(), file = file.path(output, "fetch_session_info.txt"))
message("Finished metadata download and audit: ", output)
