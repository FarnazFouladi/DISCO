# Rscript analysis/ibd_transcriptomics/code/fetch_srp042228_counts.R METADATA_DIRECTORY
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Supply metadata directory.")
raw <- file.path(args[1], "raw")
if (!dir.exists(raw)) stop("Run fetch_srp042228_metadata.R first.")
for (pkg in c("recount3", "curl")) if (!requireNamespace(pkg, quietly = TRUE)) stop("Required: ", pkg)
url <- recount3::locate_url("SRP042228", "data_sources/sra", type = "gene", organism = "human",
  annotation = "gencode_v26", recount3_url = "https://recount-opendata.s3.amazonaws.com/recount3/release")
stopifnot(length(url) == 1L)
dest <- file.path(raw, basename(url))
if (!file.exists(dest)) {
  temporary <- paste0(dest, ".partial")
  curl::curl_download(url, temporary, quiet = FALSE)
  if (!file.rename(temporary, dest)) stop("Could not finalize count download.")
}
utils::write.csv(data.frame(file = basename(dest), url = unname(url), md5 = unname(tools::md5sum(dest)),
  downloaded_at = as.character(Sys.time()), annotation = "gencode_v26"),
  file.path(args[1], "gene_counts_manifest.csv"), row.names = FALSE)
message("Count file: ", dest)
