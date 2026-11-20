# Standalone annotation audit of genes in SRP042228's GENCODE v26 count file.
# Run from the repository root with Rscript, or source() after setting options:
# ibd.exploration_dir, ibd.pathway_reference, ibd.annotation_gtf (optional local GTF).
# Requires the IBD coverage-sum download; no expression or disease labels enter flags.
script <- local({
  paths <- unlist(lapply(sys.frames(), function(x) x$ofile), use.names = FALSE)
  paths <- c(rev(paths), sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)),
    "analysis/ibd_transcriptomics/code/audit_ibd_annotation_redundancy.R")
  paths <- paths[basename(paths) == "audit_ibd_annotation_redundancy.R" & file.exists(paths)]
  if (!length(paths)) stop("Use Rscript or source() with this script's path.")
  normalizePath(paths[1])
})
root <- dirname(dirname(script))
for (pkg in c("data.table", "GenomicRanges", "rtracklayer", "curl"))
  if (!requireNamespace(pkg, quietly = TRUE)) stop("Required audit dependency: ", pkg)
helper <- file.path(dirname(script), "ibd_annotation_redundancy_helpers.R")
source(helper)
input <- getOption("ibd.exploration_dir", file.path(root, "results", "inputs"))
reference <- getOption("ibd.pathway_reference", file.path(input, "reference"))
counts_file <- file.path(input, "raw", "sra.gene_sums.SRP042228.G026.gz")
if (!file.exists(counts_file)) stop("Download SRP042228 counts first: ", counts_file)
counts_hash <- tools::md5sum(counts_file)
connection <- gzfile(counts_file, "rt")
header <- readLines(connection, n = 1L)
close(connection)
if (!identical(header, "##annotation=G026")) stop("Expected the GENCODE v26 count annotation.")
# Read only IDs, not expression values, from the actual IBD count file.
ids <- data.table::fread(counts_file, skip = 2L, select = "gene_id",
  colClasses = "character", nThread = 1L)[[1L]]
stopifnot(length(ids) > 0L, !anyNA(ids), !anyDuplicated(ids), all(nzchar(ids)))
dir.create(reference, recursive = TRUE, showWarnings = FALSE)
gtf_url <- paste0("https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/",
  "release_26/gencode.v26.chr_patch_hapl_scaff.annotation.gtf.gz")
gtf <- getOption("ibd.annotation_gtf", file.path(reference, "gencode.v26.chr_patch_hapl_scaff.annotation.gtf.gz"))
if (!file.exists(gtf)) {
  message("Downloading public GENCODE v26 ALL annotation.")
  temporary <- paste0(gtf, ".partial")
  curl::curl_download(gtf_url, temporary, quiet = FALSE)
  if (!file.rename(temporary, gtf)) stop("Could not finalize GTF download.")
}
gtf <- normalizePath(gtf, mustWork = TRUE)
gtf_hash <- tools::md5sum(gtf)
message("Importing GENCODE v26 gene/exon coordinates for ", length(ids), " IBD genes.")
annotation <- rtracklayer::import(gtf, format = "gtf", feature.type = c("gene", "exon"))
genes <- annotation[annotation$type == "gene"]
gtf_ids <- as.character(genes$gene_id)
stopifnot(!anyNA(gtf_ids), !anyDuplicated(gtf_ids))
locations <- match(ids, gtf_ids)
if (anyNA(locations)) stop("IBD count IDs missing from GTF: ", paste(head(ids[is.na(locations)]), collapse = ", "))
genes <- genes[locations]
names(genes) <- ids
exons <- annotation[annotation$type == "exon"]
rm(annotation)
cutoff <- getOption("ibd.annotation_overlap_cutoff", .95)
audit <- ibd_annotation_redundancy(genes, exons, cutoff)
if (!all(audit$genes$exon_annotation_available)) stop("Some IBD genes have no exon annotation; audit incomplete.")
candidate <- audit$pairs[audit$pairs$redundancy_candidate, , drop = FALSE]
output <- file.path(reference, "annotation_audit")
dir.create(output, recursive = TRUE, showWarnings = FALSE)
write_table <- function(x, name) utils::write.csv(x, file.path(output, name), row.names = FALSE)
write_table(data.frame(gene_id = ids), "audited_gene_ids.csv")
write_table(audit$genes, "gene_annotation_coverage.csv")
write_table(audit$pairs, "all_overlapping_gene_pairs.csv")
write_table(as.data.frame(table(audit$pairs$classification)), "classification_counts.csv")
write_table(data.frame(project = "SRP042228", annotation = "GENCODE v26 ALL",
  n_input_genes = length(ids), n_genes_with_exons = sum(audit$genes$exon_annotation_available),
  n_overlapping_pairs = nrow(audit$pairs), reciprocal_overlap_cutoff = cutoff,
  n_redundancy_candidates = nrow(candidate)), "summary.csv")
candidate_file <- file.path(reference, "redundancy_candidates_all_genes.csv")
temporary <- paste0(candidate_file, ".partial")
utils::write.csv(candidate, temporary, row.names = FALSE)
stopifnot(identical(counts_hash, tools::md5sum(counts_file)), identical(gtf_hash, tools::md5sum(gtf)))
if (!file.rename(temporary, candidate_file)) stop("Could not replace candidate CSV.")
source_files <- c(script, helper)
saveRDS(list(project = "SRP042228", annotation_release = "GENCODE v26 ALL",
  reciprocal_overlap_cutoff = cutoff, n_input_genes = length(ids), n_candidates = nrow(candidate),
  counts_file = normalizePath(counts_file), counts_md5 = unname(counts_hash),
  source_url = gtf_url, gtf = gtf, gtf_md5 = unname(gtf_hash),
  candidate_md5 = unname(tools::md5sum(candidate_file)), source_md5 = tools::md5sum(source_files),
  audited_at = Sys.time()), file.path(output, "provenance.rds"))
checks <- c(counts_file, gtf, source_files, candidate_file)
write_table(data.frame(file = checks, md5 = unname(tools::md5sum(checks))), "source_checksums.csv")
capture.output(sessionInfo(), file = file.path(output, "session_info.txt"))
message("IBD audit complete: ", length(ids), " genes, ", nrow(candidate), " candidate pairs. Saved: ", candidate_file)
