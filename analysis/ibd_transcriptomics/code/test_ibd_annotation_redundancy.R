# Offline interval-logic checks; Rscript this file from the repository root.
library(GenomicRanges)
source("analysis/ibd_transcriptomics/code/ibd_annotation_redundancy_helpers.R")
genes <- GRanges(rep("chr1", 6), IRanges(c(1,1,1,1,201,201), c(100,100,200,100,300,300)),
  strand = c("+", "+", "+", "-", "+", "+"))
names(genes) <- LETTERS[1:6]
# Repeated/overlapping transcript exons must not inflate the exon union.
exons <- GRanges(rep("chr1", 8), IRanges(c(1,41,1,1,1,201,201,201), c(60,100,100,200,100,300,295,295)),
  strand = c("+", "+", "+", "+", "-", "+", "+", "+"),
  gene_id = c("A","A","B","C","D","E","F","F"))
audit <- ibd_annotation_redundancy(genes, exons)
pair <- function(a,b) audit$pairs[audit$pairs$gene1 == a & audit$pairs$gene2 == b, ]
stopifnot(pair("A","B")$redundancy_candidate,
  pair("A","B")$exon_bp1 == 100, pair("A","B")$shared_exon_bp == 100,
  !pair("A","C")$redundancy_candidate, pair("A","C")$overlap_fraction2 == .5,
  !pair("A","D")$redundancy_candidate,
  pair("E","F")$redundancy_candidate, pair("E","F")$overlap_fraction1 == .95)
stopifnot(!ibd_annotation_redundancy(genes, exons, .96)$pairs$redundancy_candidate[
  audit$pairs$gene1 == "E" & audit$pairs$gene2 == "F"])
bad <- exons
start(bad)[1] <- 0
stopifnot(inherits(try(ibd_annotation_redundancy(genes, bad), silent = TRUE), "try-error"))
message("PASS: exon unions, exact and 95% reciprocal overlap, one-sided containment, opposite strands, and coordinate validation.")
