# Annotation-only audit. Neither expression nor disease labels determine flags.
# Coordinates are 1-based closed genomic intervals (GRanges/GTF convention).
ibd_annotation_redundancy <- function(genes, exons, reciprocal_cutoff = 0.95) {
  stopifnot(length(reciprocal_cutoff) == 1L, is.finite(reciprocal_cutoff),
            reciprocal_cutoff > 0, reciprocal_cutoff <= 1)
  ids <- names(genes)
  if (is.null(ids) || anyNA(ids) || anyDuplicated(ids)) stop("Unique gene IDs required.")
  exon_ids <- as.character(S4Vectors::mcols(exons)$gene_id)
  if (length(exon_ids) != length(exons)) stop("Exons require gene_id metadata.")
  keep <- !is.na(exon_ids) & exon_ids %in% ids
  exons <- exons[keep]
  exon_ids <- exon_ids[keep]
  # Check that exon IDs refer to the same chromosome, strand and gene span.
  matched <- genes[match(exon_ids, ids)]
  consistent <- as.character(GenomicRanges::seqnames(exons)) ==
    as.character(GenomicRanges::seqnames(matched)) &
    as.character(GenomicRanges::strand(exons)) == as.character(GenomicRanges::strand(matched)) &
    BiocGenerics::start(exons) >= BiocGenerics::start(matched) &
    BiocGenerics::end(exons) <= BiocGenerics::end(matched)
  if (!all(consistent)) stop("Exon coordinates disagree with saved gene annotation.")
  # Merge exons across ALL transcripts of each gene to avoid counting shared
  # transcript exons more than once. Matching gene spans alone is insufficient.
  unions <- GenomicRanges::reduce(S4Vectors::splitAsList(exons, exon_ids), ignore.strand = FALSE)
  union_bp <- sum(BiocGenerics::width(unions))
  lengths <- as.numeric(union_bp[match(ids, names(union_bp))])
  names(lengths) <- ids
  hits <- GenomicRanges::findOverlaps(genes, genes, ignore.strand = TRUE)
  i <- S4Vectors::queryHits(hits)
  j <- S4Vectors::subjectHits(hits)
  keep <- i < j
  i <- i[keep]; j <- j[keep]
  keys <- paste(i, j, sep = ":")
  shared <- numeric(length(i))
  flat <- unlist(unions, use.names = FALSE)
  gene_index <- rep(match(names(unions), ids), S4Vectors::elementNROWS(unions))
  exon_hits <- GenomicRanges::findOverlaps(flat, flat, ignore.strand = TRUE)
  a <- S4Vectors::queryHits(exon_hits)
  b <- S4Vectors::subjectHits(exon_hits)
  keep <- gene_index[a] < gene_index[b]
  a <- a[keep]; b <- b[keep]
  if (length(a)) {
    overlap_bp <- pmin(BiocGenerics::end(flat)[a], BiocGenerics::end(flat)[b]) -
      pmax(BiocGenerics::start(flat)[a], BiocGenerics::start(flat)[b]) + 1
    overlap_keys <- paste(gene_index[a], gene_index[b], sep = ":")
    summed <- rowsum(overlap_bp, overlap_keys, reorder = FALSE)
    location <- match(rownames(summed), keys)
    stopifnot(!anyNA(location))
    shared[location] <- summed[, 1]
  }
  assessed <- is.finite(lengths[i]) & is.finite(lengths[j])
  shared[!assessed] <- NA_real_
  strand <- as.character(GenomicRanges::strand(genes))
  same_strand <- strand[i] == strand[j] & strand[i] %in% c("+", "-")
  fraction1 <- shared / lengths[i]
  fraction2 <- shared / lengths[j]
  exact <- assessed & same_strand & shared == lengths[i] & shared == lengths[j]
  high <- assessed & same_strand & fraction1 >= reciprocal_cutoff & fraction2 >= reciprocal_cutoff
  exact[is.na(exact)] <- FALSE
  high[is.na(high)] <- FALSE
  classification <- rep("partial_exonic_overlap", length(i))
  classification[assessed & shared == 0] <- "gene_span_only_no_shared_exons"
  containment <- assessed & same_strand &
    pmax(fraction1, fraction2) >= reciprocal_cutoff & !high
  classification[which(containment)] <- "one_sided_exonic_containment"
  classification[!same_strand] <- "opposite_or_unknown_strand_overlap"
  classification[!assessed] <- "exon_annotation_missing"
  classification[high] <- "high_reciprocal_exonic_overlap"
  classification[exact] <- "identical_exon_union"
  annotation <- S4Vectors::mcols(genes)
  field <- function(name) {
    if (name %in% names(annotation)) as.character(annotation[[name]]) else rep(NA_character_, length(ids))
  }
  pairs <- data.frame(
    gene1 = ids[i], gene2 = ids[j], name1 = field("gene_name")[i],
    name2 = field("gene_name")[j], type1 = field("gene_type")[i], type2 = field("gene_type")[j],
    chromosome = as.character(GenomicRanges::seqnames(genes))[i],
    strand1 = strand[i], strand2 = strand[j], same_strand,
    start1 = BiocGenerics::start(genes)[i], end1 = BiocGenerics::end(genes)[i],
    start2 = BiocGenerics::start(genes)[j], end2 = BiocGenerics::end(genes)[j],
    exon_bp1 = lengths[i], exon_bp2 = lengths[j], shared_exon_bp = shared,
    overlap_fraction1 = fraction1, overlap_fraction2 = fraction2,
    identical_gene_span = BiocGenerics::start(genes)[i] == BiocGenerics::start(genes)[j] &
      BiocGenerics::end(genes)[i] == BiocGenerics::end(genes)[j],
    classification, redundancy_candidate = high, row.names = NULL
  )
  gene_summary <- data.frame(gene_id = ids, gene_name = field("gene_name"),
    gene_type = field("gene_type"), exon_union_bp = unname(lengths),
    exon_annotation_available = is.finite(lengths),
    n_candidate_partners = tabulate(c(i[high], j[high]), nbins = length(ids)))
  list(pairs = pairs, genes = gene_summary, exon_unions = unions)
}
