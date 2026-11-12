# Keep full biological labels in tables; use short IDs only in filesystem paths.
ibd_validate_pathway_lookup <- function(lookup) {
  if (!all(c("pathway", "pathway_id") %in% names(lookup)) || !nrow(lookup) ||
      anyNA(lookup[c("pathway", "pathway_id")]) || any(!nzchar(lookup$pathway)) ||
      anyDuplicated(lookup$pathway) || anyDuplicated(lookup$pathway_id) ||
      any(!grepl("^[A-Za-z0-9][A-Za-z0-9_-]{0,39}$", lookup$pathway_id)))
    stop("Pathway lookup requires unique pathway names and unique safe IDs of at most 40 characters.")
  lookup
}

ibd_pathway_lookup <- function(membership, pathways) {
  if (!all(c("gs_name", "gs_id") %in% names(membership))) stop("MSigDB membership must include gs_id.")
  lookup <- unique(as.data.frame(membership[membership$gs_name %in% pathways, c("gs_name", "gs_id")]))
  names(lookup) <- c("pathway", "pathway_id")
  lookup <- ibd_validate_pathway_lookup(lookup)
  if (!setequal(pathways, lookup$pathway)) stop("Missing pathway IDs in membership.")
  lookup[match(pathways, lookup$pathway), , drop = FALSE]
}

ibd_read_pathway_lookup <- function(input, pathways) {
  file <- file.path(input, "pathway_lookup.csv")
  if (file.exists(file)) {
    lookup <- ibd_validate_pathway_lookup(utils::read.csv(file, stringsAsFactors = FALSE))
    if (!all(pathways %in% lookup$pathway)) stop("Pathway lookup does not cover the requested panels.")
    return(lookup[match(pathways, lookup$pathway), , drop = FALSE])
  }
  # Older runs keep their original input paths; their new plot files still get short IDs.
  snapshot <- file.path(input, "pathway_membership.rds")
  if (file.exists(snapshot)) return(ibd_pathway_lookup(readRDS(snapshot), pathways))
  ordered <- sort(unique(pathways))
  lookup <- data.frame(pathway = ordered, pathway_id = sprintf("P%05d", seq_along(ordered)))
  lookup[match(pathways, lookup$pathway), , drop = FALSE]
}

ibd_pair_result_path <- function(input, comparison, pathway, lookup) {
  comparison_dir <- file.path(input, paste0("comparison_", comparison, "_vs_control"))
  if (file.exists(file.path(input, "pathway_lookup.csv"))) {
    id <- lookup$pathway_id[match(pathway, lookup$pathway)]
    if (length(id) != 1L || is.na(id)) stop("Missing pathway in lookup: ", pathway)
    return(file.path(comparison_dir, id, "pathway_cor_diff.csv"))
  }
  file.path(comparison_dir, pathway, paste0(pathway, "_cor_diff.csv"))
}
