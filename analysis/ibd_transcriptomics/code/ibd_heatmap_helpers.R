# Selection uses full pathway names as keys; shortened labels are display-only.
ibd_prepare_heatmap_rows <- function(summary, lookup, comparisons, min_percent = 0, label_max_chars = 60L) {
  if (!is.numeric(min_percent) || length(min_percent) != 1L || !is.finite(min_percent) ||
      min_percent < 0 || min_percent > 100)
    stop("ibd.pathway_min_percent must be a number between 0 and 100.")
  if (!is.numeric(label_max_chars) || length(label_max_chars) != 1L || !is.finite(label_max_chars) ||
      label_max_chars < 10 || label_max_chars != floor(label_max_chars))
    stop("ibd.pathway_label_max_chars must be an integer of at least 10.")
  pathways <- sort(unique(summary$pathway))
  key <- function(x) paste(x$pathway, x$comparison, sep = "::")
  expected <- expand.grid(pathway = pathways, comparison = comparisons)
  if (anyDuplicated(key(summary)) || !setequal(key(summary), key(expected)))
    stop("Heatmap requires one row per pathway and comparison.")
  percentage <- tested <- matrix(NA_real_, nrow = length(pathways), ncol = length(comparisons),
    dimnames = list(pathways, comparisons))
  idx <- cbind(match(summary$pathway, pathways), match(as.character(summary$comparison), comparisons))
  percentage[idx] <- summary$percent_significant
  tested[idx] <- summary$n_tested
  if (any(!is.finite(tested)) || any(tested < 0) || any(tested != floor(tested)))
    stop("Invalid tested-gene counts.")
  maximum <- apply(percentage, 1, function(x) if (any(is.finite(x))) max(x[is.finite(x)]) else NA_real_)
  selected <- is.finite(maximum) & maximum >= min_percent
  if (!any(selected)) stop("No pathways meet ibd.pathway_min_percent.")
  ids <- lookup$pathway_id[match(pathways, lookup$pathway)]
  if (anyNA(ids)) stop("Missing pathway IDs for heatmap labels.")
  labels <- gsub("_", " ", sub("^(HALLMARK|REACTOME)_", "", pathways))
  too_long <- nchar(labels) > label_max_chars
  labels[too_long] <- paste0(substr(labels[too_long], 1, label_max_chars - 3L), "...")
  # Retain separate rows and disambiguate names that truncate to the same label.
  duplicated_label <- duplicated(labels) | duplicated(labels, fromLast = TRUE)
  labels[duplicated_label] <- paste0(labels[duplicated_label], " [", ids[duplicated_label], "]")
  rows <- data.frame(pathway = pathways, pathway_id = ids, pathway_label = labels,
    maximum_percent_significant = unname(maximum), selected = unname(selected),
    tested_genes = apply(tested, 1, function(x) paste0("(", paste(x, collapse = ", "), ")")),
    row.names = NULL)
  list(rows = rows, percentages = percentage[selected, , drop = FALSE],
    tested = tested[selected, , drop = FALSE])
}
