# Validate the minimum co-occurrence threshold
validate_cooccur <- function(cooccur) {
  if (!is.numeric(cooccur) || length(cooccur) != 1 || is.na(cooccur) ||
      !is.finite(cooccur) || cooccur < 10 || cooccur != floor(cooccur)) {
    stop("cooccur must be a single integer greater than or equal to 10",
         call. = FALSE)
  }

  cooccur
}

# Report the sample overlap between two modalities within a group
report_sample_overlap <- function(modality1, modality2, group) {
  n_common <- length(intersect(rownames(modality1), rownames(modality2)))

  if (n_common == 0) {
    stop(sprintf("%s modalities have no common samples", group),
         call. = FALSE)
  }

  message(sprintf("%s modalities have %d common sample(s)", group, n_common))
  invisible(n_common)
}

# Align two groups to the same features within a modality
align_common_features <- function(data1, data2, modality) {
  features1 <- colnames(data1)
  features2 <- colnames(data2)

  if (is.null(features1) || is.null(features2) ||
      anyNA(features1) || anyNA(features2) ||
      any(!nzchar(features1)) || any(!nzchar(features2))) {
    stop(
      sprintf("Modality %d must have non-missing, non-empty feature names", modality),
      call. = FALSE
    )
  }

  if (anyDuplicated(features1)) {
    stop(sprintf("Feature names in data1 modality %d must be unique", modality),
         call. = FALSE)
  }

  if (anyDuplicated(features2)) {
    stop(sprintf("Feature names in data2 modality %d must be unique", modality),
         call. = FALSE)
  }

  # Preserve the feature order supplied in data1 and apply it to both groups.
  common_features <- features1[features1 %in% features2]

  if (length(common_features) == 0) {
    stop(sprintf("No common features were found for modality %d", modality),
         call. = FALSE)
  }

  dropped1 <- sum(!features1 %in% common_features)
  dropped2 <- sum(!features2 %in% common_features)
  features2_common <- features2[features2 %in% common_features]
  reordered2 <- !identical(features2_common, common_features)

  if (dropped1 > 0 || dropped2 > 0 || reordered2) {
    warning(
      sprintf(
        paste0(
          "Modality %d was aligned to %d common feature(s): ",
          "dropped %d feature(s) from data1 and %d from data2%s"
        ),
        modality,
        length(common_features),
        dropped1,
        dropped2,
        if (reordered2) "; reordered data2 to match data1" else ""
      ),
      call. = FALSE
    )
  }

  data1_aligned <- data1[, common_features, drop = FALSE]
  data2_aligned <- data2[, common_features, drop = FALSE]

  if (!identical(colnames(data1_aligned), colnames(data2_aligned))) {
    stop(sprintf("Failed to align feature order for modality %d", modality),
         call. = FALSE)
  }

  list(data1 = data1_aligned, data2 = data2_aligned)
}
