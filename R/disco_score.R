# Standard error for a Fisher-z-transformed Spearman correlation
spearman_fisher_se <- function(R, N,
                               correction = c("bonett_wright", "fieller_1.06")) {
  correction <- match.arg(correction)

  variance_multiplier <- switch(
    correction,
    bonett_wright = 1 + R^2 / 2,
    fieller_1.06 = 1.06
  )

  sqrt(variance_multiplier / pmax(N - 3, 1))
}

# Fisher's Z test for a single Spearman correlation
fisher_z_p <- function(R, N,
                       correction = c("bonett_wright", "fieller_1.06")) {
  Z <- 0.5 * log((1 + R) / (1 - R))

  # Correct for the bias related to sample size
  Z_bc <- Z - R / (2 * (N - 1))
  se <- spearman_fisher_se(R, N, correction = correction)

  2 * pnorm(-abs(Z_bc / se))
}

# Fisher's Z test for two independent Spearman correlations
fisher_z_p_2sample <- function(R1, N1, R2, N2,
                               correction = c("bonett_wright", "fieller_1.06")) {
  Z1 <- 0.5 * log((1 + R1) / (1 - R1))
  Z2 <- 0.5 * log((1 + R2) / (1 - R2))

  # Correct for different sample sizes between the two groups
  Z1_bc <- Z1 - R1 / (2 * (N1 - 1))
  Z2_bc <- Z2 - R2 / (2 * (N2 - 1))

  se1 <- spearman_fisher_se(R1, N1, correction = correction)
  se2 <- spearman_fisher_se(R2, N2, correction = correction)
  se <- sqrt(se1^2 + se2^2)
  Z <- (Z1_bc - Z2_bc) / se

  2 * pnorm(-abs(Z))
}

# Wrapper for Spearman correlation
get_corr <- function(D) {
  Hmisc::rcorr(as.matrix(D), type = "spearman")$r
}

# Highly optimized pairwise complete cases counter
get_npairs <- function(D) {
  crossprod(!is.na(D))
}

# Compare correlations only when both groups meet the co-occurrence threshold
compare_eligible_correlations <- function(R1, C1, R2, C2,
                                          cooccur, filt_threshold,
                                          correction = c("bonett_wright",
                                                         "fieller_1.06")) {
  correction <- match.arg(correction)
  eligible <- C1 >= cooccur & C2 >= cooccur
  diag(eligible) <- FALSE

  empty_matrix <- function(template) {
    matrix(
      NA_real_,
      nrow = nrow(template),
      ncol = ncol(template),
      dimnames = dimnames(template)
    )
  }

  P1 <- empty_matrix(R1)
  P2 <- empty_matrix(R2)
  pvals <- empty_matrix(R1)

  if (any(eligible)) {
    P1[eligible] <- fisher_z_p(
      R1[eligible], C1[eligible], correction = correction
    )
    P2[eligible] <- fisher_z_p(
      R2[eligible], C2[eligible], correction = correction
    )

    pvals[eligible] <- fisher_z_p_2sample(
      R1[eligible], C1[eligible], R2[eligible], C2[eligible],
      correction = correction
    )
  }

  # Retain eligible pairs correlated in at least one group.
  retained <- eligible & (P1 < filt_threshold | P2 < filt_threshold)
  retained[is.na(retained)] <- FALSE
  pvals[!retained] <- NA_real_

  list(
    pvals = pvals,
    eligible = eligible,
    retained = retained
  )
}

# Count significant correlations for the features represented by each row
add_significant_correlation_counts <- function(data, two_modal) {
  if (nrow(data) == 0) {
    data$num_sig_corr_t1 <- integer()
    data$num_sig_corr_t2 <- integer()
    return(data)
  }

  t1 <- as.character(data$t1)
  t2 <- as.character(data$t2)

  if (two_modal) {
    counts_t1 <- table(t1)
    counts_t2 <- table(t2)
    data$num_sig_corr_t1 <- unname(as.integer(counts_t1[t1]))
    data$num_sig_corr_t2 <- unname(as.integer(counts_t2[t2]))
  } else {
    # In one-modality mode, a feature can occur on either side of an edge.
    counts <- table(c(t1, t2))
    data$num_sig_corr_t1 <- unname(as.integer(counts[t1]))
    data$num_sig_corr_t2 <- unname(as.integer(counts[t2]))
  }

  data
}

# Validate the minimum co-occurrence threshold
validate_cooccur <- function(cooccur) {
  if (!is.numeric(cooccur) || length(cooccur) != 1 || is.na(cooccur) ||
      !is.finite(cooccur) || cooccur < 4 || cooccur != floor(cooccur)) {
    stop("cooccur must be a single integer greater than or equal to 4",
         call. = FALSE)
  }

  cooccur
}

# Validate the numerical boundary used in the Cauchy transformation
validate_cauchy_tolerance <- function(tol) {
  if (!is.numeric(tol) || length(tol) != 1 || is.na(tol) ||
      !is.finite(tol) || tol <= 0 || tol >= 0.5) {
    stop("tol must be a single finite number greater than 0 and less than 0.5",
         call. = FALSE)
  }

  tol
}

# Bound p-values before applying the existing Cauchy transformation
safe_tan_cauchy <- function(p_vect, tol = 1e-16) {
  p_vect <- pmax(pmin(p_vect, 1 - tol), tol)
  1 / length(p_vect) * sum(tan((0.5 - p_vect) * pi))
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

#' Calculate DISruption in COrrelations (DISCO)
#'
#' Quantifies feature-specific changes in Spearman correlation patterns between
#' two groups. Inputs may contain one data modality or two data modalities.
#' Corresponding modalities are aligned to their common features before the
#' correlations are compared.
#'
#' @param data1 A list containing one or two numeric matrices or data frames for
#'   the first group. Rows represent samples and columns represent features.
#' @param data2 A list with the same number of modalities as `data1`, containing
#'   data for the second group. Rows represent samples and columns represent
#'   features.
#' @param cooccur A single integer greater than or equal to 4. A feature pair is
#'   tested only when both groups have at least this many pairwise non-missing
#'   observations.
#' @param filt_threshold Single-group correlation p-value threshold. A pair is
#'   retained when its correlation is significant in at least one group at this
#'   threshold.
#' @param diff_threshold Minimum absolute difference between the two raw
#'   Spearman correlations used to identify significant disrupted pairs.
#' @param p_alpha Two-sample correlation-test p-value threshold.
#' @param q_alpha Adjusted two-sample correlation-test p-value threshold.
#' @param adj_method Multiple-testing adjustment method passed to
#'   [stats::p.adjust()].
#' @param output_dir Directory in which result CSV files are written. The
#'   directory is created when it does not exist.
#' @param prefix Character prefix used for result filenames.
#' @param tol Numerical boundary used to constrain p-values before the Cauchy
#'   transformation. It must be greater than 0 and less than 0.5.
#' @param spearman_correction Standard-error correction for Fisher-transformed
#'   Spearman correlations. Either `"bonett_wright"` or `"fieller_1.06"`.
#'
#' @return A list with three data frames: `cor_diff` contains retained pairwise
#'   comparisons, `cor_diff_sig` contains significant disrupted pairs, and
#'   `disco` contains feature-level DISCO results. The same tables are written
#'   to CSV files in `output_dir`.
#' @export
disco <- function(data1,
                  data2,
                  cooccur = 10,
                  filt_threshold = 0.001,
                  diff_threshold = 0.3,
                  p_alpha = 0.01,
                  q_alpha = 0.2,
                  adj_method = "BY",
                  output_dir,
                  prefix,
                  tol = 1e-16,
                  spearman_correction = c("bonett_wright",
                                          "fieller_1.06")) {

  cooccur <- validate_cooccur(cooccur)
  tol <- validate_cauchy_tolerance(tol)
  spearman_correction <- match.arg(spearman_correction)

  if (length(data1) != length(data2) || !(length(data1) %in% c(1, 2))) {
    stop("Number of data modalities should be one or two for each group")
  }

  for (i in seq_along(data1)) {
    aligned <- align_common_features(data1[[i]], data2[[i]], modality = i)
    data1[[i]] <- aligned$data1
    data2[[i]] <- aligned$data2
  }

  if (length(data1) == 2) {
    two_modal <- TRUE
    n_feat1 <- ncol(data1[[1]])
    n_feat2 <- ncol(data1[[2]])

    report_sample_overlap(data1[[1]], data1[[2]], group = "data1")
    report_sample_overlap(data2[[1]], data2[[2]], group = "data2")

    # Concatenating data modality 1 and data modality 2 for group 1
    d_merged1 <- merge(data1[[1]], data1[[2]], by = "row.names", all = TRUE) %>%
      column_to_rownames("Row.names")
    # Concatenating data modality 1 and data modality 2 for group 2
    d_merged2 <- merge(data2[[1]], data2[[2]], by = "row.names", all = TRUE) %>%
      column_to_rownames("Row.names")

  } else {
    two_modal <- FALSE
    n_feat1 <- ncol(data1[[1]])
    d_merged1 <- data1[[1]]
    d_merged2 <- data2[[1]]
  }

  # Estimate correlations and cooccurrences
  C1 <- get_npairs(d_merged1)
  C2 <- get_npairs(d_merged2)
  R1 <- get_corr(d_merged1)
  R2 <- get_corr(d_merged2)

  # Test only pairs meeting the co-occurrence threshold in both groups.
  comparison <- compare_eligible_correlations(
    R1, C1, R2, C2,
    cooccur = cooccur,
    filt_threshold = filt_threshold,
    correction = spearman_correction
  )
  pvals_tmp <- comparison$pvals

  # Create a long dataframe
  if (two_modal) {
    inds <- as.matrix(expand.grid(
      row = 1:n_feat1,
      col = (n_feat1 + 1):(n_feat1 + n_feat2)
    ))
  } else {
    inds <- which(upper.tri(R1), arr.ind = TRUE)
  }

  t1 <- rownames(R1)[inds[, 1]]
  t2 <- rownames(R1)[inds[, 2]]
  r1 <- R1[inds]
  r2 <- R2[inds]
  n1 <- C1[inds]
  n2 <- C2[inds]

  pvals_filtered <- pvals_tmp[inds]
  adjp <- p.adjust(pvals_filtered, method = adj_method)

  cat("Number of all pairwise correlations:", length(pvals_filtered), "\n")
  cat("Number of retained pairwise correlations:", sum(!is.na(pvals_filtered)), "\n")

  df_tmp <- data.frame(t1, t2, r1, r2, n1, n2, p = pvals_filtered, q = adjp) %>%
    filter(!is.na(p))

  # Find the number of significant correlations
  df_tmp_sig <- df_tmp %>%
    mutate(abs_diff_cor = abs(r1 - r2)) %>%
    filter(p < p_alpha & abs_diff_cor > diff_threshold & q < q_alpha)

  df_tmp_sig <- add_significant_correlation_counts(df_tmp_sig, two_modal)

  cat("Number of significant pairwise correlations:", nrow(df_tmp_sig), "\n")

  # Calculating disco score
  if (two_modal) {
    disco_rows <- sapply(1:n_feat1, function(i) {
      p_vect <- as.vector(na.omit(
        pvals_tmp[(n_feat1 + 1):(n_feat1 + n_feat2), i]
      ))
      if(length(p_vect) == 0) return(NULL)
      t <- safe_tan_cauchy(p_vect, tol = tol)
      pval <- pcauchy(t, location = 0, scale = 1, lower.tail = FALSE)
      data.frame(species = colnames(pvals_tmp)[i], t = t, p = pval, num_corr = length(p_vect))
    }, simplify = FALSE)
  } else {
    disco_rows <- sapply(1:n_feat1, function(i) {
      p_vect <- as.vector(na.omit(pvals_tmp[-i, i]))
      if(length(p_vect) == 0) return(NULL)
      t <- safe_tan_cauchy(p_vect, tol = tol)
      pval <- pcauchy(t, location = 0, scale = 1, lower.tail = FALSE)
      data.frame(species = colnames(pvals_tmp)[i], t = t, p = pval, num_corr = length(p_vect))
    }, simplify = FALSE)
  }

  disco_results <- bind_rows(
    data.frame(
      species = character(),
      t = numeric(),
      p = numeric(),
      num_corr = integer()
    ),
    disco_rows
  )

  # Adjustment of p-values
  disco_results <- disco_results %>%
    mutate(q = p.adjust(p, method = adj_method)) %>%
    arrange(q)

  # Create output directory if missing
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

  # Save the results
  write.csv(df_tmp, file.path(output_dir, paste0(prefix, "_cor_diff.csv")), quote = TRUE, row.names = FALSE)
  write.csv(df_tmp_sig, file.path(output_dir, paste0(prefix, "_cor_diff_sig.csv")), quote = TRUE, row.names = FALSE)
  write.csv(disco_results, file.path(output_dir, paste0(prefix, "_disco.csv")), quote = TRUE, row.names = FALSE)

  return(list(
    cor_diff = df_tmp,
    cor_diff_sig = df_tmp_sig,
    disco = disco_results
  ))
}
