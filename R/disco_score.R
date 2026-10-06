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
