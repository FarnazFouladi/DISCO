library(dplyr)
library(tibble)
library(Hmisc)

# Fisher's Z test for a single sample
fisher_z_p <- function(R, N) {
  # Bound R to avoid Inf/-Inf when R is exactly 1 or -1
  #R <- pmax(pmin(R, 0.9999), -0.9999)

  Z <- 0.5 * log((1 + R) / (1 - R))

  # Correct for the bias related to sample size
  Z_bc <- Z - R / (2 * (N - 1))
  se <- 1 / sqrt(pmax(N - 3, 1))

  2 * pnorm(-abs(Z_bc / se))
}

# Fisher's Z test for two samples
fisher_z_p_2sample <- function(R1, N1, R2, N2) {
  # Bound R to prevent Inf
  #R1 <- pmax(pmin(R1, 0.9999), -0.9999)
  #R2 <- pmax(pmin(R2, 0.9999), -0.9999)

  Z1 <- 0.5 * log((1 + R1) / (1 - R1))
  Z2 <- 0.5 * log((1 + R2) / (1 - R2))

  # Correct for different sample sizes between the two groups
  Z1_bc <- Z1 - R1 / (2 * (N1 - 1))
  Z2_bc <- Z2 - R2 / (2 * (N2 - 1))

  se <- sqrt(1 / pmax(N1 - 3, 1) + 1 / pmax(N2 - 3, 1))
  Z <- (Z1_bc - Z2_bc) / se

  pvals <- 2 * pnorm(-abs(Z))

  # Returning Z2_bc - Z1_bc (delta Z)
  return(list(pvals, Z2_bc - Z1_bc))
}

# Wrapper for Spearman correlation
get_corr <- function(D) {
  Hmisc::rcorr(as.matrix(D), type = "spearman")$r
}

# Highly optimized pairwise complete cases counter
get_npairs <- function(D) {
  crossprod(!is.na(D))
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

# Main DISCO Score Function
disco_score <- function(data1, # [List] Data for group 1
                        data2, # [List] Data for group 2
                        cooccur = 10,
                        filt_threshold = 0.001,
                        diff_threshold = 0.3,
                        p_alpha = 0.01,
                        q_alpha = 0.2,
                        adj_method = "BY",
                        output_dir,
                        prefix) {

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
    p <- ncol(data1[[1]])
    q <- ncol(data1[[2]])

    # Concatenating data modality 1 and data modality 2 for group 1
    d_merged1 <- merge(data1[[1]], data1[[2]], by = "row.names", all = TRUE) %>%
      column_to_rownames("Row.names")
    # Concatenating data modality 1 and data modality 2 for group 2
    d_merged2 <- merge(data2[[1]], data2[[2]], by = "row.names", all = TRUE) %>%
      column_to_rownames("Row.names")

  } else {
    two_modal <- FALSE
    p <- ncol(data1[[1]])
    q <- 0
    d_merged1 <- data1[[1]]
    d_merged2 <- data2[[1]]
  }

  # Enforce rule: Zeros are NAs (Working with complete data)
  #d_merged1[d_merged1 == 0] <- NA
  #d_merged2[d_merged2 == 0] <- NA

  # Estimate correlations and cooccurrences
  C1 <- get_npairs(d_merged1)
  C2 <- get_npairs(d_merged2)
  R1 <- get_corr(d_merged1)
  R2 <- get_corr(d_merged2)

  # Apply cooccurrence filter
  R1[C1 < cooccur] <- 0
  R2[C2 < cooccur] <- 0

  # Fisher's Z test for single correlations
  P1 <- fisher_z_p(R1, C1)
  P2 <- fisher_z_p(R2, C2)

  # Filter correlations with p >= threshold in both groups
  filt <- (P1 >= filt_threshold) & (P2 >= filt_threshold)
  diag(filt) <- TRUE

  # Fisher Z test - two sample test
  fisher_res <- fisher_z_p_2sample(R1, C1, R2, C2)
  pvals <- fisher_res[[1]]
  delta_R <- fisher_res[[2]]
  delta_R[is.na(delta_R)] <- 0

  # Replace p-values with NAs for correlations that are filtered out
  pvals_tmp <- pvals
  pvals_tmp[filt] <- NA

  # Create a long dataframe
  if (two_modal) {
    inds <- as.matrix(expand.grid(row = 1:p, col = (p + 1):(p + q)))
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

  # Create a matrix for adjusted p-values
  adj_pvals <- matrix(NA, nrow = nrow(pvals_tmp), ncol = ncol(pvals_tmp))
  adj_pvals[inds] <- adjp

  df_tmp <- data.frame(t1, t2, r1, r2, n1, n2, p = pvals_filtered, q = adjp) %>%
    filter(!is.na(p))

  # Find the number of significant correlations
  df_tmp_sig <- df_tmp %>%
    mutate(abs_diff_cor = abs(r1 - r2)) %>%
    filter(p < p_alpha & abs_diff_cor > diff_threshold & q < q_alpha) %>%
    group_by(t1) %>%
    mutate(num_sig_corr_t1 = n()) %>%
    ungroup() %>%
    group_by(t2) %>%
    mutate(num_sig_corr_t2 = n()) %>%
    ungroup() # Added final ungroup

  cat("Number of significant pairwise correlations:", nrow(df_tmp_sig), "\n")

  # Calculating disco score
  # Helper to prevent tan(pi/2) resulting in Inf
  safe_tan_cauchy <- function(p_vect) {
    p_vect <- pmax(pmin(p_vect, 1 - 1e-7), 1e-7)
    1 / length(p_vect) * sum(tan((0.5 - p_vect) * pi))
  }

  if (two_modal) {
    disco <- bind_rows(sapply(1:p, function(i) {
      p_vect <- as.vector(na.omit(pvals_tmp[(p + 1):(p + q), i]))
      if(length(p_vect) == 0) return(NULL)
      t <- safe_tan_cauchy(p_vect)
      pval <- pcauchy(t, location = 0, scale = 1, lower.tail = FALSE)
      data.frame(species = colnames(pvals_tmp)[i], t = t, p = pval, num_corr = length(p_vect))
    }, simplify = FALSE))
  } else {
    disco <- bind_rows(sapply(1:p, function(i) {
      p_vect <- as.vector(na.omit(pvals_tmp[-i, i]))
      if(length(p_vect) == 0) return(NULL)
      t <- safe_tan_cauchy(p_vect)
      pval <- pcauchy(t, location = 0, scale = 1, lower.tail = FALSE)
      data.frame(species = colnames(pvals_tmp)[i], t = t, p = pval, num_corr = length(p_vect))
    }, simplify = FALSE))
  }

  # Adjustment of p-values
  disco <- disco %>%
    mutate(q = p.adjust(p, method = adj_method)) %>%
    arrange(q)

  # Create output directory if missing
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

  # Save the results
  write.csv(df_tmp, file.path(output_dir, paste0(prefix, "_cor_diff.csv")), quote = TRUE, row.names = FALSE)
  write.csv(df_tmp_sig, file.path(output_dir, paste0(prefix, "_cor_diff_sig.csv")), quote = TRUE, row.names = FALSE)
  write.csv(disco, file.path(output_dir, paste0(prefix, "_disco.csv")), quote = TRUE, row.names = FALSE)
  write.csv(delta_R, file.path(output_dir, paste0(prefix, "_delta_cor.csv")), quote = TRUE, row.names = TRUE)
  write.csv(pvals_tmp, file.path(output_dir, paste0(prefix, "_p.csv")), quote = TRUE, row.names = TRUE)
  write.csv(adj_pvals, file.path(output_dir, paste0(prefix, "_q.csv")), quote = TRUE, row.names = TRUE)

  return(list(df_tmp = df_tmp, df_tmp_sig = df_tmp_sig, disco = disco, delta_R = delta_R, pvals_tmp = pvals_tmp, adj_pvals = adj_pvals))
}
