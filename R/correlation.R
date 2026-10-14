# Bonett-Wright standard error for a Fisher-z-transformed Spearman correlation
spearman_fisher_se <- function(R, N) {
  sqrt((1 + R^2 / 2) / pmax(N - 3, 1))
}

# Fisher's Z test for a single Spearman correlation
fisher_z_p <- function(R, N) {
  Z <- 0.5 * log((1 + R) / (1 - R))

  # Correct for the bias related to sample size
  Z_bc <- Z - R / (2 * (N - 1))
  se <- spearman_fisher_se(R, N)

  2 * pnorm(-abs(Z_bc / se))
}

# Fisher's Z test for two independent Spearman correlations
fisher_z_p_2sample <- function(R1, N1, R2, N2) {
  Z1 <- 0.5 * log((1 + R1) / (1 - R1))
  Z2 <- 0.5 * log((1 + R2) / (1 - R2))

  # Correct for different sample sizes between the two groups
  Z1_bc <- Z1 - R1 / (2 * (N1 - 1))
  Z2_bc <- Z2 - R2 / (2 * (N2 - 1))

  se1 <- spearman_fisher_se(R1, N1)
  se2 <- spearman_fisher_se(R2, N2)
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
                                          cooccur, filt_threshold) {
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
    P1[eligible] <- fisher_z_p(R1[eligible], C1[eligible])
    P2[eligible] <- fisher_z_p(R2[eligible], C2[eligible])

    pvals[eligible] <- fisher_z_p_2sample(
      R1[eligible], C1[eligible], R2[eligible], C2[eligible]
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
