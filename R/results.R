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
