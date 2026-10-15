# Prepare aligned inputs for repeated DISCO calculations
prepare_disco_inputs <- function(data1, data2, report_overlap = TRUE) {
  if (length(data1) != length(data2) || !(length(data1) %in% c(1, 2))) {
    stop("Number of data modalities should be one or two for each group",
         call. = FALSE)
  }

  for (i in seq_along(data1)) {
    aligned <- align_common_features(data1[[i]], data2[[i]], modality = i)
    data1[[i]] <- aligned$data1
    data2[[i]] <- aligned$data2
  }

  two_modal <- length(data1) == 2L
  n_feat1 <- ncol(data1[[1]])

  if (two_modal) {
    n_feat2 <- ncol(data1[[2]])
    if (report_overlap) {
      report_sample_overlap(data1[[1]], data1[[2]], group = "data1")
      report_sample_overlap(data2[[1]], data2[[2]], group = "data2")
    }

    d_merged1 <- merge(data1[[1]], data1[[2]],
                       by = "row.names", all = TRUE) %>%
      column_to_rownames("Row.names")
    d_merged2 <- merge(data2[[1]], data2[[2]],
                       by = "row.names", all = TRUE) %>%
      column_to_rownames("Row.names")
  } else {
    n_feat2 <- 0L
    d_merged1 <- data1[[1]]
    d_merged2 <- data2[[1]]
  }

  list(
    data1 = data1,
    data2 = data2,
    two_modal = two_modal,
    n_feat1 = n_feat1,
    n_feat2 = n_feat2,
    d_merged1 = d_merged1,
    d_merged2 = d_merged2
  )
}

# Calculate quantities that do not depend on the filtering threshold
calculate_disco_base <- function(prepared, cooccur) {
  C1 <- get_npairs(prepared$d_merged1)
  C2 <- get_npairs(prepared$d_merged2)
  R1 <- get_corr(prepared$d_merged1)
  R2 <- get_corr(prepared$d_merged2)

  correlation_tests <- calculate_correlation_tests(
    R1, C1, R2, C2, cooccur = cooccur
  )

  if (prepared$two_modal) {
    inds <- as.matrix(expand.grid(
      row = seq_len(prepared$n_feat1),
      col = prepared$n_feat1 + seq_len(prepared$n_feat2)
    ))
  } else {
    inds <- which(upper.tri(R1), arr.ind = TRUE)
  }

  list(
    R1 = R1,
    R2 = R2,
    C1 = C1,
    C2 = C2,
    correlation_tests = correlation_tests,
    inds = inds,
    two_modal = prepared$two_modal,
    n_feat1 = prepared$n_feat1,
    n_feat2 = prepared$n_feat2
  )
}

# Build pair- and feature-level results at one filtering threshold
calculate_disco_at_threshold <- function(base, filt_threshold, adj_method) {
  comparison <- filter_correlation_tests(
    base$correlation_tests,
    filt_threshold = filt_threshold
  )
  pvals_tmp <- comparison$pvals
  inds <- base$inds

  t1 <- rownames(base$R1)[inds[, 1]]
  t2 <- rownames(base$R1)[inds[, 2]]
  r1 <- base$R1[inds]
  r2 <- base$R2[inds]
  n1 <- base$C1[inds]
  n2 <- base$C2[inds]
  pvals_filtered <- pvals_tmp[inds]
  adjp <- stats::p.adjust(pvals_filtered, method = adj_method)

  cor_diff <- data.frame(
    t1 = t1,
    t2 = t2,
    r1 = r1,
    r2 = r2,
    n1 = n1,
    n2 = n2,
    p = pvals_filtered,
    q = adjp,
    stringsAsFactors = FALSE
  )
  cor_diff <- cor_diff[!is.na(cor_diff$p), , drop = FALSE]
  cor_diff <- cor_diff[order(cor_diff$q, na.last = TRUE), , drop = FALSE]
  rownames(cor_diff) <- NULL

  disco_rows <- lapply(seq_len(base$n_feat1), function(i) {
    if (base$two_modal) {
      modality2_rows <- base$n_feat1 + seq_len(base$n_feat2)
      p_vect <- as.vector(stats::na.omit(pvals_tmp[modality2_rows, i]))
    } else {
      p_vect <- as.vector(stats::na.omit(pvals_tmp[-i, i]))
    }

    if (length(p_vect) == 0L) {
      return(NULL)
    }

    statistic <- safe_tan_cauchy(p_vect)
    data.frame(
      species = colnames(pvals_tmp)[i],
      t = statistic,
      p = stats::pcauchy(
        statistic,
        location = 0,
        scale = 1,
        lower.tail = FALSE
      ),
      num_corr = length(p_vect),
      stringsAsFactors = FALSE
    )
  })

  disco_results <- dplyr::bind_rows(
    data.frame(
      species = character(),
      t = numeric(),
      p = numeric(),
      num_corr = integer(),
      stringsAsFactors = FALSE
    ),
    disco_rows
  )
  disco_results$q <- stats::p.adjust(disco_results$p, method = adj_method)
  disco_results <- disco_results[order(disco_results$q), , drop = FALSE]

  list(cor_diff = cor_diff, disco = disco_results)
}

# Run one fixed threshold for legacy sensitivity-analysis scripts
calculate_fixed_threshold_disco <- function(data1, data2, cooccur = 10,
                                            filt_threshold,
                                            adj_method = "BY",
                                            output_dir, prefix) {
  cooccur <- validate_cooccur(cooccur)
  if (!is.numeric(filt_threshold) || length(filt_threshold) != 1L ||
      is.na(filt_threshold) || !is.finite(filt_threshold) ||
      filt_threshold <= 0 || filt_threshold > 1) {
    stop(
      "filter_threshold must be a single number greater than 0 and less than or equal to 1",
      call. = FALSE
    )
  }

  prepared <- prepare_disco_inputs(data1, data2, report_overlap = TRUE)
  base <- calculate_disco_base(prepared, cooccur = cooccur)
  result <- calculate_disco_at_threshold(
    base,
    filt_threshold = filt_threshold,
    adj_method = adj_method
  )

  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  utils::write.csv(
    result$cor_diff,
    file.path(output_dir, paste0(prefix, "_cor_diff.csv")),
    quote = TRUE,
    row.names = FALSE
  )
  utils::write.csv(
    result$disco,
    file.path(output_dir, paste0(prefix, "_disco.csv")),
    quote = TRUE,
    row.names = FALSE
  )

  result
}
