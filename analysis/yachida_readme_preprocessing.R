## ---- readme-preprocess-helper ----
# Samples are rows; features are columns. Apply independently per modality.
preprocess_abundance <- function(x, prevalence_cutoff = 0.20) {
  stopifnot(
    is.matrix(x), is.numeric(x), nrow(x) > 0L, ncol(x) > 0L,
    all(is.finite(x)), all(x >= 0), all(rowSums(x) > 0),
    !is.null(rownames(x)), !anyDuplicated(rownames(x)),
    !is.null(colnames(x)), !anyDuplicated(colnames(x)),
    length(prevalence_cutoff) == 1L, is.finite(prevalence_cutoff),
    prevalence_cutoff > 0, prevalence_cutoff <= 1
  )

  # Convert the full profile to proportions, then filter across both groups.
  x <- sweep(x, 1, rowSums(x), "/")
  prevalence <- colMeans(x > 0)
  keep <- prevalence >= prevalence_cutoff
  if (sum(keep) < 2L) stop("Fewer than two features pass the prevalence filter.")
  retained <- x[, keep, drop = FALSE]

  # One half-minimum per feature, estimated across the combined cohort.
  half_min <- apply(retained, 2, function(z) min(z[z > 0]) / 2)
  imputed <- retained
  zero_positions <- which(imputed == 0, arr.ind = TRUE)
  if (nrow(zero_positions) > 0L) {
    imputed[zero_positions] <- half_min[zero_positions[, "col"]]
  }
  log_x <- log(imputed)
  clr <- sweep(log_x, 1, rowMeans(log_x), "-")

  list(
    clr = clr,
    retained_relative_abundance = retained,
    imputed_abundance = imputed,
    prevalence = data.frame(
      feature = names(prevalence), prevalence = unname(prevalence), retained = keep
    ),
    half_minimum = data.frame(feature = names(half_min), replacement = unname(half_min))
  )
}

## ---- readme-residual-helper ----
# Fit one model per feature over the combined cohort, without a disease term.
residualize_abundance <- function(clr, metadata) {
  stopifnot(identical(rownames(clr), rownames(metadata)), all(is.finite(clr)),
            is.numeric(metadata$age), is.numeric(metadata$BMI),
            all(is.finite(metadata$age)), all(is.finite(metadata$BMI)))
  model_data <- metadata
  model_data$gender <- factor(trimws(as.character(model_data$gender)))
  if (anyNA(model_data$gender) || any(!nzchar(as.character(model_data$gender))) ||
      nlevels(model_data$gender) < 2L) {
    stop("At least two non-missing gender levels are required.")
  }
  design <- stats::model.matrix(~ age + gender + BMI, data = model_data)
  fit <- stats::lm.fit(design, clr)
  if (fit$rank < ncol(design) || nrow(design) <= fit$rank) {
    stop("Covariate model must have full rank and positive residual degrees of freedom.")
  }
  residuals <- fit$residuals
  dimnames(residuals) <- dimnames(clr)
  stopifnot(all(is.finite(residuals)))
  list(residuals = residuals, coefficients = fit$coefficients,
       design = design, residual_df = nrow(design) - fit$rank)
}

## ---- internal-readme-helpers ----
# Hash relative paths so moving the checkout does not invalidate provenance.
readme_source_hashes <- function(repo_dir) {
  paths <- c("analysis/run_yachida_readme_examples.R",
             "analysis/yachida_readme_preprocessing.R",
             list.files(file.path(repo_dir, "R"), pattern = "[.]R$",
                        full.names = FALSE))
  paths[-c(1L, 2L)] <- file.path("R", paths[-c(1L, 2L)])
  stats::setNames(unname(tools::md5sum(file.path(repo_dir, paths))), paths)
}

readme_cached_result <- function(output_dir, input, recompute) {
  result_file <- file.path(output_dir, "disco_analysis.rds")
  input_file <- file.path(output_dir, "disco_inputs.rds")
  if (recompute || !file.exists(result_file) || !file.exists(input_file) ||
      !identical(readRDS(input_file), input)) return(NULL)
  result <- readRDS(result_file)
  if (!all(file.exists(result$output_files))) return(NULL)
  message("Reusing matching results: ", output_dir)
  result
}
