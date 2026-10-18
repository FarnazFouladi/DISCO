## ---- readme-preprocess-helper ----
# Samples are rows; features are columns. Apply independently per modality.
preprocess_abundance <- function(x, group, prevalence_cutoff = 0.10) {
  stopifnot(
    is.matrix(x), is.numeric(x), nrow(x) > 0L, ncol(x) > 0L,
    all(is.finite(x)), all(x >= 0), all(rowSums(x) > 0),
    !is.null(rownames(x)), !anyDuplicated(rownames(x)),
    !is.null(colnames(x)), !anyDuplicated(colnames(x)),
    length(prevalence_cutoff) == 1L, is.finite(prevalence_cutoff),
    prevalence_cutoff > 0, prevalence_cutoff <= 1,
    identical(rownames(x), names(group)), !anyNA(group),
    setequal(as.character(group), c("healthy", "CRC"))
  )

  # Require the prevalence cutoff separately in BOTH original groups.
  x <- sweep(x, 1, rowSums(x), "/")
  prevalence <- cbind(
    healthy = colMeans(x[group == "healthy", , drop = FALSE] > 0),
    CRC = colMeans(x[group == "CRC", , drop = FALSE] > 0)
  )
  keep <- apply(prevalence >= prevalence_cutoff, 1, all)
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
      feature = rownames(prevalence), prevalence, retained = unname(keep),
      row.names = NULL
    ),
    half_minimum = data.frame(feature = names(half_min), replacement = unname(half_min))
  )
}

## ---- readme-residual-helper ----
# Fit independent age/gender/BMI models within the original disease groups.
residualize_abundance <- function(clr, metadata, group) {
  stopifnot(identical(rownames(clr), rownames(metadata)), all(is.finite(clr)),
            is.numeric(metadata$age), is.numeric(metadata$BMI),
            all(is.finite(metadata$age)), all(is.finite(metadata$BMI)),
            identical(rownames(clr), names(group)), !anyNA(group),
            setequal(as.character(group), c("healthy", "CRC")))
  residuals <- matrix(NA_real_, nrow(clr), ncol(clr), dimnames = dimnames(clr))
  models <- list()
  for (group_name in c("healthy", "CRC")) {
    rows <- which(group == group_name)
    model_data <- metadata[rows, , drop = FALSE]
    model_data$gender <- factor(trimws(as.character(model_data$gender)))
    if (anyNA(model_data$gender) || any(!nzchar(as.character(model_data$gender))) ||
        nlevels(model_data$gender) < 2L) {
      stop("At least two non-missing gender levels are required in ", group_name, ".")
    }
    design <- stats::model.matrix(~ age + gender + BMI, data = model_data)
    fit <- stats::lm.fit(design, clr[rows, , drop = FALSE])
    if (fit$rank < ncol(design) || nrow(design) <= fit$rank) {
      stop("Covariate model requires full rank and positive residual degrees of freedom in ",
           group_name, ".")
    }
    residuals[rows, ] <- fit$residuals
    models[[group_name]] <- list(
      coefficients = fit$coefficients, design = design,
      sample_ids = rownames(model_data), residual_df = nrow(design) - fit$rank
    )
  }
  stopifnot(all(is.finite(residuals)))
  list(residuals = residuals, models = models)
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
