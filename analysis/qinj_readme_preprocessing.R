## ---- readme-select-helper ----
# Select one independent, gender/BMI-complete sample per QinJ participant.
select_qinj_metadata <- function(metadata) {
  metadata <- as.data.frame(metadata)
  required <- c("study_name", "sample_id", "subject_id", "disease", "gender", "BMI")
  if (!all(required %in% names(metadata))) stop("Required QinJ metadata columns are missing.")
  metadata <- metadata[!is.na(metadata$study_name) & metadata$study_name == "QinJ_2012" &
    metadata$disease %in% c("healthy", "T2D"), , drop = FALSE]
  if (!nrow(metadata) || anyNA(metadata$sample_id) || any(!nzchar(metadata$sample_id)) ||
      anyDuplicated(metadata$sample_id)) stop("Sample IDs must be non-missing and unique.")
  metadata$subject_id <- trimws(as.character(metadata$subject_id))
  if (anyNA(metadata$subject_id) || any(!nzchar(metadata$subject_id))) stop("Subject IDs are required.")
  if (!is.numeric(metadata$BMI)) stop("BMI must be numeric.")
  conflicting <- names(which(vapply(split(metadata$disease, metadata$subject_id),
    function(x) length(unique(x)) > 1L, logical(1))))
  conflict <- metadata[metadata$subject_id %in% conflicting, , drop = FALSE]
  conflict$exclusion_reason <- rep("Subject has conflicting disease groups", nrow(conflict))
  metadata <- metadata[!metadata$subject_id %in% conflicting, , drop = FALSE]
  complete <- stats::complete.cases(metadata[, c("gender", "BMI")]) &
    is.finite(metadata$BMI) & nzchar(trimws(as.character(metadata$gender)))
  excluded_covariates <- metadata[!complete, , drop = FALSE]
  excluded_covariates$exclusion_reason <- rep("Missing or invalid gender or BMI", nrow(excluded_covariates))
  metadata <- metadata[complete, , drop = FALSE]
  day <- if ("days_from_first_collection" %in% names(metadata))
    metadata$days_from_first_collection else rep(NA_real_, nrow(metadata))
  if (!is.numeric(day)) stop("Collection day must be numeric.")
  metadata <- metadata[order(metadata$subject_id, day, metadata$sample_id, na.last = TRUE), , drop = FALSE]
  repeated <- duplicated(metadata$subject_id)
  excluded_repeats <- metadata[repeated, , drop = FALSE]
  excluded_repeats$exclusion_reason <- rep("Repeated subject: retained earliest covariate-complete sample",
    nrow(excluded_repeats))
  metadata <- metadata[!repeated, , drop = FALSE]
  rownames(metadata) <- NULL # returnSamples() requires ordinary row indices.
  list(metadata = metadata, excluded_covariates = excluded_covariates,
    excluded_subjects = rbind(conflict, excluded_repeats))
}

## ---- readme-preprocess-helper ----
# Samples are rows; features are columns. Apply independently per modality.
preprocess_abundance <- function(x, group, prevalence_cutoff = 0.10) {
  stopifnot(
    is.matrix(x), is.numeric(x), nrow(x) > 0L, ncol(x) > 0L,
    all(is.na(x) | is.finite(x)), all(x >= 0, na.rm = TRUE),
    all(rowSums(x, na.rm = TRUE) > 0),
    !is.null(rownames(x)), !anyDuplicated(rownames(x)),
    !is.null(colnames(x)), !anyDuplicated(colnames(x)),
    length(prevalence_cutoff) == 1L, is.finite(prevalence_cutoff),
    prevalence_cutoff > 0, prevalence_cutoff <= 1,
    identical(rownames(x), names(group)), !anyNA(group),
    setequal(as.character(group), c("healthy", "T2D"))
  )

  # Require the prevalence cutoff separately in BOTH original groups.
  x <- sweep(x, 1, rowSums(x, na.rm = TRUE), "/")
  prevalence <- cbind(
    healthy = colMeans(x[group == "healthy", , drop = FALSE] > 0, na.rm = TRUE),
    T2D = colMeans(x[group == "T2D", , drop = FALSE] > 0, na.rm = TRUE)
  )
  keep <- rowSums(is.finite(prevalence) & prevalence >= prevalence_cutoff) == 2L
  if (sum(keep) < 2L) stop("Fewer than two features pass the prevalence filter.")
  retained <- x[, keep, drop = FALSE]

  # Robust CLR: centre each sample using only its observed positive features.
  # Zeros and existing missing values stay missing; do not impute them.
  log_x <- retained
  log_x[!is.na(log_x) & log_x == 0] <- NA_real_
  log_x <- log(log_x)
  centre <- rowMeans(log_x, na.rm = TRUE)
  centre[!is.finite(centre)] <- NA_real_
  clr <- sweep(log_x, 1, centre, "-")

  list(
    clr = clr,
    retained_relative_abundance = retained,
    transformation = "robust CLR; zeros remain NA",
    prevalence = data.frame(
      feature = rownames(prevalence), prevalence, retained = unname(keep),
      row.names = NULL
    )
  )
}

## ---- readme-residual-helper ----
# Fit independent gender/BMI models within the original disease groups.
residualize_abundance <- function(clr, metadata, group) {
  stopifnot(identical(rownames(clr), rownames(metadata)),
            all(is.na(clr) | is.finite(clr)),
            is.numeric(metadata$BMI),
            all(is.finite(metadata$BMI)),
            identical(rownames(clr), names(group)), !anyNA(group),
            setequal(as.character(group), c("healthy", "T2D")))
  residuals <- matrix(NA_real_, nrow(clr), ncol(clr), dimnames = dimnames(clr))
  models <- list()
  for (group_name in c("healthy", "T2D")) {
    rows <- which(group == group_name)
    model_data <- metadata[rows, , drop = FALSE]
    model_data$gender <- factor(trimws(as.character(model_data$gender)))
    if (anyNA(model_data$gender) || any(!nzchar(as.character(model_data$gender))) ||
        nlevels(model_data$gender) < 2L) {
      stop("At least two non-missing gender levels are required in ", group_name, ".")
    }
    design <- stats::model.matrix(~ gender + BMI, data = model_data)
    if (qr(design)$rank < ncol(design) || nrow(design) <= ncol(design)) {
      stop("Covariate model requires full rank and positive residual degrees of freedom in ",
           group_name, ".")
    }
    coefficients <- matrix(NA_real_, ncol(design), ncol(clr),
                           dimnames = list(colnames(design), colnames(clr)))
    n_observed <- residual_df <- stats::setNames(integer(ncol(clr)), colnames(clr))
    status <- stats::setNames(rep("ok", ncol(clr)), colnames(clr))
    for (i in seq_len(ncol(clr))) {
      observed <- !is.na(clr[rows, i])
      n_observed[i] <- sum(observed)
      feature_design <- design[observed, , drop = FALSE]
      if (sum(observed) <= ncol(design)) {
        status[i] <- "insufficient_observations"
        residual_df[i] <- NA_integer_
        next
      }
      if (qr(feature_design)$rank < ncol(design)) {
        status[i] <- "rank_deficient"
        residual_df[i] <- NA_integer_
        next
      }
      fit <- stats::lm.fit(feature_design, clr[rows[observed], i])
      coefficients[, i] <- fit$coefficients
      residuals[rows[observed], i] <- fit$residuals
      residual_df[i] <- fit$df.residual
    }
    models[[group_name]] <- list(
      coefficients = coefficients, design = design,
      sample_ids = rownames(model_data), residual_df = residual_df,
      n_observed = n_observed, status = status
    )
    if (any(status != "ok")) {
      warning(sum(status != "ok"), " feature model(s) could not be fitted in ",
              group_name, "; residuals left NA. See models[[group]]$status.", call. = FALSE)
    }
  }
  stopifnot(all(is.na(residuals) | is.finite(residuals)),
            all(is.na(residuals[is.na(clr)])))

  list(residuals = residuals, models = models)
}

## ---- internal-readme-helpers ----
# Track example scripts and the installed package actually used for analysis.
# Relative keys keep provenance independent of checkout/library locations.
readme_source_hashes <- function(repo_dir) {
  paths <- c("analysis/run_qinj_readme_examples.R",
             "analysis/qinj_readme_preprocessing.R")
  example_hashes <- stats::setNames(
    unname(tools::md5sum(file.path(repo_dir, paths))), paths)
  package_dir <- find.package("DISCO")
  package_paths <- c("DESCRIPTION", "NAMESPACE",
    file.path("R", list.files(file.path(package_dir, "R"))))
  package_hashes <- stats::setNames(
    unname(tools::md5sum(file.path(package_dir, package_paths))),
    file.path("installed_DISCO", package_paths))
  c(example_hashes, package_hashes)
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
