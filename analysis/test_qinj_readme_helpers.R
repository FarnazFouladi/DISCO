# Lightweight checks without downloading data or running permutations.
# Run: Rscript analysis/test_qinj_readme_helpers.R
file_argument <- grep("^--file=", commandArgs(FALSE), value = TRUE)
if (length(file_argument) != 1L) stop("Run this file with Rscript.")
analysis_dir <- dirname(normalizePath(sub("^--file=", "", file_argument[[1L]])))
source(file.path(analysis_dir, "qinj_readme_preprocessing.R"))

set.seed(19)
metadata <- data.frame(
  age = stats::rnorm(24, 50, 10), gender = rep(c("female", "male"), 12),
  BMI = stats::rnorm(24, 25, 4),
  age_category = rep(c("adult", "adult", "senior"), 8), row.names = paste0("sample", seq_len(24))
)
abundance <- matrix(stats::runif(120), nrow = 24,
                    dimnames = list(rownames(metadata), paste0("feature", seq_len(5))))
abundance[seq_len(4), 1] <- 0
group <- stats::setNames(rep(c("healthy", "T2D"), each = 12), rownames(metadata))
prepared <- preprocess_abundance(abundance, group)
zeros <- which(prepared$retained_relative_abundance == 0, arr.ind = TRUE)
stopifnot(
  all(is.na(prepared$clr[zeros])),
  identical(is.na(prepared$clr), prepared$retained_relative_abundance == 0),
  max(abs(rowMeans(prepared$clr, na.rm = TRUE))) < 1e-12
)
# Check the actual positive-only log-ratio calculation, not just centring.
for (i in seq_len(nrow(abundance))) {
  positive <- prepared$retained_relative_abundance[i, ] > 0
  log_values <- log(prepared$retained_relative_abundance[i, positive])
  stopifnot(isTRUE(all.equal(unname(prepared$clr[i, positive]),
                            unname(log_values - mean(log_values)))))
}
adjusted <- residualize_abundance(prepared$clr, metadata, group)
for (group_name in c("healthy", "T2D")) {
  rows <- which(group == group_name)
  design <- adjusted$models[[group_name]]$design
  for (i in seq_len(ncol(prepared$clr))) {
    observed <- !is.na(prepared$clr[rows, i])
    reference <- stats::lm(prepared$clr[rows[observed], i] ~ gender + BMI,
                           data = metadata[rows[observed], , drop = FALSE])
    stopifnot(isTRUE(all.equal(unname(adjusted$residuals[rows[observed], i]),
                              unname(stats::residuals(reference)))),
              max(abs(crossprod(design[observed, , drop = FALSE],
                                adjusted$residuals[rows[observed], i]))) < 1e-10)
  }
  stopifnot(!any(grepl("disease", colnames(design))))
}
stopifnot(identical(dimnames(adjusted$residuals), dimnames(prepared$clr)),
          identical(is.na(adjusted$residuals), is.na(prepared$clr)))

# With no missing values, retain the previous multivariate lm.fit behaviour.
complete <- matrix(stats::rnorm(120), nrow = 24, dimnames = dimnames(abundance))
complete_adjusted <- residualize_abundance(complete, metadata, group)
for (group_name in c("healthy", "T2D")) {
  rows <- which(group == group_name)
  design <- stats::model.matrix(~ gender + BMI, metadata[rows, ])
  reference <- stats::lm.fit(design, complete[rows, ])$residuals
  stopifnot(isTRUE(all.equal(unname(complete_adjusted$residuals[rows, ]), unname(reference))))
}

# Untestable feature-specific models stay missing, with explicit diagnostics.
sparse <- prepared$clr
sparse[4:12, 2] <- NA_real_
sparse[which(group == "healthy" & metadata$gender == "male"), 3] <- NA_real_
sparse_adjusted <- suppressWarnings(residualize_abundance(sparse, metadata, group))
stopifnot(sparse_adjusted$models$healthy$status[2] == "insufficient_observations",
          sparse_adjusted$models$healthy$status[3] == "rank_deficient",
          all(is.na(sparse_adjusted$residuals[1:12, 2:3])),
          all(sparse_adjusted$models$T2D$status == "ok"))

misaligned <- try(residualize_abundance(prepared$clr, metadata[24:1, ], group), silent = TRUE)
single_gender <- metadata
single_gender$gender <- "female"
invalid_gender <- try(residualize_abundance(prepared$clr, single_gender, group), silent = TRUE)
rank_deficient <- metadata
rank_deficient$BMI <- 25
invalid_rank <- try(residualize_abundance(prepared$clr, rank_deficient, group), silent = TRUE)
stopifnot(inherits(misaligned, "try-error"), inherits(invalid_gender, "try-error"),
          inherits(invalid_rank, "try-error"))

# A taxon common in only one group must fail, even with high pooled prevalence.
boundary_group <- stats::setNames(rep(c("healthy", "T2D"), each = 10), paste0("b", 1:20))
boundary <- cbind(
  common = rep(1, 20), exclusive = c(rep(1, 10), rep(0, 10)),
  exactly_10pct = rep(c(1, rep(0, 9)), 2)
)
rownames(boundary) <- names(boundary_group)
filtered <- preprocess_abundance(boundary, boundary_group)
stopifnot(identical(colnames(filtered$clr), c("common", "exactly_10pct")),
          all(filtered$prevalence$healthy[filtered$prevalence$retained] >= 0.1),
          all(filtered$prevalence$T2D[filtered$prevalence$retained] >= 0.1))
# A sample with no positive retained features remains an all-NA row.
tiny <- rbind(c(0, 0, 1), c(1, 2, 0), c(1, 3, 0), c(2, 1, 0))
dimnames(tiny) <- list(paste0("t", 1:4), c("a", "b", "excluded"))
tiny_group <- stats::setNames(c("healthy", "healthy", "T2D", "T2D"), rownames(tiny))
tiny_result <- preprocess_abundance(tiny, tiny_group, 0.5)
stopifnot(identical(colnames(tiny_result$clr), c("a", "b")),
          all(is.na(tiny_result$clr[1, ])), !any(is.nan(tiny_result$clr)))

# Existing missing measurements remain missing as well as zeros.
with_na <- abundance
with_na[2, 2] <- NA_real_
na_result <- preprocess_abundance(with_na, group)
stopifnot(is.na(na_result$clr[2, 2]), all(is.na(na_result$clr[zeros])))
message("Robust CLR, missing-value preservation, prevalence and feature-wise covariate-model checks passed.")

# Repeated visits, covariate missingness, and disease conflicts are audited.
selection_metadata <- data.frame(study_name = "QinJ_2012", sample_id = paste0("s",1:6),
  subject_id = c("a","a","b","c","d","d"),
  disease = c("healthy","healthy","T2D","T2D","healthy","T2D"),
  gender = c("female","female","male","female","male","male"),
  BMI = c(NA,25,26,27,28,29), age_category = rep("adult", 6), days_from_first_collection = c(0,10,0,0,0,10))
selected <- select_qinj_metadata(selection_metadata)
stopifnot(identical(selected$metadata$sample_id, c("s2","s3","s4")),
  identical(selected$excluded_covariates$sample_id, "s1"),
  setequal(selected$excluded_subjects$sample_id,c("s5","s6")))
# Age data are not required, and missing ages must not exclude participants.
without_age <- selection_metadata
without_age$age_category <- NULL
without_age_selection <- select_qinj_metadata(without_age)
stopifnot(identical(selected$metadata$sample_id, without_age_selection$metadata$sample_id))
message("QinJ participant selection checks passed.")

# Check the installed metadata, source() setup, and numerical settings without
# launching the full data download or 1000-permutation analyses.
real <- select_qinj_metadata(curatedMetagenomicData::sampleMetadata)
stopifnot(nrow(real$metadata) == 344L,
  sum(real$metadata$disease == "healthy") == 174L,
  sum(real$metadata$disease == "T2D") == 170L,
  nrow(real$excluded_covariates) == 19L)
runner <- file.path(analysis_dir, "run_qinj_readme_examples.R")
setup <- new.env(parent = globalenv())
setup$commandArgs <- function(...) c("R", paste0("--file=", runner))
for (expr in parse(runner)) {
  if (is.call(expr) && identical(expr[[1]], as.name("<-")) &&
      identical(expr[[2]], as.name("results_dir"))) break
  eval(expr, setup)
}
# Snapshot of the Nielsen settings used to create this workflow; the QinJ
# checks remain runnable after the separate Nielsen analysis is archived.
reference_settings <- list(workers = 12L, n_resamples = 1000L, random_seed = 1L,
  cooccur = 20L, adj_method = "BY", q_value = .2, global_p_cutoff = .05,
  transformation = "robust CLR; zeros remain NA")
stopifnot(identical(setup$settings[names(reference_settings)], reference_settings),
  is.numeric(setup$settings$prevalence_cutoff), length(setup$settings$prevalence_cutoff)==1L,
  is.finite(setup$settings$prevalence_cutoff), setup$settings$prevalence_cutoff > 0,
  setup$settings$prevalence_cutoff <= 1)
bootstrap <- new.env(parent = globalenv())
bootstrap$source <- function(file, ...) {
  stopifnot(identical(normalizePath(file),file.path(analysis_dir,"qinj_readme_preprocessing.R")))
  stop(structure(list(message="Path verified",call=NULL),
    class=c("qinj_bootstrap_verified","error","condition")))
}
stopifnot(tryCatch({base::source(runner,local=bootstrap); FALSE},
  qinj_bootstrap_verified=function(e) TRUE))

# Small integration run, using correlated sparse features to make empirical
# thresholds evaluable in two permutations. No real-data downloads.
set.seed(91)
n <- 96L
m <- data.frame(gender = rep(c("female","male"), n/2),
  BMI = stats::rnorm(n,25,3), age_category = rep(c("adult","adult","senior","senior"),n/4),
  row.names = paste0("s",seq_len(n)))
g <- stats::setNames(rep(c("healthy","T2D"),each=n/2),rownames(m))
x <- exp(outer(stats::rnorm(n),c(-3,-2,-1,1,2,3)) + matrix(stats::rnorm(n*6,sd=.02),n))
dimnames(x) <- list(rownames(m),paste0("f",1:6))
x[cbind(seq_len(n),rep(1:6,length.out=n))] <- 0
p <- preprocess_abundance(x,g,.2)
a <- residualize_abundance(p$clr,m,g)
path <- x[,1:4]; colnames(path) <- paste0("path",1:4)
pa <- residualize_abundance(preprocess_abundance(path,g,.2)$clr,m,g)
for (mode in c("taxa","taxa_pathways")) {
  values <- if (mode=="taxa") list(a$residuals) else list(a$residuals,pa$residuals)
  r <- suppressMessages(DISCO::disco(
    data1=lapply(values,function(z) z[g=="healthy",,drop=FALSE]),
    data2=lapply(values,function(z) z[g=="T2D",,drop=FALSE]),
    empirically_derive_threshold=TRUE, n_resamples=2L,workers=1L,
    random_seed=1L,cooccur=20L,q_value=.2,global_p_cutoff=.05,adj_method="BY",
    output_dir=tempfile(paste0("qinj-",mode)),prefix="test"))
  stopifnot(r$global_test$n_valid_resamples==2L,
    all(r$cor_diff$n1 >= 20),all(r$cor_diff$n2 >= 20),all(file.exists(r$output_files)))
}
message("QinJ real metadata, settings validation, source() startup, and both DISCO modalities passed.")
