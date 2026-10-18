# Lightweight checks without downloading data or running permutations.
# Run: Rscript analysis/test_yachida_readme_helpers.R
file_argument <- grep("^--file=", commandArgs(FALSE), value = TRUE)
if (length(file_argument) != 1L) stop("Run this file with Rscript.")
analysis_dir <- dirname(normalizePath(sub("^--file=", "", file_argument[[1L]])))
source(file.path(analysis_dir, "yachida_readme_preprocessing.R"))

set.seed(19)
metadata <- data.frame(
  age = stats::rnorm(24, 50, 10), gender = rep(c("female", "male"), 12),
  BMI = stats::rnorm(24, 25, 4), row.names = paste0("sample", seq_len(24))
)
abundance <- matrix(stats::runif(120), nrow = 24,
                    dimnames = list(rownames(metadata), paste0("feature", seq_len(5))))
abundance[seq_len(8), 1] <- 0
group <- stats::setNames(rep(c("healthy", "CRC"), each = 12), rownames(metadata))
prepared <- preprocess_abundance(abundance, group)
zeros <- which(prepared$retained_relative_abundance == 0, arr.ind = TRUE)
stopifnot(
  identical(unname(prepared$imputed_abundance[zeros]),
            unname(prepared$half_minimum$replacement[zeros[, "col"]])),
  max(abs(rowMeans(prepared$clr))) < 1e-12
)
adjusted <- residualize_abundance(prepared$clr, metadata, group)
for (group_name in c("healthy", "CRC")) {
  rows <- which(group == group_name)
  for (i in seq_len(ncol(prepared$clr))) {
    reference <- stats::lm(prepared$clr[rows, i] ~ age + gender + BMI,
                           data = metadata[rows, , drop = FALSE])
    stopifnot(isTRUE(all.equal(unname(adjusted$residuals[rows, i]),
                              unname(stats::residuals(reference)))))
  }
  design <- adjusted$models[[group_name]]$design
  stopifnot(max(abs(crossprod(design, adjusted$residuals[rows, ]))) < 1e-10,
            !any(grepl("disease", colnames(design))))
}
stopifnot(
  identical(dimnames(adjusted$residuals), dimnames(prepared$clr))
)
misaligned <- try(residualize_abundance(prepared$clr, metadata[24:1, ], group), silent = TRUE)
single_gender <- metadata
single_gender$gender <- "female"
invalid_gender <- try(residualize_abundance(prepared$clr, single_gender, group), silent = TRUE)
rank_deficient <- metadata
rank_deficient$BMI <- rank_deficient$age
invalid_rank <- try(residualize_abundance(prepared$clr, rank_deficient, group), silent = TRUE)
stopifnot(inherits(misaligned, "try-error"), inherits(invalid_gender, "try-error"),
          inherits(invalid_rank, "try-error"))

# A taxon common in only one group must fail, even with high pooled prevalence.
boundary_group <- stats::setNames(rep(c("healthy", "CRC"), each = 10), paste0("b", 1:20))
boundary <- cbind(
  common = rep(1, 20), exclusive = c(rep(1, 10), rep(0, 10)),
  exactly_10pct = rep(c(1, rep(0, 9)), 2)
)
rownames(boundary) <- names(boundary_group)
filtered <- preprocess_abundance(boundary, boundary_group)
stopifnot(identical(colnames(filtered$clr), c("common", "exactly_10pct")),
          all(filtered$prevalence$healthy[filtered$prevalence$retained] >= 0.1),
          all(filtered$prevalence$CRC[filtered$prevalence$retained] >= 0.1))
message("Yachida per-group prevalence and separate-group residualization checks passed.")
