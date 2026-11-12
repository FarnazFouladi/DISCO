# Helpers for exploratory independent-group DISCO on SRP042228 ileal biopsies.
# Choose representatives using annotation only, independently of disease/results.
# Disjoint pairs are required: overlapping clusters need an explicit policy.
ibd_load_pathway_membership <- function(reference, collection = "H", subcollection = NULL) {
  scalar_string <- function(x) is.character(x) && length(x) == 1L && !is.na(x) && nzchar(x)
  if (!scalar_string(collection)) stop("ibd.msigdb_collection must be a nonempty string.")
  if (!is.null(subcollection) && !scalar_string(subcollection))
    stop("ibd.msigdb_subcollection must be NULL or a nonempty string.")
  # Reuse the established human Hallmark snapshot from the CRC workflows.
  if (collection == "H" && is.null(subcollection)) {
    file <- file.path(reference, "hallmark_membership.rds")
    if (!file.exists(file)) stop("Missing Hallmark membership: ", file)
    return(list(membership = readRDS(file), input_file = file, package_version = NA_character_))
  }
  if (!requireNamespace("msigdbr", quietly = TRUE)) stop("msigdbr is required for other collections.")
  membership <- msigdbr::msigdbr(db_species = "HS", species = "Homo sapiens",
    collection = collection, subcollection = subcollection)
  list(membership = membership, input_file = character(),
    package_version = as.character(utils::packageVersion("msigdbr")))
}

ibd_select_pathways <- function(membership, min_genes = 0L, pathway_names = NULL) {
  if (!is.numeric(min_genes) || length(min_genes) != 1L || is.na(min_genes) ||
      !is.finite(min_genes) || min_genes < 0 || min_genes != floor(min_genes))
    stop("ibd.min_pathway_genes must be a nonnegative integer.")
  required <- c("gs_name", "gene_symbol", "ensembl_gene", "db_version")
  if (!all(required %in% names(membership)) || !nrow(membership)) stop("Incomplete pathway membership.")
  if (anyNA(membership$gs_name) || any(!nzchar(membership$gs_name))) stop("Missing pathway names.")
  # Count unique native human symbols, not repeated Ensembl mappings, before expression filtering.
  symbol <- if ("db_gene_symbol" %in% names(membership)) membership$db_gene_symbol else membership$gene_symbol
  if (anyNA(symbol) || any(!nzchar(symbol))) stop("Missing gene symbols in pathway membership.")
  unique_genes <- unique(data.frame(pathway = membership$gs_name, gene = symbol))
  sizes <- aggregate(list(n_database_genes = rep(1L, nrow(unique_genes))),
    list(pathway = unique_genes$pathway), sum)
  sizes <- sizes[order(sizes$pathway), ]
  sizes$minimum_genes <- min_genes
  sizes$passes_minimum <- sizes$n_database_genes >= min_genes
  if (!is.null(pathway_names)) {
    if (!is.character(pathway_names) || !length(pathway_names) || anyNA(pathway_names) ||
        anyDuplicated(pathway_names) || !all(pathway_names %in% sizes$pathway))
      stop("ibd.pathway_names must contain unique pathway names from the selected collection.")
    too_small <- sizes$pathway[sizes$pathway %in% pathway_names & !sizes$passes_minimum]
    if (length(too_small)) stop("Requested pathways below ibd.min_pathway_genes: ", paste(too_small, collapse = ", "))
  } else pathway_names <- sizes$pathway[sizes$passes_minimum]
  if (!length(pathway_names)) stop("No pathways meet ibd.min_pathway_genes.")
  sizes$selected <- sizes$pathway %in% pathway_names
  list(pathways = pathway_names, sizes = sizes)
}

ibd_pathway_exclusions <- function(pairs, pathway_ids) {
  out <- ibd_hallmark_exclusions(pairs, pathway_ids)
  names(out) <- sub("_in_hallmark$", "_in_collection", names(out))
  out$representative_rule <- gsub("Hallmark", "selected collection", out$representative_rule)
  out
}

ibd_redundancy_exclusions <- function(pairs) {
  required <- c("gene1", "gene2", "exon_bp1", "exon_bp2", "same_strand",
                "overlap_fraction1", "overlap_fraction2", "redundancy_candidate")
  if (!all(required %in% names(pairs))) stop("Incomplete redundancy audit table.")
  ids <- c(pairs$gene1, pairs$gene2)
  if (anyNA(ids) || any(!nzchar(ids)) || anyDuplicated(ids))
    stop("Redundancy audit must contain disjoint gene pairs with unique IDs.")
  if (anyNA(pairs[, required]) ||
      any(!pairs$same_strand | !pairs$redundancy_candidate) ||
      any(pairs$overlap_fraction1 < .95 | pairs$overlap_fraction2 < .95) ||
      any(!is.finite(pairs$exon_bp1) | !is.finite(pairs$exon_bp2)) ||
      any(pairs$exon_bp1 <= 0 | pairs$exon_bp2 <= 0))
    stop("Expected same-strand pairs with >=95% reciprocal exonic overlap.")
  keep_first <- pairs$exon_bp1 > pairs$exon_bp2 |
    (pairs$exon_bp1 == pairs$exon_bp2 & pairs$gene1 < pairs$gene2)
  pairs$retained_gene <- ifelse(keep_first, pairs$gene1, pairs$gene2)
  pairs$excluded_gene <- ifelse(keep_first, pairs$gene2, pairs$gene1)
  pairs$representative_rule <- rep("longer GTF exon union; ties: smaller gene ID", nrow(pairs))
  pairs
}

ibd_hallmark_exclusions <- function(pairs, hallmark_ids) {
  # The existing validator enforces disjoint, same-strand, >=95% overlap pairs.
  out <- ibd_redundancy_exclusions(pairs)
  stable_id <- function(x) sub("\\.[0-9]+$", "", x)
  out$gene1_in_hallmark <- stable_id(out$gene1) %in% hallmark_ids
  out$gene2_in_hallmark <- stable_id(out$gene2) %in% hallmark_ids
  prefer1 <- out$gene1_in_hallmark & !out$gene2_in_hallmark
  prefer2 <- out$gene2_in_hallmark & !out$gene1_in_hallmark
  out$retained_gene[prefer1] <- out$gene1[prefer1]
  out$excluded_gene[prefer1] <- out$gene2[prefer1]
  out$retained_gene[prefer2] <- out$gene2[prefer2]
  out$excluded_gene[prefer2] <- out$gene1[prefer2]
  out$representative_rule <- ifelse(prefer1 | prefer2, "only Hallmark member in pair",
    "both or neither in Hallmark: longer exon union; ties by gene ID")
  out
}

ibd_limma_age_sex <- function(log_cpm, metadata, batch_variable = NULL) {
  stopifnot(identical(colnames(log_cpm), metadata$sample_id),
            all(is.finite(log_cpm)), all(is.finite(metadata$age_at_diagnosis)),
            all(tolower(as.character(metadata$sex)) %in% c("female", "male")))
  design <- stats::model.matrix(~0 + analysis_group, metadata)
  age <- matrix(metadata$age_at_diagnosis - mean(metadata$age_at_diagnosis), ncol = 1,
                dimnames = list(metadata$sample_id, "age_centered"))
  sex <- factor(tolower(as.character(metadata$sex)), levels = c("female", "male"))
  covariates <- age
  if (length(unique(sex)) > 1L) {
    sex_design <- stats::model.matrix(~sex)[, -1, drop = FALSE]
    covariates <- cbind(age, sex_design)
  }
  batch <- NULL
  full <- cbind(design, covariates)
  if (!is.null(batch_variable)) {
    if (!batch_variable %in% names(metadata)) stop("Unknown batch variable.")
    value <- metadata[[batch_variable]]
    if (anyNA(value) || any(!nzchar(value)) || any(value == "Unknown")) stop("Missing batch values; specify an explicit missing-data policy.")
    batch <- factor(value)
    if (nlevels(batch) < 2L) stop("Batch variable has only one level.")
    full <- cbind(full, stats::model.matrix(~batch)[, -1, drop = FALSE])
  }
  if (qr(full)$rank < ncol(full) || nrow(full) <= ncol(full)) {
    stop("Age/sex/batch model is confounded or lacks residual degrees of freedom; correction was not applied.")
  }
  adjusted <- limma::removeBatchEffect(log_cpm, batch = batch, covariates = covariates, design = design)
  list(expression = adjusted, group_design = design, age_covariate = age, covariates = covariates,
       design_matrix = full,
       batch_variable = batch_variable, model_rank = qr(full)$rank,
       residual_df = nrow(full) - qr(full)$rank)
}


ibd_validate_comparisons <- function(comparisons) {
  if (!is.character(comparisons) || !length(comparisons) || anyNA(comparisons) ||
      anyDuplicated(comparisons) || !all(comparisons %in% c("CD", "iCD", "cCD", "UC")))
    stop("ibd.comparisons must contain unique values from CD, iCD, cCD, UC.")
  comparisons
}

ibd_load_count_inputs <- function(input) {
  files <- file.path(input, c("sample_metadata.csv", "raw/recount_metadata.rds",
    "raw/sra.gene_sums.SRP042228.G026.gz"))
  if (!all(file.exists(files))) stop("Missing input files; run the metadata audit and fetch_srp042228_counts.R: ",
    paste(files[!file.exists(files)], collapse = ", "))
  m <- utils::read.csv(files[1], stringsAsFactors = FALSE, na.strings = c("", "NA"))
  rc <- readRDS(files[2])
  required <- c("gsm", "subject_label", "diagnosis", "subtype", "age_at_diagnosis", "sex")
  if (!all(required %in% names(m))) stop("Incomplete audited sample metadata.")
  stopifnot(!anyNA(m$gsm), !anyNA(m$subject_label), !anyDuplicated(m$gsm),
    !anyDuplicated(m$subject_label), !anyDuplicated(rc$external_id),
    !anyDuplicated(rc$sra.sample_name), setequal(m$gsm, rc$sra.sample_name),
    all(m$diagnosis %in% c("CD", "UC", "Non-IBD control")),
    all(m$subtype[m$diagnosis == "CD"] %in% c("iCD", "cCD")),
    all(m$subtype[m$diagnosis != "CD"] == m$diagnosis[m$diagnosis != "CD"]))
  m <- m[match(rc$sra.sample_name, m$gsm), , drop = FALSE]
  rownames(m) <- NULL
  m$sample_id <- as.character(rc$external_id)
  m$patient_id <- m$subject_label
  m$sex <- tolower(trimws(m$sex))
  coverage <- recount3::read_counts(files[3], samples = m$sample_id)
  # Same calculation as recount3::compute_read_counts(round = TRUE).
  lengths <- as.numeric(rc[["recount_qc.star.average_mapped_length"]])
  stopifnot(identical(colnames(coverage), m$sample_id), !anyDuplicated(rownames(coverage)),
    all(is.finite(coverage)), all(coverage >= 0), length(lengths) == ncol(coverage),
    all(is.finite(lengths)), all(lengths > 0))
  counts <- round(sweep(coverage, 2, lengths, "/"))
  stopifnot(all(colSums(counts) > 0))
  list(counts = counts, metadata = m, files = files)
}

ibd_comparison_cohort <- function(metadata, comparison) {
  ibd_validate_comparisons(comparison)
  if (length(comparison) != 1L) stop("Select one comparison at a time.")
  disease <- if (comparison == "CD") metadata$diagnosis == "CD" else metadata$subtype == comparison
  control <- metadata$diagnosis == "Non-IBD control"
  complete_age <- is.finite(metadata$age_at_diagnosis) & metadata$age_at_diagnosis >= 0
  complete_sex <- !is.na(metadata$sex) & metadata$sex %in% c("female", "male")
  keep <- (disease | control) & complete_age & complete_sex
  audit <- data.frame(sample_id = metadata$sample_id, comparison = comparison,
    diagnosis = metadata$diagnosis, subtype = metadata$subtype, eligible_group = disease | control,
    complete_age = complete_age, complete_sex = complete_sex, retained = keep,
    reason = ifelse(!(disease | control), "other_disease_group",
      ifelse(!complete_age, "missing_age", ifelse(!complete_sex, "missing_sex", "retained"))))
  m <- metadata[keep, , drop = FALSE]
  m$analysis_group <- factor(ifelse(m$diagnosis == "Non-IBD control", "Control", "Disease"),
    levels = c("Control", "Disease"))
  list(metadata = m, audit = audit)
}

ibd_prepare_pathway_comparison <- function(counts, metadata, comparison, cooccur = 10L) {
  stopifnot(!anyDuplicated(metadata$sample_id), !anyDuplicated(metadata$patient_id),
    all(metadata$sample_id %in% colnames(counts)))
  cohort <- ibd_comparison_cohort(metadata, comparison)
  m <- cohort$metadata
  if (any(table(m$analysis_group) < cooccur)) stop("Too few complete-case samples in comparison: ", comparison)
  y <- edgeR::DGEList(counts[, m$sample_id, drop = FALSE])
  keep <- edgeR::filterByExpr(y, group = m$analysis_group, min.count = 10,
    min.total.count = 15, large.n = 10, min.prop = .7)
  if (sum(keep) < 2L) stop("Too few expressed genes.")
  y <- edgeR::calcNormFactors(y[keep, , keep.lib.sizes = FALSE], method = "TMM")
  log_cpm <- edgeR::cpm(y, log = TRUE, prior.count = 2)
  adjusted <- ibd_limma_age_sex(log_cpm, m)
  variances <- vapply(levels(m$analysis_group), function(g)
    apply(adjusted$expression[, m$analysis_group == g, drop = FALSE], 1, stats::var),
    numeric(nrow(adjusted$expression)))
  valid <- apply(variances, 1, function(x) all(is.finite(x) & x > 1e-10))
  list(metadata = m, cohort_audit = cohort$audit,
    filter = data.frame(gene_id = rownames(counts), retained = keep),
    variances = data.frame(gene_id = rownames(log_cpm), variances, valid = valid),
    log_cpm = log_cpm, adjustment = adjusted,
    expression = adjusted$expression[valid, , drop = FALSE],
    libraries = data.frame(sample_id = colnames(y), y$samples))
}

ibd_joint_gene_variation <- function(expression, groups, threshold = 1e-10) {
  stopifnot(ncol(expression) == length(groups), !anyNA(groups), all(table(groups) >= 2L))
  variances <- vapply(levels(groups), function(g)
    apply(expression[, groups == g, drop = FALSE], 1, stats::var), numeric(nrow(expression)))
  list(variances = variances, valid = apply(variances, 1, function(x) all(is.finite(x) & x > threshold)))
}

ibd_prepare_joint_cohort <- function(counts, metadata, cooccur = 10L) {
  groups <- c("Non-IBD control", "iCD", "cCD", "UC")
  stopifnot(!anyDuplicated(metadata$sample_id), !anyDuplicated(metadata$patient_id),
    all(metadata$sample_id %in% colnames(counts)))
  complete_age <- is.finite(metadata$age_at_diagnosis) & metadata$age_at_diagnosis >= 0
  complete_sex <- !is.na(metadata$sex) & metadata$sex %in% c("female", "male")
  eligible <- metadata$subtype %in% groups
  selected <- eligible & complete_age & complete_sex
  audit <- data.frame(sample_id = metadata$sample_id, subtype = metadata$subtype,
    complete_age = complete_age, complete_sex = complete_sex, retained = selected,
    reason = ifelse(!eligible, "other_group", ifelse(!complete_age, "missing_age",
      ifelse(!complete_sex, "missing_sex", "retained"))))
  m <- metadata[selected, , drop = FALSE]
  m$joint_group <- factor(m$subtype, levels = groups)
  if (any(table(m$joint_group) < cooccur))
    stop("Joint preprocessing requires at least cooccur complete-case samples in each of control, iCD, cCD, UC.")
  # Preserve the four disjoint disease-group effects, including when pooled CD is requested later.
  m$analysis_group <- m$joint_group
  y <- edgeR::DGEList(counts[, m$sample_id, drop = FALSE])
  keep <- edgeR::filterByExpr(y, group = m$joint_group, min.count = 10,
    min.total.count = 15, large.n = 10, min.prop = .7)
  if (sum(keep) < 2L) stop("Too few expressed genes after joint filtering.")
  y <- edgeR::calcNormFactors(y[keep, , keep.lib.sizes = FALSE], method = "TMM")
  log_cpm <- edgeR::cpm(y, log = TRUE, prior.count = 2)
  adjusted <- ibd_limma_age_sex(log_cpm, m)
  variation <- ibd_joint_gene_variation(adjusted$expression, m$joint_group)
  variances <- variation$variances
  valid <- variation$valid
  if (sum(valid) < 2L) stop("Too few genes with valid adjusted variation in every group.")
  adjusted$scope <- "joint_four_groups"
  adjusted$fitted_sample_ids <- m$sample_id
  list(metadata = m, cohort_audit = audit, preprocessing_mode = "joint_four_groups",
    filter = data.frame(gene_id = rownames(counts), retained = keep),
    variances = data.frame(gene_id = rownames(log_cpm), variances, valid = valid, check.names = FALSE),
    log_cpm = log_cpm, adjustment = adjusted, expression = adjusted$expression[valid, , drop = FALSE],
    libraries = data.frame(sample_id = colnames(y), y$samples))
}

ibd_subset_joint_comparison <- function(joint, comparison, cooccur = 10L) {
  stopifnot(identical(joint$preprocessing_mode, "joint_four_groups"))
  cohort <- ibd_comparison_cohort(joint$metadata, comparison)
  m <- cohort$metadata
  if (any(table(m$analysis_group) < cooccur)) stop("Too few samples for comparison: ", comparison)
  idx <- match(m$sample_id, joint$metadata$sample_id)
  # Subset the already fitted values, never refit/re-normalize within the comparison.
  adjusted <- joint$adjustment
  adjusted$expression <- adjusted$expression[, idx, drop = FALSE]
  for (nm in c("group_design", "age_covariate", "covariates", "design_matrix"))
    adjusted[[nm]] <- adjusted[[nm]][idx, , drop = FALSE]
  list(metadata = m, cohort_audit = cohort$audit, preprocessing_mode = "joint_four_groups",
    filter = joint$filter, variances = joint$variances,
    log_cpm = joint$log_cpm[, idx, drop = FALSE], adjustment = adjusted,
    expression = joint$expression[, idx, drop = FALSE],
    libraries = joint$libraries[match(m$sample_id, joint$libraries$sample_id), , drop = FALSE])
}

ibd_pathway_gene_audit <- function(hallmark, pathway, all_gene_ids, exclusions, prepared) {
  h <- hallmark[hallmark$gs_name == pathway, c("ensembl_gene", "gene_symbol")]
  h <- h[!duplicated(h$ensembl_gene), , drop = FALSE]
  lookup <- data.frame(ensembl_gene = sub("\\.[0-9]+$", "", all_gene_ids), gene_id = all_gene_ids)
  z <- merge(h, lookup, by = "ensembl_gene", all.x = TRUE, sort = FALSE)
  z$reason <- "eligible"
  z$reason[is.na(z$gene_id)] <- "not_mapped_to_expression_annotation"
  z$reason[z$gene_id %in% exclusions$excluded_gene] <- "annotation_redundant"
  filtered <- prepared$filter$gene_id[prepared$filter$retained]
  z$reason[z$reason == "eligible" & !z$gene_id %in% filtered] <- "low_expression"
  z$reason[z$reason == "eligible" & !z$gene_id %in% rownames(prepared$expression)] <- "no_adjusted_variation"
  z
}

ibd_run_pathway_panel <- function(prepared, audit, comparison, pathway, settings, directory,
                                  disco_fun = DISCO::disco) {
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(audit, file.path(directory, "gene_eligibility.csv"), row.names = FALSE)
  genes <- unique(audit$gene_id[audit$reason == "eligible"])
  result <- NULL
  error <- ""
  status <- "completed"
  start <- Sys.time()
  if (length(genes) < 2) {
    status <- "skipped_insufficient_genes"
  } else if (settings$preprocess_only) {
    status <- "prepared_only"
  } else {
    m <- prepared$metadata
    expr <- prepared$expression[genes, , drop = FALSE]
    result <- tryCatch(disco_fun(
      data1 = list(t(expr[, m$analysis_group == "Control", drop = FALSE])),
      data2 = list(t(expr[, m$analysis_group == "Disease", drop = FALSE])),
      cooccur = settings$cooccur, empirically_derive_threshold = settings$empirical,
      filter_threshold = settings$filter_threshold, q_value = settings$q_value,
      n_resamples = settings$n_resamples, random_seed = settings$seed,
      workers = settings$workers, adj_method = "BY", global_p_cutoff = .05,
      output_dir = directory, prefix = "pathway"), error = identity)
    if (inherits(result, "error")) {
      error <- conditionMessage(result)
      status <- "failed"
      result <- NULL
    } else saveRDS(result, file.path(directory, "disco_result.rds"))
  }
  score <- audit
  score$comparison <- comparison
  score$pathway <- pathway
  score$p <- score$q_within_pathway <- NA_real_
  score$n_partner_tests <- NA_integer_
  if (!is.null(result)) {
    idx <- match(score$gene_id, result$disco$species)
    score$p <- result$disco$p[idx]
    score$q_within_pathway <- result$disco$q[idx]
    score$n_partner_tests <- result$disco$num_corr[idx]
    score$reason[score$reason == "eligible"] <- ifelse(
      is.finite(score$p[score$reason == "eligible"]), "tested", "no_valid_partner_tests")
    if (!any(is.finite(score$p))) status <- "completed_no_testable_genes"
  } else score$reason[score$reason == "eligible"] <- status
  score <- score[order(score$q_within_pathway, na.last = TRUE), ]
  utils::write.csv(score, file.path(directory, "all_pathway_gene_scores.csv"), row.names = FALSE)
  summary <- data.frame(comparison, pathway, status, n_control = sum(prepared$metadata$analysis_group == "Control"),
    n_disease = sum(prepared$metadata$analysis_group == "Disease"), n_eligible_genes = length(genes),
    n_tested_genes = sum(is.finite(score$p)), n_significant_genes = sum(score$q_within_pathway < settings$q_value, na.rm = TRUE),
    selected_threshold = if (is.null(result)) NA_real_ else result$selected_threshold,
    global_p = if (is.null(result$global_test)) NA_real_ else result$global_test$global_p_value,
    elapsed_seconds = as.numeric(difftime(Sys.time(), start, units = "secs")), error = error)
  list(scores = score, summary = summary)
}
