# Rscript analysis/ibd_transcriptomics/code/explore_srp042228_metadata.R METADATA_DIRECTORY
# Metadata-only investigation; no expression download or DISCO run.
parse_geo_samples <- function(lines) {
  starts <- grep("^\\^SAMPLE = ", lines)
  if (!length(starts)) stop("No GEO sample records found.")
  out <- lapply(seq_along(starts), function(i) {
    end <- if (i < length(starts)) starts[i + 1L] - 1L else length(lines)
    block <- lines[seq.int(starts[i], end)]
    block <- block[grepl("^!Sample_.* = ", block)]
    data.frame(gsm = sub("^\\^SAMPLE = ", "", lines[starts[i]]),
      field = sub(" = .*", "", block), value = sub("^[^=]+ = ", "", block))
  })
  do.call(rbind, out)
}

audit_srp042228 <- function(output) {
  write_table <- function(d, name) utils::write.csv(d, file.path(output, name), row.names = FALSE, na = "")
  raw <- file.path(output, "raw")
  con <- gzfile(file.path(raw, "GSE57945_family.soft.gz"), "rt")
  lines <- tryCatch(readLines(con, warn = FALSE), finally = close(con))
  fields <- parse_geo_samples(lines)
  write_table(fields, "geo_sample_fields.csv")
  series <- lines[grepl("^!Series_.* = ", lines)]
  write_table(data.frame(field = sub(" = .*", "", series), value = sub("^[^=]+ = ", "", series)), "geo_series_fields.csv")
  ids <- unique(fields$gsm)
  values <- function(id, field) fields$value[fields$gsm == id & fields$field == field]
  characteristic <- function(id, key) {
    z <- values(id, "!Sample_characteristics_ch1")
    z <- z[startsWith(z, paste0(key, ": "))]
    if (length(z) != 1L) stop("Missing or repeated characteristic: ", id, " / ", key)
    substring(z, nchar(key) + 3L)
  }
  samples <- data.frame(gsm = ids, title = vapply(ids, function(id) values(id, "!Sample_title"), ""))
  samples$subject_label <- vapply(ids, function(id) {
    z <- values(id, "!Sample_description")
    z <- z[grepl("^CCFA_Risk_[0-9]+$", z)]
    if (length(z) != 1L) stop("Unrecognized subject label for ", id)
    z
  }, "")
  mapping <- c(tissue = "tissue", sex = "Sex", age_at_diagnosis = "age at diagnosis",
    paris_age = "paris age", diagnosis = "diagnosis", subtype = "l2 type",
    histopathology = "histopathology", deep_ulcer = "deep ulcer")
  for (nm in names(mapping)) samples[[nm]] <- vapply(ids, characteristic, "", key = mapping[[nm]])
  samples$sex <- tolower(samples$sex)
  samples$diagnosis[tolower(samples$diagnosis) == "not ibd"] <- "Non-IBD control"
  samples$subtype[tolower(samples$subtype) == "not ibd"] <- "Non-IBD control"
  samples$deep_ulcer <- tolower(samples$deep_ulcer)
  for (nm in c("histopathology", "deep_ulcer")) samples[[nm]][tolower(samples[[nm]]) == "na"] <- NA_character_
  samples$age_at_diagnosis <- as.numeric(samples$age_at_diagnosis)
  samples$complete_age_sex <- is.finite(samples$age_at_diagnosis) & samples$sex %in% c("male", "female")
  stopifnot(!anyDuplicated(samples$gsm), !anyDuplicated(samples$subject_label),
    all(samples$diagnosis %in% c("CD", "UC", "Non-IBD control")))
  ena <- utils::read.delim(file.path(raw, "ena_runs.tsv"), stringsAsFactors = FALSE)
  recount <- readRDS(file.path(raw, "recount_metadata.rds"))
  stopifnot(!anyDuplicated(ena$run_accession), !anyDuplicated(recount$external_id),
    setequal(ids, ena$sample_alias), all(recount$external_id %in% ena$run_accession),
    all(recount$sra.sample_name == ena$sample_alias[match(recount$external_id, ena$run_accession)]))
  samples$n_ena_runs <- as.integer(table(factor(ena$sample_alias, levels = ids)))
  samples$n_recount_runs <- as.integer(table(factor(recount$sra.sample_name, levels = ids)))
  write_table(samples, "sample_metadata.csv")
  runs <- cbind(ena, samples[match(ena$sample_alias, samples$gsm), setdiff(names(samples), c("gsm", "title")), drop = FALSE])
  rownames(runs) <- NULL
  runs$in_recount3 <- ena$run_accession %in% recount$external_id
  for (nm in grep("^recount_(qc|seq_qc)\\.", names(recount), value = TRUE))
    runs[[nm]] <- recount[[nm]][match(runs$run_accession, recount$external_id)]
  write_table(runs, "run_metadata.csv")
  write_table(samples[samples$n_ena_runs > 1L, ], "technical_replicates.csv")
  write_table(samples[samples$n_recount_runs == 0L, ], "missing_recount_samples.csv")
  summarize <- function(group) do.call(rbind, lapply(unique(group), function(g) {
    s <- samples[group == g, ]
    data.frame(group = g, samples = nrow(s), complete_age_sex = sum(s$complete_age_sex),
      recount3_samples = sum(s$n_recount_runs > 0), sequencing_runs = sum(s$n_ena_runs),
      female = sum(s$sex == "female"), male = sum(s$sex == "male"),
      median_age_at_diagnosis = median(s$age_at_diagnosis),
      minimum_age = min(s$age_at_diagnosis), maximum_age = max(s$age_at_diagnosis))
  }))
  counts <- summarize(samples$diagnosis)
  subtype_counts <- summarize(samples$subtype)
  write_table(counts, "diagnosis_sample_counts.csv")
  write_table(subtype_counts, "subtype_sample_counts.csv")
  comparisons <- do.call(rbind, lapply(c("CD", "UC", "iCD", "cCD"), function(g) {
    case <- if (g %in% c("CD", "UC")) samples$diagnosis == g else samples$subtype == g
    control <- samples$diagnosis == "Non-IBD control"
    ncase <- sum(case & samples$complete_age_sex)
    ncontrol <- sum(control & samples$complete_age_sex)
    data.frame(comparison = paste(g, "versus non-IBD"), disease_samples = ncase,
      control_samples = ncontrol, total_complete_age_sex = ncase + ncontrol,
      shared_batches = NA_integer_, annotated_disease_control_pairs = 0L)
  }))
  write_table(comparisons, "candidate_comparisons.csv")
  covariates <- data.frame(variable = c("Age at diagnosis", "Sex", "Diagnosis", "iCD/cCD/UC/non-IBD subtype",
    "Paris age class", "Histopathology", "Deep ulcer", "Tissue", "Treatment-naive status",
    "Recruitment site", "Sequencing batch / flowcell / lane", "Ancestry", "PCDAI", "Antibiotic exposure", "RIN"),
    n_annotated_samples = c(322, 322, 322, 322, 322, 218, 218, 322, NA, 0, 0, 0, 0, 0, 0),
    scope = c("GEO labels this age at diagnosis, including controls", "Male/female",
      "CD/UC/non-IBD", "iCD/cCD within Crohn disease", "Derived age category; avoid redundancy with continuous age",
      "CD only; 3 Undetermined, 215 definite classifications", "CD only; 76 yes, 142 no",
      "Constant: ileal biopsy", "Study-level design; no varying individual treatment field",
      rep("Not supplied in downloaded GEO characteristics", 6)))
  write_table(covariates, "covariate_availability.csv")
  for (nm in c("histopathology", "deep_ulcer", "paris_age")) {
    tab <- as.data.frame(table(diagnosis = samples$diagnosis, subtype = samples$subtype,
      value = samples[[nm]], useNA = "ifany"))
    write_table(tab[tab$Freq > 0, ], paste0(nm, "_counts.csv"))
  }
  completeness <- function(d) do.call(rbind, lapply(names(d), function(nm) {
    z <- as.character(d[[nm]])
    ok <- !is.na(z) & nzchar(trimws(z)) & !tolower(z) %in% c("na", "not provided", "missing")
    data.frame(field = nm, n_available = sum(ok), total = length(z), n_unique = length(unique(z[ok])))
  }))
  write_table(completeness(samples), "sample_covariate_completeness.csv")
  write_table(completeness(recount), "recount_field_completeness.csv")
  write_table(completeness(ena), "ena_field_completeness.csv")
  # Compare clinical metadata independently of accession matching. Ages allow the
  # two-decimal rounding used by GEO; categorical comparisons ignore letter case.
  discrepancies <- list()
  for (i in seq_len(nrow(recount))) {
    si <- match(recount$sra.sample_name[i], samples$gsm)
    pieces <- strsplit(recount$sra.sample_attributes[i], "|", fixed = TRUE)[[1]]
    attrs <- setNames(sub("^.*?;;", "", pieces), sub(";;.*$", "", pieces))
    for (nm in names(mapping)) {
      a <- attrs[mapping[[nm]]]
      b <- as.character(samples[[nm]][si])
      if (nm %in% c("diagnosis", "subtype") && !is.na(a) && tolower(a) == "not ibd") a <- "Non-IBD control"
      if (!is.na(a) && tolower(a) == "na") a <- NA_character_
      same <- if (is.na(a) || is.na(b)) is.na(a) && is.na(b) else if (nm == "age_at_diagnosis")
        abs(as.numeric(a) - as.numeric(b)) <= .011 else tolower(a) == tolower(b)
      if (!same) discrepancies[[length(discrepancies) + 1L]] <- data.frame(
        gsm = samples$gsm[si], run = recount$external_id[i], field = nm, geo = b, recount3 = unname(a))
    }
  }
  discrep <- if (length(discrepancies)) do.call(rbind, discrepancies) else
    data.frame(gsm = character(), run = character(), field = character(), geo = character(), recount3 = character())
  write_table(discrep, "geo_recount_annotation_discrepancies.csv")
  write_table(data.frame(documented_disease_control_pairs = 0L,
    repeated_subject_labels = sum(duplicated(samples$subject_label)),
    repeated_run_samples = sum(samples$n_ena_runs > 1L),
    interpretation = "Cross-sectional diagnostic ileal biopsies; normal ileum in CD is not a matched non-IBD control"), "pairing_audit.csv")
  qc <- do.call(rbind, lapply(counts$group, function(g) {
    s <- runs[runs$diagnosis == g, ]
    do.call(rbind, lapply(c("recount_qc.star.number_of_input_reads", "recount_qc.star.uniquely_mapped_reads_%"), function(nm) {
      x <- as.numeric(sub("%$", "", s[[nm]])); x <- x[is.finite(x)]
      data.frame(group = g, metric = nm, n = length(x), min = min(x), median = median(x), max = max(x))
    }))
  }))
  write_table(qc, "run_qc_summary.csv")
  grDevices::pdf(file.path(output, "metadata_exploration.pdf"), width = 11, height = 8.5)
  tryCatch({
    par(mfrow = c(2, 2), mar = c(5, 5, 4, 1))
    cols <- c(CD = "#D9534F", UC = "#8B60AD", "Non-IBD control" = "#2878B5")
    barplot(setNames(counts$samples, counts$group), col = cols[counts$group],
      main = "Current GEO samples: GSE57945", ylab = "Biopsies")
    boxplot(age_at_diagnosis ~ diagnosis, data = samples, main = "Age at diagnosis (GEO field)",
      col = cols[sort(unique(samples$diagnosis))], ylab = "Years", xlab = "")
    barplot(setNames(subtype_counts$samples, subtype_counts$group), main = "Clinical subtypes", ylab = "Biopsies")
    y <- as.numeric(runs[["recount_qc.star.number_of_input_reads"]]) / 1e6
    boxplot(y ~ runs$diagnosis, main = "recount3 read depth", ylab = "Million reads", xlab = "",
      col = cols[sort(unique(samples$diagnosis))])
  }, finally = grDevices::dev.off())
  writeLines(c("SRP042228 / GSE57945: metadata-only exploration.",
    "Current GEO contains 322 samples; its series summary still quotes the original 359-person cohort.",
    "Use current sample-level labels. 2015 corrigendum DOI:10.1172/JCI79657 reports relabeling and exclusions.",
    "The correction's historical exclusions were not reapplied to the current GEO release.",
    "No validated disease/control patient pairs; subject labels are unique, but this does not genetically verify identity.",
    "Non-IBD controls were symptomatic children investigated for IBD with normal diagnostic findings, not population healthy volunteers.",
    "All tissue is ileal biopsy. cCD denotes colon-only Crohn disease with clinically uninvolved ileum.",
    "Age is explicitly labeled age at diagnosis, including controls; no separate biopsy-age field is provided.",
    "Sex, diagnosis, Paris age class, and iCD/cCD subtype are available for all samples.",
    "Deep ulcer and histopathology annotations apply to CD; NA in UC/controls must not be recoded as no inflammation.",
    "Three CD histopathology values are Undetermined. Recorded availability is not the same as definite classification.",
    "No per-sample sequencing batch, flowcell, lane, recruitment site, ancestry, PCDAI, antibiotic exposure or RIN in GEO characteristics.",
    "Sequencing release dates are administrative dates and must not be used as sequencing batches.",
    "All samples use the same platform/layout; these constants cannot identify batch effects.",
    "GEO provides processed baselined/log-transformed normalized expression; do not blindly log-transform it again.",
    "For new DISCO analyses, uniform recount3 counts with explicit normalization are preferable to mixing processed sources.",
    "Potential comparisons: CD vs non-IBD; UC vs non-IBD; iCD/cCD separately vs non-IBD; ulceration within CD.",
    "Age/sex adjustment is feasible. Inflammation should be a planned stratification or sensitivity analysis, because it can be part of the disease effect.",
    "Shared-batch counts cannot be established without batch identifiers. No samples excluded in this metadata audit.",
    paste0("GEO/recount3 clinical discrepancies (rounding tolerance 0.011 years): ", nrow(discrep)),
    "Original article: https://www.jci.org/articles/view/75436",
    "Correction: https://www.jci.org/articles/view/79657"), file.path(output, "exploration_notes.txt"))
  capture.output(sessionInfo(), file = file.path(output, "exploration_session_info.txt"))
  print(counts)
  print(subtype_counts)
  message("Clinical annotation discrepancies: ", nrow(discrep))
  invisible(list(samples = samples, counts = counts, subtypes = subtype_counts, discrepancies = discrep))
}

if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) != 1L) stop("Supply the output directory from fetch_srp042228_metadata.R.")
  audit_srp042228(args[1])
}
