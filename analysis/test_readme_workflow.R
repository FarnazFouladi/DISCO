# Offline integration checks: Rscript analysis/test_readme_workflow.R
# Uses temporary results and synthetic counts; does not run real-data permutations.
repo <- normalizePath('.')
source('analysis/ibd_transcriptomics/code/ibd_result_files.R')
work <- tempfile('disco-readme-tests-')
dir.create(work)

# Exercise the actual IBD runner twice against a small synthetic cohort.
output_dir <- file.path(work, 'runner')
input_dir <- file.path(output_dir, 'inputs')
reference_dir <- file.path(input_dir, 'reference')
dir.create(file.path(input_dir, 'raw'), recursive = TRUE)
dir.create(reference_dir)
write.csv(data.frame(note = 'synthetic fixture'), file.path(input_dir, 'sample_metadata.csv'), row.names = FALSE)
saveRDS(list(note = 'synthetic fixture'), file.path(input_dir, 'raw/recount_metadata.rds'))
writeBin(as.raw(0), file.path(input_dir, 'raw/sra.gene_sums.SRP042228.G026.gz'))
reference <- 'analysis/ibd_transcriptomics/results/inputs/reference'
membership <- readRDS(file.path(reference, 'hallmark_membership.rds'))
pathway <- 'HALLMARK_TNFA_SIGNALING_VIA_NFKB'
membership <- membership[membership$gs_name == pathway, ]
membership <- membership[!duplicated(membership$ensembl_gene), ][1:5, ]
saveRDS(membership, file.path(reference_dir, 'hallmark_membership.rds'))
redundancy <- read.csv(file.path(reference, 'redundancy_candidates_all_genes.csv'))
write.csv(redundancy[FALSE, ], file.path(reference_dir, 'redundancy_candidates_all_genes.csv'), row.names = FALSE)
set.seed(207)
metadata <- data.frame(sample_id = paste0('sample', 1:48), patient_id = paste0('patient', 1:48),
  subtype = rep(c('Non-IBD control', 'iCD', 'cCD', 'UC'), each = 12),
  diagnosis = rep(c('Non-IBD control', 'CD', 'CD', 'UC'), each = 12),
  age_at_diagnosis = runif(48, 5, 16), sex = rep(c('female', 'male'), 24))
counts <- matrix(rnbinom(5 * 48, mu = 100, size = 3), 5,
  dimnames = list(paste0(membership$ensembl_gene, '.1'), metadata$sample_id))
input_hashes <- tools::md5sum(list.files(input_dir, recursive = TRUE, full.names = TRUE))
old_options <- options(ibd.exploration_dir = input_dir, ibd.pathway_reference = reference_dir,
  ibd.output_dir = output_dir, ibd.pathway_names = pathway, ibd.min_pathway_genes = 0L,
  ibd.pathway_preprocess_only = FALSE, ibd.comparisons = c('iCD', 'cCD', 'UC'),
  ibd.preprocessing_mode = 'joint_four_groups', ibd.msigdb_collection = 'H',
  ibd.msigdb_subcollection = NULL, ibd.disco_workers = 1L)
runner <- new.env(parent = globalenv())
runner$source <- function(file, ...) {
  sys.source(file, envir = runner)
  if (basename(file) == 'ibd_hallmark_disco_helpers.R')
    runner$ibd_load_count_inputs <- function(input) list(counts = counts, metadata = metadata)
}
run_ibd <- function() suppressWarnings(sys.source('analysis/ibd_transcriptomics/code/run_ibd_hallmark_disco.R', runner))
run_ibd()
first <- read.csv(file.path(output_dir, 'all_pathway_gene_scores.csv'))
stopifnot(all(first$reason == 'tested'), readRDS(file.path(output_dir, 'completion.rds'))$n_failed == 0)
dir.create(file.path(output_dir, 'comparison_CD_vs_control'))
writeLines('old result', file.path(output_dir, 'comparison_CD_vs_control/stale.txt'))
dir.create(file.path(output_dir, 'figure'))
writeLines('old plot', file.path(output_dir, 'figure/stale.txt'))
run_ibd()
stopifnot(identical(first, read.csv(file.path(output_dir, 'all_pathway_gene_scores.csv'))),
  !dir.exists(file.path(output_dir, 'comparison_CD_vs_control')),
  !dir.exists(file.path(output_dir, 'figure')),
  identical(input_hashes, tools::md5sum(names(input_hashes))))
options(old_options)
message('Repeated IBD runs replace generated results, preserve inputs, and retain identical numerical results.')

# Verify the maintainer command runs the complete pipeline and stops on failure.
commands <- list()
maintainer <- new.env(parent = globalenv())
maintainer$system2 <- function(command, args) { commands[[length(commands) + 1L]] <<- args; 0L }
sys.source('analysis/rebuild_readme.R', maintainer)
stopifnot(length(commands) == 4L, grepl('render_readme.R', commands[[4L]][2], fixed = TRUE))
n_commands <- 0L
maintainer$system2 <- function(command, args) { n_commands <<- n_commands + 1L; 1L }
error <- tryCatch(sys.source('analysis/rebuild_readme.R', maintainer), error = identity)
stopifnot(inherits(error, 'error'), n_commands == 1L)
message('Maintainer pipeline runs in order and stops before rendering after a failure.')

# Render an isolated README after changing saved counts and significance values.
fixture <- file.path(work, 'render')
dir.create(file.path(fixture, 'man', 'figures'), recursive = TRUE)
copy_file <- function(relative) {
  dest <- file.path(fixture, relative)
  dir.create(dirname(dest), recursive = TRUE, showWarnings = FALSE)
  stopifnot(file.copy(relative, dest, overwrite = TRUE))
}
files <- c('README.Rmd', 'analysis/render_readme.R',
  list.files('analysis/diabetes_microbiome/code', full.names = TRUE),
  'analysis/diabetes_microbiome/results/README-results.rds',
  'analysis/diabetes_microbiome/results/taxa/disco_analysis.rds',
  'analysis/diabetes_microbiome/results/taxa_pathways/disco_analysis.rds',
  'analysis/diabetes_microbiome/results/taxa/disco_inputs.rds',
  'analysis/diabetes_microbiome/results/taxa_pathways/disco_inputs.rds',
  'analysis/diabetes_microbiome/results/differential_scatterplots/top10_pairs.csv',
  'analysis/diabetes_microbiome/results/taxa_pathway_scatterplots/top1_pairs.csv',
  list.files('analysis/ibd_transcriptomics/code', full.names = TRUE),
  list.files('analysis/ibd_transcriptomics/results', full.names = TRUE),
  list.files('analysis/ibd_transcriptomics/results/figure', full.names = TRUE))
# The renderer's image sources, plus the README chunk sources above.
files <- c(files, list.files('analysis/diabetes_microbiome/results', pattern = '[.]png$', full.names = TRUE, recursive = TRUE))
for (file in unique(files[file.exists(files) & !dir.exists(files)])) copy_file(file)
ibd <- file.path(fixture, 'analysis/ibd_transcriptomics/results')
scores <- read.csv(file.path(ibd, 'all_pathway_gene_scores.csv'))
rows <- read.csv(file.path(ibd, 'figure/heatmap_row_order.csv'))
extra <- setdiff(unique(scores$pathway), rows$pathway)[1]
changed <- scores$comparison == 'UC' & scores$reason == 'tested' &
  scores$pathway %in% c(extra, 'HALLMARK_INFLAMMATORY_RESPONSE')
scores$q_within_pathway[changed] <- 0
write.csv(scores, file.path(ibd, 'all_pathway_gene_scores.csv'), row.names = FALSE)
panels <- read.csv(file.path(ibd, 'pathway_analysis_summary.csv'))
panels$n_control <- 41L
for (i in seq_len(nrow(panels))) {
  z <- scores[scores$comparison == panels$comparison[i] & scores$pathway == panels$pathway[i], ]
  panels$n_significant_genes[i] <- sum(z$reason == 'tested' & is.finite(z$q_within_pathway) & z$q_within_pathway < .05)
}
write.csv(panels, file.path(ibd, 'pathway_analysis_summary.csv'), row.names = FALSE)
samples <- read.csv(file.path(ibd, 'comparison_sample_counts.csv'))
samples$n_control <- 41L
write.csv(samples, file.path(ibd, 'comparison_sample_counts.csv'), row.names = FALSE)
# Corrupt the previous display to ensure rendering actually rebuilds it.
writeLines('stale figure', file.path(ibd, 'figure/pathway_significance_heatmap.png'))
log <- system2(file.path(R.home('bin'), 'Rscript'),
  c('--vanilla', shQuote(file.path(fixture, 'analysis/render_readme.R'))), stdout = TRUE, stderr = TRUE)
if (!is.null(attr(log, 'status'))) stop(paste(log, collapse = '\n'))
markdown <- paste(readLines(file.path(fixture, 'README.md')), collapse = '\n')
stopifnot(grepl('24 pathways', markdown, fixed = TRUE),
  grepl('100.0% in UC', markdown, fixed = TRUE),
  grepl('\\|[[:space:]]*41[[:space:]]*\\|', markdown),
  all(read.csv(file.path(ibd, 'comparison_sample_counts.csv'))$n_control == 41L),
  file.info(file.path(fixture, 'man/figures/README-ibd-hallmark-significance.png'))$size > 1000,
  !grepl('November 2, 2026', markdown, fixed = TRUE))
message('Rendering rebuilds the IBD figure and reflects changed pathway counts and percentages in README.md.')
