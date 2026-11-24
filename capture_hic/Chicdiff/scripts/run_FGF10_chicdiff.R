#!/usr/bin/env Rscript
# Differential Capture Hi-C: IMR-90 (2 replicates) vs lymphoblasts (2 replicates).
# Arima nominal ~3-kb capture design; cis <=2 Mb peak matrix already prepared.
# This script runs differential TESTING; it does NOT regenerate CHiCAGO calls,
# .chinput counts, the peak matrix, or genome-wide Hi-C contact matrices.
# Usage: Rscript run_FGF10_chicdiff.R --root /path/to/Hi-C_3kb [--analysis-dir /path/to/results] [--force]

suppressPackageStartupMessages({
  library(Chicdiff)
  library(Chicago)
})

args <- commandArgs(trailingOnly = TRUE)
if (any(args %in% c('--help', '-h'))) {
  cat('Usage: Rscript run_FGF10_chicdiff.R --root <capture-hic-project-dir> [--analysis-dir <dir>] [--force]\n')
  cat('Example: --root /path/to/Hi-C_3kb\n')
  cat('Default: refuse to overwrite existing IMR_vs_Lymph_results.Rds.\n')
  quit(status = 0)
}
allowed <- c('--root', '--analysis-dir', '--force')
options <- args[grepl('^--', args)]
unknown <- setdiff(options, allowed)
if (length(unknown)) stop('Unknown option(s): ', paste(unknown, collapse=', '))
get_option <- function(flag, default = NULL) {
  idx <- which(args == flag)
  if (!length(idx)) return(default)
  if (length(idx) != 1L || idx >= length(args) || startsWith(args[idx + 1L], '--')) {
    stop('Expected one value after ', flag)
  }
  args[idx + 1L]
}
root <- get_option('--root')
if (is.null(root)) stop('Pass --root /path/to/Hi-C_3kb to locate original inputs')
root <- normalizePath(root, mustWork = TRUE)
force <- '--force' %in% args
reference <- file.path(root, 'reference')
analysis <- normalizePath(get_option('--analysis-dir', file.path(root, 'chicdiff_results')), mustWork = TRUE)
outdir <- file.path(analysis, 'differential')
peakfile <- file.path(analysis, 'IMR_Lymph.txt')
prefix <- 'IMR_vs_Lymph'

samples <- c('IMR1', 'IMR2', 'Lymph1', 'Lymph2')
chicago_path <- function(s) file.path(root, 'output', s, 'chicago', 'data', paste0(s, '.Rds'))
count_path <- function(s) file.path(root, 'output', s, 'bam2chicago', s, paste0(s, '.chinput'))

# Condition order is IMR then Lymph. DO NOT infer coefficient direction solely
# from ordering: verify the Chicdiff/DESeq2 contrast before describing log2FC.
chicagoData <- list(
  IMR = setNames(vapply(samples[1:2], chicago_path, character(1)), samples[1:2]),
  Lymph = setNames(vapply(samples[3:4], chicago_path, character(1)), samples[3:4])
)
countData <- list(
  IMR = setNames(vapply(samples[1:2], count_path, character(1)), samples[1:2]),
  Lymph = setNames(vapply(samples[3:4], count_path, character(1)), samples[3:4])
)

needed <- c(
  peakfile, unlist(chicagoData, use.names = FALSE),
  unlist(countData, use.names = FALSE),
  file.path(reference, 'hg38_chicago_input_3kb_1.baitmap'),
  file.path(reference, 'hg38_chicago_input_3kb_1.rmap')
)
missing <- needed[!file.exists(needed)]
if (length(missing)) stop('Missing input(s):\n', paste(missing, collapse = '\n'))

# The matrix header must refer to the same sample names as chicagoData/countData.
header <- names(read.delim(peakfile, sep = '\t', nrows = 1,
                           check.names = FALSE, stringsAsFactors = FALSE))
if (!all(samples %in% header)) {
  stop('Missing columns in peak matrix: ', paste(setdiff(samples, header), collapse = ', '))
}

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
results_path <- file.path(outdir, paste0(prefix, '_results.Rds'))
if (file.exists(results_path) && !force) {
  stop('Existing Chicdiff results found: ', results_path,
       '\nNothing rerun. Use --force ONLY if you deliberately want to overwrite results.')
}

cat('Capture Hi-C: Arima nominal ~3 kb; peaks selected upstream (cis <=2 Mb).\n')
cat('Groups: IMR1+IMR2 versus Lymph1+Lymph2\n')
cat('Peak matrix: ', peakfile, '\n', sep = '')
cat('Results directory: ', outdir, '\n', sep = '')

setwd(outdir)
settings <- setChicdiffExperiment(
  designDir = reference,
  chicagoData = chicagoData,
  countData = countData,
  peakfiles = peakfile,
  outprefix = prefix
)
results <- chicdiffPipeline(settings)

# Chicdiff itself saves prefix_results.Rds, prefix_countput.Rds,
# diagnostic plots and related files. The following TSV is portable.
dat <- as.data.frame(results)
if (!all(c('baitID','OEstart','OEend','log2FoldChange','padj') %in% names(dat))) {
  warning('Unexpected result columns: ', paste(names(dat), collapse = ', '))
}
write.table(dat, file = paste0(prefix, '_full_results.tsv'), sep = '\t',
            row.names = FALSE, quote = FALSE, na = 'NA')
saveRDS(results, file = paste0(prefix, '_full_results.Rds'))
writeLines(capture.output(sessionInfo()), con = 'sessionInfo.txt')
cat('\nCompleted: ', nrow(dat), ' tested differential interaction regions.\n', sep = '')
cat('Wrote ', file.path(outdir, paste0(prefix, '_full_results.tsv')), '\n', sep = '')
cat('Note: differential testing was genome-wide across the captured peak set;\n')
cat('extract FGF10-region/bait rows in a separate script.\n')
