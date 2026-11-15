#!/usr/bin/env Rscript
# FGF10 targeted plots from existing Chicdiff tests and exact fragment-pair tables.
# Arima nominal ~3-kb Capture Hi-C design. NO differential tests rerun.
# Run extract_FGF10_bait_TSVs.py first to generate portable source TSVs.
# All custom plots can be recreated from the exported TSVs alone.
# Native Chicdiff mirrored profiles are optional and need local result RDS files.

suppressPackageStartupMessages(library(ggplot2))

root <- 'file.path'
tables <- file.path(root, 'FGF10_tables')
output <- file.path(root, 'FGF10_clean_plots')
dir.create(output, recursive = TRUE, showWarnings = FALSE)

prefix <- 'IMR_vs_Lymph'
results_rds <- file.path(root, paste0(prefix, '_results.Rds'))
counts_rds <- file.path(root, paste0(prefix, '_countput.Rds'))
baitmap <- 'file.path/hg38_chicago_input_3kb_1.baitmap'
exact_tsv <- file.path(tables, 'FGF10_exact_SNP_fragment_pairs_counts_scores.tsv')
regions_tsv <- file.path(tables, 'FGF10_all_3_baits_SNP_overlap.tsv')
snp <- 43619485L

bait_labels <- c(
  '282949' = '282949 | FGF10 promoter-overlapping',
  '282950' = '282950 | FGF10 promoter-overlapping',
  '282953' = '282953 | FGF10 promoter-proximal'
)
native_names <- c(
  '282949' = 'FGF10_promoter_bait282949_Chicdiff.png',
  '282950' = 'FGF10_promoter_bait282950_Chicdiff.png',
  '282953' = 'FGF10_proximal_bait282953_Chicdiff.png'
)

save_plot <- function(plot, file, width, height) {
  ggsave(file.path(output, file), plot = plot, width = width, height = height,
         dpi = 250, bg = 'white', limitsize = FALSE)
  message('Saved ', file.path(output, file))
}

# A: Native Chicdiff mirrored mean profiles (one FGF10-associated bait per file).
# The statistical model uses IMR1/IMR2 and Lymph1/Lymph2 independently.
# The native profile is a group-aggregated display, NOT four separate profiles.
if (all(file.exists(c(results_rds, counts_rds, baitmap))) &&
    requireNamespace('Chicdiff', quietly = TRUE)) {
  for (id in names(native_names)) {
    path <- file.path(output, native_names[[id]])
    # Open the device, explicitly draw plotDiffBaits() return value,
    # close device, and only then report successful output.
    png(path, width = 1900, height = 1250, res = 180, bg = 'white')
    tryCatch({
      p <- Chicdiff::plotDiffBaits(results_rds, counts_rds, baitmap,
                                  baits = as.numeric(id))
      if (is.null(p)) stop('plotDiffBaits returned NULL for bait ', id)
      if (inherits(p, 'grob') || inherits(p, 'gTree') || inherits(p, 'gtable')) {
        grid::grid.draw(p)
      } else {
        print(p)  # essential inside a loop
      }
    }, error = function(e) {
      stop('Failed to render bait ', id, ': ', conditionMessage(e))
    }, finally = dev.off())
    if (!file.exists(path) || is.na(file.info(path)$size) || file.info(path)$size == 0) {
      stop('PNG was not generated: ', path)
    }
    message('Saved ', path)
  }
} else {
  message('Skipping native Chicdiff mirrored profiles (need Chicdiff, RDS results, and baitmap).')
}

# B: Individual fragment-pair CHiCAGO scores and raw counts for all four replicates.
if (file.exists(exact_tsv)) {
  d <- read.delim(exact_tsv, sep = '\t', check.names = FALSE, stringsAsFactors = FALSE)
  required <- c('sample', 'baitID', 'otherEndID', 'N', 'score')
  if (!all(required %in% names(d))) {
    stop('Exact-pair TSV missing: ', paste(setdiff(required, names(d)), collapse = ', '))
  }
  d <- d[as.character(d$baitID) %in% names(bait_labels) &
           d$otherEndID == 282709, , drop = FALSE]
  if (nrow(d) == 0) stop('No rows for SNP-containing fragment 282709.')
  d$sample <- factor(d$sample, levels = c('IMR1','IMR2','Lymph1','Lymph2'))
  d$bait <- factor(as.character(d$baitID), levels = names(bait_labels),
                   labels = unname(bait_labels))
  d$N <- as.numeric(d$N)
  d$score <- as.numeric(d$score)
  sample_colors <- c(IMR1='#31578B', IMR2='#7FA4D6',
                     Lymph1='#C46B48', Lymph2='#E8AF81')
  theme_fig <- theme_bw(base_size = 12) + theme(
    legend.position='none', panel.grid.minor = element_blank(),
    strip.text = element_text(face='bold', size=10),
    axis.text.x=element_text(angle=30, hjust=1),
    plot.title=element_text(face='bold', size=14),
    plot.subtitle=element_text(size=10)
  )
  subtitle <- 'Each pair includes SNP-containing restriction fragment 282709 (chr5:43,616,691–43,619,914; GRCh38)'
  p_counts <- ggplot(d, aes(x=sample,y=N,fill=sample)) +
    geom_col(width=.67) +
    geom_text(aes(label=N), vjust=-.35, size=3.6) +
    facet_wrap(~bait,nrow=1) + scale_fill_manual(values=sample_colors) +
    scale_y_continuous(expand=expansion(mult=c(0,.20))) +
    labs(title='FGF10–SNP fragment: observed contacts in each replicate',
         subtitle=subtitle, x=NULL, y='Unnormalized interaction read count (N)') + theme_fig
  p_scores <- ggplot(d, aes(x=sample,y=score,fill=sample)) +
    geom_col(width=.67) +
    geom_hline(yintercept=5,linetype='dashed',color='#444444',linewidth=.5) +
    geom_text(aes(label=sprintf('%.2f',score)),vjust=-.35,size=3.6) +
    facet_wrap(~bait,nrow=1) + scale_fill_manual(values=sample_colors) +
    scale_y_continuous(expand=expansion(mult=c(0,.20))) +
    labs(title='FGF10–SNP fragment: CHiCAGO score in each replicate',
         subtitle=paste0(subtitle,' | Dashed line: CHiCAGO score 5, NOT a differential p-value'),
         x=NULL,y='CHiCAGO interaction score') + theme_fig
  save_plot(p_counts,'FGF10_SNP_fragment_four_replicates_counts.png',13.7,5.6)
  save_plot(p_scores,'FGF10_SNP_fragment_four_replicates_scores.png',13.7,5.6)
} else {
  message('Skipping replicate plots: ', exact_tsv, ' not found (run TSV extractor with --exact).')
}

# C: SNP-overlapping Chicdiff region-level statistics for all 3 baits.
# Multiple overlapping windows per bait share counts; NOT independent loops.
if (file.exists(regions_tsv)) {
  z <- read.delim(regions_tsv, sep='\t', check.names=FALSE, stringsAsFactors=FALSE)
  req <- c('baitID','regionID','OEchr','OEstart','OEend','log2FoldChange','weighted_padj')
  if (!all(req %in% names(z))) stop('Regional TSV missing: ',paste(setdiff(req,names(z)),collapse=', '))
  z <- z[as.character(z$baitID) %in% names(bait_labels) & z$OEchr == 'chr5' &
           z$OEstart <= snp & z$OEend >= snp, , drop=FALSE]
  if (nrow(z)) {
    z <- z[order(z$baitID,z$OEstart,z$regionID), , drop=FALSE]
    z$bait <- factor(as.character(z$baitID),levels=names(bait_labels),labels=unname(bait_labels))
    z$region <- factor(paste0('Region ',z$regionID),levels=rev(unique(paste0('Region ',z$regionID))))
    z$annotation <- sprintf('log2FC %.2f; FDR %.2g',z$log2FoldChange,z$weighted_padj)
    xleft <- min(z$OEstart)-5000
    xright <- max(z$OEend)+27000
    p_regions <- ggplot(z, aes(y=region)) +
      geom_segment(aes(x=OEstart,xend=OEend,yend=region),
                   color='#7798BE',linewidth=6,lineend='butt') +
      geom_vline(xintercept=snp,linetype='dashed',color='#AD343B',linewidth=.6) +
      geom_text(aes(x=OEend+2600,label=annotation),hjust=0,size=3.4) +
      facet_grid(bait~.,scales='free_y',space='free_y') +
      scale_x_continuous(labels=function(x) format(x,big.mark=',',scientific=FALSE)) +
      coord_cartesian(xlim=c(xleft,xright),clip='off') +
      theme_bw(base_size=12) +
      theme(panel.grid.minor=element_blank(),strip.text.y=element_text(angle=0),
            plot.title=element_text(face='bold'),axis.text.x=element_text(angle=20,hjust=1)) +
      labs(title='Differential FGF10 contacts in regions containing rs78368856',
           subtitle='Red dashed line: chr5:43,619,485; overlapping test windows are not independent loops',
           x='Genomic position on chr5 (GRCh38)',y=NULL)
    save_plot(p_regions,'FGF10_SNP_Chicdiff_all_three_baits_regions.png',13.7,7.5)
    # Preserve the familiar single-bait 282953 plot.
    single <- z[z$baitID==282953,,drop=FALSE]
    if (nrow(single)) {
      single$region <- factor(as.character(single$region),
                              levels=rev(unique(as.character(single$region))))
      p_single <- ggplot(single,aes(y=region)) +
        geom_segment(aes(x=OEstart,xend=OEend,yend=region),
                     color='#7798BE',linewidth=6,lineend='butt') +
        geom_vline(xintercept=snp,linetype='dashed',color='#AD343B',linewidth=.6) +
        geom_text(aes(x=OEend+2600,label=annotation),hjust=0,size=3.6) +
        scale_x_continuous(labels=function(x) format(x,big.mark=',',scientific=FALSE)) +
        coord_cartesian(xlim=c(xleft,xright),clip='off') +
        labs(title='Chicdiff regions containing rs78368856 (FGF10-proximal bait 282953)',
             subtitle='SNP marked by dashed line; overlapping windows are not independent loops',
             x='Genomic position on chr5 (GRCh38)',y=NULL) +
        theme_bw(base_size=12) + theme(panel.grid.minor=element_blank(),
                                       plot.title=element_text(face='bold'))
      save_plot(p_single,'FGF10_SNP_Chicdiff_overlapping_test_regions.png',13,5)
    }
  } else {
    message('No SNP-overlapping Chicdiff regions for the selected baits.')
  }
} else {
  message('Skipping region plots: ', regions_tsv, ' not found (run TSV extractor).')
}
message('Finished. Output directory: ',output)
