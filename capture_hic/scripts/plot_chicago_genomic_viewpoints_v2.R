#!/usr/bin/env Rscript
# Standalone plotting of existing Arima/R CHiCAGO analysis objects.
# DOES NOT rerun HiCUP, bam2chicago, or CHiCAGO.
# Requires: R + Bioconductor Chicago (to load chicagoData .Rds objects).

suppressPackageStartupMessages(library(Chicago))

# Configure paths at runtime. Required files are NOT part of this repository.
args <- commandArgs(trailingOnly = TRUE)
if (any(args %in% c("-h", "--help"))) {
  cat("Usage: Rscript plot_chicago_genomic_viewpoints_v2.R ",
      "--capture-root <dir> --reference-dir <dir> --outdir <dir> ",
      "[--bait-id 282953] [--position 44400000] [--chrom chr5] ",
      "[--window 2000000] [--samples IMR1,IMR2,Lymph1,Lymph2]\n",
      "Default: select the nearest bait to chr5:44,400,000.\n",
      sep = "")
  quit(status = 0)
}
allowed <- c("--capture-root", "--reference-dir", "--outdir", "--bait-id", "--position", "--chrom", "--window", "--samples")
unknown <- setdiff(args[startsWith(args, "--")], allowed)
if (length(unknown)) stop("Unknown option(s): ", paste(unknown, collapse=", "))
get_option <- function(flag, default=NULL) {
  idx <- which(args == flag)
  if (!length(idx)) return(default)
  if (length(idx) != 1L || idx == length(args) || startsWith(args[idx + 1L], "--"))
    stop("Provide exactly one value for ", flag)
  args[idx + 1L]
}
CAPTURE_ROOT <- get_option("--capture-root")
REFERENCE_DIR <- get_option("--reference-dir")
FIGURE_ROOT <- get_option("--outdir")
if (any(vapply(list(CAPTURE_ROOT, REFERENCE_DIR, FIGURE_ROOT), is.null, logical(1)))) {
  stop("Required: --capture-root, --reference-dir and --outdir (see --help)")
}
CAPTURE_ROOT <- normalizePath(CAPTURE_ROOT, mustWork=TRUE)
REFERENCE_DIR <- normalizePath(REFERENCE_DIR, mustWork=TRUE)
FIGURE_ROOT <- normalizePath(FIGURE_ROOT, mustWork=FALSE)
BAITMAP <- file.path(REFERENCE_DIR, "hg38_chicago_input_3kb_1.baitmap")
RMAP <- file.path(REFERENCE_DIR, "hg38_chicago_input_3kb_1.rmap")
SAMPLES <- trimws(strsplit(get_option("--samples", "IMR1,IMR2,Lymph1,Lymph2"), ",", fixed=TRUE)[[1]])
if (any(!nzchar(SAMPLES))) stop("Empty sample name in --samples")
VIEW_CHROM <- get_option("--chrom", "chr5")
VIEW_POS <- as.numeric(get_option("--position", "44400000"))
VIEW_WINDOW <- as.numeric(get_option("--window", "2000000"))
if (!is.finite(VIEW_POS) || !is.finite(VIEW_WINDOW) || VIEW_WINDOW <= 0)
  stop("Invalid --position or --window")
SELECT_BAIT_ID <- get_option("--bait-id")
if (!is.null(SELECT_BAIT_ID)) {
  SELECT_BAIT_ID <- as.integer(SELECT_BAIT_ID)
  if (is.na(SELECT_BAIT_ID)) stop("Invalid --bait-id")
}
VIEWPOINTS <- list(list(chr=VIEW_CHROM, pos=VIEW_POS,
                        upstream=VIEW_WINDOW, downstream=VIEW_WINDOW))
MAX_DISTANCE_TO_BAIT <- 100000
REMOVE_BAIT_TO_BAIT <- TRUE
MAX_ARCS <- 250

if (!file.exists(BAITMAP) || !file.exists(RMAP)) {
  stop("Reference file missing. Check BAITMAP and RMAP paths.")
}

baitmap <- read.delim(BAITMAP, header = FALSE, stringsAsFactors = FALSE,
                      sep = "\t", quote = "", comment.char = "")
rmap <- read.delim(RMAP, header = FALSE, stringsAsFactors = FALSE,
                   sep = "\t", quote = "", comment.char = "")
if (ncol(baitmap) < 4 || ncol(rmap) < 4) stop("Unexpected baitmap/rmap structure")
names(baitmap)[1:4] <- c("chr", "start", "end", "baitID")
names(rmap)[1:4] <- c("chr", "start", "end", "otherEndID")
baitmap$baitID <- as.integer(baitmap$baitID)
rmap$otherEndID <- as.integer(rmap$otherEndID)

# Score categories are mutually exclusive.
score_group <- function(s) {
  ifelse(s >= 7, ">=7",
    ifelse(s >= 5, "5-<7",
      ifelse(s >= 3, "3-<5", "<3")))
}
COLORS <- c("<3" = "#9AA3AE", "3-<5" = "#396ACA",
            "5-<7" = "#E39135", ">=7" = "#D83A44")

# Change labels here without changing the actual score cutoffs.
SCORE_LABELS <- c(
  "<3" = "Score <3",
  "3-<5" = "Score 3 to <5 (below cutoff)",
  "5-<7" = "Score 5 to <7 (called)",
  ">=7" = "Score >=7 (higher-scoring)"
)

find_chicago_rds <- function(sample) {
  folder <- file.path(CAPTURE_ROOT, sample, "chicago")
  if (!dir.exists(folder)) stop("No CHiCAGO results directory: ", folder)
  hits <- list.files(folder, pattern = "\\.[Rr][Dd][Ss]$", full.names = TRUE,
                     recursive = TRUE, ignore.case = TRUE)
  exact <- hits[tolower(basename(hits)) == tolower(paste0(sample, ".Rds"))]
  if (length(exact) == 1L) return(exact)
  if (length(hits) == 1L) return(hits)
  stop("Cannot identify a unique CHiCAGO .Rds for ", sample, ". Found: ",
       paste(hits, collapse = "; "), "\nUse 'find ", folder,
       " -iname \"*.rds\"' to inspect existing objects.")
}

nearest_bait <- function(vp) {
  choices <- baitmap[baitmap$chr == vp$chr, , drop = FALSE]
  if (nrow(choices) == 0L) stop("No baits on ", vp$chr)
  centers <- (choices$start + choices$end) / 2
  if (is.null(SELECT_BAIT_ID)) {
    j <- which.min(abs(centers - vp$pos))
  } else {
    matches <- which(choices$baitID == SELECT_BAIT_ID)
    if (length(matches) != 1L) stop("Bait ID not present on chromosome ", vp$chr, ": ", SELECT_BAIT_ID)
    j <- matches[[1]]
  }
  dist <- abs(centers[j] - vp$pos)
  if (is.null(SELECT_BAIT_ID) && dist > MAX_DISTANCE_TO_BAIT) {
    stop("Nearest captured bait is ", format(dist, big.mark=","),
         " bp from requested coordinate ", vp$chr, ":", vp$pos)
  }
  bait <- choices[j, , drop=FALSE]
  bait$center <- (bait$start + bait$end) / 2
  if (ncol(baitmap) >= 5L) bait$annotation <- as.character(choices[j, 5])
  message("Selected bait ID ", bait$baitID, " at ", bait$chr, ":",
          bait$start, "-", bait$end, " (", round(dist), " bp from input)")
  bait
}

# Join original CHiCAGO full scores/model data to genomic positions.
get_viewpoint <- function(cd, bait, vp) {
  x <- cd@x
  required <- c("baitID", "otherEndID", "N", "score", "Bmean", "distSign")
  missing <- setdiff(required, names(x))
  if (length(missing)) stop("Missing columns in saved chicagoData: ",
                            paste(missing, collapse=", "))
  # Subset BEFORE converting large data.table objects to data.frame.
  keep <- x$baitID == bait$baitID & !is.na(x$distSign)
  if (REMOVE_BAIT_TO_BAIT && "isBait2bait" %in% names(x)) {
    keep <- keep & !(x$isBait2bait %in% TRUE)
  }
  d <- as.data.frame(x[which(keep), ])
  if (!nrow(d)) stop("No interactions for selected bait ", bait$baitID)
  loc <- match(as.integer(d$otherEndID), rmap$otherEndID)
  d$position <- (rmap$start[loc] + rmap$end[loc]) / 2
  d$otherChr <- rmap$chr[loc]
  left <- bait$center - vp$upstream
  right <- bait$center + vp$downstream
  d <- d[!is.na(d$position) & d$otherChr == bait$chr &
         d$position >= left & d$position <= right &
         is.finite(d$N) & is.finite(d$score), , drop=FALSE]
  if (!nrow(d)) stop("No cis other ends in requested viewpoint window")
  d <- d[order(d$position), , drop=FALSE]
  d$group <- factor(score_group(d$score),
                    levels=c("<3", "3-<5", "5-<7", ">=7"))
  d
}

# Use absolute coordinates in Mb on x-axis (not bait-relative distances).
coord_axis <- function(left, right) {
  ats <- pretty(c(left, right), n = 8)
  ats <- ats[ats >= left & ats <= right]
  axis(1, at = ats, labels = sprintf("%.2f", ats / 1e6), cex.axis=0.88)
}

scatter_plot <- function(d, bait, vp, dispersion, sample) {
  left <- bait$center - vp$upstream
  right <- bait$center + vp$downstream
  model_ok <- any(is.finite(d$Bmean) & d$Bmean >= 0)
  # CHiCAGO stores the dispersion in cd@params, not cd@settings.
  # The dashed curve is the optional normal approximation used in
  # HiCExplorer's example, not a CHiCAGO significance threshold.
  upper_ok <- model_ok && length(dispersion) == 1L &&
              is.numeric(dispersion) && is.finite(dispersion) && dispersion > 0
  upper <- if (upper_ok) {
    d$Bmean + 1.96 * sqrt(pmax(0, d$Bmean + d$Bmean^2 / dispersion))
  } else numeric(0)
  ymax <- max(c(d$N, if (model_ok) d$Bmean, upper), na.rm=TRUE) * 1.10
  plot(NA_real_, NA_real_, xlim=c(left,right), ylim=c(0,max(1,ymax)),
       xaxt="n", xlab=paste0(bait$chr," position (Mb, hg38)"),
       ylab="Observed interaction counts (N)",
       main=paste0(sample, " | CHiCAGO viewpoint | bait ", bait$baitID))
  coord_axis(left, right)
  abline(v=bait$center, col="#14889B", lwd=1.4, lty=3)
  if (model_ok) lines(d$position, d$Bmean, col="#424D58", lwd=1.2)
  if (upper_ok) lines(d$position, upper, col="#424D58", lwd=1, lty=2)
  for (grp in levels(d$group)) {
    k <- which(d$group == grp)
    if (length(k)) points(d$position[k], d$N[k], pch=16,
                          col=adjustcolor(COLORS[[grp]], alpha.f=if(grp=="<3") .55 else .9),
                          cex=if(grp=="<3") .46 else .77)
  }
  groups <- levels(d$group)
  legend("topright", bty="n", cex=.79,
         legend=c(unname(SCORE_LABELS[groups]),
                  if (model_ok) "Fitted Brownian mean (Bmean)",
                  if (upper_ok) "Approx. 95% Brownian upper",
                  "Captured bait"),
         col=c(unname(COLORS[groups]),
               if (model_ok) "#424D58",
               if (upper_ok) "#424D58", "#14889B"),
         pch=c(rep(16, length(groups)), if (model_ok) NA,
               if (upper_ok) NA, NA),
         lty=c(rep(NA, length(groups)), if (model_ok) 1,
               if (upper_ok) 2, 3),
         lwd=c(rep(NA, length(groups)), if (model_ok) 1.5,
               if (upper_ok) 1, 1.4))
  mtext("N = observed raw or replicate-merged counts; Bmean = fitted Brownian mean; dashed line (if shown) = approximate upper reference.",
        side=1, line=3.5, cex=.65)
}

arcs_plot <- function(d, bait, vp, sample) {
  left <- bait$center - vp$upstream
  right <- bait$center + vp$downstream
  hits <- d[d$score >= 3, , drop=FALSE]
  if (nrow(hits) > MAX_ARCS) {
    hits <- hits[order(hits$score, decreasing=TRUE)[seq_len(MAX_ARCS)],,drop=FALSE]
  }
  hits <- hits[order(hits$score),,drop=FALSE]  # strongest arcs on top
  heights <- log1p(hits$N)
  ymax <- max(1, heights, na.rm=TRUE) * 1.15
  plot(NA_real_, NA_real_, xlim=c(left,right), ylim=c(0,ymax),
       xaxt="n", yaxt="n", xlab=paste0(bait$chr," position (Mb, hg38)"),
       ylab="Arc height: log(1 + N)",
       main=paste0(sample, " | CHiCAGO interaction arcs (score >=3) | bait ", bait$baitID))
  coord_axis(left,right)
  axis(2, las=1)
  abline(h=0, col="gray65")
  abline(v=bait$center, col="#14889B", lwd=1.4, lty=3)
  t <- seq(0,pi,length.out=80)
  for (i in seq_len(nrow(hits))) {
    a <- min(bait$center, hits$position[i])
    b <- max(bait$center, hits$position[i])
    px <- (a+b)/2 - (b-a)/2*cos(t)
    py <- heights[i]*sin(t)
    lines(px, py, col=adjustcolor(COLORS[[as.character(hits$group[i])]],alpha.f=.7),
          lwd=if (hits$score[i]>=7) 1.8 else 1.25)
  }
  legend("topright",bty="n",cex=.85,
         legend=c(unname(SCORE_LABELS[c("3-<5", "5-<7", ">=7")]), "Captured bait"),
         col=c(COLORS[["3-<5"]],COLORS[["5-<7"]],COLORS[[">=7"]],"#14889B"),
         lty=c(1,1,1,3), lwd=c(1.5,1.5,2,1.5))
  if (!nrow(hits)) mtext("No interactions with CHiCAGO score >=3 in this window", side=3, line=-2)
  mtext(paste0("Each arc joins the actual bait midpoint and other-end coordinate. ",
               nrow(hits), " interactions drawn."), side=1, line=3.5, cex=.7)
}

save_both <- function(plot_fun, prefix) {
  pdf(paste0(prefix,".pdf"),width=12,height=6.5,useDingbats=FALSE)
  tryCatch(plot_fun(), finally=dev.off())
  png(paste0(prefix,".png"),width=12,height=6.5,units="in",res=220)
  tryCatch(plot_fun(), finally=dev.off())
}

for (sample in SAMPLES) {
  message("\n### ", sample)
  tryCatch({
    in_file <- find_chicago_rds(sample)
    message("Reading CHiCAGO results: ", basename(in_file))
    cd <- readRDS(in_file)
    if (!is(cd, "chicagoData")) stop("RDS does not contain a chicagoData object")
    dispersion <- cd@params$dispersion
    if (length(dispersion) == 1 && is.finite(dispersion) && dispersion > 0) {
      message("Fitted Brownian dispersion: ", dispersion)
    } else {
      message("Dispersion not available: dashed upper line will be omitted")
    }
    dest <- file.path(FIGURE_ROOT, sample)
    dir.create(dest, recursive=TRUE, showWarnings=FALSE)
    for (vp in VIEWPOINTS) {
      tryCatch({
        bait <- nearest_bait(vp)
        d <- get_viewpoint(cd, bait, vp)
        prefix <- file.path(dest, paste0(sample,"_",bait$chr,"_bait",bait$baitID))
        save_both(function() scatter_plot(d, bait, vp, dispersion, sample), paste0(prefix,"_scatter_genomic"))
        save_both(function() arcs_plot(d, bait, vp, sample), paste0(prefix,"_arcs_genomic"))
        message("Saved viewpoint plots for ", sample)
      }, error=function(e) message("Viewpoint skipped: ", conditionMessage(e)))
    }
    rm(cd)
    invisible(gc())
  }, error=function(e) message("Sample skipped: ", conditionMessage(e)))
}
