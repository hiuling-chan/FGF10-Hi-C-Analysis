# FGF10 Capture Hi-C differential-interaction analysis

This repository folder contains scripts to reproduce targeted Chicdiff summaries and plots for **rs78368856** at **chr5:43,619,485 (GRCh38/hg38)**, near *FGF10*. The analyses compare **two independently cultured IMR-90 fibroblast samples (IMR1/IMR2)** with **two lymphoblast samples (Lymph1/Lymph2)**.

**Assay:** Arima Capture Hi-C with a nominal **~3-kb capture design**, processed using Arima-CHiC v1.5 and CHiCAGO. This is **not** the binning resolution of separate genome-wide Hi-C heatmaps (20 or 50 kb).

**Interaction set:** cis contacts **≤2 Mb**, with union peaks defined by **CHiCAGO score ≥5 in at least one sample**, via `makePeakMatrix.R --twopass --notrans --maxdist 2000000 --cutoff 5`. Chicdiff v0.6 then performs differential testing using the original four `.chinput` files and replicate-specific `.Rds` objects. Individual lower-scoring contacts can contribute to larger regional differential tests. Chicdiff tests **extended genomic regions**, so overlapping tests are **not independent loops**.

## What the fragment IDs mean

| Arima fragment ID | GRCh38 coordinates | Role | Interpretation |
|---|---|---|---|
| **282709** | chr5:43,616,691–43,619,914 | **Other-end restriction fragment** | Contains rs78368856. **Not a captured bait** in these FGF10 comparisons. |
| **282949** | chr5:44,384,622–44,388,806 | **FGF10 promoter-overlapping bait** | Overlaps the left edge of the assayed FGF10 promoter interval (chr5:44,388,434–44,390,361). |
| **282950** | chr5:44,388,807–44,392,024 | **FGF10 promoter-overlapping bait** | Overlaps the assayed FGF10 promoter interval. |
| **282953** | chr5:44,398,564–44,401,751 | **FGF10 promoter-proximal bait** | Near the FGF10 promoter but **outside** the specified assayed promoter interval. |

Coordinates are from the hg38 Arima `.rmap` and `.baitmap` files and the assayed promoter coordinates. The word “bait” refers to the captured restriction fragment; the “other end” is its interacting fragment. Here rs78368856 is on the other-end fragment **282709**.

## Scripts

1. `scripts/run_FGF10_chicdiff.R` — optional reproducible **full differential test**, with safety check that refuses to overwrite completed results unless explicitly passed `--force`. Requires the **existing** peak matrix, CHiCAGO RDS, CHiCAGO `.chinput`, and Arima design references. **You do not need to rerun this to extract or plot your current results.**
2. `scripts/extract_FGF10_bait_TSVs.py` — from the completed full Chicdiff results, exports all regions involving bait 282949, 282950, and 282953; SNP-overlapping subsets; annotations; and exact per-replicate CHiCAGO fragment-pair counts/scores (if that original table is available). Python standard library only.
3. `scripts/plot_FGF10_chicdiff_clean.R` — plots dedicated bait profiles using existing Chicdiff RDS objects (when available); sample-by-sample exact fragment counts/scores; and SNP-overlapping regional differential tests from exported TSVs. **Does not rerun statistics.**

## Reproduce on the laboratory server

Set the Chicdiff environment once:

```bash
export MAMBA_ROOT_PREFIX="$HOME/micromamba"
MM="$HOME/micromamba-bootstrap/bin/micromamba"
```

The following assumes these scripts are in the working directory. The code contains the original lab path defaults and can be edited for portability.

**A. Only if rerunning differential testing is necessary** (skip if `IMR_vs_Lymph_results.Rds` already exists):

```bash
"$MM" run -n chicdiff Rscript scripts/run_FGF10_chicdiff.R
```

If the results already exist, the script intentionally stops instead of overwriting them. `--force` explicitly permits a deliberate reanalysis. The peak matrix is read from `/nfs/volume2/Hi-C_3kb/chicdiff_results/IMR_Lymph.txt` (it must first be generated using the original ~3-kb Arima capture design and cis ≤2 Mb criteria). The script writes results under `.../chicdiff_results/differential/`.

**B. Extract FGF10 tables from completed results (no reanalysis):**

```bash
python3 scripts/extract_FGF10_bait_TSVs.py
```

Default source is `.../differential/IMR_vs_Lymph_full_results.tsv`. The script optionally reads `.../annotations/fgf10_anchor_audit/06_all_candidate_fragment_chicago_scores_counts.tsv`. Outputs are written to `.../differential/FGF10_tables/`.

**C. Generate figures from completed outputs (no reanalysis):**

```bash
"$MM" run -n chicdiff Rscript scripts/plot_FGF10_chicdiff_clean.R
```

Outputs are in `.../differential/FGF10_clean_plots/`. Individual sample counts and scores come from the exact fragment-level CHiCAGO TSV, while statistical bars/intervals come from **Chicdiff’s broader region-level results**. For Native `plotDiffBaits()` profiles, the script also reads the local results and count RDS objects and the original baitmap. The mirrored `IMR`/`Lymph` profiles aggregate the two biological replicates per group for display; they do not show the four replicates individually.

## Exported TSVs

| File | What it contains |
|---|---|
| `FGF10_all_3_baits_all.tsv` | Every tested region for baits 282949, 282950, 282953 |
| `FGF10_promoter_baits_all.tsv` | Every tested region for promoter-overlapping baits 282949 and 282950 |
| `FGF10_promoter_baits_SNP_overlap.tsv` | Tested promoter-bait regions whose other-end genomic intervals contain the SNP |
| `FGF10_proximal_bait282953_all.tsv` | Every tested region for promoter-proximal bait 282953 |
| `FGF10_proximal_bait282953_SNP_overlap.tsv` | Tested regions for 282953 containing the SNP |
| `FGF10_all_3_baits_SNP_overlap.tsv` | Combined SNP-overlapping regions for all three baits |
| `FGF10_exact_SNP_fragment_pairs_counts_scores.tsv` | **Single ~3-kb other-end fragment pair** scores and raw counts (four replicates × three baits); only created if the original CHiCAGO exact-pair table is available |
| `FGF10_bait_and_SNP_fragment_annotations.tsv` | Coordinates and explanations of 282709, 282949, 282950, 282953 |

The statistical result columns include `log2FoldChange`, `pvalue`, `padj`, and `weighted_padj`. Before attributing the sign of `log2FoldChange` to IMR or Lymph, **verify the coefficient/contrast orientation from the saved Chicdiff model**. A small adjusted p-value supports a **between-cell-type difference in a genomic interaction region**. It **does not** demonstrate that rs78368856 itself alters chromatin looping or *FGF10* expression.

## Recommended GitHub contents

**Include:** these scripts, this README, the FGF10 tables, the relevant figures, `sessionInfo.txt` (after checking for system/private paths) and optionally `IMR_vs_Lymph_full_results.tsv` (~1.2 MB) with laboratory approval. If publishing the peak list or peak matrix, document the exact selection criteria and ensure you have permission for the underlying data.

**Do not upload by default:** `chicdiff_run.log` (can contain local paths and errors), raw FASTQ/BAM files, `.chinput` files, CHiCAGO `.Rds` objects, `IMR_vs_Lymph_countput.Rds`, raw patient genomes, or Arima reference maps without redistribution permission. Keep the **full log privately** for troubleshooting; a short README plus version information usually serves public reproducibility better. A sanitized excerpt can be shared if a specific installation or analysis issue must be documented.

## References

- Cairns J et al. Chicdiff: a computational pipeline for detecting differential chromosomal interactions in Capture Hi-C data. *Bioinformatics* (2019). DOI: [10.1093/bioinformatics/btz450](https://doi.org/10.1093/bioinformatics/btz450).
- Chicdiff source and documentation: https://github.com/RegulatoryGenomicsGroup/chicdiff
- CHiCAGO/Chicdiff workflow protocol: https://pmc.ncbi.nlm.nih.gov/articles/PMC7612634/
