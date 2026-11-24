# Genome-wide Hi-C processing and FGF10 locus visualization

This directory documents the genome-wide Hi-C processing and regional visualization used to examine chromatin organization around **FGF10** and the candidate modifier **rs78368856** (GRCh38/hg38, chr5:43,619,485). It is intended as a **code-and-methods record**, not a complete sequencing-data deposit. Original experimental datasets and most processed outputs are not distributed here.

## Samples and experimental context

Four independently cultured samples were analyzed:

| Sample identifier | Cell type |
| --- | --- |
| IMR1 | IMR-90 human lung fibroblasts, replicate 1 |
| IMR2 | IMR-90 human lung fibroblasts, replicate 2 |
| Lymph1 | Lymphoblasts, replicate 1 |
| Lymph2 | Lymphoblasts, replicate 2 |

The genome-wide Hi-C results were visualized together with separately processed **Capture Hi-C** interactions from the same samples. The genome-wide and Capture Hi-C experiments have distinct data-processing workflows.

## Genome-wide Hi-C analysis

Sequencing reads were analyzed against **GRCh38/hg38** using **nf-core/hic v2.1.0** executed through **Nextflow** with the Docker profile. The pipeline was configured for the Arima restriction design, using an hg38 restriction-fragment BED file containing **GATC/GANTC** sites (`hg38_GATC_GANTC.bed`). The pipeline used a custom Nextflow configuration file and a samplesheet referring to trimmed FASTQ inputs.

The following command is transcribed from the `.nextflow.log` launch line:

```bash
nextflow run nf-core/hic -r 2.1.0 -profile docker \
  -c <hic_long_align.config> \
  --input <samplesheet_trimmed.csv> \
  --outdir <hic_results> \
  --fasta <hg38.fa> \
  --restriction_fragments <hg38_GATC_GANTC.bed> \
  --digestion arima \
  --min_restriction_fragment_size 10 \
  --max_restriction_fragment_size 100000 \
  --min_insert_size 100 \
  --max_insert_size 1000 \
  --split_fastq \
  --fastq_chunks_size 20000000 \
  --bin_size 5000,10000,20000,40000,50000,100000,250000,500000,1000000 \
  --res_tads 20000,40000 \
  --save_raw_maps \
  --save_reference \
  --max_cpus 20 \
  -resume
```

The original run also used local Docker/Nextflow software and reference resources, which are not bundled with this repository.

The downstream regional figures used **balanced `.cool` matrices** at **20-kb** and **50-kb** resolutions. Chromatin insulation was represented using the **300-kb-window `log2_insulation_score_300000`** column from the sample-specific **20-kb** insulation TSVs. The figure script extracts the chr5 values and writes a temporary bedGraph track; it does not calculate insulation scores from raw Hi-C pairs.

## Figure generation: `plot_chr5_hic_all_replicates.py`

The original manuscript figure script was **`plot_fgf10_hic_all_replicates_v7.py`**. It generated integrated chromosome 5 regional views for the four samples:

| View | GRCh38 region | Genome-wide Hi-C bin size | Plotted Hi-C depth |
| --- | --- | --- | --- |
| Close | chr5:42,796,905–44,441,081 | 20 kb | 1 Mb |
| Wide | chr5:42,300,000–46,000,000 | 50 kb | 2 Mb |

The script assembled these layers using the pyGenomeTracks Python API:

1. **Genome-wide Hi-C:** sample-specific balanced Cooler matrices; triangular contact map displays used `log1p` transformation and the `Reds` colormap.
2. **Capture Hi-C:** sample-specific `*.cis_le_2Mb.bedpe` files from the separate Capture Hi-C analysis, limited to cis contacts with both endpoints in the plotted view. The script **does not filter to a single bait or to FGF10-specific baits**. Arcs are divided into **5 to <7 (orange)**, and **≥7 (red)**. A score of **≥5** is the conventional CHiCAGO high-confidence threshold.
3. **Insulation:** a regional bedGraph derived from the sample-specific 20-kb insulation TSV (`log2_insulation_score_300000`).
4. **Gene annotations:** simplified gene-span BED tracks derived from GENCODE human release 50 (GRCh38).
5. **Additional annotations:** public CTCF tracks for lung fibroblasts (ENCFF777ODE) and IMR-90 cells (ENCFF307XFM), and publicly accessible ATAC-seq signal tracks for IMR-90, HUVECs and placenta. External tracks are included only when accessible.
6. **Variant marker:** a dotted line at rs78368856 was displayed in the IMR1/IMR2 panels by the v7 plotting code. Its selective display is a plotting choice.

### Example figure-generation command

The original script accepts configurable input and output directories:

```bash
python plot_fgf10_hic_all_replicates.py \
  --project-root <genome-wide_Hi-C_project> \
  --capture-root <capture_Hi-C_output> \
  --gtf <gencode.v50.basic.annotation.gtf.gz> \
  --outdir <figures/chr5_all_samples> \
  --png
```

**Required local inputs:** balanced sample-specific 20-kb and 50-kb Cooler files; sample-specific 20-kb insulation TSVs containing the 300-kb window column; sample-specific Capture Hi-C `cis_le_2Mb.bedpe` files; GENCODE v50 annotation; and accessible external annotation tracks if those tracks are to appear in the figure. The script is provided as analysis provenance and may require path or dependency adjustments in a new environment.

## Software and references

- nf-core/hic **v2.1.0**: https://nf-co.re/hic/2.1.0
- Nextflow: https://www.nextflow.io/
- HiC-Pro and Bowtie2: used within the nf-core/hic processing workflow
- Cooler and pyGenomeTracks: used for regional contact-matrix visualization and track assembly
- CHiCAGO: used to score interactions in the **separate Capture Hi-C** workflow
- GENCODE human release **50**, GRCh38: https://www.gencodegenes.org/
