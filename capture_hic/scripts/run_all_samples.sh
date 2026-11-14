#!/bin/bash

export TMPDIR=/tmp

# Config 
BOWTIE2_PATH="/usr/bin/bowtie2"
BOWTIE2_INDEX="/data/reference/hg38"
DIGEST="/data/reference/Digest_hg38_DpnII_Arima_None_20-31-36_15-07-2025.txt"
HICUP_DIR="/opt/hicup"
CHICAGO_DIR="/opt/arima/chicagoTools"
BED_PROBES="/data/reference/3536031_Covered_good.bed"
BAITMAP="/data/reference/hg38_chicago_input_3kb_1.baitmap"
RMAP="/data/reference/hg38_chicago_input_3kb_1.rmap"
DESIGN_DIR="/data/reference"
GENOME="hg38"
RESOLUTION="3kb"
THREADS=24

SAMPLES_DIR="/data/samples"
OUTPUT_DIR="/data/output"

# List of sample names
SAMPLES=("IMR1" "IMR2" "Lymph1" "Lymph2")

echo "Starting Arima-HiC pipeline on all samples..."
echo "Samples to process: ${SAMPLES[*]}"
echo "=============================="

for sample in "${SAMPLES[@]}"; do
  echo "Processing sample: $sample"

  sample_in="${SAMPLES_DIR}/${sample}"
  fastq1="${sample_in}/${sample}_S1_R1_001.fastq.gz"
  fastq2="${sample_in}/${sample}_S1_R2_001.fastq.gz"
  sample_out="${OUTPUT_DIR}/${sample}"
  hicup_out="${sample_out}/hicup"
  bam2chicago_out="${sample_out}/bam2chicago/${sample}"
  chicago_out="${sample_out}/chicago"

  hicup_bam="${hicup_out}/${sample}_S1_R1_2_001.hicup.bam"

  chinput_file="${bam2chicago_out}/${sample}.chinput"

  chicago_washu_text="${chicago_out}/data/${sample}_washU_text.txt"
  chicago_washu_track="${chicago_out}/data/${sample}_washU_track.txt"
  chicago_done_flag="${chicago_out}/CHICAGO_DONE.flag"

  if [[ ! -f "$fastq1" || ! -f "$fastq2" ]]; then
    echo "ERROR: FASTQ files for $sample missing! Skipping sample."
    continue
  fi

  mkdir -p "$hicup_out" "$bam2chicago_out" "$chicago_out"

  # Step 1: HiCUP
  if [[ ! -f "$hicup_bam" ]]; then
    echo "Running HiCUP..."
    /opt/arima/Arima-CHiC-v1.5.sh -W 1 -Y 0 -Z 0 -P 0 \
      -A "$BOWTIE2_PATH" \
      -X "$BOWTIE2_INDEX" \
      -d "$DIGEST" \
      -H "$HICUP_DIR" \
      -C "$CHICAGO_DIR" \
      -I "${fastq1},${fastq2}" \
      -o "$sample_out" \
      -p "$sample" \
      -b "$BED_PROBES" \
      -R "$RMAP" \
      -B "$BAITMAP" \
      -D "$DESIGN_DIR" \
      -O "$GENOME" \
      -r "$RESOLUTION" \
      -t "$THREADS"
    if [[ $? -ne 0 ]]; then
      echo "HiCUP step failed for $sample. Exiting."
      exit 1
    fi

    echo "Cleaning large HiCUP intermediate files"
    rm -f "$hicup_out"/*.sam "$hicup_out"/*.pair.bam "$hicup_out"/*.filt.bam
  else
    echo "HiCUP already done for $sample, skipping."
  fi

  # Step 2: bam2chicago (.chinput generator)
  if [[ ! -f "$chinput_file" ]]; then
    echo "Running bam2chicago to generate .chinput..."
    /opt/arima/Arima-CHiC-v1.5.sh -W 0 -Y 1 -Z 0 -P 0 \
      -A "$BOWTIE2_PATH" \
      -X "$BOWTIE2_INDEX" \
      -d "$DIGEST" \
      -H "$HICUP_DIR" \
      -C "$CHICAGO_DIR" \
      -I "${fastq1},${fastq2}" \
      -o "$sample_out" \
      -p "$sample" \
      -b "$BED_PROBES" \
      -R "$RMAP" \
      -B "$BAITMAP" \
      -D "$DESIGN_DIR" \
      -O "$GENOME" \
      -r "$RESOLUTION" \
      -t "$THREADS"
    if [[ $? -ne 0 ]]; then
      echo "bam2chicago failed for $sample. Exiting."
      exit 1
    fi
  else
    echo ".chinput already exists for $sample, skipping bam2chicago."
  fi

  # Step 3: CHiCAGO
  if [[ -f "$chinput_file" && ! -f "$chicago_washu_text" && ! -f "$chicago_washu_track" ]]; then
    echo "Running CHiCAGO..."
    /opt/arima/Arima-CHiC-v1.5.sh -W 0 -Y 0 -Z 1 -P 1 \
      -A "$BOWTIE2_PATH" \
      -X "$BOWTIE2_INDEX" \
      -d "$DIGEST" \
      -H "$HICUP_DIR" \
      -C "$CHICAGO_DIR" \
      -I "${fastq1},${fastq2}" \
      -o "$sample_out" \
      -p "$sample" \
      -b "$BED_PROBES" \
      -R "$RMAP" \
      -B "$BAITMAP" \
      -D "$DESIGN_DIR" \
      -O "$GENOME" \
      -r "$RESOLUTION" \
      -t "$THREADS"
    if [[ $? -ne 0 ]]; then
      echo "CHiCAGO failed for $sample. Exiting."
      exit 1
    fi

    if [[ -f "$chicago_washu_text" || -f "$chicago_washu_track" ]]; then
      touch "$chicago_done_flag"
    else
      echo "WARNING: CHiCAGO command finished but expected washU output files were not found for $sample."
      exit 1
    fi

  elif [[ ! -f "$chinput_file" ]]; then
    echo "No .chinput file for $sample; skipping CHiCAGO."
  else
    echo "CHiCAGO already done for $sample, skipping."
  fi

  echo "Finished processing $sample"
  echo "------------------------------"
done

echo "All samples complete."
