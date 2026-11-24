#!/usr/bin/env bash
# Portable wrapper for the original three Arima-CHiC v1.5 stages.
# Local sequencing data, design references, and installation are not distributed.
set -euo pipefail

usage() {
  cat <<'EOF'
Usage: run_all_samples.sh --stage hicup|bam2chicago|chicago|all [--execute]

Configuration is supplied by environment variables (see ../config.example.env):
  ARIMA_SCRIPT, BOWTIE2_BIN, BOWTIE2_INDEX, DIGEST_FILE,
  HICUP_DIR, CHICAGO_TOOLS_DIR, PROBE_BED, BAITMAP, RMAP,
  DESIGN_DIR, SAMPLES_DIR, OUTPUT_DIR.
Optional: SAMPLES (comma separated; defaults to IMR1,IMR2,Lymph1,Lymph2),
          THREADS (default 24), GENOME (default hg38), RESOLUTION (default 3kb).

Default is a dry run that prints stage/sample summaries, not a full command.
Provide --execute to process data. This is NOT an automatic resume script.
Do not re-run completed stages without confirming output and input state.
EOF
}

stage=''
execute=0
while (($#)); do
  case "$1" in
    --stage) [[ $# -ge 2 ]] || { usage >&2; exit 2; }; stage=$2; shift 2 ;;
    --execute) execute=1; shift ;;
    -h|--help) usage; exit 0 ;;
    *) echo "Unknown argument: $1" >&2; usage >&2; exit 2 ;;
  esac
done
case "$stage" in
  hicup|bam2chicago|chicago|all) ;;
  *) echo 'Missing or invalid --stage' >&2; usage >&2; exit 2 ;;
esac

required=(ARIMA_SCRIPT BOWTIE2_BIN BOWTIE2_INDEX DIGEST_FILE
          HICUP_DIR CHICAGO_TOOLS_DIR PROBE_BED BAITMAP RMAP
          DESIGN_DIR SAMPLES_DIR OUTPUT_DIR)
for key in "${required[@]}"; do
  if [[ -z "${!key:-}" ]]; then
    echo "Missing environment setting: $key" >&2
    exit 2
  fi
done

: "${GENOME:=hg38}"
: "${RESOLUTION:=3kb}"
: "${THREADS:=24}"
: "${SAMPLES:=IMR1,IMR2,Lymph1,Lymph2}"
IFS=',' read -r -a sample_names <<< "$SAMPLES"

# These options match the original Arima-CHiC v1.5 wrapper configuration.
run_stage() {
  local sample=$1 stage_name=$2 fastq1=$3 fastq2=$4 sample_out=$5
  local -a mode=()
  case "$stage_name" in
    hicup)      mode=(-W 1 -Y 0 -Z 0 -P 0) ;;
    bam2chicago) mode=(-W 0 -Y 1 -Z 0 -P 0) ;;
    chicago)    mode=(-W 0 -Y 0 -Z 1 -P 1) ;;
  esac
  echo "Sample: $sample | Stage: $stage_name"
  if ((execute == 0)); then
    echo '  Dry run; no pipeline commands executed.'
    return
  fi
  bash "$ARIMA_SCRIPT" "${mode[@]}" \
    -A "$BOWTIE2_BIN" -X "$BOWTIE2_INDEX" -d "$DIGEST_FILE" \
    -H "$HICUP_DIR" -C "$CHICAGO_TOOLS_DIR" \
    -I "${fastq1},${fastq2}" -o "$sample_out" -p "$sample" \
    -b "$PROBE_BED" -R "$RMAP" -B "$BAITMAP" \
    -D "$DESIGN_DIR" -O "$GENOME" -r "$RESOLUTION" -t "$THREADS"
}

if ((execute)); then
  for file in "$ARIMA_SCRIPT" "$BOWTIE2_BIN" "$DIGEST_FILE" "$PROBE_BED" "$BAITMAP" "$RMAP"; do
    [[ -f "$file" ]] || { echo "Missing input: $file" >&2; exit 2; }
  done
  for dir in "$HICUP_DIR" "$CHICAGO_TOOLS_DIR" "$DESIGN_DIR" "$SAMPLES_DIR"; do
    [[ -d "$dir" ]] || { echo "Missing directory: $dir" >&2; exit 2; }
  done
fi

for sample in "${sample_names[@]}"; do
  [[ "$sample" =~ ^[A-Za-z0-9_-]+$ ]] || { echo "Invalid sample name: $sample" >&2; exit 2; }
  fastq1="${SAMPLES_DIR}/${sample}/${sample}_S1_R1_001.fastq.gz"
  fastq2="${SAMPLES_DIR}/${sample}/${sample}_S1_R2_001.fastq.gz"
  sample_out="${OUTPUT_DIR}/${sample}"
  if ((execute)); then
    [[ -s "$fastq1" && -s "$fastq2" ]] || {
      echo "Missing FASTQ(s) for sample ${sample}; stopping." >&2
      exit 1
    }
    mkdir -p "$sample_out"
  fi
  if [[ "$stage" == all ]]; then
    for phase in hicup bam2chicago chicago; do
      run_stage "$sample" "$phase" "$fastq1" "$fastq2" "$sample_out"
    done
  else
    run_stage "$sample" "$stage" "$fastq1" "$fastq2" "$sample_out"
  fi
done

if ((execute)); then echo 'Finished requested stages.'; fi
