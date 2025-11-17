#!/usr/bin/env bash
set -euo pipefail


# Check for mandatory arguments
if [ "$#" -ne 4 ]; then
    echo "Usage: $0 <BASE_DIRECTORY> <GENOME_PATH> <START_TYPE (pod5|bam|archive)> <INPUT_PATH>" >&2
    exit 1
fi


MAIN_DIR="$1"
GENOME="$2"
START_TYPE="$3"
INPUT_PATH="$4"
RUN_PREFIX="DS98"


# =========================
# Configuration (All under data/)
# =========================
DATA_DIR="${MAIN_DIR}/data"
BASECALLED_DIR="${DATA_DIR}/basecalled"
DEMUX_DIR="${DATA_DIR}/demultiplexed"
FASTQ_DIR="${DATA_DIR}/fastq"
PER_SAMPLE_DIR="${DATA_DIR}/per_sample"
POD5_TMP_DIR="${MAIN_DIR}/pod5_tmp"


BASECALLED_BAM="${BASECALLED_DIR}/${RUN_PREFIX}_TDAC.bam"
DORADO_KIT="SQK-NBD114-24"
SAMPLE_CSV="${MAIN_DIR}/samples.csv"
RENAME_TABLE="${DATA_DIR}/rename_table.tsv"


# =========================
# Sanity checks and setup
# =========================
echo "[TDAC.sh] Using project root: ${MAIN_DIR}"
if [[ ! -s "${SAMPLE_CSV}" ]]; then
  echo "[ERROR] Missing samples CSV at ${SAMPLE_CSV}. Expected header: 'Sample,Barcode'." >&2
  exit 1
fi


mkdir -p "${BASECALLED_DIR}" "${DEMUX_DIR}" "${FASTQ_DIR}" "${PER_SAMPLE_DIR}" "${POD5_TMP_DIR}"


# Helper: NBxx -> xx
nb_to_num() { echo "${1#NB}"; }


# =========================
# 1) Basecalling (Conditional)
# =========================
POD5_INPUT_DIR=""


if [[ "${START_TYPE}" == "pod5" ]]; then
    echo "[*] Starting from pod5 folder."
    if [[ ! -d "${INPUT_PATH}" ]]; then
        echo "[ERROR] pod5 path not found: ${INPUT_PATH}" >&2
        exit 1
    fi
    POD5_INPUT_DIR="${INPUT_PATH}"
    echo "[*] Basecalling with dorado from: ${POD5_INPUT_DIR}"
    mamba run -n nano dorado basecaller hac --kit-name "${DORADO_KIT}" "${POD5_INPUT_DIR}" --no-trim > "${BASECALLED_BAM}"


elif [[ "${START_TYPE}" == "archive" ]]; then
    echo "[*] Starting from pod5 archive."
    if [[ ! -s "${INPUT_PATH}" ]]; then
        echo "[ERROR] pod5 archive not found or empty: ${INPUT_PATH}" >&2
        exit 1
    fi
    echo "[*] Cleaning temporary extraction directory: ${POD5_TMP_DIR}"
    rm -rf "${POD5_TMP_DIR:?}"/*
    
    echo "[*] Extracting ${INPUT_PATH} (multithreaded with pigz)..."
    if mamba run -n tdac pigz --version &> /dev/null; then
        echo "[*] Using pigz for multithreaded extraction..."
        mamba run -n tdac pigz -dc "${INPUT_PATH}" | tar -xf - -C "${POD5_TMP_DIR}" --strip-components=1 2>&1 | grep -v "Removing leading" || {
            echo "[!] Extraction with strip-components failed, trying without stripping..."
            mamba run -n tdac pigz -dc "${INPUT_PATH}" | tar -xf - -C "${POD5_TMP_DIR}"
        }
    else
        echo "[!] pigz not found in tdac environment, using standard tar..."
        echo "[TIP] Install pigz: mamba install -n tdac pigz -y"
        tar -xzf "${INPUT_PATH}" -C "${POD5_TMP_DIR}"
    fi
    
    POD5_SEARCH=$(find "${POD5_TMP_DIR}" -type f -name "*.pod5" -print -quit)
    
    if [[ -z "${POD5_SEARCH}" ]]; then
        echo "[ERROR] No pod5 files found in extracted archive!" >&2
        ls -lR "${POD5_TMP_DIR}" >&2
        exit 1
    fi
    
    POD5_INPUT_DIR=$(dirname "${POD5_SEARCH}")
    POD5_COUNT=$(find "${POD5_INPUT_DIR}" -name "*.pod5" | wc -l)
    echo "[*] Found ${POD5_COUNT} pod5 files in: ${POD5_INPUT_DIR}"
    
    echo "[*] Basecalling with dorado..."
    mamba run -n nano dorado basecaller hac --kit-name "${DORADO_KIT}" "${POD5_INPUT_DIR}" --no-trim > "${BASECALLED_BAM}"
    
    echo "[*] Cleaning up temporary files..."
    rm -rf "${POD5_TMP_DIR:?}"/*


elif [[ "${START_TYPE}" == "bam" ]]; then
    echo "[*] Skipping basecalling, starting from existing BAM."
    
    if [[ ! -s "${INPUT_PATH}" ]]; then
        echo "[ERROR] BAM file not found or empty: ${INPUT_PATH}" >&2
        exit 1
    fi
    
    INPUT_ABS=$(realpath "${INPUT_PATH}" 2>/dev/null)
    if [[ ! -f "${INPUT_ABS}" ]]; then
        echo "[ERROR] Failed to resolve input path: ${INPUT_PATH}" >&2
        exit 1
    fi
    
    mkdir -p "${BASECALLED_DIR}"
    CANONICAL_DIR=$(realpath "${BASECALLED_DIR}")
    CANONICAL_PATH="${CANONICAL_DIR}/$(basename "${BASECALLED_BAM}")"
    
    if [[ "${INPUT_ABS}" == "${CANONICAL_PATH}" ]]; then
        echo "[*] ✓ Input BAM already at canonical location"
    else
        echo "[*] Creating symlink: $(basename ${INPUT_PATH}) -> ${BASECALLED_BAM}"
        rm -f "${BASECALLED_BAM}"
        ln -sf "${INPUT_ABS}" "${BASECALLED_BAM}"
        
        if [[ $? -eq 0 ]]; then
            echo "[*] ✓ Symlink created successfully"
        else
            echo "[ERROR] Failed to create symlink" >&2
            exit 1
        fi
    fi
    
    echo "[*] Using: ${BASECALLED_BAM}"
else
    echo "[ERROR] Invalid start type: '${START_TYPE}'. Must be 'pod5', 'bam', or 'archive'." >&2
    exit 1
fi


# =========================
# 2) Demultiplexing
# =========================
echo "[*] Demultiplexing with ${DORADO_KIT}..."
mamba run -n nano dorado demux --kit-name "${DORADO_KIT}" --output-dir "${DEMUX_DIR}" "${BASECALLED_BAM}"


# =========================
# 3) Build index of demux file -> barcode number
# =========================
declare -A BARCODE_TO_FILE
while IFS= read -r f; do
  bn="$(basename "$f")"
  if [[ "$bn" =~ _barcode([0-9]+)\.bam$ ]]; then
    bnum="${BASH_REMATCH[1]}"
    BARCODE_TO_FILE["$bnum"]="$bn"
  fi
done < <(find "${DEMUX_DIR}" -maxdepth 1 -type f -name "*.bam" | sort)


# =========================
# 4) Convert demux BAM -> FASTQ.gz (ONLY samples in CSV) - DIRECT FIND
# =========================
echo "[*] Converting demultiplexed BAMs to FASTQ.gz (filtered by samples.csv)..."


# Get unique barcodes from CSV
awk -F',' 'NR>1 {print $2}' "${SAMPLE_CSV}" | tr -d '\r' | sort -u | while read nb; do
  bnum="${nb#NB}"
  
  # Find matching BAM directly in demux dir
  bam_file=$(find "${DEMUX_DIR}" -maxdepth 1 -name "*_barcode${bnum}.bam" -print -quit)
  
  if [[ -n "${bam_file}" ]]; then
    bn=$(basename "${bam_file}")
    out_fastq="${FASTQ_DIR}/${bn%.bam}.fastq.gz"
    echo "Converting barcode $bnum: $bn"
    mamba run -n nano samtools fastq -c 6 -@ 4 -0 "${out_fastq}" "${bam_file}"
  else
    echo "[!] Warning: Barcode ${nb} (num ${bnum}) not found in ${DEMUX_DIR}"
  fi
done


# =========================
# 5) Merge per-sample across barcodes (UNIQUE samples only)
# =========================
echo "[*] Merging FASTQ files per sample..."


# Get unique samples and process each once
awk -F',' 'NR>1 {print $1}' "${SAMPLE_CSV}" | tr -d '\r' | sort -u | while read -r sample; do
  [[ -z "${sample}" ]] && continue
  
  # Get all barcodes for this sample from CSV
  mapfile -t barcode_list < <(awk -F',' -v s="${sample}" 'NR>1 && $1==s {print $2}' "${SAMPLE_CSV}" | tr -d '\r')
  
  out_fq="${PER_SAMPLE_DIR}/${RUN_PREFIX}_${sample}.fastq.gz"
  
  if (( ${#barcode_list[@]} == 1 )); then
    # Single barcode - create symlink
    nb="${barcode_list[0]}"
    bnum="$(nb_to_num "${nb}")"
    demux_bam="${BARCODE_TO_FILE[$bnum]:-}"
    
    if [[ -n "${demux_bam}" ]]; then
      base="${demux_bam%.bam}"
      src_fastq="${FASTQ_DIR}/${base}.fastq.gz"
      
      if [[ -f "${src_fastq}" ]]; then
        echo "[*] Linking ${sample} (1 barcode: ${nb}) -> ${out_fq}"
        ln -sf "$(realpath "${src_fastq}")" "${out_fq}"
      fi
    fi
    
  else
    # Multiple barcodes - merge them
    fq_files=()
    echo "[*] Merging ${#barcode_list[@]} barcodes for ${sample}: ${barcode_list[*]}"
    
    for nb in "${barcode_list[@]}"; do
      bnum="$(nb_to_num "${nb}")"
      demux_bam="${BARCODE_TO_FILE[$bnum]:-}"
      
      if [[ -n "${demux_bam}" ]]; then
        base="${demux_bam%.bam}"
        src_fastq="${FASTQ_DIR}/${base}.fastq.gz"
        if [[ -f "${src_fastq}" ]]; then
          fq_files+=("${src_fastq}")
        fi
      fi
    done
    
    if (( ${#fq_files[@]} > 0 )); then
        echo "[*] Merging ${#fq_files[@]} FASTQ files for ${sample} -> ${out_fq}"  
        rm -f "${out_fq}"
        cat "${fq_files[@]}" > "${out_fq}"
        
        # Validate merged gzip integrity
        if gzip -t "${out_fq}" 2>/dev/null; then
          echo "[✓] Merged FASTQ gzip integrity verified"
        else
          echo "[!] ERROR: Merged FASTQ gzip check failed" >&2
          rm -f "${out_fq}"
          exit 1
        fi
    fi
  fi
done


# =========================
# 6) Create rename table
# =========================
echo -e "Sample\tBarcode\tBarcodeNum\tMergedFASTQ" > "${RENAME_TABLE}"


echo "[*] Building rename table with merged FASTQ paths..."
while IFS=',' read -r sample nb; do
  [[ -z "${sample}" || -z "${nb}" ]] && continue
  
  merged_fastq="${PER_SAMPLE_DIR}/${RUN_PREFIX}_${sample}.fastq.gz"
  
  if [[ ! -f "${merged_fastq}" ]]; then
    echo "[!] Warning: Merged FASTQ not found: ${merged_fastq}" >&2
    continue
  fi
  
  bnum="$(nb_to_num "${nb}")"
  echo -e "${sample}\t${nb}\t${bnum}\t${merged_fastq}" >> "${RENAME_TABLE}"
done < <(tail -n +2 "${SAMPLE_CSV}" | tr -d '\r')


# =========================
# 7) Final Cleanup
# =========================
if [[ "${START_TYPE}" == "archive" ]]; then
    echo "[*] Final cleanup of pod5_tmp..."
    rm -rf "${POD5_TMP_DIR:?}"
fi


echo ""
echo "============================================================"
echo "✓ TDAC file preprocessing complete"
echo "============================================================"
echo "Rename table: ${RENAME_TABLE}"
echo "Per-sample FASTQs: ${PER_SAMPLE_DIR}"
SAMPLE_COUNT=$(tail -n +2 "${RENAME_TABLE}" | cut -f1 | sort | uniq | wc -l)
echo "Total samples processed: ${SAMPLE_COUNT}"
echo "============================================================"
