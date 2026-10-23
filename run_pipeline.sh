#!/usr/bin/env bash
# ==============================================================================
# Biokart WES Pipeline - One-Click Execution Wrapper
# Runs the full Next Generation Sequencing Whole Exome Sequencing pipeline
# ==============================================================================

set -euo pipefail

SOURCE="${BASH_SOURCE[0]}"
while [ -h "$SOURCE" ]; do
    DIR="$( cd -P "$( dirname "$SOURCE" )" >/dev/null 2>&1 && pwd )"
    SOURCE="$(readlink "$SOURCE")"
    [[ $SOURCE != /* ]] && SOURCE="$DIR/$SOURCE"
done
SCRIPT_DIR="$( cd -P "$( dirname "$SOURCE" )" >/dev/null 2>&1 && pwd )"
cd "${SCRIPT_DIR}"

# Set local Nextflow home directory to prevent root/home permission issues
export NXF_HOME="${SCRIPT_DIR}/.nextflow"
mkdir -p "${NXF_HOME}"

# Determine Nextflow binary
if [ -x "${SCRIPT_DIR}/bin/nextflow" ]; then
    NEXTFLOW_BIN="${SCRIPT_DIR}/bin/nextflow"
elif command -v nextflow >/dev/null 2>&1; then
    NEXTFLOW_BIN="$(command -v nextflow)"
else
    echo "==> Nextflow not found. Downloading Nextflow..."
    curl -fsSL https://get.nextflow.io | bash
    mv nextflow "${SCRIPT_DIR}/bin/nextflow"
    chmod +x "${SCRIPT_DIR}/bin/nextflow"
    NEXTFLOW_BIN="${SCRIPT_DIR}/bin/nextflow"
fi

# Default parameters
INPUT_DIR="${SCRIPT_DIR}/data"
INPUT_PATTERN=""
SAMPLE_SHEET=""
GENOME="hg19"
FASTA="${SCRIPT_DIR}/resources/genomes/hg19/hg19.fa"
BED="${SCRIPT_DIR}/resources/exome_targets.bed"
OUTDIR="${SCRIPT_DIR}/results"
TRIMMER="cutadapt"
QUALITY_CUTOFF="30"
ANNOTATOR="annovar"
PROFILE="standard"
TEST_MODE=false

# Auto-detect container engine
if ! command -v docker >/dev/null 2>&1; then
    PROFILE="nocontainer"
fi

# Help message
usage() {
    cat <<EOF
Usage: ./run_pipeline.sh [OPTIONS]

Options:
  -i, --input <pattern>       Input FASTQ glob pattern (default: 'data/*_{R1,R2}*.fastq.gz')
  -s, --sample-sheet <file>   Path to sample sheet CSV (overrides --input)
  -g, --genome <hg19>         Reference genome identifier (default: hg19)
  -f, --fasta <path>          Path to reference FASTA (default: resources/genomes/hg19/hg19.fa)
  -b, --bed <path>            Path to exome target BED file
  -o, --outdir <dir>          Output directory (default: ./results)
  -p, --profile <profile>     Nextflow profile: standard (Docker), nocontainer (default: auto-detected)
  -t, --test                  Run self-contained pipeline test with mock data
  -h, --help                  Show this help message

Examples:
  # 1. Run quick test mode
  ./run_pipeline.sh --test

  # 2. Run with FASTQ files in data/
  ./run_pipeline.sh -i 'data/*_{R1,R2}*.fastq.gz'

  # 3. Run with a sample sheet and custom exome BED
  ./run_pipeline.sh -s resources/sample_sheet.csv -b resources/my_targets.bed
EOF
    exit 0
}

# Parse CLI options
while [[ $# -gt 0 ]]; do
    case "$1" in
        -i|--input)
            INPUT_PATTERN="$2"
            shift 2
            ;;
        -s|--sample-sheet)
            SAMPLE_SHEET="$2"
            shift 2
            ;;
        -g|--genome)
            GENOME="$2"
            shift 2
            ;;
        -f|--fasta)
            FASTA="$2"
            shift 2
            ;;
        -b|--bed)
            BED="$2"
            shift 2
            ;;
        -o|--outdir)
            OUTDIR="$2"
            shift 2
            ;;
        -p|--profile)
            PROFILE="$2"
            shift 2
            ;;
        -t|--test)
            TEST_MODE=true
            shift
            ;;
        -h|--help)
            usage
            ;;
        *)
            echo "Unknown option: $1"
            usage
            ;;
    esac
done

echo "========================================================================"
echo "          BIOKART WHOLE EXOME SEQUENCING (WES) PIPELINE                "
echo "========================================================================"
echo " Nextflow Binary : ${NEXTFLOW_BIN}"
echo " Execution Profile: ${PROFILE}"
echo " Reference Genome : ${GENOME}"
echo " Output Directory : ${OUTDIR}"
echo "========================================================================"

if [ "${TEST_MODE}" = true ]; then
    echo "==> Running Self-Contained Test Suite..."
    ./tests/setup_test_data.sh
    python3 tests/test_report_generation.py
    echo ""
    echo "==> Test Completed Successfully!"
    echo "Generated Sample Tables:"
    echo " - Table 1: ${SCRIPT_DIR}/tests/mock_run/reports/Table1_QC_stats.tsv"
    echo " - Table 2: ${SCRIPT_DIR}/tests/mock_run/reports/Table2_Variant_stats.tsv"
    exit 0
fi

# Build Nextflow CLI arguments
NF_ARGS=(
    run "${SCRIPT_DIR}/main.nf"
    -profile "${PROFILE}"
    --genome "${GENOME}"
    --outdir "${OUTDIR}"
    --trimmer "${TRIMMER}"
    --quality_cutoff "${QUALITY_CUTOFF}"
    --annotator "${ANNOTATOR}"
)

if [ -n "${SAMPLE_SHEET}" ]; then
    NF_ARGS+=(--sample_sheet "${SAMPLE_SHEET}")
elif [ -n "${INPUT_PATTERN}" ]; then
    NF_ARGS+=(--input "${INPUT_PATTERN}")
else
    mkdir -p "${INPUT_DIR}"
    NF_ARGS+=(--input "${INPUT_DIR}/*_{R1,R2}*.fastq.gz")
fi

if [ -f "${FASTA}" ]; then
    NF_ARGS+=(--fasta "${FASTA}")
fi

if [ -f "${BED}" ]; then
    NF_ARGS+=(--bed "${BED}")
fi

# Execute Nextflow
echo "==> Launching Pipeline Execution..."
"${NEXTFLOW_BIN}" "${NF_ARGS[@]}"

echo ""
echo "========================================================================"
echo "                       PIPELINE EXECUTION COMPLETE                      "
echo "========================================================================"
echo " Deliverables Generated:"
echo "  - Table 1 (QC Before/After Trimming) : ${OUTDIR}/reports/Table1_QC_stats.html"
echo "  - Table 2 (Total SNVs & InDels)      : ${OUTDIR}/reports/Table2_Variant_stats.html"
echo "  - ANNOVAR Functional Annotations     : ${OUTDIR}/annotation/"
echo "  - MultiQC Aggregated Report          : ${OUTDIR}/reports/multiqc/multiqc_report.html"
echo "  - Recalibrated BAMs                  : ${OUTDIR}/aligned/"
echo "  - Filtered VCFs                      : ${OUTDIR}/variants/"
echo "========================================================================"
