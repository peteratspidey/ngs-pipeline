#!/usr/bin/env bash

set -euo pipefail

SOURCE="${BASH_SOURCE[0]}"
while [ -h "$SOURCE" ]; do
    DIR="$( cd -P "$( dirname "$SOURCE" )" >/dev/null 2>&1 && pwd )"
    SOURCE="$(readlink "$SOURCE")"
    [[ $SOURCE != /* ]] && SOURCE="$DIR/$SOURCE"
done
SCRIPT_DIR="$( cd -P "$( dirname "$SOURCE" )" >/dev/null 2>&1 && pwd )"
cd "${SCRIPT_DIR}"

export NXF_HOME="${SCRIPT_DIR}/.nextflow"
mkdir -p "${NXF_HOME}"

# Ensure bin directory exists
mkdir -p "${SCRIPT_DIR}/bin"

USE_DOCKER_NEXTFLOW=false
NEXTFLOW_BIN=""

ensure_nextflow() {
    # 1. Check local bin/nextflow
    if [ -x "${SCRIPT_DIR}/bin/nextflow" ]; then
        NEXTFLOW_BIN="${SCRIPT_DIR}/bin/nextflow"
        return 0
    fi

    # 2. Check system PATH
    if command -v nextflow >/dev/null 2>&1; then
        NEXTFLOW_BIN="$(command -v nextflow)"
        return 0
    fi

    echo "========================================================================"
    echo " Nextflow not found. Initializing automatic download..."
    echo "========================================================================"

    # Method A: Direct release binary download (fast, reliable, no GitHub API limits)
    echo "--> Downloading Nextflow executable from GitHub Releases..."
    if curl -fsSL -o "${SCRIPT_DIR}/bin/nextflow" "https://github.com/nextflow-io/nextflow/releases/download/v23.10.1/nextflow" 2>/dev/null; then
        chmod +x "${SCRIPT_DIR}/bin/nextflow"
        NEXTFLOW_BIN="${SCRIPT_DIR}/bin/nextflow"
        echo " Nextflow downloaded successfully to ${NEXTFLOW_BIN}"
        return 0
    fi

    # Method B: Official get.nextflow.io installer
    echo "--> Trying get.nextflow.io installer..."
    if curl -fsSL https://get.nextflow.io 2>/dev/null | bash 2>/dev/null; then
        if [ -f nextflow ]; then
            mv nextflow "${SCRIPT_DIR}/bin/nextflow"
            chmod +x "${SCRIPT_DIR}/bin/nextflow"
            NEXTFLOW_BIN="${SCRIPT_DIR}/bin/nextflow"
            echo " Nextflow downloaded successfully."
            return 0
        fi
    fi

    # Method C: wget fallback
    if command -v wget >/dev/null 2>&1; then
        echo "--> Trying wget download..."
        if wget -q -O "${SCRIPT_DIR}/bin/nextflow" "https://github.com/nextflow-io/nextflow/releases/download/v23.10.1/nextflow"; then
            chmod +x "${SCRIPT_DIR}/bin/nextflow"
            NEXTFLOW_BIN="${SCRIPT_DIR}/bin/nextflow"
            echo " Nextflow downloaded successfully via wget."
            return 0
        fi
    fi

    # Method D: Docker fallback
    if command -v docker >/dev/null 2>&1; then
        echo "==> Host Nextflow binary download blocked. Falling back to containerized Nextflow runner via Docker..."
        USE_DOCKER_NEXTFLOW=true
        return 0
    fi

    echo "ERROR: Could not download Nextflow automatically and Docker is not available."
    echo "Please run: curl -fsSL https://get.nextflow.io | bash && mv nextflow bin/"
    exit 1
}

ensure_java() {
    local has_java=false
    local java_major=0

    if command -v java >/dev/null 2>&1; then
        has_java=true
        local ver_str
        ver_str=$(java -version 2>&1 | head -n1)
        java_major=$(echo "$ver_str" | sed -n 's/.*version "\([0-9]*\)\..*/\1/p')
        if [ "$java_major" = "1" ]; then
            java_major=$(echo "$ver_str" | sed -n 's/.*version "1\.\([0-9]*\)\..*/\1/p')
        fi
        java_major="${java_major:-0}"
    fi

    # Nextflow requires Java >= 17 (up to 26)
    if [ "${has_java}" = true ] && [ "${java_major}" -ge 17 ] && [ "${java_major}" -le 26 ]; then
        return 0
    fi

    if [ "${has_java}" = true ]; then
        echo "==> Host Java version (${java_major}) is below the required Java 17."
    else
        echo "==> Java runtime not found on host."
    fi

    # If Docker is available, run Nextflow containerized (bundled Java 17)
    if command -v docker >/dev/null 2>&1; then
        echo "==> Docker detected! Nextflow will execute inside a container (bundled with Java 17)."
        USE_DOCKER_NEXTFLOW=true
        return 0
    fi

    echo "==> Attempting auto-install of OpenJDK 17..."
    if command -v apt-get >/dev/null 2>&1 && [ "$(id -u)" -eq 0 ]; then
        apt-get update && apt-get install -y openjdk-17-jre-headless
    elif command -v conda >/dev/null 2>&1; then
        conda install -y -c conda-forge openjdk=17
    else
        echo "========================================================================"
        echo " ERROR: Nextflow requires Java 17 or later."
        echo " Run this command on your machine to install Java 17:"
        echo "   sudo apt update && sudo apt install -y openjdk-17-jre-headless"
        echo "========================================================================"
        exit 1
    fi
}

ensure_nextflow
ensure_java

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
if command -v docker >/dev/null 2>&1; then
    PROFILE="standard"
elif command -v apptainer >/dev/null 2>&1; then
    PROFILE="apptainer"
elif command -v singularity >/dev/null 2>&1; then
    PROFILE="singularity"
else
    PROFILE="nocontainer"
fi

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
  -p, --profile <profile>     Nextflow profile: standard (Docker), quadro (high-perf workstation), apptainer, singularity, nocontainer
  -t, --test                  Run self-contained pipeline test with mock data
  -h, --help                  Show this help message

Examples:
  ./run_pipeline.sh --test
  ./run_pipeline.sh -i 'data/*_{R1,R2}*.fastq.gz'
  ./run_pipeline.sh -s resources/sample_sheet.csv -b resources/my_targets.bed
EOF
    exit 0
}

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
echo "          NGS WHOLE EXOME SEQUENCING (WES) PIPELINE                "
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

if [ "${USE_DOCKER_NEXTFLOW}" = true ]; then
    if ! docker image inspect nextflow/nextflow:23.10.1 >/dev/null 2>&1; then
        echo "--> Pulling Nextflow container image (with auto-retry on connection drops)..."
        for attempt in 1 2 3 4; do
            if docker pull nextflow/nextflow:23.10.1; then
                break
            fi
            echo "--> Download attempt ${attempt} interrupted by network. Retrying in 3 seconds..."
            sleep 3
        done
    fi
    echo "==> Launching Pipeline via Dockerized Nextflow..."
    docker run -i --rm \
        -u "$(id -u):$(id -g)" \
        -v /var/run/docker.sock:/var/run/docker.sock \
        -v "${SCRIPT_DIR}:${SCRIPT_DIR}" \
        -w "${SCRIPT_DIR}" \
        -e NXF_HOME="${NXF_HOME}" \
        nextflow/nextflow:23.10.1 "${NF_ARGS[@]}"
else
    echo "==> Launching Pipeline Execution..."
    "${NEXTFLOW_BIN}" "${NF_ARGS[@]}"
fi

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
