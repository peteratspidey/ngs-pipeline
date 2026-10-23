#!/usr/bin/env bash
# ==============================================================================
# Download small human test FASTQ paired-end dataset for pipeline testing
# ==============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
DATA_DIR="${SCRIPT_DIR}/data"

mkdir -p "${DATA_DIR}"

SAMPLE_NAME="${1:-TEST_SAMPLE}"

echo "======================================================================"
echo " Downloading test paired-end FASTQ dataset (~16 MB each)"
echo " Target: ${DATA_DIR}/${SAMPLE_NAME}_R1.fastq.gz"
echo "         ${DATA_DIR}/${SAMPLE_NAME}_R2.fastq.gz"
echo "======================================================================"

URL_R1="https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/genomics/homo_sapiens/illumina/fastq/test_1.fastq.gz"
URL_R2="https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/genomics/homo_sapiens/illumina/fastq/test_2.fastq.gz"

echo "--> Fetching Read 1 (R1)..."
curl -L -o "${DATA_DIR}/${SAMPLE_NAME}_R1.fastq.gz" "${URL_R1}"

echo "--> Fetching Read 2 (R2)..."
curl -L -o "${DATA_DIR}/${SAMPLE_NAME}_R2.fastq.gz" "${URL_R2}"

echo ""
echo " Successfully downloaded test files to:"
ls -lh "${DATA_DIR}/${SAMPLE_NAME}_R1.fastq.gz" "${DATA_DIR}/${SAMPLE_NAME}_R2.fastq.gz"
echo ""
echo "To run the pipeline with these test files, run:"
echo "  ./run_pipeline.sh"
echo "======================================================================"

