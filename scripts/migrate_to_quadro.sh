#!/usr/bin/env bash
# ==============================================================================
# Biokart WES Pipeline - Quadro Workstation Migration Helper
# Facilitates moving and deploying the pipeline onto a Quadro workstation
# ==============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "${SCRIPT_DIR}"

GITHUB_REPO="https://github.com/peteratspidey/ngs-pipeline.git"

usage() {
    cat <<EOF
Usage: ./scripts/migrate_to_quadro.sh [OPTIONS]

Options:
  -b, --bundle               Create a standalone portable .tar.gz bundle for USB/SCP transfer
  -s, --sync <user@host:dir> Sync directly to Quadro workstation via rsync over SSH
  -p, --pre-pull             Pre-pull all container images required by the pipeline
  -g, --git-info             Show Git clone commands to run on the Quadro system
  -h, --help                 Display this help message

Examples:
  # 1. Create a transfer bundle
  ./scripts/migrate_to_quadro.sh --bundle

  # 2. Directly sync to Quadro machine over SSH
  ./scripts/migrate_to_quadro.sh --sync user@192.168.1.100:/home/user/ngs-pipeline

  # 3. Pre-pull all Docker container images
  ./scripts/migrate_to_quadro.sh --pre-pull
EOF
    exit 0
}

# Pre-pull container images
pre_pull_images() {
    echo "==> Pre-pulling container images for offline / fast execution..."
    IMAGES=(
        "quay.io/biocontainers/fastqc:0.12.1--hdfd78af_0"
        "quay.io/biocontainers/cutadapt:4.4--py39hf95cd2a_1"
        "quay.io/biocontainers/mulled-v2-fe8faa35dbf6dc65a0f7f5d4ea12e31a79f73e40:21960a4acb7e3e414f43403865b6e3952e46b080-0"
        "broadinstitute/gatk:4.5.0.0"
        "quay.io/biocontainers/bcftools:1.19--h8b25389_0"
        "quay.io/biocontainers/multiqc:1.19--pyhdfd78af_0"
        "python:3.11-slim"
    )

    if command -v docker >/dev/null 2>&1; then
        for img in "${IMAGES[@]}"; do
            echo "--> Pulling ${img}..."
            docker pull "${img}"
        done
        echo "==> All Docker container images successfully pulled!"
    elif command -v apptainer >/dev/null 2>&1; then
        mkdir -p apptainer_cache
        for img in "${IMAGES[@]}"; do
            echo "--> Pulling ${img} into Apptainer cache..."
            apptainer pull --dir apptainer_cache "docker://${img}"
        done
        echo "==> All Apptainer container images successfully cached!"
    else
        echo "Neither Docker nor Apptainer found. Please install a container engine first."
        exit 1
    fi
}

# Create a clean portable bundle
create_bundle() {
    BUNDLE_NAME="ngs-pipeline-quadro-$(date +%Y%m%d).tar.gz"
    echo "==> Creating clean portable archive: ${BUNDLE_NAME}..."
    tar -czf "${SCRIPT_DIR}/${BUNDLE_NAME}" \
        --exclude=".nextflow" \
        --exclude="work" \
        --exclude="results" \
        --exclude="*.log" \
        --exclude="singularity_cache" \
        --exclude="apptainer_cache" \
        --exclude="*.tar.gz" \
        -C "$(dirname "${SCRIPT_DIR}")" "$(basename "${SCRIPT_DIR}")"

    echo "==> Bundle created successfully: ${SCRIPT_DIR}/${BUNDLE_NAME}"
    echo "    Size: $(du -h "${SCRIPT_DIR}/${BUNDLE_NAME}" | cut -f1)"
    echo ""
    echo "To unpack on your Quadro system:"
    echo "  tar -xzf ${BUNDLE_NAME}"
    echo "  cd ngs-pipeline"
    echo "  ./run_pipeline.sh -p quadro"
}

# Sync directly via rsync over SSH
sync_to_host() {
    local target="$1"
    echo "==> Syncing ngs-pipeline directly to ${target} via rsync..."
    rsync -avzP \
        --exclude=".nextflow" \
        --exclude="work" \
        --exclude="results" \
        --exclude="*.log" \
        --exclude="singularity_cache" \
        --exclude="apptainer_cache" \
        "${SCRIPT_DIR}/" "${target}/"

    echo "==> Direct sync complete!"
    echo "On your Quadro system, run:"
    echo "  ./run_pipeline.sh -p quadro"
}

# Git info
show_git_info() {
    echo "========================================================================"
    echo "           MIGRATING VIA GITHUB TO YOUR QUADRO WORKSTATION              "
    echo "========================================================================"
    echo "On your Quadro machine terminal, run:"
    echo ""
    echo "  git clone ${GITHUB_REPO}"
    echo "  cd ngs-pipeline"
    echo ""
    echo "To run with Quadro high-performance configuration:"
    echo "  ./run_pipeline.sh -p quadro"
    echo "========================================================================"
}

if [ $# -eq 0 ]; then
    show_git_info
    exit 0
fi

while [ $# -gt 0 ]; do
    case "$1" in
        -b|--bundle)
            create_bundle
            shift
            ;;
        -s|--sync)
            if [ -z "${2:-}" ]; then
                echo "Error: --sync requires a remote target (e.g. user@hostname:/path)"
                exit 1
            fi
            sync_to_host "$2"
            shift 2
            ;;
        -p|--pre-pull)
            pre_pull_images
            shift
            ;;
        -g|--git-info)
            show_git_info
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
