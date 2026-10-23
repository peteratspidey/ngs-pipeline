#!/usr/bin/env bash
# ==============================================================================
# NGS WES Pipeline - Reference Data Download & Preparation Script
# Downloads and prepares human reference genome hg19 (GRCh37), dbSNP,
# Mills/1000G Indels, and ANNOVAR databases for reproducible analysis.
# ==============================================================================

set -euo pipefail

BASE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
GENOMES_DIR="${BASE_DIR}/resources/genomes/hg19"
ANNOVAR_DIR="${BASE_DIR}/resources/annovar"

mkdir -p "${GENOMES_DIR}"
mkdir -p "${ANNOVAR_DIR}"

echo "======================================================================"
echo " Preparing Reference Files for NGS WES Pipeline (hg19 / GRCh37)"
echo " Base Directory: ${BASE_DIR}"
echo "======================================================================"

# 1. Download and Index hg19 Reference Genome
echo "==> 1. Setting up hg19 FASTA reference..."
FASTA="${GENOMES_DIR}/hg19.fa"
if [ ! -f "${FASTA}" ]; then
    echo "Downloading hg19 FASTA from UCSC..."
    wget -c -P "${GENOMES_DIR}" http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz
    gunzip -f "${GENOMES_DIR}/hg19.fa.gz"
else
    echo "hg19 FASTA already exists at ${FASTA}"
fi

# Index FASTA with samtools
if [ ! -f "${FASTA}.fai" ]; then
    if command -v samtools >/dev/null 2>&1; then
        echo "Indexing hg19 with samtools faidx..."
        samtools faidx "${FASTA}"
    elif command -v docker >/dev/null 2>&1; then
        echo "Indexing hg19 with samtools via Docker..."
        docker run --rm -v "${GENOMES_DIR}:/data" quay.io/biocontainers/mulled-v2-fe8faa35dbf6dc65a0f7f5d4ea12e31a79f73e40:f45ad9036aa41bb10f875a330fa877d8869018a1-0 samtools faidx /data/hg19.fa
    fi
fi

# Create sequence dictionary with GATK or Picard
DICT="${GENOMES_DIR}/hg19.dict"
if [ ! -f "${DICT}" ]; then
    if command -v gatk >/dev/null 2>&1; then
        echo "Creating sequence dictionary with GATK..."
        gatk CreateSequenceDictionary -R "${FASTA}" -O "${DICT}"
    elif command -v picard >/dev/null 2>&1; then
        echo "Creating sequence dictionary with Picard..."
        picard CreateSequenceDictionary R="${FASTA}" O="${DICT}"
    elif command -v docker >/dev/null 2>&1; then
        echo "Creating sequence dictionary with GATK via Docker..."
        docker run --rm -v "${GENOMES_DIR}:/data" broadinstitute/gatk:4.5.0.0 gatk CreateSequenceDictionary -R /data/hg19.fa -O /data/hg19.dict
    fi
fi

# Index for BWA-MEM
if [ ! -f "${GENOMES_DIR}/hg19.bwt" ]; then
    if command -v bwa >/dev/null 2>&1; then
        echo "Building BWA index for hg19..."
        bwa index -p "${GENOMES_DIR}/hg19" "${FASTA}"
    elif command -v docker >/dev/null 2>&1; then
        echo "Building BWA index for hg19 via Docker..."
        docker run --rm -v "${GENOMES_DIR}:/data" quay.io/biocontainers/mulled-v2-fe8faa35dbf6dc65a0f7f5d4ea12e31a79f73e40:f45ad9036aa41bb10f875a330fa877d8869018a1-0 bwa index -p /data/hg19 /data/hg19.fa
    fi
fi

# 2. Download Known Sites for BQSR and Variant Calling
echo "==> 2. Setting up dbSNP & Known InDels for BQSR..."
DBSNP="${GENOMES_DIR}/dbsnp_138.b37.vcf.gz"
if [ ! -f "${DBSNP}" ]; then
    echo "Downloading dbSNP 138 (b37/hg19)..."
    wget -c -P "${GENOMES_DIR}" https://ftp.ncbi.nih.gov/snp/organisms/human_9606_b150_GRCh37p13/VCF/00-All.vcf.gz -O "${DBSNP}" || true
fi

# 3. ANNOVAR Human Databases
echo "==> 3. Setting up ANNOVAR databases for hg19..."
if command -v annotate_variation.pl >/dev/null 2>&1; then
    echo "Downloading ANNOVAR refGene database..."
    annotate_variation.pl -buildver hg19 -downdb -webfrom annovar refGene "${ANNOVAR_DIR}" || true
    echo "Downloading ANNOVAR avsnp150 database..."
    annotate_variation.pl -buildver hg19 -downdb -webfrom annovar avsnp150 "${ANNOVAR_DIR}" || true
    echo "Downloading ANNOVAR clinvar database..."
    annotate_variation.pl -buildver hg19 -downdb -webfrom annovar clinvar_20221231 "${ANNOVAR_DIR}" || true
    echo "Downloading ANNOVAR 1000G database..."
    annotate_variation.pl -buildver hg19 -downdb -webfrom annovar 1000g2015aug "${ANNOVAR_DIR}" || true
    echo "Downloading ANNOVAR exac03 database..."
    annotate_variation.pl -buildver hg19 -downdb -webfrom annovar exac03 "${ANNOVAR_DIR}" || true
    echo "Downloading ANNOVAR dbnsfp database (SIFT, PolyPhen, CADD)..."
    annotate_variation.pl -buildver hg19 -downdb -webfrom annovar dbnsfp42a "${ANNOVAR_DIR}" || true
else
    echo "Notice: 'annotate_variation.pl' not in PATH. ANNOVAR fallback annotator will format output using 'bin/annotate_vcf.py'."
fi

echo "======================================================================"
echo " Reference setup completed successfully."
echo "======================================================================"

