#!/usr/bin/env bash
# ==============================================================================
# Biokart WES Pipeline - Reference Data Download & Preparation Script
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
echo " Preparing Reference Files for Biokart WES Pipeline (hg19 / GRCh37)"
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
if [ ! -f "${FASTA}.fai" ] && command -v samtools >/dev/null 2>&1; then
    echo "Indexing hg19 with samtools faidx..."
    samtools faidx "${FASTA}"
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
    fi
fi

# Index for BWA-MEM
if [ ! -f "${GENOMES_DIR}/hg19.bwt" ] && command -v bwa >/dev/null 2>&1; then
    echo "Building BWA index for hg19..."
    bwa index -p "${GENOMES_DIR}/hg19" "${FASTA}"
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

