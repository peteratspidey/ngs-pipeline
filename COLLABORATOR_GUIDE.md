# Biokart Whole Exome Sequencing (WES) Pipeline
## Collaborator & Deployment Handover Guide

> **Repository**: [https://github.com/peteratspidey/ngs-pipeline](https://github.com/peteratspidey/ngs-pipeline)  
> **Workflow Specification**: Biokart India WES Protocol (SNV & small InDel discovery and annotation)  
> **Target Environment**: Workstation / HPC / NVIDIA Quadro Server (Docker or Apptainer)

---

## 1. Project Overview

This repository contains an end-to-end, reproducible **Post-NGS Analysis Pipeline** implemented in **Nextflow DSL2** with containerized processes. It automates secondary and tertiary whole exome analysis from raw demultiplexed FASTQ reads to clinical QC metrics and annotated variant tables.

```
                    ┌─────────────────────────────────────────┐
                    │ Raw Paired-End FASTQs (*_R1 / *_R2)     │
                    └───────────────────┬─────────────────────┘
                                        ▼
                    ┌─────────────────────────────────────────┐
                    │ FastQC (Pre-trimming QC)                │
                    └───────────────────┬─────────────────────┘
                                        ▼
                    ┌─────────────────────────────────────────┐
                    │ Cutadapt Trimming (Q30 cutoff, 5'/3')   │
                    └───────────────────┬─────────────────────┘
                                        ▼
                    ┌─────────────────────────────────────────┐
                    │ FastQC (Post-trimming HQ reads)         │
                    └───────────────────┬─────────────────────┘
                                        ▼
                    ┌─────────────────────────────────────────┐
                    │ BWA-MEM Alignment (hg19 / GRCh37)       │
                    └───────────────────┬─────────────────────┘
                                        ▼
                    ┌─────────────────────────────────────────┐
                    │ Picard MarkDuplicates + GATK BQSR       │
                    └───────────────────┬─────────────────────┘
                                        ▼
                    ┌─────────────────────────────────────────┐
                    │ GATK HaplotypeCaller (-L Exome Targets) │
                    └───────────────────┬─────────────────────┘
                                        ▼
                    ┌─────────────────────────────────────────┐
                    │ ANNOVAR Functional Annotation           │
                    │ (refGene, avsnp150, 1000g, ClinVar, ...)│
                    └───────────────────┬─────────────────────┘
                                        ▼
  ┌─────────────────────────┬─────────────────────────┬─────────────────────────┐
  ▼                         ▼                         ▼                         ▼
Table 1 (QC Statistics)   Table 2 (SNV/InDel Counts) MultiQC Dashboard       Annotated VCFs
```

---

## 2. System Requirements

| Component | Minimum | Recommended (Quadro Workstation) |
| :--- | :--- | :--- |
| **Operating System** | Linux (Ubuntu 20.04/22.04+, RHEL/Rocky 8/9) | Ubuntu 22.04 LTS / 24.04 LTS |
| **CPU** | 8 cores / 16 threads | 16–32+ physical cores (AMD Threadripper / Xeon) |
| **RAM** | 32 GB | 64 GB – 128+ GB |
| **Storage** | 100 GB free space (SSD/NVMe) | 500 GB+ NVMe scratch space |
| **Container Engine** | Docker 20.10+ OR Apptainer 1.1+ | Docker or Apptainer (Singularity) |
| **Runtime** | Java 11+ (OpenJRE) | Java 17 LTS / OpenJDK |

---

## 3. Quick Setup on Your Workstation

### Step 1: Clone the Repository
```bash
git clone https://github.com/peteratspidey/ngs-pipeline.git
cd ngs-pipeline
```

### Step 2: (Optional) Pre-Pull Container Images
Avoid download delays during live runs by caching the container images in advance:
```bash
./scripts/migrate_to_quadro.sh --pre-pull
```

---

## 4. How to Run the Pipeline

The project includes a unified runner script ([`run_pipeline.sh`](file:///home/peter/Desktop/NGS/ngs-pipeline/run_pipeline.sh)) that handles Nextflow execution, container configuration, and directory binding.

### A. Run Self-Contained Sanity Test (~2 Seconds)
Verifies the reporting engine, FastQC parsing, Table 1 & Table 2 generator, and annotation fallbacks:
```bash
./run_pipeline.sh --test
```

### B. Download Mini Test FASTQs & Run
Download a paired-end human test dataset (~16 MB each) and execute the full pipeline:
```bash
# 1. Download test reads into data/
./scripts/download_test_data.sh TEST_SAMPLE

# 2. Run with standard Docker
./run_pipeline.sh

# 3. OR Run with high-performance Quadro profile
./run_pipeline.sh -p quadro
```

### C. Run on Production Real Exome Sequencing Data
Place your raw paired-end sequencing reads into `data/` formatted as:
`SAMPLE_R1.fastq.gz` and `SAMPLE_R2.fastq.gz`.

```bash
# Standard run with Docker:
./run_pipeline.sh -i 'data/*_{R1,R2}*.fastq.gz' -b resources/exome_targets.bed

# Quadro workstation run (32 cores, 128GB RAM, --ipc=host):
./run_pipeline.sh -i 'data/*_{R1,R2}*.fastq.gz' -b resources/exome_targets.bed -p quadro

# Institutional / Rootless HPC run with Apptainer:
./run_pipeline.sh -i 'data/*_{R1,R2}*.fastq.gz' -b resources/exome_targets.bed -p quadro_apptainer
```

---

## 5. Command-Line Reference

| Flag | Description | Default |
| :--- | :--- | :--- |
| `-i, --input <glob>` | Input FASTQ paired-end pattern | `'data/*_{R1,R2}*.fastq.gz'` |
| `-s, --sample-sheet <csv>` | Path to sample sheet CSV (overrides `--input`) | None |
| `-b, --bed <path>` | Exome capture target BED file | `resources/exome_targets.bed` |
| `-f, --fasta <path>` | Reference genome FASTA (hg19) | `resources/genomes/hg19/hg19.fa` |
| `-o, --outdir <dir>` | Results destination folder | `./results` |
| `-p, --profile <name>` | Execution profile (`standard`, `quadro`, `quadro_apptainer`, `nocontainer`) | Auto-detected |
| `-t, --test` | Run self-contained test suite | `false` |
| `-h, --help` | Display CLI help menu | |

---

## 6. Downloading Full Reference Genomes & Databases

For production whole-exome variant calling on `hg19`, download and index the reference FASTA, dbSNP VCF, and ANNOVAR databases:
```bash
./scripts/download_references.sh
```
This fetches:
- **Reference Genome**: UCSC / 1000G `hg19.fa` + BWA index + Samtools `.fai` + Picard `.dict`.
- **Known Sites for BQSR**: dbSNP (`dbsnp_138.hg19.vcf.gz`) + Mills and 1000G gold standard InDels.
- **ANNOVAR Databases**: `refGene`, `avsnp150`, `1000g2015aug_all`, `exac03`, `clinvar_20221231`, `dbnsfp42a`.

---

## 7. Deliverables & Output Directory Structure

Upon completion, all outputs are deposited in `results/`:

```
results/
├── reports/
│   ├── Table1_QC_stats.html          <-- Exact Table 1: Raw vs High-Quality QC metrics
│   ├── Table1_QC_stats.tsv           <-- Machine-readable TSV format
│   ├── Table2_Variant_stats.html     <-- Exact Table 2: Total SNVs & small InDels
│   └── Table2_Variant_stats.tsv      <-- Machine-readable TSV format
├── qc/
│   ├── fastqc_raw/                   <-- FastQC HTML & ZIP reports for raw reads
│   ├── fastqc_trimmed/               <-- FastQC HTML & ZIP reports for Q30-trimmed reads
│   └── multiqc/
│       └── multiqc_report.html       <-- Interactive MultiQC dashboard
├── trimmed/                          <-- Cutadapt adapter-free, Q30-filtered FASTQs
├── aligned/
│   ├── *.markdup.bam / *.bai         <-- Deduplicated sorted alignments
│   └── *.bqsr.bam / *.bai            <-- Recalibrated analysis-ready BAMs
├── variants/
│   ├── *.raw.vcf.gz                  <-- Raw GATK HaplotypeCaller variants
│   └── *.filtered.vcf.gz             <-- Hard-filtered analysis-ready VCF
└── annotation/
    ├── *.hg19_multianno.txt          <-- ANNOVAR tab-delimited annotation table
    └── *.hg19_multianno.vcf          <-- Fully annotated VCF file
```

---

## 8. Troubleshooting & FAQ

### Q1: "Got permission denied while trying to connect to the Docker daemon"
Ensure your user belongs to the `docker` group on the workstation:
```bash
sudo usermod -aG docker $USER
newgrp docker
```

### Q2: Running on an institutional server without Docker sudo access
Use the Apptainer/Singularity profile:
```bash
./run_pipeline.sh -p quadro_apptainer
```

### Q3: How to resume a pipeline after interruption
Nextflow caches intermediate files in `work/`. Simply rerun your previous command with `-resume`:
```bash
./bin/nextflow run main.nf -profile quadro -resume
```

---

## 9. Contact & Support
- **Repository Issues**: [https://github.com/peteratspidey/ngs-pipeline/issues](https://github.com/peteratspidey/ngs-pipeline/issues)
- **Author**: Peter (`peteratspidey`)

