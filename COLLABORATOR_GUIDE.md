# Whole Exome Sequencing (WES) Pipeline
## Quadro Workstation Quickstart Guide

This guide walks you through setting up and running the Whole Exome Sequencing post-NGS analysis pipeline on your **NVIDIA Quadro workstation**.

### Step 0: One-Time System Prerequisite (Java 17)

Nextflow requires **Java 17 or later**. Install it with:

```bash
sudo apt update && sudo apt install -y openjdk-17-jre-headless
```

---

### Step 1: Clone the Repository

Open a terminal on your Quadro workstation and run:

```bash
git clone https://github.com/peteratspidey/ngs-pipeline.git
cd ngs-pipeline
```

---

### Step 2: Run the Quick Test (Takes ~2 Seconds)

Before running real data, run the built-in test to verify that the report generator, Table 1, and Table 2 logic work properly:

```bash
./run_pipeline.sh --test
```

If you see:
```text
ALL REPORT & ANNOTATION TESTS PASSED SUCCESSFULLY!
```
Everything is set up correctly and ready to go.

---

### Step 3: Working with FASTQ Files

#### Option A: Download Mini Test FASTQ Reads (~16 MB)
To test the pipeline on real Illumina paired-end sequencing reads:

```bash
./scripts/download_test_data.sh TEST_SAMPLE
```
*(This automatically downloads `TEST_SAMPLE_R1.fastq.gz` and `TEST_SAMPLE_R2.fastq.gz` into the `data/` directory)*

#### Option B: Use Your Own Sequencing Data
Copy your paired-end FASTQ files into the `data/` folder:
- `data/YOURSAMPLE_R1.fastq.gz`
- `data/YOURSAMPLE_R2.fastq.gz`

---

### Step 4: Run the Pipeline on Quadro

Run the pipeline using the dedicated **Quadro profile** (automatically utilizes high CPU cores, 128 GB RAM, and Docker):

```bash
./run_pipeline.sh -p quadro
```

> **Note**: If your FASTQ files are in a custom directory outside `data/`, simply pass the path with `-i`:
> ```bash
> ./run_pipeline.sh -p quadro -i '/path/to/my_fastqs/*_{R1,R2}*.fastq.gz'
> ```

---

### Step 5: Check Your Results

All results will be saved in the `results/` folder:

| File / Folder | What It Contains |
| :--- | :--- |
| **`results/reports/Table1_QC_stats.html`** | **Table 1**: Pre- & post-trimming QC statistics (Reads, GC%, Duplications, GB throughput) |
| **`results/reports/Table2_Variant_stats.html`** | **Table 2**: Total SNVs & small InDels per sample |
| **`results/reports/multiqc/multiqc_report.html`** | Interactive MultiQC quality report |
| **`results/aligned/`** | Recalibrated, analysis-ready BAM files |
| **`results/variants/`** | Filtered VCF variant files |
| **`results/annotation/`** | ANNOVAR functional annotations |

---

### Quadro Workstation Specs & Settings

The `-p quadro` profile is pre-configured specifically for your workstation hardware:
- **CPUs**: Automatically allocates up to 32 parallel threads for BWA-MEM alignment and GATK variant calling.
- **Memory**: Allocates up to 128 GB RAM with shared memory access (`--ipc=host`).
- **Containers**: Uses Docker automatically—no need to install bioinformatics tools manually.
- **Nextflow**: Bootstraps automatically from `bin/nextflow` (or runs inside Docker if Java is not installed on the host).
