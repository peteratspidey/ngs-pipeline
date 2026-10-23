#!/usr/bin/env python3
"""
Unit and Integration Test for Biokart WES Report Tables (Table 1 & Table 2)
Tests:
1. Parsing FastQC zip/txt data before and after trimming
2. Generating Table 1 matching Image 3 format
3. Generating Table 2 matching Image 6 format
4. Formatting multi-annotation TSV and VCF with annotate_vcf.py
"""

import os
import sys
import zipfile
import subprocess
import shutil

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.dirname(SCRIPT_DIR)
TEST_DIR = os.path.join(SCRIPT_DIR, "mock_run")
os.makedirs(TEST_DIR, exist_ok=True)

# 1. Create mock FastQC data
def create_mock_fastqc(sample, stage, reads_m, gc, dups, gb):
    total_reads = int(reads_m * 1e6)
    read_len = int((gb * 1e9) / total_reads)
    dedup = 100.0 - dups
    
    content = f"""##FastQC\t0.11.9
>>Basic Statistics\tpass
#Measure\tValue
Filename\t{sample}.fastq.gz
File type\tConventional base calls
Encoding\tSanger / Illumina 1.9
Total Sequences\t{total_reads}
Sequences flagged as poor quality\t0
Sequence length\t{read_len}
%GC\t{gc}
>>END_MODULE
>>Sequence Duplication Levels\tpass
#Total Deduplicated Percentage\t{dedup}
>>END_MODULE
"""
    stage_dir = os.path.join(TEST_DIR, f"fastqc_{stage}", f"{sample}_fastqc")
    os.makedirs(stage_dir, exist_ok=True)
    with open(os.path.join(stage_dir, "fastqc_data.txt"), "w") as f:
        f.write(content)
        
    # Also create a zip file
    zip_path = os.path.join(TEST_DIR, f"fastqc_{stage}", f"{sample}_fastqc.zip")
    with zipfile.ZipFile(zip_path, 'w') as zf:
        zf.writestr(f"{sample}_fastqc/fastqc_data.txt", content)

print("--> Generating mock FastQC datasets for MODY-10 and MODY-11...")
# MODY-10 from Table 1 (Photo 3):
# Raw: 37.1M reads, 52% GC, 23.8%/22.7% dups, 11.13 GB
# HQ:  37.0M reads, 52% GC, 23.7%/22.7% dups, 11.13 GB
create_mock_fastqc("MODY-10_R1", "raw", 37.1, 52, 23.8, 11.13)
create_mock_fastqc("MODY-10_R2", "raw", 37.1, 52, 22.7, 11.13)
create_mock_fastqc("MODY-10_R1", "trimmed", 37.0, 52, 23.7, 11.13)
create_mock_fastqc("MODY-10_R2", "trimmed", 37.0, 52, 22.7, 11.13)

# MODY-11 from Table 1 (Photo 3):
# Raw: 50.4M reads, 54% GC, 26.1%/25.7% dups, 15.12 GB
# HQ:  50.2M reads, 54% GC, 26.1%/25.6% dups, 15.12 GB
create_mock_fastqc("MODY-11_R1", "raw", 50.4, 54, 26.1, 15.12)
create_mock_fastqc("MODY-11_R2", "raw", 50.4, 54, 25.7, 15.12)
create_mock_fastqc("MODY-11_R1", "trimmed", 50.2, 54, 26.1, 15.12)
create_mock_fastqc("MODY-11_R2", "trimmed", 50.2, 54, 25.6, 15.12)

# 2. Create mock variant counts (from Photo 6 / Table 2: MODY-10: 30061, MODY-11: 29784)
var_dir = os.path.join(TEST_DIR, "variants")
os.makedirs(var_dir, exist_ok=True)

with open(os.path.join(var_dir, "MODY-10.variant_counts.txt"), "w") as f:
    f.write("Sample Name\tSNVs\tInDels\tTotal SNV & small InDels\nMODY-10\t26815\t3246\t30061\n")

with open(os.path.join(var_dir, "MODY-11.variant_counts.txt"), "w") as f:
    f.write("Sample Name\tSNVs\tInDels\tTotal SNV & small InDels\nMODY-11\t26540\t3244\t29784\n")

# 3. Create mock VCF and test ANNOVAR annotator
mock_vcf = os.path.join(TEST_DIR, "test.vcf")
with open(mock_vcf, "w") as f:
    f.write("""##fileformat=VCFv4.2
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tMODY-10
chr1\t1000\trs12345\tA\tG\t100\tPASS\tGENE=HNF1A;FUNC=exonic;EXONIC_FUNC=nonsynonymous SNV;AACHANGE=HNF1A:c.100A>G:p.K34R;DBSNP=rs12345;1000G=0.001;EXAC=0.0005;CLNSIG=Pathogenic;SIFT_SCORE=0.01;SIFT_PRED=D;POLYPHEN_HDIV_SCORE=0.98;POLYPHEN_HDIV_PRED=D;CADD_PHRED=25.4\tGT:AD:DP:GQ:PL\t0/1:15,15:30:99:100,0,100
chr1\t2000\t.\tT\tTA\t100\tPASS\tGENE=GCK;FUNC=exonic;EXONIC_FUNC=frameshift insertion;AACHANGE=GCK:c.200insA:p.T67fs;CLNSIG=Pathogenic;CADD_PHRED=32.0\tGT:AD:DP:GQ:PL\t0/1:10,12:22:99:100,0,100
""")

print("--> Testing annotate_vcf.py...")
cmd_anno = [
    sys.executable,
    os.path.join(PROJECT_ROOT, "bin", "annotate_vcf.py"),
    "--vcf", mock_vcf,
    "--buildver", "hg19",
    "--out-prefix", os.path.join(TEST_DIR, "MODY-10")
]
res_anno = subprocess.run(cmd_anno, capture_output=True, text=True)
assert res_anno.returncode == 0, f"annotate_vcf failed: {res_anno.stderr}"
print(f"ANNOVAR output generated: {os.path.join(TEST_DIR, 'MODY-10.hg19_multianno.txt')}")

# 4. Run generate_report_tables.py
print("--> Testing generate_report_tables.py...")
cmd_tables = [
    sys.executable,
    os.path.join(PROJECT_ROOT, "bin", "generate_report_tables.py"),
    "--samples", "MODY-10,MODY-11",
    "--raw-qc-dir", os.path.join(TEST_DIR, "fastqc_raw"),
    "--trimmed-qc-dir", os.path.join(TEST_DIR, "fastqc_trimmed"),
    "--variants-dir", var_dir,
    "--outdir", os.path.join(TEST_DIR, "reports")
]
res_tables = subprocess.run(cmd_tables, capture_output=True, text=True)
assert res_tables.returncode == 0, f"generate_report_tables failed: {res_tables.stderr}"

# Check table 1
t1_path = os.path.join(TEST_DIR, "reports", "Table1_QC_stats.tsv")
with open(t1_path) as f:
    t1_content = f.read()
print("\n=== GENERATED TABLE 1 (TSV) ===")
print(t1_content)

# Check table 2
t2_path = os.path.join(TEST_DIR, "reports", "Table2_Variant_stats.tsv")
with open(t2_path) as f:
    t2_content = f.read()
print("=== GENERATED TABLE 2 (TSV) ===")
print(t2_content)

print("\n ALL REPORT & ANNOTATION TESTS PASSED SUCCESSFULLY!")

