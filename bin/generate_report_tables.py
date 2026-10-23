#!/usr/bin/env python3
"""
Generate NGS WES Pipeline QC and Variant Summary Tables
Matches:
- Table 1: QC stats of paired reads sequences before and after quality trimming
- Table 2: Number of variants called in the sample (Total SNV & small InDels)
"""

import sys
import os
import re
import glob
import json
import zipfile
import argparse

def parse_fastqc_data_file(content):
    stats = {
        'total_reads': 0,
        'gc_content': 0,
        'percent_duplicates': 0.0,
        'read_length': 0,
        'amount_gb': 0.0
    }
    
    in_dup_section = False
    dup_percent_val = 0.0
    
    for line in content.splitlines():
        line = line.strip()
        if not line:
            continue
            
        if line.startswith("Total Sequences"):
            parts = line.split('\t')
            if len(parts) > 1:
                stats['total_reads'] = int(parts[1])
        elif line.startswith("%GC"):
            parts = line.split('\t')
            if len(parts) > 1:
                stats['gc_content'] = int(round(float(parts[1])))
        elif line.startswith("Sequence length"):
            parts = line.split('\t')
            if len(parts) > 1:
                len_str = parts[1].split('-')[-1]
                stats['read_length'] = int(len_str)
        elif line.startswith("#Total Deduplicated Percentage"):
            parts = line.split('\t')
            if len(parts) > 1:
                dedup = float(parts[1])
                dup_percent_val = round(100.0 - dedup, 1)
                stats['percent_duplicates'] = dup_percent_val
        elif line.startswith(">>Sequence Duplication Levels"):
            in_dup_section = True
        elif in_dup_section and line.startswith(">>END_MODULE"):
            in_dup_section = False

    # Calculate data amount in GB (total reads * read length in bases / 1e9)
    if stats['total_reads'] > 0 and stats['read_length'] > 0:
        stats['amount_gb'] = round((stats['total_reads'] * stats['read_length']) / 1e9, 2)
        
    return stats

def extract_fastqc_stats(fastqc_path):
    if not os.path.exists(fastqc_path):
        return None
        
    if os.path.isdir(fastqc_path):
        data_txt = os.path.join(fastqc_path, 'fastqc_data.txt')
        if os.path.exists(data_txt):
            with open(data_txt, 'r', encoding='utf-8', errors='ignore') as f:
                return parse_fastqc_data_file(f.read())
                
    elif fastqc_path.endswith('.zip'):
        try:
            with zipfile.ZipFile(fastqc_path, 'r') as zf:
                for name in zf.namelist():
                    if name.endswith('fastqc_data.txt'):
                        with zf.open(name) as f:
                            content = f.read().decode('utf-8', errors='ignore')
                            return parse_fastqc_data_file(content)
        except Exception as e:
            print(f"Error reading zip {fastqc_path}: {e}", file=sys.stderr)
            
    return None

def find_fastqc_for_sample(qc_dir, sample_id, read_pair):
    patterns = [
        f"{sample_id}*{read_pair}*fastqc*",
        f"{sample_id}*_{read_pair}*",
        f"{sample_id}*"
    ]
    for pat in patterns:
        matches = glob.glob(os.path.join(qc_dir, "**", pat), recursive=True)
        if matches:
            for m in matches:
                if read_pair.lower() in os.path.basename(m).lower() or read_pair in os.path.basename(m):
                    return m
            return matches[0]
    return None

def generate_table1(samples, raw_qc_dir, trimmed_qc_dir, out_tsv, out_html):
    rows = []
    
    for sample in samples:
        for pair in ['R1', 'R2']:
            file_name = f"{sample}_{pair}"
            
            raw_path = find_fastqc_for_sample(raw_qc_dir, sample, pair)
            trim_path = find_fastqc_for_sample(trimmed_qc_dir, sample, pair)
            
            raw_stats = extract_fastqc_stats(raw_path) if raw_path else None
            trim_stats = extract_fastqc_stats(trim_path) if trim_path else None
            
            r_reads = f"{raw_stats['total_reads'] / 1e6:.1f}" if raw_stats else "N/A"
            r_gc = f"{raw_stats['gc_content']}" if raw_stats else "N/A"
            r_dups = f"{raw_stats['percent_duplicates']:.1f}" if raw_stats else "N/A"
            r_gb = f"{raw_stats['amount_gb']:.2f}" if raw_stats else "N/A"
            
            t_reads = f"{trim_stats['total_reads'] / 1e6:.1f}" if trim_stats else (r_reads if r_reads != "N/A" else "N/A")
            t_gc = f"{trim_stats['gc_content']}" if trim_stats else (r_gc if r_gc != "N/A" else "N/A")
            t_dups = f"{trim_stats['percent_duplicates']:.1f}" if trim_stats else (r_dups if r_dups != "N/A" else "N/A")
            t_gb = f"{trim_stats['amount_gb']:.2f}" if trim_stats else (r_gb if r_gb != "N/A" else "N/A")
            
            rows.append({
                'Sample Name': sample if pair == 'R1' else '',
                'Sample Key': sample,
                'Files': file_name,
                'Raw Reads (M)': r_reads,
                'Raw % GC': r_gc,
                'Raw % Dups': r_dups,
                'Raw GB': r_gb,
                'HQ Reads (M)': t_reads,
                'HQ % GC': t_gc,
                'HQ % Dups': t_dups,
                'HQ GB': t_gb
            })

    # Write TSV
    with open(out_tsv, 'w', encoding='utf-8') as f:
        f.write("Sample Name\tFiles\tRaw_Reads_M\tRaw_GC_pct\tRaw_Dups_pct\tRaw_GB\tHQ_Reads_M\tHQ_GC_pct\tHQ_Dups_pct\tHQ_GB\n")
        for r in rows:
            sample_col = r['Sample Key']
            f.write(f"{sample_col}\t{r['Files']}\t{r['Raw Reads (M)']}\t{r['Raw % GC']}\t{r['Raw % Dups']}\t{r['Raw GB']}\t{r['HQ Reads (M)']}\t{r['HQ % GC']}\t{r['HQ % Dups']}\t{r['HQ GB']}\n")

    # Write HTML
    html_lines = [
        "<!DOCTYPE html>",
        "<html><head><meta charset='utf-8'>",
        "<title>Table 1: QC Stats of Paired Reads Before and After Quality Trimming</title>",
        "<style>",
        "table { border-collapse: collapse; width: 100%; font-family: sans-serif; font-size: 13px; }",
        "th, td { border: 1px solid #ddd; padding: 8px; text-align: center; }",
        "th { background-color: #f2f4f8; font-weight: bold; }",
        ".header-group { background-color: #e5ebf5; font-size: 14px; }",
        "tr:nth-child(even) { background-color: #fafbfc; }",
        "</style></head><body>",
        "<h3>Table 1: Depicts the QC stats of paired reads sequences before and after quality trimming</h3>",
        "<table>",
        "<thead>",
        "<tr>",
        "<th rowspan='2'>Sample Name</th>",
        "<th rowspan='2'>Files</th>",
        "<th colspan='4' class='header-group'>QC stats of raw fastq reads</th>",
        "<th colspan='4' class='header-group'>QC stats of HQ reads</th>",
        "</tr>",
        "<tr>",
        "<th>No of reads (millions)</th><th>% GC</th><th>% Dups</th><th>Amount of data (GB)</th>",
        "<th>No of reads (millions)</th><th>% GC</th><th>% Dups</th><th>Amount of data (GB)</th>",
        "</tr>",
        "</thead><tbody>"
    ]
    
    for r in rows:
        html_lines.append("<tr>")
        html_lines.append(f"<td><b>{r['Sample Name']}</b></td>")
        html_lines.append(f"<td>{r['Files']}</td>")
        html_lines.append(f"<td>{r['Raw Reads (M)']}</td>")
        html_lines.append(f"<td>{r['Raw % GC']}</td>")
        html_lines.append(f"<td>{r['Raw % Dups']}</td>")
        html_lines.append(f"<td>{r['Raw GB']}</td>")
        html_lines.append(f"<td>{r['HQ Reads (M)']}</td>")
        html_lines.append(f"<td>{r['HQ % GC']}</td>")
        html_lines.append(f"<td>{r['HQ % Dups']}</td>")
        html_lines.append(f"<td>{r['HQ GB']}</td>")
        html_lines.append("</tr>")
        
    html_lines.extend(["</tbody></table></body></html>"])
    
    with open(out_html, 'w', encoding='utf-8') as f:
        f.write("\n".join(html_lines))
    print(f"[Report] Generated Table 1: {out_tsv} and {out_html}")

def generate_table2(samples, variants_dir, out_tsv, out_html):
    rows = []
    
    for sample in samples:
        total = 0
        snvs = 0
        indels = 0
        
        # Check for count file
        count_files = glob.glob(os.path.join(variants_dir, f"{sample}*.variant_counts.txt"))
        if count_files:
            with open(count_files[0], 'r', encoding='utf-8') as f:
                for line in f:
                    if line.startswith(sample):
                        parts = line.strip().split('\t')
                        if len(parts) >= 4:
                            snvs = int(parts[1])
                            indels = int(parts[2])
                            total = int(parts[3])
        else:
            # Check vcf directly
            vcf_files = glob.glob(os.path.join(variants_dir, f"{sample}*.vcf*"))
            if vcf_files:
                import gzip
                vcf_f = vcf_files[0]
                opener = gzip.open if vcf_f.endswith('.gz') else open
                try:
                    with opener(vcf_f, 'rt', encoding='utf-8', errors='ignore') as f:
                        for line in f:
                            if not line.startswith('#'):
                                parts = line.split('\t')
                                if len(parts) > 4:
                                    ref, alt = parts[3], parts[4]
                                    total += 1
                                    if len(ref) == 1 and len(alt) == 1:
                                        snvs += 1
                                    else:
                                        indels += 1
                except Exception as e:
                    print(f"Error parsing VCF {vcf_f}: {e}", file=sys.stderr)
                    
        rows.append({
            'Sample Name': sample,
            'SNVs': snvs,
            'InDels': indels,
            'Total SNV & small InDels': total if total > 0 else (snvs + indels)
        })

    # Write TSV
    with open(out_tsv, 'w', encoding='utf-8') as f:
        f.write("Sample Name\tSNVs\tInDels\tTotal SNV & small InDels\n")
        for r in rows:
            f.write(f"{r['Sample Name']}\t{r['SNVs']}\t{r['InDels']}\t{r['Total SNV & small InDels']}\n")

    # Write HTML
    html_lines = [
        "<!DOCTYPE html>",
        "<html><head><meta charset='utf-8'>",
        "<title>Table 2: Number of Variants Called in the Sample</title>",
        "<style>",
        "table { border-collapse: collapse; width: 60%; font-family: sans-serif; font-size: 13px; margin: 20px 0; }",
        "th, td { border: 1px solid #ddd; padding: 10px; text-align: center; }",
        "th { background-color: #f2f4f8; font-weight: bold; }",
        "tr:nth-child(even) { background-color: #fafbfc; }",
        "</style></head><body>",
        "<h3>Table 2: Depicts the number of variants called in the sample</h3>",
        "<table>",
        "<thead><tr><th>Sample Name</th><th>SNVs</th><th>small InDels</th><th>Total SNV & small InDels</th></tr></thead>",
        "<tbody>"
    ]
    
    for r in rows:
        html_lines.append(f"<tr><td><b>{r['Sample Name']}</b></td><td>{r['SNVs']}</td><td>{r['InDels']}</td><td><b>{r['Total SNV & small InDels']}</b></td></tr>")
        
    html_lines.extend(["</tbody></table></body></html>"])
    
    with open(out_html, 'w', encoding='utf-8') as f:
        f.write("\n".join(html_lines))
    print(f"[Report] Generated Table 2: {out_tsv} and {out_html}")

def main():
    parser = argparse.ArgumentParser(description="Generate NGS NGS QC and Variant Tables")
    parser.add_argument("--samples", required=True, help="Comma-separated sample names or sample list")
    parser.add_argument("--raw-qc-dir", default="./results/qc/fastqc_raw", help="Raw FastQC directory")
    parser.add_argument("--trimmed-qc-dir", default="./results/qc/fastqc_trimmed", help="Trimmed FastQC directory")
    parser.add_argument("--variants-dir", default="./results/variants", help="Variants directory")
    parser.add_argument("--outdir", default="./results/reports", help="Output directory for reports")
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    sample_list = [s.strip() for s in args.samples.split(',') if s.strip()]

    t1_tsv = os.path.join(args.outdir, "Table1_QC_stats.tsv")
    t1_html = os.path.join(args.outdir, "Table1_QC_stats.html")
    generate_table1(sample_list, args.raw_qc_dir, args.trimmed_qc_dir, t1_tsv, t1_html)

    t2_tsv = os.path.join(args.outdir, "Table2_Variant_stats.tsv")
    t2_html = os.path.join(args.outdir, "Table2_Variant_stats.html")
    generate_table2(sample_list, args.variants_dir, t2_tsv, t2_html)

if __name__ == "__main__":
    main()

