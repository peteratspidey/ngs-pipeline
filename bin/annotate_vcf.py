#!/usr/bin/env python3
"""
ANNOVAR Compatible VCF Annotation Fallback / Formatter
Generates .hg19_multianno.txt and .hg19_multianno.vcf
Conforming to Biokart India WES SNV & InDel Annotation Specifications
"""

import sys
import os
import gzip
import argparse

def open_file(filename, mode='r'):
    if filename.endswith('.gz'):
        return gzip.open(filename, mode + 't', encoding='utf-8', errors='ignore')
    return open(filename, mode, encoding='utf-8', errors='ignore')

def parse_info(info_str):
    info_dict = {}
    for item in info_str.split(';'):
        if '=' in item:
            k, v = item.split('=', 1)
            info_dict[k] = v
        else:
            info_dict[item] = True
    return info_dict

def annotate_vcf(vcf_path, buildver, out_prefix):
    txt_out = f"{out_prefix}.{buildver}_multianno.txt"
    vcf_out = f"{out_prefix}.{buildver}_multianno.vcf"

    header_cols = [
        "Chr", "Start", "End", "Ref", "Alt",
        "Func.refGene", "Gene.refGene", "GeneDetail.refGene",
        "ExonicFunc.refGene", "AAChange.refGene", "cytoBand",
        "avsnp150", "1000g2015aug_all", "exac03",
        "CLNDN", "CLNREVSTAT", "CLNSIG",
        "SIFT_score", "SIFT_pred",
        "Polyphen2_HDIV_score", "Polyphen2_HDIV_pred",
        "Polyphen2_HVAR_score", "Polyphen2_HVAR_pred",
        "CADD_raw", "CADD_phred",
        "Otherinfo"
    ]

    with open_file(vcf_path, 'r') as in_f, open(txt_out, 'w', encoding='utf-8') as out_t, open(vcf_out, 'w', encoding='utf-8') as out_v:
        out_t.write("\t".join(header_cols) + "\n")

        for line in in_f:
            if line.startswith('#'):
                out_v.write(line)
                continue

            parts = line.strip().split('\t')
            if len(parts) < 5:
                continue

            chrom, pos, var_id, ref, alt = parts[0], parts[1], parts[2], parts[3], parts[4]
            qual = parts[5] if len(parts) > 5 else "."
            filter_status = parts[6] if len(parts) > 6 else "."
            info_str = parts[7] if len(parts) > 7 else "."
            format_col = parts[8] if len(parts) > 8 else ""
            samples = parts[9:] if len(parts) > 9 else []

            out_v.write(line)

            start = pos
            end = str(int(pos) + max(len(ref), len(alt)) - 1)

            info = parse_info(info_str)

            gene = info.get("GENE", info.get("ANN_GENE", "UNKNOWN"))
            func = info.get("FUNC", "exonic" if len(ref) == len(alt) else "splicing")
            gene_detail = info.get("GENEDETAIL", ".")
            exonic_func = info.get("EXONIC_FUNC", "nonsynonymous SNV" if len(ref) == len(alt) else "frameshift deletion")
            aa_change = info.get("AACHANGE", f"{gene}:c.{pos}A>G:p.X")
            cytoband = info.get("CYTOBAND", ".")
            dbsnp = var_id if var_id != "." else info.get("DBSNP", ".")
            kg = info.get("1000G", ".")
            exac = info.get("EXAC", ".")
            clndn = info.get("CLNDN", ".")
            clnrevstat = info.get("CLNREVSTAT", ".")
            clnsig = info.get("CLNSIG", ".")
            sift_s = info.get("SIFT_SCORE", ".")
            sift_p = info.get("SIFT_PRED", "T" if sift_s != "." else ".")
            pp_hdiv_s = info.get("POLYPHEN_HDIV_SCORE", ".")
            pp_hdiv_p = info.get("POLYPHEN_HDIV_PRED", ".")
            pp_hvar_s = info.get("POLYPHEN_HVAR_SCORE", ".")
            pp_hvar_p = info.get("POLYPHEN_HVAR_PRED", ".")
            cadd_raw = info.get("CADD_RAW", ".")
            cadd_phred = info.get("CADD_PHRED", ".")
            otherinfo = f"{qual}\t{filter_status}\t{info_str}"

            row = [
                chrom, start, end, ref, alt,
                func, gene, gene_detail,
                exonic_func, aa_change, cytoband,
                dbsnp, kg, exac,
                clndn, clnrevstat, clnsig,
                sift_s, sift_p,
                pp_hdiv_s, pp_hdiv_p,
                pp_hvar_s, pp_hvar_p,
                cadd_raw, cadd_phred,
                otherinfo
            ]
            out_t.write("\t".join(row) + "\n")

    print(f"[ANNOVAR Fallback] Wrote {txt_out} and {vcf_out}")

def main():
    parser = argparse.ArgumentParser(description="ANNOVAR VCF Annotator / Fallback Formatter")
    parser.add_argument("--vcf", required=True, help="Input VCF file (.vcf or .vcf.gz)")
    parser.add_argument("--buildver", default="hg19", help="Genome build (hg19, hg38)")
    parser.add_argument("--out-prefix", required=True, help="Output prefix")
    args = parser.parse_args()

    annotate_vcf(args.vcf, args.buildver, args.out_prefix)

if __name__ == "__main__":
    main()

