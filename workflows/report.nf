#!/usr/bin/env nextflow

/**
 * Final Report Generation Workflow
 * Creates MultiQC report, Biokart Table 1 (QC stats), Table 2 (Variants called), and executive report
 */

nextflow.enable.dsl=2

// ============================================================================
// PROCESSES
// ============================================================================

process MULTIQC_FINAL {
    tag "Final MultiQC Report"
    label 'report'
    publishDir "${params.outdir}/reports/multiqc", mode: params.publish_mode
    
    input:
    path(qc_files)
    
    output:
    tuple path("multiqc_report.html"), path("multiqc_data"), emit: report
    
    script:
    """
    multiqc \
        --force \
        --outdir . \
        --title "NGS Exome Sequencing Analysis Report" \
        . 2>&1 | tee multiqc_final.log
    """
}

process GENERATE_SUMMARY_TABLES {
    tag "Generate Table 1 & Table 2"
    label 'report'
    publishDir "${params.outdir}/reports", mode: params.publish_mode
    
    input:
    val(sample_names)
    path(raw_qc_files)
    path(trimmed_qc_files)
    path(variant_files)
    
    output:
    path("Table1_QC_stats.tsv"), emit: table1_tsv
    path("Table1_QC_stats.html"), emit: table1_html
    path("Table2_Variant_stats.tsv"), emit: table2_tsv
    path("Table2_Variant_stats.html"), emit: table2_html
    
    script:
    def samples_arg = sample_names.join(',')
    """
    python3 ${workflow.projectDir}/bin/generate_report_tables.py \
        --samples "${samples_arg}" \
        --raw-qc-dir . \
        --trimmed-qc-dir . \
        --variants-dir . \
        --outdir .
    """
}

// ============================================================================
// WORKFLOW
// ============================================================================

workflow WORKFLOW_REPORT {
    take:
    samples_list
    raw_qc_files
    trimmed_qc_files
    aligned_files
    variant_files
    
    main:
    // 1. Generate Biokart Summary Tables (Table 1 & Table 2)
    tables_ch = GENERATE_SUMMARY_TABLES(
        samples_list,
        raw_qc_files.collect().ifEmpty([]),
        trimmed_qc_files.collect().ifEmpty([]),
        variant_files.collect().ifEmpty([])
    )
    
    // 2. Aggregate all QC outputs for MultiQC
    all_qc_ch = raw_qc_files.flatten()
        .mix(trimmed_qc_files.flatten())
        .mix(aligned_files.flatten())
        .mix(variant_files.flatten())
        .collect()
        
    multiqc_ch = MULTIQC_FINAL(all_qc_ch)
    
    emit:
    table1_tsv  = tables_ch.table1_tsv
    table1_html = tables_ch.table1_html
    table2_tsv  = tables_ch.table2_tsv
    table2_html = tables_ch.table2_html
    multiqc     = multiqc_ch.report
}