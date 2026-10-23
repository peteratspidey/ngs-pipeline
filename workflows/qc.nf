#!/usr/bin/env nextflow

/**
 * Quality Control Workflow
 * Runs FastQC and MultiQC on input reads
 */

nextflow.enable.dsl=2

// ============================================================================
// PROCESSES
// ============================================================================

process FASTQC {
    tag "${sample_id} (${stage})"
    label 'qc'
    publishDir path: { "${params.outdir}/qc/fastqc_${stage}" }, mode: params.publish_mode
    
    input:
    tuple val(sample_id), path(reads), val(stage)
    
    output:
    tuple val(sample_id), path("*_fastqc.{zip,html}"), val(stage), emit: fastqc
    
    script:
    """
    fastqc \
        --threads ${task.cpus} \
        --outdir . \
        ${reads} \
        2>&1 | tee fastqc.log
    """
}

process MULTIQC {
    tag "MultiQC (${stage})"
    label 'qc'
    publishDir path: { "${params.outdir}/qc/multiqc_${stage}" }, mode: params.publish_mode
    
    input:
    path(qc_files)
    val(stage)
    
    output:
    tuple val("multiqc_${stage}"), path("multiqc_report.html"), path("multiqc_data"), emit: multiqc
    
    script:
    """
    multiqc \
        --force \
        --outdir . \
        --title "NGS Pipeline QC Report (${stage})" \
        . \
        2>&1 | tee multiqc.log
    """
}

// ============================================================================
// WORKFLOW
// ============================================================================

workflow WORKFLOW_QC {
    take:
    reads_ch
    stage

    main:
    staged_reads_ch = reads_ch.map { sample_id, r1, r2 ->
        [sample_id, [r1, r2], stage]
    }
    fastqc_ch = FASTQC(staged_reads_ch)
    
    results_ch = reads_ch
        .combine(fastqc_ch, by: 0)
        .map { sample_id, r1, r2, fastqc_files, qc_stage ->
            [sample_id, r1, r2, [fastqc: fastqc_files, stage: qc_stage]]
        }

    emit:
    results_ch
    fastqc_files = fastqc_ch.map { it[1] }
}
