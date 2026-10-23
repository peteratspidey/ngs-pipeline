#!/usr/bin/env nextflow

/**
 * Trimming Workflow
 * Supports fastp (default) and Trimmomatic for adapter/quality trimming
 */

nextflow.enable.dsl=2

// ============================================================================
// PROCESSES
// ============================================================================

process FASTP {
    tag "${sample_id}"
    label 'trim'
    publishDir "${params.outdir}/trimmed", mode: params.publish_mode, pattern: "*.trimmed.fastq.gz"
    publishDir "${params.outdir}/qc/fastp", mode: params.publish_mode, pattern: "*.{html,json}"
    
    input:
    tuple val(sample_id), path(r1), path(r2)
    
    output:
    tuple val(sample_id), path("${sample_id}_R1.trimmed.fastq.gz"), path("${sample_id}_R2.trimmed.fastq.gz"), path("${sample_id}_fastp.html"), path("${sample_id}_fastp.json"), emit: trimmed
    
    script:
    """
    fastp \
        --in1 ${r1} \
        --in2 ${r2} \
        --out1 ${sample_id}_R1.trimmed.fastq.gz \
        --out2 ${sample_id}_R2.trimmed.fastq.gz \
        --html ${sample_id}_fastp.html \
        --json ${sample_id}_fastp.json \
        --thread ${task.cpus} \
        --detect_adapter_for_pe \
        --trim_poly_g \
        --trim_poly_x \
        --qualified_quality_phred 20 \
        --unqualified_percent_limit 40 \
        --length_required 36 \
        --compression 6 \
        2>&1 | tee fastp.log
    """
}

process CUTADAPT {
    tag "${sample_id}"
    label 'trim'
    publishDir "${params.outdir}/trimmed", mode: params.publish_mode, pattern: "*.trimmed.fastq.gz"
    publishDir "${params.outdir}/qc/cutadapt", mode: params.publish_mode, pattern: "*.log"

    input:
    tuple val(sample_id), path(r1), path(r2)

    output:
    tuple val(sample_id), path("${sample_id}_R1.trimmed.fastq.gz"), path("${sample_id}_R2.trimmed.fastq.gz"), path("${sample_id}_cutadapt.log"), emit: trimmed

    script:
    def adapter_r1 = params.adapter_r1 ?: "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA"
    def adapter_r2 = params.adapter_r2 ?: "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"
    def quality_cutoff = params.quality_cutoff ?: 30
    def min_length = params.min_length ?: 36
    """
    cutadapt \
        -j ${task.cpus} \
        -q ${quality_cutoff},${quality_cutoff} \
        -m ${min_length} \
        -a ${adapter_r1} \
        -A ${adapter_r2} \
        -o ${sample_id}_R1.trimmed.fastq.gz \
        -p ${sample_id}_R2.trimmed.fastq.gz \
        ${r1} ${r2} > ${sample_id}_cutadapt.log 2>&1
    """
}

process TRIMMOMATIC {
    tag "${sample_id}"
    label 'trim'
    publishDir "${params.outdir}/trimmed", mode: params.publish_mode, pattern: "*.trimmed.fastq.gz"
    
    input:
    tuple val(sample_id), path(r1), path(r2)
    val adapters_file
    
    output:
    tuple val(sample_id), path("${sample_id}_R1.trimmed.fastq.gz"), path("${sample_id}_R2.trimmed.fastq.gz"), path("${sample_id}_R1.unpaired.fastq.gz"), path("${sample_id}_R2.unpaired.fastq.gz"), emit: trimmed
    
    script:
    """
    trimmomatic PE -threads ${task.cpus} -phred33 \
        ${r1} ${r2} \
        ${sample_id}_R1.trimmed.fastq.gz ${sample_id}_R1.unpaired.fastq.gz \
        ${sample_id}_R2.trimmed.fastq.gz ${sample_id}_R2.unpaired.fastq.gz \
        ILLUMINACLIP:${adapters_file}:2:30:10 \
        LEADING:3 TRAILING:3 \
        SLIDINGWINDOW:4:20 \
        MINLEN:36 \
        2>&1 | tee trimmomatic.log
    
    gzip -f ${sample_id}_R1.unpaired.fastq 2>/dev/null || true
    gzip -f ${sample_id}_R2.unpaired.fastq 2>/dev/null || true
    """
}

// ============================================================================
// WORKFLOW
// ============================================================================

workflow WORKFLOW_TRIM {
    take:
    reads_ch

    main:
    // Select trimmer based on parameter (cutadapt, fastp, trimmomatic)
    if (params.trimmer == 'cutadapt') {
        trimmed_ch = CUTADAPT(reads_ch)
        results_ch = trimmed_ch.map { sample_id, r1, r2, log_file ->
            def qc = [cutadapt_log: log_file]
            return [sample_id, r1, r2, qc]
        }
    } else if (params.trimmer == 'fastp') {
        trimmed_ch = FASTP(reads_ch)
        results_ch = trimmed_ch.map { sample_id, r1, r2, html, json ->
            def qc = [fastp_html: html, fastp_json: json]
            return [sample_id, r1, r2, qc]
        }
    } else {
        adapters_ch = Channel.fromPath("${workflow.projectDir}/resources/adapters/TruSeq3-PE.fa", checkIfExists: true)
        trimmed_ch = TRIMMOMATIC(reads_ch, adapters_ch.first())
        results_ch = trimmed_ch.map { sample_id, r1, r2, up1, up2 ->
            def qc = [unpaired_r1: up1, unpaired_r2: up2]
            return [sample_id, r1, r2, qc]
        }
    }

    emit:
    results_ch
}
