#!/usr/bin/env nextflow

/**
 * Alignment Workflow
 * Supports BWA-MEM (default) and Bowtie2 for read alignment
 * Includes BAM sorting, marking duplicates, and BQSR
 */

nextflow.enable.dsl=2

// ============================================================================
// PROCESSES
// ============================================================================

process BWA_INDEX {
    tag "BWA Index: ${genome}"
    label 'index'
    
    input:
    path(fasta)
    val(genome)
    
    output:
    tuple val(genome), path("${genome}.*"), emit: index
    
    script:
    """
    bwa index -p ${genome} ${fasta} 2>&1 | tee bwa_index.log
    """
}

process BWA_MEM {
    tag "${sample_id}"
    label 'align'
    
    input:
    tuple val(sample_id), path(r1), path(r2)
    tuple val(genome), path(index_files)
    
    output:
    tuple val(sample_id), path("${sample_id}.bam"), emit: bam
    
    script:
    """
    bwa mem \
        -t ${task.cpus} \
        -M \
        -R "@RG\\tID:${sample_id}\\tSM:${sample_id}\\tPL:ILLUMINA" \
        ${genome} \
        ${r1} ${r2} | \
    samtools sort -@ ${task.cpus} -O BAM -o ${sample_id}.bam - 2>&1 | tee bwa_mem.log
    
    samtools index ${sample_id}.bam
    """
}

process BOWTIE2_INDEX {
    tag "Bowtie2 Index: ${genome}"
    label 'index'
    
    input:
    path(fasta)
    val(genome)
    
    output:
    tuple val(genome), path("${genome}.*"), emit: index
    
    script:
    """
    bowtie2-build ${fasta} ${genome} 2>&1 | tee bowtie2_build.log
    """
}

process BOWTIE2 {
    tag "${sample_id}"
    label 'align'
    
    input:
    tuple val(sample_id), path(r1), path(r2)
    tuple val(genome), path(index_files)
    
    output:
    tuple val(sample_id), path("${sample_id}.bam"), emit: bam
    
    script:
    """
    bowtie2 \
        --threads ${task.cpus} \
        -x ${genome} \
        -1 ${r1} -2 ${r2} \
        --rg-id ${sample_id} \
        --rg "SM:${sample_id}" \
        --rg "PL:ILLUMINA" | \
    samtools sort -@ ${task.cpus} -O BAM -o ${sample_id}.bam - 2>&1 | tee bowtie2.log
    
    samtools index ${sample_id}.bam
    """
}

process MARKDUPLICATES {
    tag "${sample_id}"
    label 'qc'
    publishDir "${params.outdir}/aligned", mode: params.publish_mode, pattern: "*.markdup.{bam,bam.bai,bai}"
    publishDir "${params.outdir}/reports/markduplicates", mode: params.publish_mode, pattern: "*.markdup_metrics.txt"
    
    input:
    tuple val(sample_id), path(bam)
    
    output:
    tuple val(sample_id), path("${sample_id}.markdup.bam"), path("${sample_id}.markdup.bam.bai"), path("${sample_id}.markdup_metrics.txt"), emit: markdup
    
    script:
    """
    gatk MarkDuplicates \
        -I ${bam} \
        -O ${sample_id}.markdup.bam \
        -M ${sample_id}.markdup_metrics.txt \
        --CREATE_INDEX true \
        --VALIDATION_STRINGENCY SILENT \
        2>&1 | tee markdup.log

    if [ -f "${sample_id}.markdup.bai" ] && [ ! -f "${sample_id}.markdup.bam.bai" ]; then
        ln -sf "${sample_id}.markdup.bai" "${sample_id}.markdup.bam.bai"
    fi
    """
}

process BQSR {
    tag "${sample_id}"
    label 'qc'
    publishDir "${params.outdir}/aligned", mode: params.publish_mode, pattern: "*.bqsr.{bam,bam.bai,bai}"
    publishDir "${params.outdir}/reports/bqsr", mode: params.publish_mode, pattern: "*.bqsr.table"
    
    input:
    tuple val(sample_id), path(bam), path(bai), path(metrics)
    path(dbsnp)
    path(known_indels)
    path(fasta)
    path(fasta_fai)
    path(fasta_dict)
    
    output:
    tuple val(sample_id), path("${sample_id}.bqsr.bam"), path("${sample_id}.bqsr.bam.bai"), path("${sample_id}.bqsr.table"), emit: bqsr
    
    script:
    """
    // Base Recalibration - Step 1: BaseRecalibrator
    gatk BaseRecalibrator \
        -I ${bam} \
        -R ${fasta} \
        --known-sites ${dbsnp} \
        --known-sites ${known_indels} \
        -O ${sample_id}.bqsr.table \
        2>&1 | tee bqsr1.log
    
    // Base Recalibration - Step 2: ApplyBQSR
    gatk ApplyBQSR \
        -R ${fasta} \
        -I ${bam} \
        --bqsr-recal-file ${sample_id}.bqsr.table \
        -O ${sample_id}.bqsr.bam \
        2>&1 | tee bqsr2.log
    
    if command -v samtools >/dev/null 2>&1; then
        samtools index ${sample_id}.bqsr.bam
    else
        gatk BuildBamIndex -I ${sample_id}.bqsr.bam
    fi

    if [ -f "${sample_id}.bqsr.bai" ] && [ ! -f "${sample_id}.bqsr.bam.bai" ]; then
        ln -sf "${sample_id}.bqsr.bai" "${sample_id}.bqsr.bam.bai"
    fi
    """
}

process ALIGNMENT_STATS {
    tag "${sample_id}"
    label 'qc'
    publishDir "${params.outdir}/reports/alignment_stats", mode: params.publish_mode
    
    input:
    tuple val(sample_id), path(bam), path(bai)
    path(fasta)
    path(bed_file)
    
    output:
    tuple val(sample_id), path("${sample_id}.flagstat.txt"), path("${sample_id}.idxstats.txt"), path("${sample_id}.coverage.txt"), emit: stats
    
    script:
    """
    // Flagstat
    samtools flagstat ${bam} > ${sample_id}.flagstat.txt 2>&1
    
    // Idxstats
    samtools idxstats ${bam} > ${sample_id}.idxstats.txt 2>&1
    
    // Exome Coverage (if BED provided)
    if [ -f "${bed_file}" ] && [ "${bed_file.name}" != "NO_BED" ]; then
        samtools depth -a -b ${bed_file} ${bam} | \
            awk '{sum+=\$3} END {if (NR>0) print "Mean coverage: " sum/NR; else print "Mean coverage: 0"}' \
            > ${sample_id}.coverage.txt 2>&1
    else
        samtools depth -a ${bam} | \
            awk '{sum+=\$3} END {if (NR>0) print "Mean coverage: " sum/NR; else print "Mean coverage: 0"}' \
            > ${sample_id}.coverage.txt 2>&1
    fi
    """
}

// ============================================================================
// WORKFLOW
// ============================================================================

workflow WORKFLOW_ALIGN {
    take:
    reads_ch
    genome_ch
    
    main:
    // Align reads using chosen aligner (default BWA-MEM)
    if (params.aligner == 'bowtie2') {
        index_ch = BOWTIE2_INDEX(genome_ch.map { it[1] }, genome_ch.map { it[0] })
        aligned_ch = BOWTIE2(reads_ch, index_ch)
    } else {
        index_ch = BWA_INDEX(genome_ch.map { it[1] }, genome_ch.map { it[0] })
        aligned_ch = BWA_MEM(reads_ch, index_ch)
    }
    
    // Mark duplicates with Picard / GATK
    markdup_ch = MARKDUPLICATES(aligned_ch)
    
    // BQSR (Base-recalibration) if enabled and known sites exist
    fasta_file = file(params.fasta)
    fai_file   = file("${params.fasta}.fai")
    dict_file  = file(params.fasta.replaceAll(/\.(fa|fasta)$/, '.dict'))
    bed_file   = params.bed ? file(params.bed) : file("${workflow.projectDir}/resources/NO_BED")

    if (params.run_bqsr && params.dbsnp && params.known_indels && file(params.dbsnp).exists() && file(params.known_indels).exists()) {
        dbsnp_file = file(params.dbsnp)
        indels_file = file(params.known_indels)
        bqsr_ch = BQSR(
            markdup_ch,
            dbsnp_file,
            indels_file,
            fasta_file,
            fai_file,
            dict_file
        )
        final_bam_ch = bqsr_ch.map { sample_id, bam, bai, table ->
            [sample_id, bam, bai, [bqsr_table: table]]
        }
    } else {
        final_bam_ch = markdup_ch.map { sample_id, bam, bai, metrics ->
            [sample_id, bam, bai, [markdup_metrics: metrics]]
        }
    }
    
    // Alignment stats
    stats_in_ch = final_bam_ch.map { sample_id, bam, bai, extra -> [sample_id, bam, bai] }
    stats_ch = ALIGNMENT_STATS(stats_in_ch, fasta_file, bed_file)
    
    // Combine BAMs, metrics, and stats
    results_ch = final_bam_ch
        .join(markdup_ch.map { sample_id, bam, bai, metrics -> [sample_id, metrics] })
        .join(stats_ch)
        .map { sample_id, bam, bai, extra, metrics, flagstat, idxstats, coverage ->
            def res = [:]
            res.bam = bam
            res.bai = bai
            res.markdup_metrics = metrics
            res.flagstat = flagstat
            res.idxstats = idxstats
            res.coverage = coverage
            if (extra.bqsr_table) res.bqsr_table = extra.bqsr_table
            return [sample_id, res]
        }

    emit:
    results_ch
}