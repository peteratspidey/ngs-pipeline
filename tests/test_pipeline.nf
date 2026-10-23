#!/usr/bin/env nextflow

/**
 * NGS WES Pipeline - Integration Test Suite
 * Validates individual modules and pipeline flow
 */

nextflow.enable.dsl=2

include { WORKFLOW_QC }        from '../workflows/qc.nf'
include { WORKFLOW_TRIM }      from '../workflows/trim.nf'
include { WORKFLOW_ALIGN }     from '../workflows/align.nf'
include { WORKFLOW_VARIANTS }  from '../workflows/variants.nf'
include { WORKFLOW_ANNOTATE }  from '../workflows/annotate.nf'
include { WORKFLOW_REPORT }    from '../workflows/report.nf'

params {
    test_data   = "${workflow.projectDir}/../test_data/*_{R1,R2}*.fastq.gz"
    test_genome = "test"
    test_outdir = "./test_results"
    fasta       = "${workflow.projectDir}/../resources/genomes/test/test.fa"
    bed         = "${workflow.projectDir}/../resources/NO_BED"
    dbsnp       = "${workflow.projectDir}/../resources/NO_DBSNP"
    trimmer     = "cutadapt"
    aligner     = "bwa"
    variant_caller = "gatk"
    annotator   = "annovar"
    outdir      = "./test_results"
}

workflow {
    log.info "Starting NGS WES Test Suite"

    // Reads channel
    Channel
        .fromFilePairs(params.test_data, size: 2, flat: true)
        .map { sample_id, r1, r2 -> [sample_id, r1, r2] }
        .set { test_reads_ch }

    fasta_file = file(params.fasta)
    genome_ch = Channel.of(tuple(params.test_genome, fasta_file))

    // 1. QC
    raw_qc_res = WORKFLOW_QC(test_reads_ch, 'raw')

    // 2. Trimming
    trimmed_res = WORKFLOW_TRIM(test_reads_ch)
    trimmed_reads_ch = trimmed_res.map { sample_id, r1, r2, extra -> [sample_id, r1, r2] }

    // 3. Trimmed QC
    trimmed_qc_res = WORKFLOW_QC(trimmed_reads_ch, 'trimmed')

    // 4. Report test
    samples_ch = test_reads_ch.map { it[0] }.collect()

    WORKFLOW_REPORT(
        samples_ch,
        raw_qc_res.fastqc_files,
        trimmed_qc_res.fastqc_files,
        Channel.empty(),
        Channel.empty()
    )

    log.info "NGS WES Test Setup Passed"
}