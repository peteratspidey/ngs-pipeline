#!/usr/bin/env nextflow

/**
 * NGS WES Pipeline - Next Generation Sequencing Analysis Pipeline
 * 
 * Whole Exome Sequencing (WES) SNV and InDel Analysis Pipeline:
 * - Quality control: FastQC & MultiQC (Raw and Trimmed)
 * - Filtering & Adapter Trimming: cutadapt (Q30 quality score cutoff)
 * - Alignment: BWA-MEM against hg19 / GRCh37
 * - Post-Alignment: Picard/GATK MarkDuplicates & Base Quality Score Recalibration (BQSR)
 * - Variant Calling: GATK HaplotypeCaller with exome capture BED target
 * - Annotation: ANNOVAR (refGene, 1000G, dbSNP, ExAC, ClinVar, SIFT, PolyPhen-2, CADD)
 * - Reporting: Table 1 (QC stats before & after trimming) & Table 2 (Variants called) + MultiQC
 */

nextflow.enable.dsl=2

// ============================================================================
// WORKFLOW INCLUDES
// ============================================================================

include { WORKFLOW_QC }           from './workflows/qc.nf'
include { WORKFLOW_TRIM }         from './workflows/trim.nf'
include { WORKFLOW_ALIGN }        from './workflows/align.nf'
include { WORKFLOW_VARIANTS }     from './workflows/variants.nf'
include { WORKFLOW_ANNOTATE }     from './workflows/annotate.nf'
include { WORKFLOW_REPORT }       from './workflows/report.nf'

// ============================================================================
// MAIN WORKFLOW
// ============================================================================

workflow {
    log.info """
    ===================================================================
     NGS WES ANALYSIS PIPELINE (hg19 / GRCh37)
    ===================================================================
     Genome          : ${params.genome}
     Trimmer         : ${params.trimmer} (Quality cutoff: Q${params.quality_cutoff ?: 30})
     Aligner         : ${params.aligner}
     Variant Caller  : ${params.variant_caller}
     Annotator       : ${params.annotator}
     Capture BED     : ${params.bed ?: 'None (whole genome mode)'}
     Output Directory: ${params.outdir}
    ===================================================================
    """

    // --------------------------------------------------------------------------
    // INPUT CHANNELS
    // --------------------------------------------------------------------------
    if (params.sample_sheet && file(params.sample_sheet).exists()) {
        Channel
            .fromPath(params.sample_sheet)
            .splitCsv(sep: ',', header: true)
            .map { row -> [row.sample_id, file(row.fastq_r1), file(row.fastq_r2)] }
            .set { grouped_reads_ch }
    } else {
        Channel
            .fromFilePairs(params.input, size: 2, flat: true)
            .ifEmpty { error "No paired-end FASTQ files found matching: ${params.input}" }
            .map { sample_id, r1, r2 -> [sample_id, r1, r2] }
            .set { grouped_reads_ch }
    }

    // Reference genome channel & validation
    if (!params.fasta) {
        def default_hg19 = "${workflow.projectDir}/resources/genomes/hg19/hg19.fa"
        def test_fa      = "${workflow.projectDir}/resources/genomes/test/test.fa"
        if (file(default_hg19).exists()) {
            params.fasta = default_hg19
        } else if (file(test_fa).exists()) {
            params.fasta = test_fa
        } else {
            error "\n" +
                  "========================================================================\n" +
                  " ERROR: Reference genome FASTA not found!\n" +
                  " \n" +
                  " Please do one of the following:\n" +
                  "   1. Download hg19 reference: ./scripts/download_references.sh\n" +
                  "   2. Specify your FASTA file:  --fasta /path/to/hg19.fa\n" +
                  "   3. Run quick self-contained test: ./run_pipeline.sh --test\n" +
                  "========================================================================"
        }
    }
    fasta_file = file(params.fasta)
    genome_ch = Channel.of(tuple(params.genome, fasta_file))

    // --------------------------------------------------------------------------
    // 1. QUALITY CONTROL (Raw reads)
    // --------------------------------------------------------------------------
    if (!params.skip_qc) {
        raw_qc_res = WORKFLOW_QC(grouped_reads_ch, 'raw')
        raw_qc_files = raw_qc_res.fastqc_files
    } else {
        raw_qc_files = Channel.empty()
    }

    // --------------------------------------------------------------------------
    // 2. TRIMMING (Cutadapt with Q30 cutoff & adapter removal)
    // --------------------------------------------------------------------------
    if (!params.skip_trim) {
        trimmed_res = WORKFLOW_TRIM(grouped_reads_ch)
        trimmed_reads_ch = trimmed_res.map { sample_id, r1, r2, extra -> [sample_id, r1, r2] }
        trimmed_qc_logs  = trimmed_res.map { sample_id, r1, r2, extra -> extra.values() }
    } else {
        trimmed_reads_ch = grouped_reads_ch
        trimmed_qc_logs  = Channel.empty()
    }

    // --------------------------------------------------------------------------
    // 3. QUALITY CONTROL (Trimmed HQ reads)
    // --------------------------------------------------------------------------
    if (!params.skip_qc) {
        trimmed_qc_res = WORKFLOW_QC(trimmed_reads_ch, 'trimmed')
        trimmed_qc_files = trimmed_qc_res.fastqc_files
    } else {
        trimmed_qc_files = Channel.empty()
    }

    // --------------------------------------------------------------------------
    // 4. ALIGNMENT (BWA-MEM + MarkDuplicates + BQSR)
    // --------------------------------------------------------------------------
    if (!params.skip_align) {
        aligned_ch = WORKFLOW_ALIGN(trimmed_reads_ch, genome_ch)
        aligned_qc_files = aligned_ch.map { sample_id, res ->
            [res.flagstat, res.idxstats, res.coverage, res.markdup_metrics]
        }
    } else {
        aligned_ch = Channel.empty()
        aligned_qc_files = Channel.empty()
    }

    // --------------------------------------------------------------------------
    // 5. VARIANT CALLING (GATK HaplotypeCaller with Exome BED)
    // --------------------------------------------------------------------------
    if (!params.skip_variants) {
        variants_ch = WORKFLOW_VARIANTS(aligned_ch, genome_ch)
        variant_qc_files = variants_ch.map { sample_id, res -> res.variant_counts }
    } else {
        variants_ch = Channel.empty()
        variant_qc_files = Channel.empty()
    }

    // --------------------------------------------------------------------------
    // 6. ANNOTATION (ANNOVAR multi-annotation tab-delimited output)
    // --------------------------------------------------------------------------
    if (!params.skip_annotate) {
        annotated_ch = WORKFLOW_ANNOTATE(variants_ch, genome_ch)
    } else {
        annotated_ch = Channel.empty()
    }

    // --------------------------------------------------------------------------
    // 7. FINAL REPORTS & NGS TABLES 1 & 2
    // --------------------------------------------------------------------------
    samples_list_ch = grouped_reads_ch.map { it[0] }.collect()

    WORKFLOW_REPORT(
        samples_list_ch,
        raw_qc_files,
        trimmed_qc_files,
        aligned_qc_files,
        variant_qc_files
    )
}
