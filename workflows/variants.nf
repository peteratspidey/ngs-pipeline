#!/usr/bin/env nextflow

/**
 * Variant Calling Workflow
 * Supports GATK HaplotypeCaller, FreeBayes, and joint genotyping
 */

nextflow.enable.dsl=2

// ============================================================================
// PROCESSES
// ============================================================================

process GATK_HAPLOTYPECALLER {
    tag "${sample_id}"
    label 'variants'
    publishDir "${params.outdir}/variants/raw", mode: params.publish_mode, pattern: "*.raw.vcf.gz*"
    
    input:
    tuple val(sample_id), path(bam), path(bai)
    path(fasta)
    path(fasta_fai)
    path(fasta_dict)
    path(dbsnp)
    path(bed_file)
    
    output:
    tuple val(sample_id), path("${sample_id}.raw.vcf.gz"), path("${sample_id}.raw.vcf.gz.tbi"), emit: vcf
    
    script:
    def dbsnp_opt = (dbsnp && dbsnp.name != "NO_DBSNP" && dbsnp.exists()) ? "--dbsnp ${dbsnp}" : ""
    def bed_opt = (bed_file && bed_file.name != "NO_BED" && bed_file.exists()) ? "-L ${bed_file}" : ""
    """
    gatk HaplotypeCaller \
        -R ${fasta} \
        -I ${bam} \
        -O ${sample_id}.raw.vcf.gz \
        ${dbsnp_opt} \
        ${bed_opt} \
        --native-pair-hmm-threads ${task.cpus} \
        2>&1 | tee haplotypecaller.log
    
    gatk IndexFeatureFile -I ${sample_id}.raw.vcf.gz 2>&1 | tee index.log
    """
}

process GATK_FILTER_VARIANTS {
    tag "${sample_id}"
    label 'variants'
    publishDir "${params.outdir}/variants", mode: params.publish_mode, pattern: "*.filtered.vcf.gz*"
    publishDir "${params.outdir}/reports/variants", mode: params.publish_mode, pattern: "*.stats.txt"
    
    input:
    tuple val(sample_id), path(raw_vcf), path(raw_tbi)
    path(fasta)
    path(fasta_fai)
    path(fasta_dict)
    
    output:
    tuple val(sample_id), path("${sample_id}.filtered.vcf.gz"), path("${sample_id}.filtered.vcf.gz.tbi"), path("${sample_id}.variant_counts.txt"), emit: filtered
    
    script:
    """
    # 1. Select SNPs
    gatk SelectVariants \
        -R ${fasta} \
        -V ${raw_vcf} \
        --select-type-to-include SNP \
        -O ${sample_id}.snps.raw.vcf.gz 2>&1 | tee select_snp.log

    # 2. Filter SNPs
    gatk VariantFiltration \
        -R ${fasta} \
        -V ${sample_id}.snps.raw.vcf.gz \
        -O ${sample_id}.snps.filtered.vcf.gz \
        --filter-name "QD_filter" --filter-expression "QD < 2.0" \
        --filter-name "FS_filter" --filter-expression "FS > 60.0" \
        --filter-name "MQ_filter" --filter-expression "MQ < 40.0" \
        --filter-name "SOR_filter" --filter-expression "SOR > 3.0" \
        --filter-name "MQRankSum_filter" --filter-expression "MQRankSum < -12.5" \
        --filter-name "ReadPosRankSum_filter" --filter-expression "ReadPosRankSum < -8.0" \
        2>&1 | tee filter_snp.log

    # 3. Select InDels
    gatk SelectVariants \
        -R ${fasta} \
        -V ${raw_vcf} \
        --select-type-to-include INDEL \
        -O ${sample_id}.indels.raw.vcf.gz 2>&1 | tee select_indel.log

    # 4. Filter InDels
    gatk VariantFiltration \
        -R ${fasta} \
        -V ${sample_id}.indels.raw.vcf.gz \
        -O ${sample_id}.indels.filtered.vcf.gz \
        --filter-name "QD_filter" --filter-expression "QD < 2.0" \
        --filter-name "FS_filter" --filter-expression "FS > 200.0" \
        --filter-name "ReadPosRankSum_filter" --filter-expression "ReadPosRankSum < -20.0" \
        2>&1 | tee filter_indel.log

    # 5. Merge filtered variants
    gatk MergeVcfs \
        -I ${sample_id}.snps.filtered.vcf.gz \
        -I ${sample_id}.indels.filtered.vcf.gz \
        -O ${sample_id}.filtered.vcf.gz 2>&1 | tee merge.log

    # 6. Extract counts for Table 2
    SNP_COUNT=\$(bcftools view -v snps -f PASS ${sample_id}.filtered.vcf.gz 2>/dev/null | grep -v "^#" | wc -l || echo "0")
    INDEL_COUNT=\$(bcftools view -v indels -f PASS ${sample_id}.filtered.vcf.gz 2>/dev/null | grep -v "^#" | wc -l || echo "0")
    TOTAL_COUNT=\$((SNP_COUNT + INDEL_COUNT))

    cat > ${sample_id}.variant_counts.txt <<EOF
Sample Name\tSNVs\tInDels\tTotal SNV & small InDels
${sample_id}\t\${SNP_COUNT}\t\${INDEL_COUNT}\t\${TOTAL_COUNT}
EOF

    bcftools stats ${sample_id}.filtered.vcf.gz > ${sample_id}.stats.txt 2>&1 || true
    """
}

// ============================================================================
// WORKFLOW
// ============================================================================

workflow WORKFLOW_VARIANTS {
    take:
    aligned_ch
    genome_ch

    main:
    fasta_file = file(params.fasta)
    fai_file   = file("${params.fasta}.fai")
    dict_file  = file(params.fasta.replaceAll(/\.(fa|fasta)$/, '.dict'))
    dbsnp_file = (params.dbsnp && file(params.dbsnp).exists()) ? file(params.dbsnp) : file("${workflow.projectDir}/resources/NO_DBSNP")
    bed_file   = (params.bed && file(params.bed).exists()) ? file(params.bed) : file("${workflow.projectDir}/resources/NO_BED")

    in_ch = aligned_ch.map { sample_id, results ->
        [sample_id, results.bam, results.bai]
    }

    raw_vcf_ch = GATK_HAPLOTYPECALLER(
        in_ch,
        fasta_file,
        fai_file,
        dict_file,
        dbsnp_file,
        bed_file
    )

    filtered_ch = GATK_FILTER_VARIANTS(
        raw_vcf_ch,
        fasta_file,
        fai_file,
        dict_file
    )

    results_ch = filtered_ch.map { sample_id, vcf, tbi, counts ->
        def res = [:]
        res.vcf = vcf
        res.tbi = tbi
        res.variant_counts = counts
        return [sample_id, res]
    }

    emit:
    results_ch
}
