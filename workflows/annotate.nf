#!/usr/bin/env nextflow

/**
 * Variant Annotation Workflow
 * Supports ANNOVAR (NGS WES specification), SnpEff, and VEP
 */

nextflow.enable.dsl=2

// ============================================================================
// PROCESSES
// ============================================================================

process ANNOVAR {
    tag "${sample_id}"
    label 'annotate'
    publishDir "${params.outdir}/annotation", mode: params.publish_mode
    
    input:
    tuple val(sample_id), path(vcf), path(tbi)
    val(genome)
    path(annovar_dir)
    
    output:
    tuple val(sample_id), path("${sample_id}.${genome}_multianno.txt"), path("${sample_id}.${genome}_multianno.vcf"), emit: annotated
    
    script:
    def protocols = params.annovar_protocols ?: "refGene,avsnp150,1000g2015aug_all,exac03,clinvar_20221231,dbnsfp42a"
    def operations = params.annovar_operations ?: "g,f,f,f,f,f"
    """
    if command -v table_annovar.pl >/dev/null 2>&1; then
        table_annovar.pl ${vcf} ${annovar_dir} \
            -buildver ${genome} \
            -out ${sample_id} \
            -remove \
            -protocol ${protocols} \
            -operation ${operations} \
            -nastring . \
            -vcfinput \
            -polish 2>&1 | tee annovar.log
    elif [ -f "${annovar_dir}/table_annovar.pl" ]; then
        perl ${annovar_dir}/table_annovar.pl ${vcf} ${annovar_dir} \
            -buildver ${genome} \
            -out ${sample_id} \
            -remove \
            -protocol ${protocols} \
            -operation ${operations} \
            -nastring . \
            -vcfinput \
            -polish 2>&1 | tee annovar.log
    else
        python3 ${workflow.projectDir}/bin/annotate_vcf.py \
            --vcf ${vcf} \
            --buildver ${genome} \
            --out-prefix ${sample_id}
    fi
    """
}

process SNPEFF_BUILD {
    tag "SnpEff Build: ${genome}"
    label 'annotate'
    
    input:
    val(genome)
    path(fasta)
    path(gtf)
    
    output:
    path("snpEff.config"), emit: config
    
    script:
    """
    cat > snpEff.config <<EOF
    ${genome}.genome : ${genome}
    ${genome}.fasta : ${fasta}
    ${genome}.gtf : ${gtf}
EOF
    snpEff build -c snpEff.config -gtf22 -v ${genome} 2>&1 | tee snpeff_build.log
    """
}

process SNPEFF {
    tag "${sample_id}"
    label 'annotate'
    publishDir "${params.outdir}/annotation", mode: params.publish_mode
    
    input:
    tuple val(sample_id), path(vcf), path(tbi)
    val(genome)
    path(snpeff_config)
    
    output:
    tuple val(sample_id), path("${sample_id}.snpeff.vcf.gz"), path("${sample_id}.snpeff.vcf.gz.tbi"), path("${sample_id}.snpeff_summary.html"), emit: annotated
    
    script:
    """
    cp ${snpeff_config} snpEff.config
    snpEff -c snpEff.config \
        -v ${genome} \
        -stats ${sample_id}.snpeff_summary.html \
        -i vcf -o vcf \
        ${vcf} > ${sample_id}.snpeff.vcf 2>&1 | tee snpeff.log
    
    bgzip -c ${sample_id}.snpeff.vcf > ${sample_id}.snpeff.vcf.gz
    tabix -p vcf ${sample_id}.snpeff.vcf.gz
    """
}

process VEP {
    tag "${sample_id}"
    label 'annotate'
    publishDir "${params.outdir}/annotation", mode: params.publish_mode
    
    input:
    tuple val(sample_id), path(vcf), path(tbi)
    path(vep_cache)
    path(fasta)
    
    output:
    tuple val(sample_id), path("${sample_id}.vep.vcf.gz"), path("${sample_id}.vep.vcf.gz.tbi"), path("${sample_id}.vep_summary.html"), emit: annotated
    
    script:
    """
    vep \
        --cache \
        --dir ${vep_cache} \
        --fasta ${fasta} \
        --input_file ${vcf} \
        --output_file ${sample_id}.vep.vcf \
        --format vcf \
        --vcf \
        --symbol \
        --fork ${task.cpus} \
        --force_overwrite \
        --stats_file ${sample_id}.vep_summary.html \
        2>&1 | tee vep.log
    
    bgzip -c ${sample_id}.vep.vcf > ${sample_id}.vep.vcf.gz
    tabix -p vcf ${sample_id}.vep.vcf.gz
    """
}

// ============================================================================
// WORKFLOW
// ============================================================================

workflow WORKFLOW_ANNOTATE {
    take:
    variants_ch
    genome_ch

    main:
    in_vcf_ch = variants_ch.map { sample_id, results ->
        [sample_id, results.vcf, results.tbi]
    }

    if (params.annotator == 'annovar') {
        annovar_db = (params.annovar_dir && file(params.annovar_dir).exists()) ? file(params.annovar_dir) : file("${workflow.projectDir}/resources")
        annotated_ch = ANNOVAR(in_vcf_ch, params.genome, annovar_db)
        results_ch = annotated_ch.map { sample_id, multianno_txt, multianno_vcf ->
            def res = [:]
            res.multianno_txt = multianno_txt
            res.multianno_vcf = multianno_vcf
            return [sample_id, res]
        }
    } else if (params.annotator == 'snpeff') {
        fasta_file = file(params.fasta)
        gtf_file = file(params.gtf)
        snpeff_config_ch = SNPEFF_BUILD(params.genome, fasta_file, gtf_file)
        annotated_ch = SNPEFF(in_vcf_ch, params.genome, snpeff_config_ch)
        results_ch = annotated_ch.map { sample_id, vcf, tbi, summary ->
            def res = [:]
            res.snpeff_vcf = vcf
            res.snpeff_summary = summary
            return [sample_id, res]
        }
    } else {
        vep_cache = file(params.vep_cache)
        fasta_file = file(params.fasta)
        annotated_ch = VEP(in_vcf_ch, vep_cache, fasta_file)
        results_ch = annotated_ch.map { sample_id, vcf, tbi, summary ->
            def res = [:]
            res.vep_vcf = vcf
            res.vep_summary = summary
            return [sample_id, res]
        }
    }

    emit:
    results_ch
}