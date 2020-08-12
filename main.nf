params.multiqc = "$baseDir/multiqc"
params.outdir = "$baseDir/results"

log.info """\
    This is a NEXTFLOW PIPELINE

    genome: ${params.genome}
    read directory: ${params.reads}
    outdir: ${params.outdir}

    """
.stripIndent()


ch_multiqc_config = file("$baseDir/assets/multiqc_config.yaml", checkIfExists: true)

process index {

    label 'bwa'

    input:
    path genome from params.genome

    // println "cpus: $task.cpus"
    output:
    file "*.{amb,ann,bwt,pac,sa}" into bwa_index_ch

    script:
    """
    bwa index $genome
    """
}

Channel
    .fromFilePairs(params.reads, checkIfExists: true)
    .map { it -> [ it[0].split('\\.')[0], it[1] ] }
    .into{ reads_ch1; reads_ch2 }


process align {

    label 'bwa'
    tag "$sample_id"
    publishDir "$params.outdir/bam"

    input:
    path genome from params.genome
    file '*' from bwa_index_ch
    tuple sample_id, path(reads) from reads_ch1

    output:
    tuple sample_id, "${sample_id}.RG.bam" into bam_ch

    script:
    """
    bwa mem -t $task.cpus $genome ${reads[0]} ${reads[1]} | samblaster --addMateTags --removeDups | samtools sort - | samtools view -Sb - > ${sample_id}.bam
    picard AddOrReplaceReadGroups -INPUT ${sample_id}.bam -OUTPUT ${sample_id}.RG.bam -VALIDATION_STRINGENCY LENIENT -RGID ${sample_id} -RGLB HUM -RGPL illumina -RGPU 1 -RGSM ${sample_id}
    samtools index ${sample_id}.RG.bam
    """
}


process fastqc {

    label 'fastqc'
    tag "FASTQC on $sample_id"

    input:
    tuple sample_id, path(reads) from reads_ch2

    output:
    path "fastqc_${sample_id}_logs" into fastqc_ch

    script:
    """
    mkdir fastqc_${sample_id}_logs
    fastqc -o fastqc_${sample_id}_logs -f fastq -q ${reads}
    """
}

process bamstats {

  label 'bamtools'
  tag "bamtools stats on $sample_id"
  echo true

  input:
  tuple sample_id, file(bam) from bam_ch

  output:
  path "${sample_id}.stats" into bamstats_ch

  script:
  """
  bamtools stats -in $bam > ${sample_id}.stats
  """

}


process multiqc {

    publishDir params.outdir, mode: 'copy'
    label 'fastqc'

    input:
    file (multiqc_config) from ch_multiqc_config
    path '*' from fastqc_ch.collect()
    path '*' from bamstats_ch.collect()

    output:
    path 'multiqc_report.html' into ch_multiqc_report

    script:
    """
    multiqc --config $multiqc_config .
    """
}

workflow.onComplete {

    c_green = "\033[0;32m";
    c_purple = "\033[0;35m";
    c_red =  "\033[0;31m";
    c_reset = "\033[0m";

    if (workflow.stats.ignoredCount > 0 && workflow.success) {
        log.info "-${c_purple}Warning, pipeline completed, but with errored process(es) ${c_reset}-"
        log.info "-${c_red}Number of ignored errored process(es) : ${workflow.stats.ignoredCount} ${c_reset}-"
        log.info "-${c_green}Number of successfully ran process(es) : ${workflow.stats.succeedCount} ${c_reset}-"
    }

    if (workflow.success) {
        log.info "-${c_purple}[nf-bwa-mem]${c_green} Pipeline completed successfully${c_reset}-"
    } else {
        log.info "-${c_purple}[nf-bwa-mem]${c_red} Pipeline completed with errors${c_reset}-"
    }

}
