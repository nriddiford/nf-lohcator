params.reads = "$baseDir/data/*{forward,reverse}.fq.gz"
params.genome = "$baseDir/data/X.fasta"
params.multiqc = "$baseDir/multiqc"
params.outdir = "$baseDir/results"

log.info """\
    This is a NEXTFLOW RNA-SEQ PIPELINE

    genome: ${params.genome}
    read directory: ${params.reads}
    outdir: ${params.outdir}

    """
.stripIndent()


process index {

    cpus 2

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
    .into{ reads_ch1; reads_ch2 }


process align {

    tag "$sample_id"
    publishDir "$params.outdir/bam"
    cpus 2

    input:
    path genome from params.genome
    file '*' from bwa_index_ch
    tuple sample_id, path(reads) from reads_ch1

    output:
    file '*.bam' into bam_ch

    """
    bwa mem $genome ${reads[0]} ${reads[1]} > ${sample_id}.bam
    """
}


process fastqc {
    tag "FASTQC on $sampl_id"

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


process multiqc {

    publishDir params.outdir, mode: 'copy'

    input:
    path '*' from fastqc_ch.collect()

    output:
    path 'multiqc_report.html' into ch_multiqc_report

    script:
    """
    multiqc .
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
        checkHostname()
        log.info "-${c_purple}[nf-bwa-mem]${c_red} Pipeline completed with errors${c_reset}-"
    }

}
