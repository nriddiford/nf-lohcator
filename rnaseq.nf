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
    echo fastqc -o fastqc_${sample_id}_logs -f fastq -q ${reads}
    """
}
