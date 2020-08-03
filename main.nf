#!/usr/bin/env nextflow

// Channel.from(1, 2, 3, 6)
//     .filter{ it > 1 }
//     .map{ it -> [it, it * it]}
//     .view{ num, sqr -> "This is > 1: $num\nThe square of $num = $sqr" }

params.genome = "$baseDir/data/dmel_6.12.fa"
println "Using genome: $params.genome"

println("GroupByTuple:")
Channel
    .from([1,'A'], [1,'B'], [2,'C'], [3,'B'], [1,'C'], [2,'A'])
    .groupTuple(by:0)
    .view()

left = Channel.from(['X', 1], ['Y', 2], ['Z', 3], ['P', 7])
right = Channel.from(['Z', 6], ['Y', 5], ['X', 4])

println("join:")
left
    .join(right, remainder: true)
    .view()

words_ch = Channel
    .from('hello', 'world', 'again')
    .map{ word -> [word, word.size()]}

process greeting {
    echo true

    input:
    tuple val(word), val(word_size) from words_ch

    script:
    """
    echo word: $word, length: $word_size
    """
}

files_ch = Channel
    .fromPath('data/*.fq.gz', checkIfExists: true)
    .map{ file -> [file.name, file ]}

process fileecho {
    echo true

    input:
    tuple val(fname), file(reads) from files_ch

    script:
    """
    echo Alignment of sample $fname with $reads
    """
}

// Split a channel into two
reads_ch = Channel
    .fromFilePairs('data/*{forward,reverse}.fq.gz', checkIfExists: true)
    .into{ reads_ch1; reads_ch2 }

process align_sample {
    tag "$sample_id"

    input:
    path genome from params.genome
    tuple val(sample_id), path(reads) from reads_ch1

    // when:
    // genome.name =~ /6.13/

    output:
    tuple val(sample_id), file("${sample_id}.bam") into bam_ch

    script:
    """
    echo bwa mem $genome $sample_id > ${sample_id}.bam
    """
}

bam_ch.view()

process alignSequences {
    tag "$sample_id"

    methods = ['small_widow', 'medium_window', 'large_window']

    input:
    tuple val(sample_id), path(reads) from reads_ch2
    each mode from methods

    """
    echo somefun -in $reads --mode $mode
    """
}

//
// process index_sample {
//
//     input:
//     file sample_bam from bam_ch
//
//     output:
//     file sample_bai into bai_ch
//
//     script:
//     """
//     echo samtools index $sample_bam
//     """
// }
