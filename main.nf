params.multiqc = "$baseDir/multiqc"

c_green = "\033[0;32m";
c_purple = "\033[0;35m";
c_red =  "\033[0;31m";
c_reset = "\033[0m";
b_green = "\033[1;32m";
// c_blue = "\e[1;34mLight Blue Text\e[0m"

// BLUE='\033[0;34m'
c_blue='\033[0;34m'

log.info """\
===
This is a NEXTFLOW PIPELINE for detecting LOH in matched tumour normal pairs.

Running with the following user-defined options:
--- ${b_green}config${c_reset} ---
  --samplePlan: ${c_blue}${params.samplePlan}${c_reset}
  --outdir: ${c_blue}${params.outputDir}${c_reset}
  --condaCacheDir ${c_blue}${params.condaCacheDir}${c_reset}

--- ${b_green}Trimmomatic (v0.39)${c_reset} ---
  --adapters [ILLUMINACLIP:]  ${c_blue}${params.adapters}${c_reset}

--- ${b_green}Bwa (v0.7.17)${c_reset} ---
  --genome  ${c_blue}${params.genome}${c_reset}

--- ${b_green}Varscan (v2.4.4)${c_reset} ---
  --varscan_normal_coverage [--min-coverage-normal] ${params.varscan_normal_coverage}
  --varscan_tumour_coverage [--min-coverage-tumor] ${params.varscan_tumour_coverage}
  --varscan_tumour_purity [--tumor-purity] ${params.varscan_tumour_purity}

--- ${b_green}lohcator ($baseDir/bin/lohcator.py)${c_reset} ---
  --sample_ids [--config] ${c_blue}${params.sample_ids}${c_reset}
  --lohcator_chromosome [--chromosome] ${params.lohcator_chromosome}

===
    """
.stripIndent()

ch_multiqc_config = file("$baseDir/assets/multiqc_config.yaml", checkIfExists: true)

Channel
    .fromPath(params.samplePlan)
    .splitCsv(header:true)
    .map{ row -> tuple(
      row.tumour_id, file(row.tr1), file(row.tr2),
      row.normal_id, file(row.nr1), file(row.nr2)
      )}
    .into{ raw_reads_normal_ch; raw_reads_tumour_ch; raw_reads_test }


process n_trimmomatic {
    label 'trimmomatic'
    tag "$normal_id"
    publishDir "$params.outputDir/processed_reads"
    echo true

    input:
    set tumour_id, _, _, normal_id, file(r1), file(r2) from raw_reads_normal_ch

    output:
    tuple sample_id, tumour_id, "${sample_id}.*.fq.gz" into trimmed_reads_normal_ch

    script:
    sample_id = normal_id
    """
    trimmomatic PE \
        -threads 8 \
        -phred33 \
        $r1 $r2 \
        ${sample_id}.forward.fq.gz \
        ${sample_id}.unpaired_1.fastq.gz \
        ${sample_id}.reverse.fq.gz \
        ${sample_id}.unpaired_2.fastq.gz \
        ILLUMINACLIP:$params.adapters:2:30:10 \
        LEADING:3 \
        TRAILING:3 \
        SLIDINGWINDOW:4:15 \
        MINLEN:30
    """
}

process t_trimmomatic {
    label 'trimmomatic'
    tag "$tumour_id"
    publishDir "$params.outputDir/processed_reads"
    echo true

    input:
    set tumour_id, file(r1), file(r2), _, _, _ from raw_reads_tumour_ch

    output:
    tuple sample_id, tumour_id, "${sample_id}.*.fq.gz" into trimmed_reads_tumour_ch

    script:
    sample_id = tumour_id
    """
    trimmomatic PE \
        -threads 8 \
        -phred33 \
        $r1 $r2 \
        ${sample_id}.forward.fq.gz \
        ${sample_id}.unpaired_1.fastq.gz \
        ${sample_id}.reverse.fq.gz \
        ${sample_id}.unpaired_2.fastq.gz \
        ILLUMINACLIP:$params.adapters:2:30:10 \
        LEADING:3 \
        TRAILING:3 \
        SLIDINGWINDOW:4:15 \
        MINLEN:30
    """
}

trimmed_reads_normal_ch
  .into{ n_reads_ch1; n_reads_ch2 }

trimmed_reads_tumour_ch
  .into{ t_reads_ch1; t_reads_ch2; t_reads_test }


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


process align_t {
    label 'medCpu'
    label 'medMem'
    label 'bwa'
    tag "$sample_id"
    publishDir "$params.outputDir/bam"
    echo true

    input:
    path genome from params.genome
    file '*' from bwa_index_ch
    tuple sample_id, tumour_id, path(reads) from t_reads_ch1

    output:
    tuple sample_id, tumour_id, "${sample_id}.RG.bam" into (t_bam_ch, t_bamstats_in_ch)

    script:
    """
    bwa mem -t $task.cpus $genome ${reads[0]} ${reads[1]} | samblaster --addMateTags --removeDups | samtools sort - | samtools view -Sb - > ${sample_id}.bam
    picard AddOrReplaceReadGroups -INPUT ${sample_id}.bam -OUTPUT ${sample_id}.RG.bam -VALIDATION_STRINGENCY LENIENT -RGID ${sample_id} -RGLB HUM -RGPL illumina -RGPU 1 -RGSM ${sample_id}
    samtools index ${sample_id}.RG.bam
    """
}

process align_n {
    label 'medCpu'
    label 'medMem'
    label 'bwa'
    tag "$sample_id"
    publishDir "$params.outputDir/bam"
    echo true

    input:
    path genome from params.genome
    file '*' from bwa_index_ch
    tuple sample_id, tumour_id, path(reads) from n_reads_ch1

    output:
    tuple sample_id, tumour_id, "${sample_id}.RG.bam" into (n_bam_ch, n_bamstats_in_ch)

    script:
    """
    bwa mem -t $task.cpus $genome ${reads[0]} ${reads[1]} | samblaster --addMateTags --removeDups | samtools sort - | samtools view -Sb - > ${sample_id}.bam
    picard AddOrReplaceReadGroups -INPUT ${sample_id}.bam -OUTPUT ${sample_id}.RG.bam -VALIDATION_STRINGENCY LENIENT -RGID ${sample_id} -RGLB HUM -RGPL illumina -RGPU 1 -RGSM ${sample_id}
    samtools index ${sample_id}.RG.bam
    """
}

t_bam_ch
  .join(n_bam_ch, by:[1])
  .into{ tn_pileup; tn_varscan; tn_freebayes }

t_bamstats_in_ch
  .concat(n_bamstats_in_ch)
  .set{ bamstats_in_ch }


// tn_freebayes.view()


process pileup {
  label 'varscan'
  tag "$tumour_id"

  input:
  path genome from params.genome
  tuple tumour_id, _, path(tumor_bam), normal_id, path(normal_bam) from tn_pileup

  output:
  tuple tumour_id, "${tumour_id}.pileup" into pileup_ch

  script:
  """
  samtools mpileup -C50 -q 1 -f $genome $normal_bam $tumor_bam > ${tumour_id}.pileup
  """
}

// process freebayes {
//     label 'freebayes'
//     tag "$tumour_id"
//     echo true
//     publishDir "$params.outputDir/vcf"
//
//     input:
//     path genome from params.genome
//     path unmappable_genome from params.unmappable_genome
//     tuple tumour_id, _, path(tumor_bam), normal_id, path(normal_bam) from tn_freebayes
//
//     output:
//     tuple tumour_id, "${tumour_id}*.vcf*" into freebayes_out_ch
//     tuple tumour_id, "${tumour_id}_snps_filt.vcf.gz", "${tumour_id}_snps_filt.vcf.gz.tbi" into freebayes_raw_out_ch
//
//     script:
//     """
//     echo "Running freebays on ${tumour_id} vs ${normal_id}"
//
//     freebayes -0 -f $genome \
//       --pooled-discrete \
//       --genotype-qualities \
//       --min-coverage 20 \
//       $tumor_bam \
//       $normal_bam \
//       -v ${tumour_id}_raw.vcf
//
//     vcfintersect ${tumour_id}_raw.vcf -b ${params.unmappable_genome} -v > ${tumour_id}_mappable.vcf
//
//     vcfallelicprimitives ${tumour_id}_mappable.vcf | \
//       vt decompose_blocksub - | \
//       vt normalize -q -r $genome - | \
//       vcfsamplediff -s VT ${normal_id} ${tumour_id} - | \
//       vcffilter -f "DP > 20" \
//         -f "QUAL > 1 & QUAL / AO > 10" \
//         -f "SAF > 0 & SAR > 0" \
//         -f "RPR > 0 & RPL > 0" \
//         -f "TYPE = snp" > ${tumour_id}_snps_filt.vcf
//
//     vcffilter -f "VT = somatic" ${tumour_id}_snps_filt.vcf > ${tumour_id}_freebayes.vcf
//
//     bgzip -f ${tumour_id}_snps_filt.vcf
//     tabix -p vcf ${tumour_id}_snps_filt.vcf.gz
//
//     bgzip -f ${tumour_id}_freebayes.vcf
//     tabix -p vcf ${tumour_id}_freebayes.vcf.gz
//     """
// }

// freebayes_raw_out_ch.view()


process varscan {
  label 'varscan'
  tag "$tumour_id"
  echo true
  publishDir "$params.outputDir/varscan"
  publishDir "$params.outputDir/vcf", pattern: "*.vcf"


  input:
  path genome from params.genome
  path unmappable_genome from params.unmappable_genome
  tuple tumour_id, "${tumour_id}.pileup" from pileup_ch

  output:
  tuple tumour_id, "${tumour_id}.snp", "${tumour_id}.indel", "${tumour_id}.*.hc", "${tumour_id}*.vcf" into varscan_out_ch


  script:
  """
  varscan somatic ${tumour_id}.pileup \
    ${tumour_id} \
    --mpileup 1 \
    --min-coverage-normal ${params.varscan_normal_coverage} \
    --min-coverage-tumor ${params.varscan_tumour_coverage} \
    --tumor-purity ${params.varscan_tumour_purity} \
    --strand-filter 1
    varscan processSomatic ${tumour_id}.snp
    varscan processSomatic ${tumour_id}.indel

    python '$baseDir/bin/VarScan2_format_converter.py' ${tumour_id}.snp.LOH.hc > ${tumour_id}_snp.LOH.hc.vcf
    vcfintersect ${tumour_id}_snp.LOH.hc.vcf -b ${params.unmappable_genome} -v > ${tumour_id}_snp.LOH.hc.filt.vcf

    #python '$baseDir/bin/VarScan2_format_converter.py' ${tumour_id}.snp.Somatic.hc > ${tumour_id}_snp.Somatic.hc.vcf
    #python '$baseDir/bin/VarScan2_format_converter.py' ${tumour_id}.indel.Somatic.hc > ${tumour_id}_indel.Somatic.hc.vcf
    #vcfintersect ${tumour_id}_snp.Somatic.hc.vcf -b ${params.unmappable_genome} -v > ${tumour_id}_snp.Somatic.hc.filt.vcf
    #vcfintersect ${tumour_id}_indel.Somatic.hc.vcf -b ${params.unmappable_genome} -v > ${tumour_id}_indel.Somatic.hc.filt.vcf

    cat ${tumour_id}.snp.*.hc > ${tumour_id}.snp.hc
  """
}


process lohcator {
    label 'lohcator'
    tag "$tumour_id"
    echo true
    publishDir "$params.outputDir/bed"
    // publishDir "$params.outputDir/vcf", pattern: "*.vcf"

    input:
    // tuple tumour_id, "${tumour_id}_snps_filt.vcf.gz", "${tumour_id}_snps_filt.vcf.gz.tbi"  from freebayes_raw_out_ch
    tuple tumour_id, "${tumour_id}.snp" from varscan_out_ch

    output:
    tuple tumour_id, "${tumour_id}_LOH_regions.bed" into lohcator_out_ch

    script:
    """
    echo Running lohcator $tumour_id for Varscan file ${tumour_id}.snp
    python $baseDir/bin/lohcator.py \
        --varscan_file ${tumour_id}.snp \
        --config ${params.sample_ids} \
        --chromosome ${params.lohcator_chromosome}
    """
}


// varscan_out_ch.view()

// Collect stats on various stages

process fastqc {

    label 'fastqc'
    tag "$sample_id"

    input:
    tuple sample_id, pair_id, path(reads) from n_reads_ch2.mix(t_reads_ch2)

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
  tag "$sample_id"

  input:
  tuple sample_id, _, file(bam) from bamstats_in_ch

  output:
  path "${sample_id}.stats" into bamstats_ch

  script:
  """
  bamtools stats -in $bam > ${sample_id}.stats
  """

}

process multiqc {

    publishDir "$params.outputDir", mode: 'copy'
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

// END processes

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
        log.info "-${c_purple}[nf-lohcator]${c_green} Pipeline completed successfully${c_reset}-"
    } else {
        log.info "-${c_purple}[nf-lohcator]${c_red} Pipeline completed with errors${c_reset}-"
    }

}
