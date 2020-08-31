params.multiqc = "$baseDir/multiqc"


log.info """\
    This is a NEXTFLOW PIPELINE

    genome: ${params.genome}
    read directory: ${params.reads}
    outdir: ${params.outputDir}

    """
.stripIndent()

ch_multiqc_config = file("$baseDir/assets/multiqc_config.yaml", checkIfExists: true)

// Get unprocessed reads from params.samplePlan
// Channel
//     .fromPath(params.samplePlan)
//     .splitCsv(header:true)
//     .map{ row -> tuple(row.sample_id, row.group, row.tissue, file(row.r1), file(row.r2)) }
//     .into{ raw_reads_ch1; raw_reads_ch2  }


Channel
    .fromPath("$baseDir/sample_plan_local2.csv")
    .splitCsv(header:true)
    .map{ row -> tuple(
      row.tumour_id, file(row.tr1), file(row.tr2),
      row.normal_id, file(row.nr1), file(row.nr2)
      )}
    .into{ raw_reads_normal_ch; raw_reads_tumour_ch; raw_reads_test }

// Channel
//     .fromPath(params.sampleDescription)
//     .splitCsv(header:true)
//     .map{ row -> tuple(row.tumour_id, row.normal_id) }
//     .into{ tumour_normal_mappings_ch; test_ch }

// tumours = Channel.create()
// normals = Channel.create()
//
// Channel
//     .fromPath(params.sampleDescription)
//     .splitCsv(header:true)
//     .map{ row -> tuple(row.tumour_id, row.normal_id) }
//     .separate( tumours, normals )

// Get tumour normal pairings from params.sampleDescription
// Channel
//     .fromPath(params.sampleDescription)
//     .splitCsv(header:true)
//     .flatMap{ row -> tuple(tumour: row.tumour_id, normal:row.normal_id) }
//     .into{ tumour_normal_mappings_ch; test_ch }

// test_ch.view()

process n_trimmomatic {
    label 'trimmomatic'
    tag "$normal_id"
    publishDir "$params.outputDir/processed_reads"
    echo true

    input:
    // set sample_id, group, tissue, file(r1), file(r2) from raw_reads_ch1
    set tumour_id, _, _, normal_id, file(r1), file(r2) from raw_reads_normal_ch

    output:
    tuple sample_id, tumour_id, "${sample_id}.*.fq.gz" into trimmed_reads_normal_ch

    script:
    // sample_id = tumour_id + '_normal'
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
    // set sample_id, group, tissue, file(r1), file(r2) from raw_reads_ch1
    set tumour_id, file(r1), file(r2), _, _, _ from raw_reads_tumour_ch

    output:
    tuple sample_id, tumour_id, "${sample_id}.*.fq.gz" into trimmed_reads_tumour_ch

    script:
    // sample_id = tumour_id + '_tumour'
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

// Channel
//     .fromFilePairs(params.reads, checkIfExists: true)
//     .map { it -> [ it[0].split('\\.')[0], it[1] ] }
//     .into{ reads_ch1; reads_ch2 }


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
    // tuple tumour_id, normal_id from tumour_normal_mappings_ch.view()

    output:
    tuple sample_id, tumour_id, "${sample_id}.bam" into (t_bam_ch, t_bamstats_in_ch)

    script:
    """
    bwa mem -t $task.cpus $genome ${reads[0]} ${reads[1]} | samblaster --addMateTags --removeDups | samtools sort - | samtools view -Sb - > ${sample_id}.bam
    #picard AddOrReplaceReadGroups -INPUT ${sample_id}.bam -OUTPUT ${sample_id}.RG.bam -VALIDATION_STRINGENCY LENIENT -RGID ${sample_id} -RGLB HUM -RGPL illumina -RGPU 1 -RGSM ${sample_id}
    samtools index ${sample_id}.bam
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
    // tuple tumour_id, normal_id from tumour_normal_mappings_ch.view()

    output:
    tuple sample_id, tumour_id, "${sample_id}.bam" into (n_bam_ch, n_bamstats_in_ch)

    script:
    """
    bwa mem -t $task.cpus $genome ${reads[0]} ${reads[1]} | samblaster --addMateTags --removeDups | samtools sort - | samtools view -Sb - > ${sample_id}.bam
    #picard AddOrReplaceReadGroups -INPUT ${sample_id}.bam -OUTPUT ${sample_id}.RG.bam -VALIDATION_STRINGENCY LENIENT -RGID ${sample_id} -RGLB HUM -RGPL illumina -RGPU 1 -RGSM ${sample_id}
    samtools index ${sample_id}.bam
    """
}

t_bam_ch
  .join(n_bam_ch, by:[1])
  .into{ tn_pileup; tn_varscan }

t_bamstats_in_ch
  .concat(n_bamstats_in_ch)
  .set{ bamstats_in_ch }

process pileup {

  label 'varscan'
  tag "$tumour_id"
  echo true

  input:
  path genome from params.genome
  tuple tumour_id, _, path(tumor_bam), normal_id, path(normal_bam) from tn_pileup

  // output:
  // tuple tumour_id, "${tumour_id}.pileup" into pileup_ch
  //

  // when:
  // sample_id.contains(tumour_id)
  script:
  println "samtools mpileup -C50 -q 1 -f $genome $normal_bam $tumor_bam > ${tumour_id}.pileup"

  """
  echo samtools mpileup -C50 -q 1 -f $genome $normal_bam $tumor_bam > ${tumour_id}.pileup
  """
}

// pileup_ch.view()

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

/********************************
 * Header log info
 */

// log.info """=======================================================
//
//  @nf-bwa-test@ workflow v${workflow.manifest.version}
// ======================================================="""
// def summary = [:]
//
// summary['Max Memory']     = params.maxMemory
// summary['Max CPUs']       = params.maxCpus
// summary['Max Time']       = params.maxTime
// summary['Container Engine'] = workflow.containerEngine
// summary['Current home']   = "$HOME"
// summary['Current user']   = "$USER"
// summary['Current path']   = "$PWD"
// summary['Working dir']    = workflow.workDir
// summary['Output dir']     = params.outputDir
// summary['Config Profile'] = workflow.profile
//
// if(params.email) summary['E-mail Address'] = params.email
// log.info summary.collect { k,v -> "${k.padRight(15)}: $v" }.join("\n")
// log.info "========================================="

// TODO - ADD YOUR NEXTFLOW PROCESS HERE

// workflow.onComplete {
//
//     // pipeline_report.html
//
//     def reportFields = [:]
//     reportFields['version'] = workflow.manifest.version
//     reportFields['runName'] = customRunName ?: workflow.runName
//     reportFields['success'] = workflow.success
//     reportFields['dateComplete'] = workflow.complete
//     reportFields['duration'] = workflow.duration
//     reportFields['exitStatus'] = workflow.exitStatus
//     reportFields['errorMessage'] = (workflow.errorMessage ?: 'None')
//     reportFields['errorReport'] = (workflow.errorReport ?: 'None')
//     reportFields['commandLine'] = workflow.commandLine
//     reportFields['projectDir'] = workflow.projectDir
//     reportFields['summary'] = summary
//     reportFields['summary']['Date Started'] = workflow.start
//     reportFields['summary']['Date Completed'] = workflow.complete
//     reportFields['summary']['Pipeline script file path'] = workflow.scriptFile
//     reportFields['summary']['Pipeline script hash ID'] = workflow.scriptId
//     if(workflow.repository) reportFields['summary']['Pipeline repository Git URL'] = workflow.repository
//     if(workflow.commitId) reportFields['summary']['Pipeline repository Git Commit'] = workflow.commitId
//     if(workflow.revision) reportFields['summary']['Pipeline Git branch/tag'] = workflow.revision
//
//
//     // Render the TXT template
//     def engine = new groovy.text.GStringTemplateEngine()
//     def tf = new File("$baseDir/assets/onCompleteTemplate.txt")
//     def txtTemplate = engine.createTemplate(tf).make(reportFields)
//     def reportTxt = txtTemplate.toString()
//
//     // Render the HTML template
//     def hf = new File("$baseDir/assets/onCompleteTemplate.html")
//     def htmlTemplate = engine.createTemplate(hf).make(reportFields)
//     def reportHtml = htmlTemplate.toString()
//
//     // Write summary e-mail HTML to a file
//     def outputSummaryDir = new File( "${params.summaryDir}/" )
//     if( !outputSummaryDir.exists() ) {
//       outputSummaryDir.mkdirs()
//     }
//     def outputHtmlFile = new File( outputSummaryDir, "pipelineReport.html" )
//     outputHtmlFile.withWriter { w -> w << reportHtml }
//     def outputTxtFile = new File( outputSummaryDir, "pipelineReport.txt" )
//     outputTxtFile.withWriter { w -> w << reportTxt }
//
//     // onComplete file
//
//     File woc = new File("${params.outputDir}/onComplete.txt")
//     Map endSummary = [:]
//     endSummary['Completed on'] = workflow.complete
//     endSummary['Duration']     = workflow.duration
//     endSummary['Success']      = workflow.success
//     endSummary['exit status']  = workflow.exitStatus
//     endSummary['Error report'] = workflow.errorReport ?: '-'
//     String endWfSummary = endSummary.collect { k,v -> "${k.padRight(30, '.')}: $v" }.join("\n")
//     println endWfSummary
//     String execInfo = "${fullSum}\nExecution summary\n${logSep}\n${endWfSummary}\n${logSep}\n"
//     woc.write(execInfo)
//
//     // final logs
//
//     c_green = "\033[0;32m";
//     c_purple = "\033[0;35m";
//     c_red =  "\033[0;31m";
//     c_reset = "\033[0m";
//
//     if (workflow.stats.ignoredCount > 0 && workflow.success) {
//         log.info "-${c_purple}Warning, pipeline completed, but with errored process(es) ${c_reset}-"
//         log.info "-${c_red}Number of ignored errored process(es) : ${workflow.stats.ignoredCount} ${c_reset}-"
//         log.info "-${c_green}Number of successfully ran process(es) : ${workflow.stats.succeedCount} ${c_reset}-"
//     }
//
//     if (workflow.success) {
//         log.info "-${c_purple}[nf-bwa-mem]${c_green} Pipeline completed successfully${c_reset}-"
//     } else {
//         log.info "-${c_purple}[nf-bwa-mem]${c_red} Pipeline completed with errors${c_reset}-"
//     }
//
//
//     if(workflow.success){
//       log.info "-${c_purple}[nf-bwa-mem]${c_green} Pipeline completed successfully${c_reset}-"
//     }else{
//       log.info "-${c_purple}[nf-bwa-mem]${c_red} Pipeline completed with errors${c_reset}-"
//         if( workflow.profile == 'test'){
//             log.error "====================================================\n" +
//                     "  WARNING! You are running with the profile 'test' only\n" +
//                     "  pipeline config profile, which runs on the head node\n" +
//                     "  and assumes all software is on the PATH.\n" +
//                     "  This is probably why everything broke.\n" +
//                     "  Please use `-profile test,conda` or `-profile test,singularity` to run on local.\n" +
//                     "  Please use `-profile test,conda,cluster` or `-profile test,singularity,cluster` to run on your cluster.\n" +
//                     "============================================================"
//         }
//     }
//
// }

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
