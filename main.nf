params.multiqc = "$baseDir/multiqc"


log.info """\
    This is a NEXTFLOW PIPELINE

    genome: ${params.genome}
    read directory: ${params.reads}
    outdir: ${params.outputDir}

    """
.stripIndent()

ch_multiqc_config = file("$baseDir/assets/multiqc_config.yaml", checkIfExists: true)

params.samplePlan = 'sample_plan.csv'

Channel
    .fromPath(params.samplePlan)
    .splitCsv(header:true)
    .map{ row -> tuple(row.sample_id, row.group, row.tissue, file(row.r1), file(row.r2)) }
    .set{ raw_reads_ch }


process trimmomatic {
    label 'trimmomatic'
    tag "$sample_id"
    publishDir "$params.outputDir/processed_reads"
    echo true

    input:
    set sample_id, group, tissue, file(r1), file(r2) from raw_reads_ch

    script:
    """
    echo your_command --sample $sample_id --group $group --tissue $tissue --reads $r1 $r2
    """
}


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
    publishDir "$params.outputDir/bam"

    input:
    path genome from params.genome
    file '*' from bwa_index_ch
    tuple sample_id, path(reads) from reads_ch1

    output:
    tuple sample_id, "${sample_id}.RG.bam" into bam_ch
    tuple sample_id, "${sample_id}.RG.bam" into bamstats_in_ch

    script:
    """
    bwa mem -t $task.cpus $genome ${reads[0]} ${reads[1]} | samblaster --addMateTags --removeDups | samtools sort - | samtools view -Sb - > ${sample_id}.bam
    picard AddOrReplaceReadGroups -INPUT ${sample_id}.bam -OUTPUT ${sample_id}.RG.bam -VALIDATION_STRINGENCY LENIENT -RGID ${sample_id} -RGLB HUM -RGPL illumina -RGPU 1 -RGSM ${sample_id}
    samtools index ${sample_id}.RG.bam
    """
}

process pileup {

  label 'varscan'
  tag "$sample_id"

  echo true

  input:
  path genome from params.genome
  tuple sample_id, file(bam) from bam_ch

  // output:
  // path "${sample_id}.stats" into bamstats_ch

  script:
  """
  echo samtools mpileup -C50 -q 1 -f $genome normal tumour > ${sample_id}.pileup
  """

}

process fastqc {

    label 'fastqc'
    tag "$sample_id"

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
  tag "$sample_id"

  input:
  tuple sample_id, file(bam) from bamstats_in_ch

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
