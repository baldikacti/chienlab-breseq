#!/usr/bin/env nextflow

nextflow.enable.dsl=2

/*
 * Input parameters: read pairs, reference and output
 * The configuration is in nextflow.config file 
 * Params are stored in the params.config file 
 */

// this prevents a warning of undefined parameter
params.help             = false

// Input parameters
params.genbank = "references/NC_011916.gbk"
params.reads = "${params.project}/*_R{1,2}_001.fastq.gz"
params.outdir = "${params.project}/breseq_analysis_results"

// this prints the input parameters
log.info"""
CHIENLAB BRESEQ PIPELINE
=============================================
reads                           : ${params.reads}
output				            : ${params.outdir}
""".stripIndent()

// this prints the help in case you use --help parameter in the command line and it stops the pipeline
if (params.help) {
    log.info 'This is the ChienLab\'s Nextflow Breseq pipeline'
    log.info 'Please define path to your project in params.config file!\n'
    log.info 'Enjoy!'
    log.info '\n'
    exit 1
}

// Input channels
read_pairs_ch =  Channel.fromFilePairs(params.reads, checkIfExists: true)
gb_ch = Channel.fromPath(params.genbank, checkIfExists: true)

//  The default workflow
workflow {

    run_breseq(read_pairs_ch, gb_ch)

}

// Workflow steps
process run_breseq {
	label 'process_high'
	conda './envs/breseq.yaml'
    tag "$pair_id"
	
	publishDir "${params.outdir}", mode: 'move'

	// ‘each’ is used to use the reference genbank file multiple times,
	// otherwise only the first read set would have been processed
	input:
	tuple val(pair_id), path(reads)
	each reference

	output:
	path "${pair_id}/*"

	script:
	"""
	breseq -j ${task.cpus} -r $reference -n $pair_id -o $pair_id $reads
	"""
}

/*
================================================================================
    Completion summary
================================================================================
*/

c_green = "\033[0;32m";
c_reset = "\033[0m"

workflow.onComplete {
    log.info"""
    Execution status: ${ workflow.success ? 'OK' : 'failed' }
    ${c_green}Results are reported here: $params.outdir${c_reset}
    ￼""".stripIndent()
}