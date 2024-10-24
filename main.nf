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
if (params.gbk)  { ch_gbk = file(params.gbk, checkIfExists: true)  } else { exit 1, 'At least one reference needs to be specified' }
if (params.gbk2) { ch_gbk2 = file(params.gbk2, checkIfExists: true) }

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
//gb_ch = Channel.fromPath(params.gbk, checkIfExists: true)

//  The default workflow
workflow {

	if(params.gbk2) {
		run_breseq_multi(read_pairs_ch, ch_gbk, ch_gbk2)
	} else {
		run_breseq_single(read_pairs_ch, ch_gbk)
	}


}

// Workflow steps
process run_breseq_single {
	label 'process_high'
    tag "$pair_id"
	conda "bioconda::breseq=0.39.0"
	
	publishDir "${params.outdir}", mode: 'move'

	// ‘each’ is used to use the reference genbank file multiple times,
	// otherwise only the first read set would have been processed
	input:
	tuple val(pair_id), path(reads)
	each ref1

	output:
	path "${pair_id}/*"

	script:
	"""
	breseq -j ${task.cpus} -r $ref1 -n $pair_id -o $pair_id $reads
	"""
}

process run_breseq_multi {
	label 'process_high'
    tag "$pair_id"
	conda "bioconda::breseq=0.39.0"
	
	publishDir "${params.outdir}", mode: 'move'

	// ‘each’ is used to use the reference genbank file multiple times,
	// otherwise only the first read set would have been processed
	input:
	tuple val(pair_id), path(reads)
	each ref1
	each ref2

	output:
	path "${pair_id}/*"

	script:
	"""
	breseq -j ${task.cpus} -r $ref1 -r $ref2 -n $pair_id -o $pair_id $reads
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