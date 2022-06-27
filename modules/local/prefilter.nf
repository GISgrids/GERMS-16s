#!/usr/bin/env nextflow


def BWAVERSION = '0.7.17'
def pysamVERSION = '0.19.1'
params.yaml="prefilter.yml"

/*
fasta_ch=Channel.from(params.fasta)
spikein_ch=Channel.from(params.spike_in_name)
inputcheck_ch=Channel.from(params.samplesheet)
*/

process pre_filter {
	tag "$meta.id"
	label "process_medium"
	conda "${params.yaml}"
	//container <docker container path>'
		
   	input:
	tuple val(meta) , path(fasta)
	val bwa_index
	val spike_in
	
   	output:
	tuple val(meta) , path ("*.flt.trimmed.fastq.gz"), emit: PF_reads
	file '*_ratios.txt', emit:ratios

   	when: 
   	task.ext.when == null || task.ext.when

	script:

	"""
	bwa mem \\
	-t $tasks.cpus \\
	-M ${bwa_index} ${R1} ${R2} \\
	$args | \\

	S16_prefilter.py \\
	-i - \\
	-1 ${PF_reads[0]} \\
	-2 ${PF_reads[1]} \\
	-s ${spike_in} \\
	> ${sample_id}_ratios.txt \\

	cat <<-END_VERSIONS >> versions.yml 
	"${task.process}":
		BWA: $BWAVERSION
		pysam: $pysamVERSION
	END_VERSIONS 




	"""

   	"""
   	bwa mem \\
	-t $task.cpus \\
	-k 30 \\
	-M ${bwa_index} ${R1} ${R2} | \\





	S16_prefilter.py \\
	-i - -1 R1.flt.trimmed.fastq.gz -2 R2.flt.trimmed.fastq -s ${spike_in} > ${sample_id}_ratios.txt 2>> SUCCESS_step3_pre_filter
	echo "COMPLETED Step3 (pre_filter) : ${sample_id}" >> SUCCESS_step3_pre_filter
	"""
}


/*
include { INPUT_CHECK } from '../subworkflows/local/input_check'


include { BuildBWAindexes} from '../module/local/BuildBWAindexes'

workflow{

	INPUT_CHECK (
	        inputcheck_ch
	    )
	    .reads
	    .map {
	        meta, fastq ->
	            meta.id = meta.id.split('_')[0..-2].join('_')
	            [ meta, fastq ] }
	    .groupTuple(by: [0])
	    .branch {
	        meta, fastq ->
	            single  : fastq.size() == 1
	                return [ meta, fastq.flatten() ]
	            multiple: fastq.size() > 1
	                return [ meta, fastq.flatten() ]
	    }
	ch_buildbwa=BuildBWAindexes.out.bwa_index
	prefilter (INPUT_CHECK.out.reads, ch_buildbwa, spikein_ch)
}
*/
