#!/usr/bin/env nextflow

def VERSION = '0.5.2'
params.yaml="emirge_vs_gg.yml"

process emirge_vs_gg {
    	
	tag "$meta.id"
	label "process_medium"
    conda "${params.yaml}"
	
	//publishDir path: "${params.outdir}/${sample_id}" , mode: 'copy' , pattern: 'greengenes-hits-graphmap.bam'
	//publishDir path: "${params.outdir}/AllSample/rawtable" , mode: 'copy' , pattern: "*_raw-table.txt"
	publishDir path: "${params.outdir}/${sample_id}/logs" , mode: 'copy' , pattern: 'SUCCESS_*'
	
    input:
    tuple val(meta) , file(emirge) 
    val src 
    val greengenes
    path taxonomy
	
    output:
	tuple val(meta) , file('greengenes-hits-graphmap.bam') , file("${sample_id}_raw-table.txt"), emit: emirge_vs_gg
    file 'SUCCESS_*', emit: success
	// the part following graphmap should actually be on a new container
	script:   
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
	graphmap \\
		align \\
		$args \\
		-@ $task.cpus \\
		-d ${emirge} \\
		--out out.bam 2 \\ 
        >> SUCCESS_step6_emirgeVSgg
		$args \\
    
    
    ident_to_bam.py out.bam greengenes-hits-graphmap.bam ${greengenes} 2>> SUCCESS_step6_emirgeVSgg
	classify_hits.py -f -q ${emirge} -i greengenes-hits-graphmap.bam -t ${taxonomy} -o ${meta}_raw-table.txt 2>> SUCCESS_step7_classify_hits
	rm out.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        emirge_vs_gg: $VERSION
    END_VERSIONS
    """

}
