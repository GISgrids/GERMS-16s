def VERSION = '0.0.7' // Version information provided by tool on CLI BUT I AM LAZY to parse it...so deal with it

process BuildBWAindexes {
    tag "${fasta}"
    label 'process_low'

    //container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    //'https://depot.galaxyproject.org/singularity/mulled-v2-e5d375990341c5aef3c9aff74f96f66f65375ef6:38aed4501da19db366dc7c8d52d31d94e760cfaf-0' :
    //'quay.io/biocontainers/mulled-v2-e5d375990341c5aef3c9aff74f96f66f65375ef6:38aed4501da19db366dc7c8d52d31d94e760cfaf-0' }"	

    input:
    path(fasta)

    output:
    path("${fasta}.*"), emit: bwa_index

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    """
    bwa index ${fasta}
    """
}