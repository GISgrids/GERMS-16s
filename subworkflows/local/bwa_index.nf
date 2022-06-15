include { BuildBWAindexes } from '../../modules/local/BuildBWAindexes.nf'

workflow BWA_INDEX_FASTA {
    take:
    fasta // file: /path/to/fasta file     
    
    main:
    ch_bwa_index = Channel.empty()
    ch_bwa_index = BuildBWAindexes ( fasta ).bwa_index

    emit:
    index       = ch_bwa_index       //    path: bwa_index
}