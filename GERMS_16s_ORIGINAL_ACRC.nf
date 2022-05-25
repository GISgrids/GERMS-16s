/*	This is a remake of the original 16S pipeline, designed for
*	Illumina shotgun sequencing of 16S rRNA amplicon sequences (Ong et
*	al., 2013, PMID 23579286). At its core it's running EMIRGE (Miller et
*	al., 2011, PMID 21595876) or EMIRGE amplicon (Miller et al., 2013;
*	PMID 23405248) for reconstructing the sequences. Sequences are then
*	mapped against a preclustered version of Greengenes for
*	classification.
*/

//includeConfig 'nextflow.config'

/* 
 * pipeline input parameters 
 */
params.reads = "fastq/*/*_R{1,2}_unaligned_*.fastq.gz"

// Output location
params.outdir = "results_16s"

// Compute and codes 
params.s3_dir = "/home/users/astar/gis/lizh/GERMS_16s_resource"
params.src_dir = "${params.s3_dir}/src"
params.FAMAS = "${params.s3_dir}/software/famas-0.0.7/bin/famas"
params.yaml = "16s.yml"

// Databases
params.ssu = "${params.s3_dir}/database/ssu/SILVA_138.1_SSURef_NR99_tax_silva_trunc.ge800bp.le3500bp.0.97.fixed.fasta"
params.ssu_db = "${params.s3_dir}/database/ssu/SILVA_138.1_SSURef_NR99_tax_silva_trunc.ge800bp.le3500bp.0.97.fixed"
params.gg_13_05 = "${params.s3_dir}/database/greengenes/gg_13_5_tax_with_hOTUs_99_reps.fasta"
params.gg_taxonomy = "${params.s3_dir}/database/greengenes/gg_13_5_tax_with_hOTUs_99_reps_nou.taxonomy"
params.spike_in_name = "Plasmodium.knowlesi.profilin"


// 	Tools parameters
//	a. emirge
params.max_read_len = 151
params.ins_size = 200
params.ins_stdev = 40
params.iter_arg = 10 

// 	b. primer_trimmer
params.minlen = 200
params.maxlen = 1400

log.info """\
         GERMS-16s - N F   P I P E L I N E    
         ===================================
         SSU 		  	: ${params.ssu}
		 Greengenes	  	: ${params.gg_13_05}
         reads        	: ${params.reads}
         outdir       	: ${params.outdir}
		 src_dir		: ${params.s3_dir}
         """
         .stripIndent()

 
/* 
 * 	Create Channel
 */
reads_ch = Channel
			.fromFilePairs( params.reads )
			.map { 	item -> 
					sampleName = item[0];
					sampleName = sampleName.split('_')[0]
					files  = item[1];
					return [ sampleName, files ]  }
		
/* 
* 	Step 1: trim - famas
*	LEGACY CODE FROM AQUILA
*	params: min3qual='3', minlen='60'
*	/mnt/software/stow/famas-0.0.7/bin/famas -i {input.fq1} -j {input.fq2} -o {output.fq1} -p {output.fq2} -q {params.min3qual} -l {params.minlen}	
*	Updated to fastp trim
*	https://github.com/OpenGene/fastp#all-options
*/
process famas {
    	
	tag "$sample_id"
	label "process_low"
	
	//publishDir path: "${params.outdir}/${sample_id}" , mode: 'copy' , pattern: '*.fastq.gz'
	publishDir path: "${params.outdir}/${sample_id}/logs" , mode: 'copy' , pattern: 'SUCCESS_*'
	
    input:
    tuple val(sample_id) , path(fastq) from reads_ch
	val famas from params.FAMAS
	
    output:
	tuple val(sample_id) , file('R1.trimmed.fastq.gz') , file('R2.trimmed.fastq.gz') into reads_trimmed_ch
	file 'SUCCESS_*' into log_ch_1
	
	script:
    """
	${famas} -i ${fastq[0]} -j ${fastq[1]} -o R1.trimmed.fastq.gz -p R2.trimmed.fastq.gz -q 3 -l 60 2>> SUCCESS_step1_famas
	echo "COMPLETED Step1 (famas) : ${sample_id}" >> SUCCESS_step1_famas
	"""
}

// Step2 : Ignore concat for now - FIXME later

/* 
* 	Step 3: pre_filter
*	LEGACY CODE FROM AQUILA
*	params:   ssu_fa=config['SSU_FA'], spikein_name=config['SPIKEIN-NAME']
*	{config[BWA]} mem -t {threads} -k 30 -M {params.ssu_fa} {input.fq1} {input.fq2} | {config[PREFILTER]} -i - -1 {output.fq1} -2 {output.fq2} -s {params.spikein_name} > {output.ratios}	
*/
process pre_filter {
    	
	tag "$sample_id"
	label "process_medium"
	conda "${params.yaml}"
		
	publishDir path: "${params.outdir}/${sample_id}/results" , mode: 'copy' , pattern: '*_ratios.txt'
	publishDir path: "${params.outdir}/${sample_id}/logs" , mode: 'copy' , pattern: 'SUCCESS_*'
		
    input:
	tuple val(sample_id) , file(R1) , file(R2) from reads_trimmed_ch
	val bwa_index from params.ssu
	val s3 from params.s3_dir
	val spike_in from params.spike_in_name
	
    output:
	tuple val(sample_id) , file('R1.flt.trimmed.fastq.gz') , file('R2.flt.trimmed.fastq') into reads_filtered_ch
	file '*_ratios.txt' into ratio_file_ch
    file 'SUCCESS_*' into log_ch_3
	
	script:   
    """
    bwa mem -t $task.cpus -k 30 -M ${bwa_index} ${R1} ${R2} | ${s3}/src/S16_prefilter.py -i - -1 R1.flt.trimmed.fastq.gz -2 R2.flt.trimmed.fastq -s ${spike_in} > ${sample_id}_ratios.txt 2>> SUCCESS_step3_pre_filter
	echo "COMPLETED Step3 (pre_filter) : ${sample_id}" >> SUCCESS_step3_pre_filter
	"""
}

/* 
* 	Step 4: emirge
*	NOTE: emirge only takes in fastq (not gz) for read2
*	LEGACY CODE FROM AQUILA
*	params :    max_read_len=config['MAX_READ_LEN'], ins_size=config['INS_SIZE'], ins_stdev=config['INS_STDEV'], ssu_fa=config['SSU_FA'], ssu_db=config['SSU_DB'],
*              iter_arg="-n 1" if config['DEBUG'] else "-n 10"
*	{EMIRGE} -l {params.max_read_len} -i {params.ins_size} -s {params.ins_stdev} --phred33 {params.iter_arg} -a {threads} emirge -f {params.ssu_fa} -b {params.ssu_db} -1 {input.fq1} -2 {input.fq2} && touch {output}'	
*/
process emirge {
    	
	tag "$sample_id"
	label "process_medium"
	conda "${params.yaml}"
	
	publishDir path: "${params.outdir}/${sample_id}" , mode: 'copy' , pattern: 'emirge_out.fa'
	// publishDir path: "${params.outdir}/AllSample/summary" , mode: 'copy' , pattern: "*.alignment_summary.txt"
	publishDir path: "${params.outdir}/${sample_id}/logs" , mode: 'copy' , pattern: 'SUCCESS_*'
	
    input:
	tuple val(sample_id) , file(R1) , file(R2) from reads_filtered_ch
	val max_read_len from params.max_read_len
	val ins_size from params.ins_size
	val ins_stdev from params.ins_stdev
	val iter_arg from params.iter_arg
	val ssu_fa from params.ssu
	val ssu_db from params.ssu_db
	val s3 from params.s3_dir
	
    output:
	tuple val(sample_id) , file('emirge_out.fa') into emirge_ch
	file '*.alignment_summary.txt' into alignment_summary_ch
	file 'SUCCESS_*' into log_ch_4
	
	shell:   
    """
    emirge_amplicon.py emirge -l !{max_read_len} -i !{ins_size} -s !{ins_stdev} --phred33 -n !{iter_arg} -a !{task.cpus} -f !{ssu_fa} -b !{ssu_db} -1 !{R1} -2 !{R2}  2>> SUCCESS_step4_emirge
	emirge_rename_fasta.py emirge/iter.!{iter_arg} > emirge_out.fa 2>>SUCCESS_step4_emirge
	zgrep '^# ' emirge/iter.10/bowtie.iter.10.log.gz > alignment_text.csv
	echo "${sample_id}\t\$(grep -oP ' alignment: \\K(\\d+)\\s' alignment_text.csv)" >> ${sample_id}.alignment_summary.txt
	echo "COMPLETED Step4 (EMIRGE) : !{sample_id}" >> SUCCESS_step4_emirge
	"""
}

/* 
* 	Step 5: primer_trimmer
*	LEGACY CODE FROM AQUILA
*	params: minlen='200', maxlen='1400'
*	{config[PRIMER_TRIMMER]} -f -i {input} -o {output} --minlen {params.minlen} --maxlen {params.maxlen}	
*/
process primer_trimmer {
    	
	tag "$sample_id"
	label "process_low"
	
	publishDir path: "${params.outdir}/${sample_id}" , mode: 'copy' , pattern: 'emirge_outprimer_trimmed.fa'
	publishDir path: "${params.outdir}/${sample_id}/logs" , mode: 'copy' , pattern: 'SUCCESS_*'
    
	input:
    tuple val(sample_id) , file(emirge) from emirge_ch
	val src from params.src_dir
	val minlen from params.minlen
	val maxlen from params.maxlen
	
    output:
	tuple val(sample_id) , file('emirge_outprimer_trimmed.fa') into emirge_primerTrimmed_ch
    file 'SUCCESS_*' into log_ch_5
	
	script:   
    """
    ${src}/primer_trimmer.py -f -i $emirge -o emirge_outprimer_trimmed.fa --minlen $minlen --maxlen $maxlen 2>> SUCCESS_step5_primerTrim
	echo "COMPLETED Step5 (primer_trim) : ${sample_id}" >> SUCCESS_step5_primerTrim
	"""
}

/* 
* 	Step 6: Graphmap2 to Greengenes 13_05
*	LEGACY CODE FROM AQUILA
*	{config[GRAPHMAP]} align -t {threads} -r {input.ref} -d {input.fa} | {config[IDENT_TO_BAM]}  - {output} {input.ref}	
*/
/* 
* 	Step 7: classify hits
*	LEGACY CODE FROM AQUILA
*	{config[CLASSIFY_HITS]} -f -q {input.query} -i {input.hits} -t {input.gg_tax} -o {output}	
*/
process emirge_vs_gg_N_classifyHits {
    	
	tag "$sample_id"
	label "process_medium"
	conda "${params.yaml}"
	
	publishDir path: "${params.outdir}/${sample_id}" , mode: 'copy' , pattern: 'greengenes-hits-graphmap.bam'
	publishDir path: "${params.outdir}/AllSample/rawtable" , mode: 'copy' , pattern: "*_raw-table.txt"
	publishDir path: "${params.outdir}/${sample_id}/logs" , mode: 'copy' , pattern: 'SUCCESS_*'
	
    input:
    tuple val(sample_id) , file(emirge) from emirge_primerTrimmed_ch
	val src from params.src_dir
	val greengenes from params.gg_13_05
	path taxonomy from params.gg_taxonomy
	val s3 from params.s3_dir
	
    output:
	tuple val(sample_id) , file('greengenes-hits-graphmap.bam') , file("${sample_id}_raw-table.txt") into emirge_vs_gg_ch
    file 'SUCCESS_*' into log_ch_6
	
	script:   
    """
    graphmap align -t $task.cpus -r ${greengenes} -d ${emirge} --out out.bam 2>> SUCCESS_step6_emirgeVSgg
	${src}/ident_to_bam.py out.bam greengenes-hits-graphmap.bam ${greengenes} 2>> SUCCESS_step6_emirgeVSgg
	${src}/classify_hits.py -f -q ${emirge} -i greengenes-hits-graphmap.bam -t ${taxonomy} -o ${sample_id}_raw-table.txt 2>> SUCCESS_step7_classify_hits
	rm out.bam
	echo "COMPLETED Step6 (emirge_vs_gg) : ${sample_id}" >> SUCCESS_step6_emirgeVSgg
	echo "COMPLETED Step7 (classify_hits) : ${sample_id}" >> SUCCESS_step7_classify_hits
	"""
}

/* 
* 	Step 8: rarefaction
*	LEGACY CODE FROM AQUILA
*	{config[PLOT_RAREFACTION]} -o {output} -t {input} -f	
*/
/* 
* 	Step 9: convert_table
*	LEGACY CODE FROM AQUILA
*	{config[CONVERT_TABLE]} {input} results/ && touch {output}	
*/
process rarefaction_convertTable {
    	
	tag "$sample_id"
	label "process_medium"
	
	publishDir path: "${params.outdir}/AllSample" , mode: 'copy' , pattern: 'results/*'
	publishDir path: "${params.outdir}/${sample_id}/logs" , mode: 'copy' , pattern: 'SUCCESS_*'
	
    input:
    tuple val(sample_id) , file(greengene_bam) , file(raw_table) from emirge_vs_gg_ch
	val src from params.src_dir
	
    output:
	tuple val(sample_id) , file('results/*') into rarefaction_ch
    file 'SUCCESS_*' into log_ch_8
	
	shell:   
    """
	mkdir -p results
	!{src}/plot_rarefaction.py -o results/rarefaction.pdf -t !{raw_table} -f 2>> SUCCESS_step8_rarefaction
	echo "COMPLETED Step8 (rarefaction) : !{sample_id}" >> SUCCESS_step8_rarefaction
	!{src}/convert_table.py !{raw_table} results/ && touch 'convert_tables.succeeded' 2>> SUCCESS_step9_convert_tables
	echo "COMPLETED Step9 (convert_table) : !{sample_id}" >> SUCCESS_step9_convert_tables
	find results/ -type f -name '*' -printf "mv '%h/%f' '%h/!{sample_id}_%f' ; " | bash -
	"""
}

alignment_summary_ch
	.collectFile(name: 'coverage.txt', newLine: false , storeDir: "${params.outdir}/AllSample")
