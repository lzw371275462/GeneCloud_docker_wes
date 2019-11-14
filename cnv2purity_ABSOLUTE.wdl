####input
####a list contain cnv+vcfmaf in format:
####		<taskID	normal_name	tumor_name	cnv_acs_seg_path	vcf_maf_path>

workflow cnv2purity{
	File cnv_lists
	String out_dir_name

	String mnt_db_dir                       ####e.g.:/mnt/GeneCloud/luozw/oncoWES_db/
	String mnt_input_dir
	String mnt_out_dir
	String task_queue

	String out_dir=mnt_out_dir+"/"+out_dir_name
	
	####
	Array [Array[String]] cnv_lists_infos=read_tsv(cnv_lists)
	
	scatter(cnv_lists_info in cnv_lists_infos){
		call cnv2purity_do{
			input:
				taskID=cnv_lists_info[0],
				normalName=cnv_lists_info[1],
				sampleName=cnv_lists_info[2],
				cnv_file=cnv_lists_info[3],
				maf_file=cnv_lists_info[4],
				out_dir=out_dir,
				mnt_db_dir=mnt_db_dir,
				mnt_input_dir=mnt_input_dir,
				mnt_out_dir=mnt_out_dir,
				task_queue=task_queue
		}
	}

	####
	output {
		Array[String] purity_outs=cnv2purity_do.log_out
		Array[String] work_stat=cnv2purity_do.log_file
	}
}

task cnv2purity_do {
	String taskID
	String sampleName
	String normalName
	String cnv_file
	String maf_file
	String out_dir
	String mnt_db_dir
	String mnt_input_dir
	String mnt_out_dir
	String task_queue
	
	####
	command <<<
		set -e
		####
		if [ ! -d ${out_dir}/${taskID}/3_cnv2purity ] ; then 
			mkdir -p  ${out_dir}/${taskID}/3_cnv2purity
		fi
		
		cd ${out_dir}/${taskID}/3_cnv2purity

		echo -e "${taskID} 3_cnv2purity_work_start@"$(date +%Y-%m-%d_%H:%M:%S)"@"$(hostname) > log
		
		####step_1: revise maf_file header
		if [ ! -f ./${sampleName}.maf ] ; then
			less ${maf_file}|sed -e 's/Amino_acids/protein_change/' -e 's/HGNC_ID/UniProt_AAxform_aapos/'  > ./${sampleName}.maf
		fi
		
		####step_2:run ABSOLUTE
		if [ ! -f cancer.PP-calls_tab.txt ] ; then
			/usr/local/bin/miniconda3/bin/Rscript /usr/local/bin/miniconda3/bin/run_ABSOLUTE.R \
			${cnv_file} \
			${out_dir}/${taskID}/3_cnv2purity/${sampleName}.maf \
			${out_dir}/${taskID}/3_cnv2purity/ \
			1>run_ABSOLUTE\.stdout 2>run_ABSOLUTE\.stderr
			if [[ $? -ne 0 && -f cancer.PP-calls_tab.txt ]] ; then rm cancer.PP-calls_tab.txt ; fi
		fi

		####
		echo -e "${taskID} 3_cnv2purity_work_end@"$(date +%Y-%m-%d_%H:%M:%S)"@"$(hostname) >> log
		
		echo -e "${taskID}\t${normalName}\t${sampleName}\t${out_dir}/${taskID}/3_cnv2purity/cancer.PP-calls_tab.txt" >log.out
	>>>
	
	runtime {
		docker:"oncowes/absolutev1.2:v1"
		memory:"5G"
		num_proc:"4"
		maxRetries:3
		mnt_db_dir:mnt_db_dir
		mnt_input_dir:mnt_input_dir
		mnt_out_dir:mnt_out_dir
		task_queue:task_queue
	}
	
	output{
		String log_out="${out_dir}/${taskID}/3_cnv2purity/log.out"
		String log_file="${out_dir}/${taskID}/3_cnv2purity/log"
	}
}
