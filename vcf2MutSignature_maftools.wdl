####input
####a list contain maf files in format:
####		<taskID	maf_paths>

workflow maf2MutSignature{
	File maf_file_lists
	String out_dir_name
	
	String mnt_db_dir			####e.g.:/mnt/GeneCloud/luozw/oncoWES_db/
	String mnt_input_dir
	String mnt_out_dir
	String task_queue
	String out_dir=mnt_out_dir+"/"+out_dir_name
	
	Array [Array[String]] maf_infos=read_tsv(maf_file_lists)
	
	scatter(maf_info in maf_infos){
		call vcf2MutSignature_do{
			input:
				taskName=maf_info[0],
				maf_file=maf_info[1],
				mnt_db_dir=mnt_db_dir,
				mnt_input_dir=mnt_input_dir,
				mnt_out_dir=mnt_out_dir,
				task_queue=task_queue,
				out_dir=out_dir
		}
	}
	####
	output {
		Array[String] signature_outs=vcf2MutSignature_do.log_out
		Array[String] work_stats=vcf2MutSignature_do.log_file
	}
}

##########
task vcf2MutSignature_do {
	String taskName
	String maf_file
	String out_dir
	String mnt_db_dir
	String mnt_input_dir
	String mnt_out_dir
	String task_queue
	
	####
	command <<<
		set -e
		
		if [ ! -d ${out_dir}/${taskName}/ ] ; then
			mkdir -p  ${out_dir}/${taskName}/
		fi

		cd ${out_dir}/${taskName}/
		echo -e "${taskName} vcf2MutSignature_work_start@"$(date +%Y-%m-%d_%H:%M:%S)"@"$(hostname) > log
		
		####
		export PATH=/usr/local/jdk1.8.0_231/bin/:$PATH
		export LD_LIBRARY_PATH=/usr/local/jdk1.8.0_231/jre/lib/amd64/:/usr/local/jdk1.8.0_231/jre/lib/amd64/server/:$LD_LIBRARY_PATH

		if [ ! -f titv_table.txt ] ; then
			/usr/local/bin/miniconda3/bin/Rscript \
			/usr/local/bin/maf2signature.R \
			--maf ${maf_file} \
			1>titv_table.txt.stdout 2>titv_table.txt.stderr
			if [[ $? -ne 0 && -f titv_table.txt ]] ; then rm titv_table.txt ; fi
		fi
		
		echo -e "${taskName} vcf2MutSignature_work_end@"$(date +%Y-%m-%d_%H:%M:%S)"@"$(hostname) >> log
		echo -e "maf2MutSignature_result\t${out_dir}/${taskName}/titv_table.txt" > log.out 
	>>>
	
	runtime {
		docker:"oncowes/vcf2signature:v3"
		memory:"5G"
		num_proc:"4"
		maxRetries:3
		mnt_db_dir:mnt_db_dir
		mnt_input_dir:mnt_input_dir
		mnt_out_dir:mnt_out_dir
		task_queue:task_queue
	}
	
	output{
		String mut_signature_out="${out_dir}/${taskName}/titv_table.txt"
		String log_out="${out_dir}/${taskName}/log.out"
		String log_file="${out_dir}/${taskName}/log"
	}
}
