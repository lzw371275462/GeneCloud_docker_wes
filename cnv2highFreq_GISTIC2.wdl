####input
####a list contain cnv in format:
####		<taskName	cnv_called_seg_file>

workflow cnv2highFreq{
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
		call cnv2highFreq_do{
			input:
				taskName=cnv_lists_info[0],
				cnv_file=cnv_lists_info[1],
				out_dir=out_dir,
				mnt_db_dir=mnt_db_dir,
				mnt_input_dir=mnt_input_dir,
				mnt_out_dir=mnt_out_dir,
				task_queue=task_queue
		}
	}

	####
	output {
		Array[String] cnv_highFreq_outs=cnv2highFreq_do.log_out
		Array[String] work_stat=cnv2highFreq_do.log_file
	}
}

task cnv2highFreq_do {
	String taskName
	String cnv_file
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
		echo -e "${taskName} cnv2highFreq_work_start@"$(date +%Y-%m-%d_%H:%M:%S)"@"$(hostname) > log
		
		export PATH=/usr/local/Polyspace/R2019a/bin/:$PATH
		
		####step_1:run GISTIC2
		if [ ! -f scores.gistic ] ; then
			/usr/local/bin/gistic2 \
			-b ${out_dir}/${taskName}/ \
			-seg ${cnv_file} \
			-refgene /usr/local/bin/refgenefiles/hg19.mat \
			-ta 0.15 \
			-td 0.15 \
			-js 8 \
			-qvt 0.1 \
			-rx 1 \
			-broad 1 \
			-brlen 0.8 \
			-res 0.01 \
			-conf 0.9 \
			-genegistic 1 \
			-arb 1 \
			-savegene 1 \
			-smallmem 0 \
			-gcm extreme \
			1>run_GISTIC2.stdout 2>run_GISTIC2.stderr
			if [[ $? -ne 0 && -f scores.gistic ]] ; then rm scores.gistic ; fi
		fi

		####
		echo -e "${taskName} cnv2highFreq_work_end@"$(date +%Y-%m-%d_%H:%M:%S)"@"$(hostname) >> log
		echo -e "${taskName}\t${out_dir}/${taskName}/scores.gistic" >log.out
	>>>
	
	runtime {
		docker:"oncowes/gistic2:v1"
		cpu: "1"
		memory:"6G"
		num_proc:"4"
		maxRetries:3
		mnt_db_dir:mnt_db_dir
		mnt_input_dir:mnt_input_dir
		mnt_out_dir:mnt_out_dir
		task_queue:task_queue
	}
	
	output{
		String log_out="${out_dir}/${taskName}/log.out"
		String log_file="${out_dir}/${taskName}/log"
	}
}
