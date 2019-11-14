####input
####bam list in format:
####            <taskID	normalName	normalBam	tumorName	tumorBam>


workflow bam2hlaType{
	File bam_lists
	String out_dir_name
	
	String mnt_db_dir			####e.g.:/mnt/GeneCloud/luozw/oncoWES_db/
	String mnt_input_dir
	String mnt_out_dir
	String task_queue
	
	String out_dir=mnt_out_dir+"/"+out_dir_name
	
	String reference_path_2=mnt_db_dir+"/hlaType_db/hla_reference_dna.fasta"
	String optitype_config_path_2=mnt_db_dir+"/hlaType_db/config.ini.example"
	String all_hla_type_path_2=mnt_db_dir+"/hlaType_db/hlaAllType.txt"

	####
	Array [Array[String]] bam_lists_infos=read_tsv(bam_lists)
	
	scatter(bam_lists_info in bam_lists_infos){
		call sam2fq_do{
			input:
				task_ID=bam_lists_info[0],
				normal_name=bam_lists_info[1],
				normal_bam=bam_lists_info[2],			
				tumor_name=bam_lists_info[3],
				tumor_bam=bam_lists_info[4],
				out_dir_path=out_dir,
				mnt_db_dir=mnt_db_dir,
				mnt_input_dir=mnt_input_dir,
				mnt_out_dir=mnt_out_dir,
				task_queue=task_queue
		}
		call fq2bam_do{
			input:
				task_ID=bam_lists_info[0],
				normal_name=bam_lists_info[1],
				normal_bam=bam_lists_info[2],			
				tumor_name=bam_lists_info[3],
				tumor_bam=bam_lists_info[4],
				out_dir_path=out_dir,
				fq_1=sam2fq_do.fq_1,
				fq_2=sam2fq_do.fq_2,
				reference_path=reference_path_2,
				mnt_db_dir=mnt_db_dir,
				mnt_input_dir=mnt_input_dir,
				mnt_out_dir=mnt_out_dir,
				task_queue=task_queue
		}
		call bam2fq_do{
			input:
				task_ID=bam_lists_info[0],
				normal_name=bam_lists_info[1],
				normal_bam=bam_lists_info[2],			
				tumor_name=bam_lists_info[3],
				tumor_bam=bam_lists_info[4],
				out_dir_path=out_dir,
				bam_1=fq2bam_do.bam1,
				bam_2=fq2bam_do.bam2,
				mnt_db_dir=mnt_db_dir,
				mnt_input_dir=mnt_input_dir,
				mnt_out_dir=mnt_out_dir,
				task_queue=task_queue
		}
		call bam2hlaType_do{
			input:
				task_ID=bam_lists_info[0],
				normal_name=bam_lists_info[1],
				normal_bam=bam_lists_info[2],			
				tumor_name=bam_lists_info[3],
				tumor_bam=bam_lists_info[4],
				out_dir_path=out_dir,
				fq_1=bam2fq_do.fq_1,
				fq_2=bam2fq_do.fq_2,
				optitype_config_path=optitype_config_path_2,
				all_hla_type_path=all_hla_type_path_2,				
				mnt_db_dir=mnt_db_dir,
				mnt_input_dir=mnt_input_dir,
				mnt_out_dir=mnt_out_dir,
				task_queue=task_queue
		}
	}
	####
	output {
		Array[String] hlaType_out_infos=bam2hlaType_do.log_out
		Array[String] work_stat=bam2hlaType_do.log_file
	}
}


task sam2fq_do {
	String task_ID
	String normal_name
	String normal_bam
	String tumor_name
	String tumor_bam
	String out_dir_path
	String mnt_db_dir
	String mnt_input_dir
	String mnt_out_dir
	String task_queue
	
	command <<<
		set -e
		if [ ! -d ${out_dir_path}/${task_ID}/2_bam2hlaType/class_I_optiType ] ; then 
			mkdir -p ${out_dir_path}/${task_ID}/2_bam2hlaType/class_I_optiType
		fi
		cd ${out_dir_path}/${task_ID}/2_bam2hlaType/class_I_optiType
		echo -e "${task_ID} bam2hlaType_work_start@"$(date +%Y-%m-%d_%H:%M:%S)"@"$(hostname) > log
		
		####step_1 bam2fq
		if [[ ! -f origin_1.fq || ! -f origin_2.fq ]] ; then
			gatk SamToFastq --FASTQ origin_1.fq --SECOND_END_FASTQ  origin_2.fq --INPUT ${normal_bam} 1>origin_1.fq.stdout 2>origin_1.fq.stderr
			if [[ $? -ne 0 && -f origin_1.fq ]] ; then rm origin_1.fq ; fi
		fi

		echo -e "${task_ID} bam2hlaType_work_end@"$(date +%Y-%m-%d_%H:%M:%S)"@"$(hostname) >> log
	>>>

	runtime {
		docker:"oncowes/gatk4:latest"
		cpu: "1"
		memory:"5G"
		num_proc:"8"
		maxRetries:3
		mnt_db_dir:mnt_db_dir
		mnt_input_dir:mnt_input_dir
		mnt_out_dir:mnt_out_dir
		task_queue:task_queue
	}
	
	output{
		String log_file="${out_dir_path}/${task_ID}/2_bam2hlaType/class_I_optiType/log"
		String fq_1="${out_dir_path}/${task_ID}/2_bam2hlaType/class_I_optiType/origin_1.fq"
		String fq_2="${out_dir_path}/${task_ID}/2_bam2hlaType/class_I_optiType/origin_2.fq"
	}
}

task fq2bam_do {
	String task_ID
	String normal_name
	String normal_bam
	String tumor_name
	String tumor_bam
	String out_dir_path
	String fq_1
	String fq_2
	String reference_path
	String mnt_db_dir
	String mnt_input_dir
	String mnt_out_dir
	String task_queue
	
	command <<<
		set -e
		
		if [ ! -d ${out_dir_path}/${task_ID}/2_bam2hlaType/class_I_optiType ] ; then 
			mkdir -p ${out_dir_path}/${task_ID}/2_bam2hlaType/class_I_optiType
		fi
		cd ${out_dir_path}/${task_ID}/2_bam2hlaType/class_I_optiType
		
		echo -e "${task_ID} bam2hlaType_work_start@"$(date +%Y-%m-%d_%H:%M:%S)"@"$(hostname) >> log

		####step_2 fq2bam
		if [ ! -f origin_1.fq.bam ] ; then
			/usr/local/bin/miniconda3/bin/razers3 -i 95 -m 1 -dr 0 -o origin_1.fq.bam ${reference_path} ${fq_1} 1>origin_1.fq.bam.stdout 2>origin_1.fq.bam.stderr
			if [[ $? -ne 0 && -f origin_1.fq.bam ]] ; then rm origin_1.fq.bam ; fi
		fi
		
		if [ ! -f origin_2.fq.bam ] ; then
			/usr/local/bin/miniconda3/bin/razers3 -i 95 -m 1 -dr 0 -o origin_2.fq.bam ${reference_path} ${fq_2} 1>origin_2.fq.bam.stdout 2>origin_2.fq.bam.stderr
			if [[ $? -ne 0 && -f origin_2.fq.bam ]] ; then rm origin_2.fq.bam ; fi
		fi
		
		echo -e "${task_ID} bam2hlaType_work_end@"$(date +%Y-%m-%d_%H:%M:%S)"@"$(hostname) >> log
		
	>>>

	runtime {
		docker:"oncowes/tools:v6"
		cpu: "1"
		memory:"10G"
		num_proc:"8"
		maxRetries:3
		mnt_db_dir:mnt_db_dir
		mnt_input_dir:mnt_input_dir
		mnt_out_dir:mnt_out_dir
		task_queue:task_queue
	}
	
	output{
		String log_file="${out_dir_path}/${task_ID}/2_bam2hlaType/class_I_optiType/log"
		String bam1="${out_dir_path}/${task_ID}/2_bam2hlaType/class_I_optiType/origin_1.fq.bam"
		String bam2="${out_dir_path}/${task_ID}/2_bam2hlaType/class_I_optiType/origin_2.fq.bam"
	}
}

task bam2fq_do {
	String task_ID
	String normal_name
	String normal_bam
	String tumor_name
	String tumor_bam
	String out_dir_path
	String bam_1
	String bam_2
	String mnt_db_dir
	String mnt_input_dir
	String mnt_out_dir
	String task_queue
	
	command <<<
		set -e
		
		if [ ! -d ${out_dir_path}/${task_ID}/2_bam2hlaType/class_I_optiType ] ; then 
			mkdir -p ${out_dir_path}/${task_ID}/2_bam2hlaType/class_I_optiType
		fi
		cd ${out_dir_path}/${task_ID}/2_bam2hlaType/class_I_optiType
		
		echo -e "${task_ID} bam2hlaType_work_start@"$(date +%Y-%m-%d_%H:%M:%S)"@"$(hostname) >> log
		
		####step_3 bam2fq
		if [[ ! -f origin_1.fq.bam_1.fq || ! -f origin_2.fq.bam_2.fq ]] ; then
			/usr/local/bin/miniconda3/bin/samtools bam2fq ${bam_1} > origin_1.fq.bam_1.fq
			/usr/local/bin/miniconda3/bin/samtools bam2fq ${bam_2} > origin_2.fq.bam_2.fq
			if [[ $? -ne 0 && -f origin_1.fq.bam_1.fq ]] ; then rm origin_1.fq.bam_1.fq ; fi
		fi

		echo -e "${task_ID} bam2hlaType_work_end@"$(date +%Y-%m-%d_%H:%M:%S)"@"$(hostname) >> log
	>>>

	runtime {
		docker:"oncowes/tools:v3"
		cpu: "1"
		memory:"8G"
		num_proc:"8"
		maxRetries:3
		mnt_db_dir:mnt_db_dir
		mnt_input_dir:mnt_input_dir
		mnt_out_dir:mnt_out_dir
		task_queue:task_queue
	}
	
	output{
		String log_file="${out_dir_path}/${task_ID}/2_bam2hlaType/class_I_optiType/log"
		String fq_1="${out_dir_path}/${task_ID}/2_bam2hlaType/class_I_optiType/origin_1.fq.bam_1.fq"
		String fq_2="${out_dir_path}/${task_ID}/2_bam2hlaType/class_I_optiType/origin_2.fq.bam_2.fq"
	}
}
task bam2hlaType_do {
	String task_ID
	String normal_name
	String normal_bam
	String tumor_name
	String tumor_bam
	String out_dir_path
	String fq_1
	String fq_2
	String optitype_config_path
	String all_hla_type_path
	String mnt_db_dir
	String mnt_input_dir
	String mnt_out_dir
	String task_queue
	
	command <<<
		set -e
		
		if [ ! -d ${out_dir_path}/${task_ID}/2_bam2hlaType/class_I_optiType ] ; then 
			mkdir -p ${out_dir_path}/${task_ID}/2_bam2hlaType/class_I_optiType
		fi
		cd ${out_dir_path}/${task_ID}/2_bam2hlaType/class_I_optiType
		
		echo -e "${task_ID} bam2hlaType_work_start@"$(date +%Y-%m-%d_%H:%M:%S)"@"$(hostname) >> log
	
		alias razers3="/usr/local/bin/razers3"
		mkdir -p /path/to/
		ln -s /usr/local/bin/razers3 /path/to/razers3	

		####step_4 analysis OptiTypePipeline
		if [ ! -f ${normal_name}_result.tsv ] ; then
			python /usr/local/bin/OptiType/OptiTypePipeline.py \
			--input ${fq_1} ${fq_2} \
			--dna \
			--outdir ./ \
			--prefix ${normal_name} \
			--config ${optitype_config_path} \
			1>${normal_name}_result.tsv.stdout 2>${normal_name}_result.tsv.stderr
			if [[ $? -ne 0 && -f ${normal_name}_result.tsv ]] ; then rm ${normal_name}_result.tsv ; fi
		fi

		####step_5 result transform
		less ${normal_name}_result.tsv |tail -1|sed 's/\s/\n/g'|grep -E 'A|B|C'|sed -e 's/[*:]/_/g' -e 's/A/hla_a/' -e 's/B/hla_b/' -e 's/C/hla_c/'|while read a ; do
			grep $a ${all_hla_type_path} | head -1
		done > ${normal_name}_result.tsv.transform
		
		####step_6 neoantigen transform
		less ${normal_name}_result.tsv |tail -1|cut -f 2-7|sed -e 's/\t/,/g' -e 's/A/HLA-A/g' -e 's/B/HLA-B/g' -e 's/C/HLA-C/g' >${normal_name}_result.tsv.neoantigen
		chmod 777 *
		
		echo -e "${task_ID} bam2hlaType_work_end@"$(date +%Y-%m-%d_%H:%M:%S)"@"$(hostname) >> log
		echo -e "${task_ID}\t${normal_name}\t"$(readlink -f ${normal_name}_result.tsv.transform)"\t"$(readlink -f ${normal_name}_result.tsv.neoantigen) >log.out
		
	>>>

	runtime {
		docker:"oncowes/hlatype_optitype:v3"
		cpu: "1"
		memory:"10G"
		num_proc:"8"
		maxRetries:3
		mnt_db_dir:mnt_db_dir
		mnt_input_dir:mnt_input_dir
		mnt_out_dir:mnt_out_dir
		task_queue:task_queue	
	}
	
	output{
		String log_out="${out_dir_path}/${task_ID}/2_bam2hlaType/class_I_optiType/log.out"
		String log_file="${out_dir_path}/${task_ID}/2_bam2hlaType/class_I_optiType/log"
	}
}
