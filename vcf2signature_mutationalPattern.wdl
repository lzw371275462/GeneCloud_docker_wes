####input
#####a list contain vcfs in format:
#####			 <taskID	vcf_lists_cfg>

workflow MutSigCaller {
	File vcf_lists_cfg
	
	String out_dir_name
	String mnt_db_dir
	String mnt_input_dir
	String mnt_out_dir
	String task_queue
	String out_dir=mnt_out_dir+"/"+out_dir_name

	String RefVersion
	String? Rank1
	String Rank=select_first([Rank1, 2])
	String? Cosmic1
	String Cosmic=select_first([Cosmic1, "legacy"])

	####
	Array[Array[String]] vcf_lists_infos=read_tsv(vcf_lists_cfg)
	scatter(vcf_lists_info in vcf_lists_infos){
		call MutationalPatterns2mutsig_do{
			input:
				taskName=vcf_lists_info[0],
				vcf_file=vcf_lists_info[1],
				out_dir=out_dir,
				mnt_db_dir=mnt_db_dir,
				mnt_input_dir=mnt_input_dir,
				mnt_out_dir=mnt_out_dir,
				task_queue=task_queue,
				RefVersion=RefVersion,
				Rank=Rank,
				Cosmic=Cosmic
		}
	}

	output {
		Array[String] work_stat=MutationalPatterns2mutsig_do.log_file
		Array[String] work_out=MutationalPatterns2mutsig_do.log_out
		Array[String] signature_results=MutationalPatterns2mutsig_do.signature_result
	}
}

	

task MutationalPatterns2mutsig_do {
	String taskName
	String vcf_file
	String out_dir
	String mnt_db_dir
	String mnt_input_dir
	String mnt_out_dir
	String task_queue
	String RefVersion
	String Rank
	String Cosmic
	
	command <<<
		set -e

		if [ ! -d ${out_dir}/${taskName} ] ; then
				mkdir -p ${out_dir}/${taskName}
		fi
		cd ${out_dir}/${taskName}

		if [ ! -d raw_vcf ] ; then mkdir -p raw_vcf ; fi
		cat ${vcf_file}|while read a b ; do 
			vcf_basename=$(basename $a)
			if [ ! -f raw_vcf/$vcf_basename ] ; then
				ln -s $a raw_vcf/
			fi
		done
		
		echo -e "vcf2Signature_MutationalPatterns_work_start@"$(date +%Y-%m-%d_%H:%M:%S)"@"$(hostname) > log

		export PATH=/usr/local/jdk1.8.0_231/bin/:/usr/local/bin/miniconda3/bin/:$PATH
		export LD_LIBRARY_PATH=/usr/local/jdk1.8.0_231/jre/lib/amd64/:/usr/local/jdk1.8.0_231/jre/lib/amd64/server/:$LD_LIBRARY_PATH
		if [ ! -f base_substitution_types.xls ] ; then
			/usr/local/bin/miniconda3/bin/Rscript /usr/local/bin/mutationalPatterns2signature.R \
			--ref ${RefVersion} \
			--vcf ${out_dir}/${taskName}/raw_vcf/ \
			--rank ${Rank} \
			--cosmic ${Cosmic} \
			--outdir ${out_dir}/${taskName}/  \
			1>base_substitution_types.xls.stdout 2>base_substitution_types.xls.stderr
			if [[ $? -ne 0 && -f base_substitution_types.xls ]] ; then rm base_substitution_types.xls ; fi
		fi
		
		echo -e "vcf2Signature_MutationalPatterns_work_start@"$(date +%Y-%m-%d_%H:%M:%S)"@"$(hostname) >> log
		echo -e "vcf2Signature_MutationalPatterns_result\t${out_dir}/${taskName}/base_substitution_types.xls" >log.out
	 >>>

	runtime {
		docker:"oncowes/vcf2signature:v3"
		cpu: "1"
		memory:"2G"
		num_proc:"1"
		maxRetrie:3
		mnt_db_dir:mnt_db_dir
		mnt_input_dir:mnt_input_dir
		mnt_out_dir:mnt_out_dir
		task_queue:task_queue
	}

	output {
		String signature_result = "${out_dir}/${taskName}/base_substitution_types.xls"
		String log_file="${out_dir}/${taskName}/log"
		String log_out="${out_dir}/${taskName}/log.out"
	}
}
