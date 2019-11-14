####input
####bam list in format:
####            <taskID normalName      normalBam       tumorName       tumorBam>

workflow bam2geneCov{
	File bamlists_file
	String out_dir_name
	String normalMinDepth=6
	String tumorMinDepth=8
	String minMapq=20
	String bpClassTypes='AT,CpG,CG'
	
	String mnt_db_dir			####e.g.:/mnt/GeneCloud/luozw/oncoWES_db/
	String mnt_input_dir
	String mnt_out_dir
	String task_queue
	String docker_user='$UID'
	String out_dir=mnt_out_dir+"/"+out_dir_name
	
	String all_gene_tsv=mnt_db_dir+"/maf2signatureGene_db/all_gene.tsv"
	String ref_fa=mnt_db_dir+"/hs37d5_ref/hs37d5.fa"
	
	Array [Array[String]] bam_lists_infos=read_tsv(bamlists_file)
	
	scatter(bam_lists_info in bam_lists_infos){
		call bam2geneCov_do{
			input:
				taskName=bam_lists_info[0],
				normal_name=bam_lists_info[1],
				normal_bam=bam_lists_info[2],
				tumor_name=bam_lists_info[3],
				tumor_bam=bam_lists_info[4],
				all_gene_tsv=all_gene_tsv,
				ref_fa=ref_fa,
				NormalMinDepth=normalMinDepth,
				TumorMinDepth=tumorMinDepth,
				MinMapq=minMapq,
				BpClassTypes=bpClassTypes,
				out_dir=out_dir,
				mnt_db_dir=mnt_db_dir,
				mnt_input_dir=mnt_input_dir,
				mnt_out_dir=mnt_out_dir,
				task_queue=task_queue,
				docker_user=docker_user
		}
	}
	
	####
	output {
		Array[String] music2_coverage_result=bam2geneCov_do.music2_coverage
		Array[String] music2_logs=bam2geneCov_do.log_file
		Array[String] music2_log_outs=bam2geneCov_do.log_out
	}
}

##########
task bam2geneCov_do {
	String taskName
	String normal_name
	String normal_bam
	String tumor_name
	String tumor_bam
	
	String all_gene_tsv
	String ref_fa
	String NormalMinDepth
	String TumorMinDepth
	String MinMapq
	String BpClassTypes
	String out_dir
	
	String mnt_db_dir
	String mnt_input_dir
	String mnt_out_dir
	String task_queue
	String docker_user

	command {
		set -e
		
		if [ ! -d ${out_dir}/${taskName}/2_bam2geneCov_music2 ] ; then
			mkdir -p  ${out_dir}/${taskName}/2_bam2geneCov_music2
		fi
		cd ${out_dir}/${taskName}/2_bam2geneCov_music2
		echo -e "${taskName} 2_bam2geneCov_music2_work_start@"$(date +%Y-%m-%d_%H:%M:%S)"@"$(hostname) > log
		
		export PATH=/usr/local/bin/miniconda3/bin/:$PATH
		####step1 Calculate_coverage
		if [ ! -f Calculate_coverage.SUCCESS ] ; then
			music2 bmr calc-covg-helper \
			--normal-tumor-bam-pair="${tumor_name}	${normal_bam}	${tumor_bam}" \
			--roi-file=${all_gene_tsv} \
			--reference-sequence=${ref_fa} \
			--output-file=./${tumor_name}.covg \
			--normal-min-depth=${NormalMinDepth} \
			--tumor-min-depth=${TumorMinDepth} \
			--min-mapq=${MinMapq}  \
			--bp-class-types ${BpClassTypes} \
			&& touch Calculate_coverage.SUCCESS
		fi
		
		echo -e "${taskName} 2_bam2geneCov_music2_work_end@"$(date +%Y-%m-%d_%H:%M:%S)"@"$(hostname) >> log
		echo -e "bam2geneCov_music2_result\t${out_dir}/${taskName}/2_bam2geneCov_music2/${tumor_name}.covg" > log.out 
	}

	runtime {
		docker:"oncowes/music2:v3"
		memory:"5G"
		num_proc:"4"
		maxRetries:3
		mnt_db_dir:mnt_db_dir
		mnt_input_dir:mnt_input_dir
		mnt_out_dir:mnt_out_dir
		task_queue:task_queue
		docker_user:docker_user
	}
	
	output {
		String music2_coverage = "${out_dir}/${taskName}/2_bam2geneCov_music2/${tumor_name}.covg"
		String log_file = "${out_dir}/${taskName}/2_bam2geneCov_music2/log"
		String log_out = "${out_dir}/${taskName}/2_bam2geneCov_music2/log_out"
	}
}
