####input
####fqs list in format:
####            <taskID	sampleName	tissue_type	SM	LB	PU	fq1	fq2>


import "fq2bam_base.wdl" as fq2bam_base

workflow fq2bam {
	File fq_conf_file
	String chip_bed_path
	String platform
	String out_dir_name

	String mnt_db_dir			####e.g.:/mnt/GeneCloud/luozw/oncoWES_db/
	String mnt_input_dir
	String mnt_out_dir
	String task_queue
	String out_dir=mnt_out_dir+"/"+out_dir_name

	String fasta_path_2=mnt_db_dir+"/hs37d5_ref/hs37d5.fa"
	
	####get_conf_info
	call extrackFqConfig{
		input:
			config_file=fq_conf_file,
			out_dir=out_dir
	}
	
	####bwa2bam_start
	scatter(conf_info in extrackFqConfig.conf_infos){
		String taskID=conf_info[0]
		String sampleName=conf_info[1]

		call fq2bam_base.fq2bam_base{
			input:
				fq_conf_file=conf_info[2],
				platform=platform,
				mnt_db_dir=mnt_db_dir,
				mnt_input_dir=mnt_input_dir,
				mnt_out_dir=mnt_out_dir,
				out_dir_name=out_dir_name,
				task_queue=task_queue
		}
		
		call MergeSamFiles_do{
				input:
					taskID=taskID,
					sam_name=sampleName,
					out_dir=out_dir,
					input_bams=fq2bam_base.bam_out,
					mnt_db_dir=mnt_db_dir,
					mnt_input_dir=mnt_input_dir,
					mnt_out_dir=mnt_out_dir,
					task_queue=task_queue
		}
		
		call BamQC_do{		
				input:
					taskID=taskID,
					sam_name=sampleName,
					out_dir=out_dir,
					input_bam=MergeSamFiles_do.bam_out,
					interval_bed=chip_bed_path,
					fasta=fasta_path_2,
					mnt_db_dir=mnt_db_dir,
					mnt_input_dir=mnt_input_dir,
					mnt_out_dir=mnt_out_dir,
					task_queue=task_queue
		}
	}
	
	####
	output{
		Array[String] bam_infos=MergeSamFiles_do.bam_out
		Array[String] bam_qc_infos=BamQC_do.qc_out
		Array[String] work_stat=BamQC_do.log_file
	}
}


####################################
##sub tasks                       ##
####################################
task extrackFqConfig{
	File config_file
	String out_dir

	command <<<
		set -e
		less ${config_file}|cut -f 1|sort -u|while read a ; do
			if [ ! -d ${out_dir}/$a ] ; then
				mkdir -p ${out_dir}/$a				
			fi
		done
		
		less ${config_file}|grep -v ^#|awk '$3~/tumor|normal/'|sort -u|awk '{print $0 > "${out_dir}/"$1"/cfg_sample."$2}'
		
		less ${config_file}|grep -v ^#|awk '$3~/tumor|normal/'|cut -f 1-2|sort -u|while read a b ; do
			echo -e "$a\t$b\t${out_dir}/$a/cfg_sample.$b"
		done
	>>>
	
	runtime {
		maxRetries:3
		backend:"Local"
	}
	
	output{
		Array[Array[String]] conf_infos = read_tsv(stdout())
	}
}

task MergeSamFiles_do {
	String taskID
	String sam_name
	String out_dir
	Array[String] input_bams

	String mnt_db_dir
	String mnt_input_dir
	String mnt_out_dir
	String task_queue
	
	####
	command <<<
		set -e
		if [ ! -d "${out_dir}/${taskID}/1_bwa2bam/${sam_name}/" ] ; then
			mkdir -p "${out_dir}/${taskID}/1_bwa2bam/${sam_name}/"
		fi
		cd "${out_dir}/${taskID}/1_bwa2bam/${sam_name}/"

		echo -e "${taskID} MergeSamFiles_start@"$(date +%Y-%m-%d_%H:%M:%S)"@"$(hostname) >> log

		####step_5.1:run MergeSamFiles
		if [ ! -f ${sam_name}_bqsr_merged.bam ] ; then
			gatk --java-options "-Xms10G -Xmx10G" \
			MergeSamFiles \
			--INPUT ${sep=' --INPUT ' input_bams} \
			--OUTPUT ${sam_name}_bqsr_merged.bam \
			--SORT_ORDER coordinate \
			1>${sam_name}_bqsr_merged.bam.stdout 2>${sam_name}_bqsr_merged.bam.stderr
			if [[ $? -ne 0 && -f ${sam_name}_bqsr_merged.bam ]] ; then rm ${sam_name}_bqsr_merged.bam ; fi
		fi
		
		####step_5.2:run BuildBamIndex
		if [ ! -f ${sam_name}_bqsr_merged.bai ] ; then
			gatk --java-options "-Xms5G -Xmx5G" \
			BuildBamIndex \
			--INPUT ${sam_name}_bqsr_merged.bam \
			1>${sam_name}_bqsr_merged.bai.stdout 2>${sam_name}_bqsr_merged.bai.stderr
			if [[ $? -ne 0 && -f ${sam_name}_bqsr_merged.bai ]] ; then rm ${sam_name}_bqsr_merged.bai ; fi
		fi

		if [ ! -f ${sam_name}_bqsr_merged.bam.bai ] ; then
			ln -s ${sam_name}_bqsr_merged.bai ${sam_name}_bqsr_merged.bam.bai
		fi
		
		echo -e "${taskID} MergeSamFiles_work_end@"$(date +%Y-%m-%d_%H:%M:%S)"@"$(hostname) >> log
	>>>
	
	runtime {
		docker:"oncowes/gatk4:latest"
		memory:"5G"
		num_proc:"4"
		maxRetries:3
		mnt_db_dir:mnt_db_dir
		mnt_input_dir:mnt_input_dir
		mnt_out_dir:mnt_out_dir
		task_queue:task_queue
	}
	
	output{
		String bam_out="${out_dir}/${taskID}/1_bwa2bam/${sam_name}/${sam_name}_bqsr_merged.bam"
		String log_file="${out_dir}/${taskID}/1_bwa2bam/${sam_name}/log"
	}
}

task BamQC_do {
	String taskID
	String sam_name
	String out_dir
	String input_bam
	String interval_bed
	String fasta
	
	String mnt_db_dir
	String mnt_input_dir
	String mnt_out_dir
	String task_queue
	
	####
	command <<<
		set -e
		if [ ! -d "${out_dir}/${taskID}/1_bwa2bam/${sam_name}/" ] ; then
			mkdir -p "${out_dir}/${taskID}/1_bwa2bam/${sam_name}/"
		fi
		cd "${out_dir}/${taskID}/1_bwa2bam/${sam_name}/"

		echo -e "${taskID} BamQC_start@"$(date +%Y-%m-%d_%H:%M:%S)"@"$(hostname) >> log
		
		####step_7 QC bam
		export PATH=/usr/local/bin/miniconda3/bin/:$PATH

		if [ ! -f bamqc.information.tsv ] ; then
			/usr/local/bin/NCbamInfo \
			-b ${interval_bed} \
			-r ${fasta} \
			-o ./ \
			-t \
			-p bamqc \
			${input_bam} \
			1>bamqc.information.tsv.stdout 2>bamqc.information.tsv.stderr
			if [[ $? -ne 0 && -f bamqc.information.tsv ]] ; then rm bamqc.information.tsv ; fi
		fi
		
		chmod 777 *
		
		echo -e "${taskID} BamQC_end@"$(date +%Y-%m-%d_%H:%M:%S)"@"$(hostname) >> log
	>>>
	
	runtime {
		docker:"oncowes/tools:v6"
		memory:"5G"
		num_proc:"4"
		maxRetries:3
		mnt_db_dir:mnt_db_dir
		mnt_input_dir:mnt_input_dir
		mnt_out_dir:mnt_out_dir
		task_queue:task_queue
	}
	
	output{
		String qc_out="${out_dir}/${taskID}/1_bwa2bam/${sam_name}/bamqc.information.tsv"
		String log_file="${out_dir}/${taskID}/1_bwa2bam/${sam_name}/log"
	}
}
