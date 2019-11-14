####input
####fqs list in format:
####            <taskID	sampleName	tissue_type	SM	LB	PU	fq1	fq2>

workflow fq2bam_base {
	####input_parameters
	File fq_conf_file
	String out_dir_name
	String platform
	
	String mnt_db_dir			####e.g.:/mnt/GeneCloud/luozw/oncoWES_db/
	String mnt_input_dir
	String mnt_out_dir
	String task_queue
	String out_dir=mnt_out_dir+"/"+out_dir_name
	
	String fasta_path_2=mnt_db_dir+"/hs37d5_ref/hs37d5.fa"
	String database_Mills_and_1000G_gold_standard_2=mnt_db_dir+"/bwa2bam_db/Mills_and_1000G_gold_standard.indels.b37.vcf"
	String databse_dbsnp_138_b37_2=mnt_db_dir+"/bwa2bam_db/dbsnp_138.b37.del100.vcf.gz"
	String database_known_1000G_indels_2=mnt_db_dir+"/bwa2bam_db/1000G_phase1.indels.b37.vcf.gz"

	Array [Array[String]] fq_pairs_info=read_tsv(fq_conf_file)
	
	scatter(fq_pair in fq_pairs_info){
		call bwa2bam_do{
			input:
				taskID=fq_pair[0],
				sam_name=fq_pair[1],
				tissue_type=fq_pair[2],
				SM=fq_pair[3],
				LB=fq_pair[4],
				PU=fq_pair[5],
				fq1=fq_pair[6],
				fq2=fq_pair[7],
				platform=platform,
				out_dir=out_dir,					
				fasta=fasta_path_2,
				mnt_db_dir=mnt_db_dir,
				mnt_input_dir=mnt_input_dir,
				mnt_out_dir=mnt_out_dir,
				task_queue=task_queue
		}
		
		call MarkDuplicates_do{
			input:
				taskID=fq_pair[0],
				sam_name=fq_pair[1],
				out_dir=out_dir,
				input_bam=bwa2bam_do.bam_out,
				mnt_db_dir=mnt_db_dir,
				mnt_input_dir=mnt_input_dir,
				mnt_out_dir=mnt_out_dir,
				task_queue=task_queue				
		}
		
		call BaseRecalibrator_do{
			input:
				taskID=fq_pair[0],
				sam_name=fq_pair[1],
				out_dir=out_dir,
				input_bam=MarkDuplicates_do.bam_out,
				ref_fasta=fasta_path_2,
				database_Mills_and_1000G_gold_standard=database_Mills_and_1000G_gold_standard_2,
				databse_dbsnp_138_b37=databse_dbsnp_138_b37_2,
				database_known_1000G_indels=database_known_1000G_indels_2,
				mnt_db_dir=mnt_db_dir,
				mnt_input_dir=mnt_input_dir,
				mnt_out_dir=mnt_out_dir,
				task_queue=task_queue
		}
		
		call ApplyBQSR_do{
			input:
				taskID=fq_pair[0],
				sam_name=fq_pair[1],
				out_dir=out_dir,
				input_bam=MarkDuplicates_do.bam_out,
				ref_fasta=fasta_path_2,
				input_table=BaseRecalibrator_do.table_out,
				mnt_db_dir=mnt_db_dir,
				mnt_input_dir=mnt_input_dir,
				mnt_out_dir=mnt_out_dir,
				task_queue=task_queue
		}
	}

	####
	output{
		Array[String] bam_out=ApplyBQSR_do.bam_out
		Array[String] work_stat=ApplyBQSR_do.log_file
	}
}


####################################
##sub tasks                       ##
####################################

task bwa2bam_do {
	String taskID
	String sam_name
	String tissue_type
	String SM
	String LB
	String PU
	String fq1
	String fq2
	String platform
	String out_dir
	
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

		echo -e "${taskID} fq2bam_work_start@"$(date +%Y-%m-%d_%H:%M:%S)"@"$(hostname)\(sample:${sam_name}\) >> log

		ID=${SM}_${LB}_${PU}
		####step1:run bwa
		if [ ! -f $ID"_fixmate.bam".SUCCESS ] ; then
			/usr/local/bin/miniconda3/bin/bwa mem -Y -M -R "@RG\tID:"$ID"\tSM:"${SM}"\tLB:"${LB}"\tPU:"${PU}"\tPL:"${platform} -t 16 -K 10000000 ${fasta} ${fq1} ${fq2} | \
			/usr/local/bin/miniconda3/bin/samtools fixmate - $ID"_fixmate.bam" \
			1>$ID"_fixmate.bam.stdout" 2>$ID"_fixmate.bam.stderr" \
			&& touch $ID"_fixmate.bam".SUCCESS
		fi
		
		####step2:sort bam
		if [ ! -f $ID"_sorted.bam".SUCCESS ] ; then
			/usr/local/bin/miniconda3/bin/samtools sort --reference ${fasta} -o $ID"_sorted.bam" -@ 16 $ID"_fixmate.bam" \
			1>$ID"_sorted.bam".stdout 2>$ID"_sorted.bam".stderr \
			&& touch $ID"_sorted.bam".SUCCESS \
			&& rm $ID"_fixmate.bam" 
		fi
		
		echo -e "${taskID} fq2bam_work_end@"$(date +%Y-%m-%d_%H:%M:%S)"@"$(hostname)\(sample:${sam_name}\) >> log
	>>>
	
	runtime {
		docker:"oncowes/tools:v6"
		memory:"10G"
		num_proc:"16"
		maxRetries:3
		mnt_db_dir:mnt_db_dir
		mnt_input_dir:mnt_input_dir
		mnt_out_dir:mnt_out_dir
		task_queue:task_queue
	}
	
	output{
		String bam_out="${out_dir}/${taskID}/1_bwa2bam/${sam_name}/${SM}_${LB}_${PU}_sorted.bam"
		String log_file="${out_dir}/${taskID}/1_bwa2bam/${sam_name}/log"
	}
}
task MarkDuplicates_do {
	String taskID
	String sam_name
	String out_dir
	String input_bam
	
	String mnt_db_dir
	String mnt_input_dir
	String mnt_out_dir
	String task_queue

	command <<<
		set -e

		if [ ! -d "${out_dir}/${taskID}/1_bwa2bam/${sam_name}/" ] ; then
			mkdir -p "${out_dir}/${taskID}/1_bwa2bam/${sam_name}/"
		fi
		cd "${out_dir}/${taskID}/1_bwa2bam/${sam_name}/"

		echo -e "${taskID} MarkDuplicates_work_start@"$(date +%Y-%m-%d_%H:%M:%S)"@"$(hostname) >> log

		####step_2:run MarkDuplicates
		if [ ! -f ${input_bam}_markDup.bam ] ; then
			gatk --java-options "-Xms10G -Xmx10G" \
			MarkDuplicates \
			--INPUT ${input_bam} \
			--OUTPUT ${input_bam}_markDup.bam \
			--METRICS_FILE ${input_bam}_markDup.metrics \
			--VALIDATION_STRINGENCY SILENT \
			--OPTICAL_DUPLICATE_PIXEL_DISTANCE 2500 \
			--ASSUME_SORT_ORDER "coordinate" \
			--CLEAR_DT false \
			--CREATE_MD5_FILE true \
			1>${input_bam}_markDup.bam.stdout 2>${input_bam}_markDup.bam.stderr
			if [[ $? -ne 0 && -f ${input_bam}_markDup.bam ]] ; then rm ${input_bam}_markDup.bam ; fi
		fi
		####
		
		echo -e "${taskID} MarkDuplicates_work_end@"$(date +%Y-%m-%d_%H:%M:%S)"@"$(hostname) >> log
	>>>
	
	runtime {
		docker:"oncowes/gatk4:latest"
		memory:"10G"
		num_proc:"4"
		maxRetries:3
		mnt_db_dir:mnt_db_dir
		mnt_input_dir:mnt_input_dir
		mnt_out_dir:mnt_out_dir
		task_queue:task_queue
	}
	
	output{
		String bam_out="${input_bam}_markDup.bam"
		String log_file="${out_dir}/${taskID}/1_bwa2bam/${sam_name}/log"
	}
}
	
task BaseRecalibrator_do {
	String taskID
	String sam_name
	String out_dir
	String input_bam
	String ref_fasta
	
	String database_Mills_and_1000G_gold_standard
	String databse_dbsnp_138_b37
	String database_known_1000G_indels
	
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

		echo -e "${taskID} BaseRecalibrator_work_start@"$(date +%Y-%m-%d_%H:%M:%S)"@"$(hostname) >> log

		####step_3:run BaseRecalibrator
		if [ ! -f ${input_bam}_bqsr.table ] ; then
			gatk --java-options "-Xms10G -Xmx10G" \
			BaseRecalibrator \
			-R ${ref_fasta} \
			-I ${input_bam} \
			--use-original-qualities \
			-O ${input_bam}_bqsr.table \
			--known-sites ${database_Mills_and_1000G_gold_standard} \
			--known-sites ${databse_dbsnp_138_b37} \
			--known-sites ${database_known_1000G_indels} \
			1>${input_bam}_bqsr.table.stdout 2>${input_bam}_bqsr.table.stderr
			if [[ $? -ne 0 && -f ${input_bam}_bqsr.table ]] ; then rm ${input_bam}_bqsr.table ; fi
		fi
		
		echo -e "${taskID} BaseRecalibrator_work_end@"$(date +%Y-%m-%d_%H:%M:%S)"@"$(hostname) >> log
	>>>
	
	runtime {
		docker:"oncowes/gatk4:latest"
		memory:"10G"
		num_proc:"4"
		maxRetries:3
		mnt_db_dir:mnt_db_dir
		mnt_input_dir:mnt_input_dir
		mnt_out_dir:mnt_out_dir
		task_queue:task_queue
	}
	
	output{
		String table_out="${input_bam}_bqsr.table"
		String log_file="${out_dir}/${taskID}/1_bwa2bam/${sam_name}/log"
	}
}
task ApplyBQSR_do {
	String taskID
	String sam_name
	String out_dir
	String input_bam
	String input_table
	String ref_fasta
	
	String mnt_db_dir
	String mnt_input_dir
	String mnt_out_dir
	String task_queue
	
	command <<<
		set -e
		
		if [ ! -d "${out_dir}/${taskID}/1_bwa2bam/${sam_name}/" ] ; then
			mkdir -p "${out_dir}/${taskID}/1_bwa2bam/${sam_name}/"
		fi
		cd "${out_dir}/${taskID}/1_bwa2bam/${sam_name}/"

		echo -e "${taskID} ApplyBQSR_work_start@"$(date +%Y-%m-%d_%H:%M:%S)"@"$(hostname) >> log
		
		####step_4:run ApplyBQSR
		if [ ! -f ${input_bam}_applyBQSR.bam ] ; then
			gatk --java-options "-Xms10G -Xmx10G" \
			ApplyBQSR \
			-R ${ref_fasta} \
			-I ${input_bam} \
			-bqsr ${input_table} \
			-O ${input_bam}_applyBQSR.bam \
			--add-output-sam-program-record \
			--use-original-qualities \
			1>${input_bam}_applyBQSR.bam.stdout 2>${input_bam}_applyBQSR.bam.stderr
			if [[ $? -ne 0 && -f ${input_bam}_applyBQSR.bam ]] ; then rm ${input_bam}_applyBQSR.bam ; fi
		fi

		echo -e "${taskID} ApplyBQSR_work_end@"$(date +%Y-%m-%d_%H:%M:%S)"@"$(hostname) >> log
	>>>
	
	runtime {
		docker:"oncowes/gatk4:latest"
		memory:"10G"
		num_proc:"4"
		maxRetries:3
		mnt_db_dir:mnt_db_dir
		mnt_input_dir:mnt_input_dir
		mnt_out_dir:mnt_out_dir
		task_queue:task_queue
	}
	
	output{
		String bam_out="${input_bam}_applyBQSR.bam"
		String log_file="${out_dir}/${taskID}/1_bwa2bam/${sam_name}/log"
	}
}
