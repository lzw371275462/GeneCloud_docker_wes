####input
####normal_bam list in format:
####            <taskID	normalName	normalBam	tumorName	tumorBam>

workflow bam2cnvPON{
	File bam_list_file
	String out_dir_name
	String chip_bed_path2
	
	String mnt_db_dir                       ####e.g.:/mnt/GeneCloud/luozw/oncoWES_db/
	String mnt_input_dir
	String mnt_out_dir
	String task_queue

	String out_dir=mnt_out_dir+"/"+out_dir_name
	
	Boolean? no_BamQC_value1
	Boolean no_BamQC_value2=select_first([no_BamQC_value1,true])
	
	Array [Array[String]] bam_lists_infos=read_tsv(bam_list_file)
	
	scatter(bam_lists_info in bam_lists_infos){
		call CalculateTargetCoverage{
			input:
				taskID=bam_lists_info[0],
				normal_name=bam_lists_info[1],
				normal_bam=bam_lists_info[2],
				out_dir_path=out_dir,
				chip_bed_path=chip_bed_path2,
				mnt_db_dir=mnt_db_dir,
				mnt_input_dir=mnt_input_dir,
				mnt_out_dir=mnt_out_dir,
				task_queue=task_queue
		}
	}
	
	call bam2cnvPON_do{
		input:
			cov_infos=CalculateTargetCoverage.normalBamCov,
			out_dir_path=out_dir,
			chip_bed_path=chip_bed_path2,
			bamQC=no_BamQC_value2,
			mnt_db_dir=mnt_db_dir,
			mnt_input_dir=mnt_input_dir,
			mnt_out_dir=mnt_out_dir,
			task_queue=task_queue
	}

	output {
		String pon_db_file=bam2cnvPON_do.pon_db
		String work_stat=bam2cnvPON_do.log_file
	}
}

task CalculateTargetCoverage{
	String taskID
	String normal_name
	String normal_bam
	String out_dir_path
	String chip_bed_path
	String mnt_db_dir
	String mnt_input_dir
	String mnt_out_dir
	String task_queue


	command <<<
		set -e
		
		if [ ! -d ${out_dir_path}/pon_cnv/ ] ; then 
			mkdir -p ${out_dir_path}/pon_cnv/
		fi
		cd ${out_dir_path}/pon_cnv/
		
		echo -e "${taskID}_${normal_name} bam2ponCal_work_start@"$(date +%Y-%m-%d_%H:%M:%S)"@"$(hostname) > ${taskID}_${normal_name}.log
		
		export PATH=/usr/bin/:$PATH
		
		####CalculateTargetCoverage
		if [ ! -f ${normal_name}_normal.coverage.tsv ] ; then
			/gatk/gatk-launch --javaOptions  "-Xmx10G" CalculateTargetCoverage \
			--TMP_DIR ./${normal_name}_tmp_dir \
			-T ${chip_bed_path} \
			-I ${normal_bam} \
			-transform PCOV \
			-groupBy SAMPLE \
			-targetInfo FULL \
			--disableReadFilter NotDuplicateReadFilter \
			-O ${normal_name}_normal.coverage.tsv \
			1>${normal_name}_normal.coverage.tsv.stdout 2>${normal_name}_normal.coverage.tsv.stderr
			if [[ $? -ne 0 && -f ${normal_name}_normal.coverage.tsv ]] ; then rm ${normal_name}_normal.coverage.tsv ; fi
		fi
		sed -i 's/name\tnormal/name\t${normal_name}/g' ${normal_name}_normal.coverage.tsv

		echo -e "${taskID}_${normal_name} bam2ponCal_work_end@"$(date +%Y-%m-%d_%H:%M:%S)"@"$(hostname) >> ${taskID}_${normal_name}.log
		
	>>>
	
	runtime {
		docker: "oncowes/gatk4_beta1:latest"
		cpu: "1"
		memory:"10G"
		sge_queue:"all.q"
		num_proc:"8"
		maxRetries:3
		mnt_db_dir:mnt_db_dir
		mnt_input_dir:mnt_input_dir
		mnt_out_dir:mnt_out_dir
		task_queue:task_queue
	}
	
	output{
		String normalBamCov="${out_dir_path}/pon_cnv/${normal_name}_normal.coverage.tsv"
	}
}

task bam2cnvPON_do {
	Array[String] cov_infos
	String out_dir_path
	String chip_bed_path
	Boolean bamQC
	String mnt_db_dir
	String mnt_input_dir
	String mnt_out_dir
	String task_queue


	command <<<
		set -e
		
		if [ ! -d ${out_dir_path}/pon_cnv/ ] ; then 
			mkdir -p ${out_dir_path}/pon_cnv/
		fi
		
		cd ${out_dir_path}/pon_cnv/
		echo -e "bam2ponCreat_work_start@"$(date +%Y-%m-%d_%H:%M:%S)"@"$(hostname) > log
		less ${write_lines(cov_infos)} |sort -u > ./pon_cov.list
		
		####step_1: CombineReadCounts
		if [ ! -f pon_cov.list.mergedCov.tsv ] ; then
			/gatk/gatk-launch --javaOptions "-Xmx10G" CombineReadCounts \
			--TMP_DIR ./tmp_dir \
			-inputList pon_cov.list \
			-O pon_cov.list.mergedCov.tsv \
			-MOF 200 \
			1>pon_cov.list.mergedCov.tsv.stdout 2>pon_cov.list.mergedCov.tsv.stderr
			if [[ $? -ne 0 && -f pon_cov.list.mergedCov.tsv ]] ; then rm pon_cov.list.mergedCov.tsv ; fi
		fi
		
		####step_2: CreatePanelOfNormals
		####--noQC,-noQC:Boolean Skip the QC step.  PoN creation will be substantially faster, but greater risk of bad samples being introduced into the PoN.
		####CreatePanelOfNormals会自动QC过滤样本，当样本数太少时，可能都被过滤完了
		if [ ! -f pon_cov.list.mergedCov.tsv.pon ] ; then
			/gatk/gatk-launch --javaOptions "-Xmx10G" CreatePanelOfNormals \
			--TMP_DIR ./tmp_dir \
			-I pon_cov.list.mergedCov.tsv \
			-O pon_cov.list.mergedCov.tsv.pon \
			--disableSpark \
			--noQC ${bamQC} \
			1>pon_cov.list.mergedCov.tsv.pon.stdout 2>pon_cov.list.mergedCov.tsv.pon.stderr
			if [[ $? -ne 0 && -f pon_cov.list.mergedCov.tsv.pon ]] ; then rm pon_cov.list.mergedCov.tsv.pon ; fi
		fi
		echo -e "bam2ponCreat_work_end@"$(date +%Y-%m-%d_%H:%M:%S)"@"$(hostname) >> log
		
	>>>

	runtime {
		docker: "oncowes/gatk4_beta1:latest"
		cpu: "1"
		memory:"10G"
		sge_queue:"all.q"
		num_proc:"8"
		maxRetries:3
		mnt_db_dir:mnt_db_dir
		mnt_input_dir:mnt_input_dir
		mnt_out_dir:mnt_out_dir
		task_queue:task_queue
	}
	
	output{
		String pon_db="${out_dir_path}/pon_cnv/pon_cov.list.mergedCov.tsv.pon"
		String log_file="${out_dir_path}/pon_cnv/log"
	}
}
