####input
####bam list in format:
####            <taskID	normalName	normalBam	tumorName	tumorBam>

workflow bam2cnv{
	File bamlists_file
	String pondb
	String out_dir_name
	String chip_bed_path2
	
	String mnt_db_dir                       ####e.g.:/mnt/GeneCloud/luozw/oncoWES_db/
	String mnt_input_dir
	String mnt_out_dir
	String task_queue

	String out_dir=mnt_out_dir+"/"+out_dir_name

	String ref_fa_path2=mnt_db_dir+"/hs37d5_ref/hs37d5.fa"
	String inderval_path2=mnt_db_dir+"/hs37d5_ref/Broad.human.exome.b37.interval_list"
	String ref_dict_path2=mnt_db_dir+"/hs37d5_ref/hs37d5.dict"

	####
	Array [Array[String]] bam_lists_infos=read_tsv(bamlists_file)
	
	scatter(bam_lists_info in bam_lists_infos){
		
		####tumor_cov_tsv
		call CalculateTargetCoverage_do{
			input:
				task_ID=bam_lists_info[0],
				normal_name=bam_lists_info[1],
				normal_bam=bam_lists_info[2],			
				tumor_name=bam_lists_info[3],
				tumor_bam=bam_lists_info[4],
				out_dir_path=out_dir,
				chip_bed_path=chip_bed_path2,
				mnt_db_dir=mnt_db_dir,
				mnt_input_dir=mnt_input_dir,
				mnt_out_dir=mnt_out_dir,
				task_queue=task_queue
		}
		
		####tumor_ptn_tsv+tumor_tn_tsv
		call NormalizeSomaticReadCounts_do{
			input:
				task_ID=bam_lists_info[0],
				normal_name=bam_lists_info[1],
				normal_bam=bam_lists_info[2],			
				tumor_name=bam_lists_info[3],
				tumor_bam=bam_lists_info[4],
				out_dir_path=out_dir,
				pondb_path=pondb,
				tumor_cov_tsv=CalculateTargetCoverage_do.tumor_cov_tsv,
				mnt_db_dir=mnt_db_dir,
				mnt_input_dir=mnt_input_dir,
				mnt_out_dir=mnt_out_dir,
				task_queue=task_queue
		}
		
		####tumor_tn_seg
		call PerformSegmentation_do{
			input:
				task_ID=bam_lists_info[0],
				normal_name=bam_lists_info[1],
				normal_bam=bam_lists_info[2],			
				tumor_name=bam_lists_info[3],
				tumor_bam=bam_lists_info[4],
				out_dir_path=out_dir,
				tumor_tn_tsv=NormalizeSomaticReadCounts_do.tumor_tn_tsv,
				mnt_db_dir=mnt_db_dir,
				mnt_input_dir=mnt_input_dir,
				mnt_out_dir=mnt_out_dir,
				task_queue=task_queue
		}
		
		#### QC PLOT
		call PlotSegmentedCopyRatio_do{
			input:
				task_ID=bam_lists_info[0],
				normal_name=bam_lists_info[1],
				normal_bam=bam_lists_info[2],			
				tumor_name=bam_lists_info[3],
				tumor_bam=bam_lists_info[4],
				out_dir_path=out_dir,
				tumor_tn_tsv=NormalizeSomaticReadCounts_do.tumor_tn_tsv,
				tumor_ptn_tsv=NormalizeSomaticReadCounts_do.tumor_ptn_tsv,
				tumor_tn_seg=PerformSegmentation_do.tumor_tn_seg,
				ref_dict=ref_dict_path2,
				mnt_db_dir=mnt_db_dir,
				mnt_input_dir=mnt_input_dir,
				mnt_out_dir=mnt_out_dir,
				task_queue=task_queue
		}
		
		####tumor_tn_seg_calledseg
		call CallSegments_do{
			input:
				task_ID=bam_lists_info[0],
				normal_name=bam_lists_info[1],
				normal_bam=bam_lists_info[2],			
				tumor_name=bam_lists_info[3],
				tumor_bam=bam_lists_info[4],
				out_dir_path=out_dir,
				tumor_tn_tsv=NormalizeSomaticReadCounts_do.tumor_tn_tsv,
				tumor_tn_seg=PerformSegmentation_do.tumor_tn_seg,
				mnt_db_dir=mnt_db_dir,
				mnt_input_dir=mnt_input_dir,
				mnt_out_dir=mnt_out_dir,
				task_queue=task_queue
		}
		####tumor_Hets+normal_Hets
		call GetHetCoverage_do{
			input:
				task_ID=bam_lists_info[0],
				normal_name=bam_lists_info[1],
				normal_bam=bam_lists_info[2],			
				tumor_name=bam_lists_info[3],
				tumor_bam=bam_lists_info[4],
				out_dir_path=out_dir,
				inderval_path=inderval_path2,
				ref_fa=ref_fa_path2,
				mnt_db_dir=mnt_db_dir,
				mnt_input_dir=mnt_input_dir,
				mnt_out_dir=mnt_out_dir,
				task_queue=task_queue
		}
		
		####tumor_ACNVs
		call AllelicCNV_do{
			input:
				task_ID=bam_lists_info[0],
				normal_name=bam_lists_info[1],
				normal_bam=bam_lists_info[2],			
				tumor_name=bam_lists_info[3],
				tumor_bam=bam_lists_info[4],
				out_dir_path=out_dir,
				tumor_tn_seg_calledseg=CallSegments_do.tumor_tn_seg_calledseg,
				tumor_tn_tsv=NormalizeSomaticReadCounts_do.tumor_tn_tsv,
				tumor_Hets=GetHetCoverage_do.tumor_Hets,
				mnt_db_dir=mnt_db_dir,
				mnt_input_dir=mnt_input_dir,
				mnt_out_dir=mnt_out_dir,
				task_queue=task_queue
		}
		
		call PlotACNVResults_do{
			input:
				task_ID=bam_lists_info[0],
				normal_name=bam_lists_info[1],
				normal_bam=bam_lists_info[2],			
				tumor_name=bam_lists_info[3],
				tumor_bam=bam_lists_info[4],
				out_dir_path=out_dir,
				tumor_tn_tsv=NormalizeSomaticReadCounts_do.tumor_tn_tsv,
				tumor_Hets=GetHetCoverage_do.tumor_Hets,
				tumor_ACNV=AllelicCNV_do.tumor_ACNV,
				ref_dict=ref_dict_path2,
				mnt_db_dir=mnt_db_dir,
				mnt_input_dir=mnt_input_dir,
				mnt_out_dir=mnt_out_dir,
				task_queue=task_queue
		}
		
		call ConvertACNVResults_do{
			input:
				task_ID=bam_lists_info[0],
				normal_name=bam_lists_info[1],
				normal_bam=bam_lists_info[2],			
				tumor_name=bam_lists_info[3],
				tumor_bam=bam_lists_info[4],
				out_dir_path=out_dir,
				tumor_tn_tsv=NormalizeSomaticReadCounts_do.tumor_tn_tsv,
				tumor_Hets=GetHetCoverage_do.tumor_Hets,
				tumor_ACNV=AllelicCNV_do.tumor_ACNV,
				mnt_db_dir=mnt_db_dir,
				mnt_input_dir=mnt_input_dir,
				mnt_out_dir=mnt_out_dir,
				task_queue=task_queue
		}
	}
	####
	output {
		Array[String] cnv_out_infos=ConvertACNVResults_do.tumor_final_CNA
		Array[String] work_stat=ConvertACNVResults_do.log_file
	}
}


task CalculateTargetCoverage_do {
	String task_ID
	String normal_name
	String normal_bam
	String tumor_name
	String tumor_bam
	String out_dir_path
	String chip_bed_path

	String mnt_db_dir
	String mnt_input_dir
	String mnt_out_dir
	String task_queue
	
	command <<<
		set -e
		
		if [ ! -d ${out_dir_path}/${task_ID}/2_bam2cnv ] ; then 
			mkdir -p ${out_dir_path}/${task_ID}/2_bam2cnv
		fi
		cd ${out_dir_path}/${task_ID}/2_bam2cnv
		
		echo -e "${task_ID} CalculateTargetCoverage_work_start@"$(date +%Y-%m-%d_%H:%M:%S)"@"$(hostname) > log
		####step_1: CalculateTargetCoverage
		if [ ! -f ${tumor_name}_tumorCov.tsv ] ; then
			/gatk/gatk-launch --javaOptions "-Xms10G -Xmx10G" CalculateTargetCoverage \
			--TMP_DIR ./tmp_dir \
			-T ${chip_bed_path} \
			-I ${tumor_bam} \
			-transform PCOV \
			-groupBy SAMPLE \
			-targetInfo FULL \
			--disableReadFilter NotDuplicateReadFilter \
			-O ${tumor_name}_tumorCov.tsv \
			1>${tumor_name}_tumorCov.tsv.stdout 2>${tumor_name}_tumorCov.tsv.stderr
			if [[ $? -ne 0 && -f ${tumor_name}_tumorCov.tsv ]] ; then rm ${tumor_name}_tumorCov.tsv ; fi
		fi
		
		echo -e "${task_ID} CalculateTargetCoverage_work_end@"$(date +%Y-%m-%d_%H:%M:%S)"@"$(hostname) >> log
	>>>

	runtime {
		docker: "oncowes/gatk4_beta1:latest"
		cpu: "1"
		memory:"10G"
		num_proc:"16"
		maxRetries:3
		mnt_db_dir:mnt_db_dir
		mnt_input_dir:mnt_input_dir
		mnt_out_dir:mnt_out_dir
		task_queue:task_queue
	}
	
	output{
		String tumor_cov_tsv="${out_dir_path}/${task_ID}/2_bam2cnv/${tumor_name}_tumorCov.tsv"
	}
}

task NormalizeSomaticReadCounts_do {
	String task_ID
	String normal_name
	String normal_bam
	String tumor_name
	String tumor_bam
	String out_dir_path
	String tumor_cov_tsv
	String pondb_path
	
	String mnt_db_dir
	String mnt_input_dir
	String mnt_out_dir
	String task_queue
	
	command <<<
		set -e
		
		if [ ! -d ${out_dir_path}/${task_ID}/2_bam2cnv ] ; then 
			mkdir -p ${out_dir_path}/${task_ID}/2_bam2cnv
		fi
		cd ${out_dir_path}/${task_ID}/2_bam2cnv
		
		echo -e "${task_ID} NormalizeSomaticReadCounts_work_start@"$(date +%Y-%m-%d_%H:%M:%S)"@"$(hostname) >> log
		####step_2: NormalizeSomaticReadCounts
		if [[ ! -f ${tumor_name}_tumorCov.tn.tsv || ! -f ${tumor_name}_tumorCov.ptn.tsv ]] ; then
			/gatk/gatk-launch --javaOptions "-Xms10G -Xmx10G" NormalizeSomaticReadCounts \
			--TMP_DIR ./tmp_dir \
			-I ${tumor_cov_tsv} \
			-PON ${pondb_path} \
			-PTN ${tumor_name}_tumorCov.ptn.tsv \
			-TN ${tumor_name}_tumorCov.tn.tsv \
			1>${tumor_name}_tumorCov.tn.tsv.stdout 2>${tumor_name}_tumorCov.tn.tsv.stderr
			if [[ $? -ne 0 && -f ${tumor_name}_tumorCov.tn.tsv ]] ; then rm ${tumor_name}_tumorCov.tn.tsv ; fi
		fi
		
		echo -e "${task_ID} NormalizeSomaticReadCounts_work_end@"$(date +%Y-%m-%d_%H:%M:%S)"@"$(hostname) >> log
	>>>

	runtime {
		docker: "oncowes/gatk4_beta1:latest"
		cpu: "1"
		memory:"10G"
		num_proc:"16"
		maxRetries:3
		mnt_db_dir:mnt_db_dir
		mnt_input_dir:mnt_input_dir
		mnt_out_dir:mnt_out_dir
		task_queue:task_queue
	}
	
	output{
		String tumor_ptn_tsv="${out_dir_path}/${task_ID}/2_bam2cnv/${tumor_name}_tumorCov.ptn.tsv"
		String tumor_tn_tsv="${out_dir_path}/${task_ID}/2_bam2cnv/${tumor_name}_tumorCov.tn.tsv"
	}
}

task PerformSegmentation_do {
	String task_ID
	String normal_name
	String normal_bam
	String tumor_name
	String tumor_bam
	String out_dir_path
	String tumor_tn_tsv
	
	String mnt_db_dir
	String mnt_input_dir
	String mnt_out_dir
	String task_queue
	
	command <<<
		set -e
		
		if [ ! -d ${out_dir_path}/${task_ID}/2_bam2cnv ] ; then 
			mkdir -p ${out_dir_path}/${task_ID}/2_bam2cnv
		fi
		cd ${out_dir_path}/${task_ID}/2_bam2cnv
		
		echo -e "${task_ID} PerformSegmentation_work_start@"$(date +%Y-%m-%d_%H:%M:%S)"@"$(hostname) >> log
		####step_3: PerformSegmentation
		if [ ! -f ${tumor_name}_tumorCov.tn.tsv.seg ] ; then
			/gatk/gatk-launch --javaOptions "-Xms5G -Xmx5G" PerformSegmentation \
			--TMP_DIR ./tmp_dir \
			-TN ${tumor_tn_tsv} \
			-O ${tumor_name}_tumorCov.tn.tsv.seg \
			-LOG \
			1>${tumor_name}_tumorCov.tn.tsv.seg.stdout 2>${tumor_name}_tumorCov.tn.tsv.seg.stderr
			if [[ $? -ne 0 && -f ${tumor_name}_tumorCov.tn.tsv.seg ]] ; then rm ${tumor_name}_tumorCov.tn.tsv.seg ; fi
		fi
		
		echo -e "${task_ID} PerformSegmentation_work_end@"$(date +%Y-%m-%d_%H:%M:%S)"@"$(hostname) >> log
	>>>

	runtime {
		docker: "oncowes/gatk4_beta1:latest"
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
		String tumor_tn_seg="${out_dir_path}/${task_ID}/2_bam2cnv/${tumor_name}_tumorCov.tn.tsv.seg"
	}
}

task PlotSegmentedCopyRatio_do {
	String task_ID
	String normal_name
	String normal_bam
	String tumor_name
	String tumor_bam
	String out_dir_path
	String tumor_tn_tsv
	String tumor_ptn_tsv
	String tumor_tn_seg
	String ref_dict
	
	String mnt_db_dir
	String mnt_input_dir
	String mnt_out_dir
	String task_queue
	
	command <<<
		set -e
		
		if [ ! -d ${out_dir_path}/${task_ID}/2_bam2cnv ] ; then 
			mkdir -p ${out_dir_path}/${task_ID}/2_bam2cnv
		fi
		cd ${out_dir_path}/${task_ID}/2_bam2cnv
		
		echo -e "${task_ID} PlotSegmentedCopyRatio_work_start@"$(date +%Y-%m-%d_%H:%M:%S)"@"$(hostname) >> log
		if [ ! -d QC_plot ] ; then 
			mkdir -p QC_plot
		fi
		####step_4: PlotSegmentedCopyRatio
		if [ ! -f QC_plot/${tumor_name}_Before_After.png ] ; then
			/gatk/gatk-launch --javaOptions "-Xmx10G" PlotSegmentedCopyRatio \
			--TMP_DIR ./tmp_dir \
			-TN ${tumor_tn_tsv} \
			-PTN ${tumor_ptn_tsv} \
			-S  ${tumor_tn_seg} \
			-pre ${tumor_name} \
			-LOG \
			-SD ${ref_dict} \
			--output QC_plot/ \
			1>QC_plot.stdout 2>QC_plot.stderr
			if [[ $? -ne 0 && -f QC_plot/${tumor_name}_Before_After.png ]] ; then rm QC_plot/${tumor_name}_Before_After.png ; fi
		fi
		
		echo -e "${task_ID} PlotSegmentedCopyRatio_work_end@"$(date +%Y-%m-%d_%H:%M:%S)"@"$(hostname) >> log
	>>>

	runtime {
		docker: "oncowes/gatk4_beta1:latest"
		cpu: "1"
		memory:"10G"
		num_proc:"16"
		maxRetries:3
		mnt_db_dir:mnt_db_dir
		mnt_input_dir:mnt_input_dir
		mnt_out_dir:mnt_out_dir
		task_queue:task_queue
	}
	
	output{
		String SegmentedCopyRatio_plot="${out_dir_path}/${task_ID}/2_bam2cnv/QC_plot"
	}
}


task CallSegments_do {
	String task_ID
	String normal_name
	String normal_bam
	String tumor_name
	String tumor_bam
	String out_dir_path
	String tumor_tn_tsv
	String tumor_tn_seg
	
	String mnt_db_dir
	String mnt_input_dir
	String mnt_out_dir
	String task_queue
	
	command <<<
		set -e
		
		if [ ! -d ${out_dir_path}/${task_ID}/2_bam2cnv ] ; then 
			mkdir -p ${out_dir_path}/${task_ID}/2_bam2cnv
		fi
		cd ${out_dir_path}/${task_ID}/2_bam2cnv
		
		echo -e "${task_ID} CallSegments_work_start@"$(date +%Y-%m-%d_%H:%M:%S)"@"$(hostname) >> log
		export PATH=/usr/bin/:$PATH
		####step_5: CallSegments
		if [ ! -f ${tumor_name}_tumorCov.tn.tsv.called.seg ] ; then
			/gatk/gatk-launch --javaOptions "-Xmx10G" CallSegments \
			--TMP_DIR ./tmp_dir \
			-TN ${tumor_tn_tsv} \
			-S  ${tumor_tn_seg} \
			-O  ${tumor_name}_tumorCov.tn.tsv.called.seg \
			1>${tumor_name}_tumorCov.tn.tsv.called.seg.stdout 2>${tumor_name}_tumorCov.tn.tsv.called.seg.stderr
			if [[ $? -ne 0 && -f ${tumor_name}_tumorCov.tn.tsv.called.seg ]] ; then rm ${tumor_name}_tumorCov.tn.tsv.called.seg ; fi
		fi
		
		echo -e "${task_ID} CallSegments_work_end@"$(date +%Y-%m-%d_%H:%M:%S)"@"$(hostname) >> log
	>>>

	runtime {
		docker: "oncowes/gatk4_beta1:latest"
		cpu: "1"
		memory:"10G"
		num_proc:"16"
		maxRetries:3
		mnt_db_dir:mnt_db_dir
		mnt_input_dir:mnt_input_dir
		mnt_out_dir:mnt_out_dir
		task_queue:task_queue
	}
	
	output{
		String tumor_tn_seg_calledseg="${out_dir_path}/${task_ID}/2_bam2cnv/${tumor_name}_tumorCov.tn.tsv.called.seg"
	}
}

task GetHetCoverage_do {
	String task_ID
	String normal_name
	String normal_bam
	String tumor_name
	String tumor_bam
	String out_dir_path
	String inderval_path
	String ref_fa
	
	String mnt_db_dir
	String mnt_input_dir
	String mnt_out_dir
	String task_queue
	
	command <<<
		set -e
		
		if [ ! -d ${out_dir_path}/${task_ID}/2_bam2cnv ] ; then 
			mkdir -p ${out_dir_path}/${task_ID}/2_bam2cnv
		fi
		cd ${out_dir_path}/${task_ID}/2_bam2cnv
		
		echo -e "${task_ID} GetHetCoverage_work_start@"$(date +%Y-%m-%d_%H:%M:%S)"@"$(hostname) >> log
		####step_6: GetHetCoverage
		if [[ ! -f ${normal_name}_normal.Hets.txt || ! -f ${tumor_name}_tumor.Hets.txt ]] ; then
			/gatk/gatk-launch --javaOptions "-Xms10G -Xmx10G" GetHetCoverage \
			--TMP_DIR ./tmp_dir \
			-N ${normal_bam} \
			-NHET ${normal_name}_normal.Hets.txt \
			-T ${tumor_bam} \
			-THET ${tumor_name}_tumor.Hets.txt \
			-SNP ${inderval_path} \
			-R ${ref_fa} \
			1>${tumor_name}_tumor.Hets.txt.stdout 2>${tumor_name}_tumor.Hets.txt.stderr
			if [[ $? -ne 0 && -f ${normal_name}_normal.Hets.txt ]] ; then rm ${normal_name}_normal.Hets.txt ; fi
		fi
		
		echo -e "${task_ID} GetHetCoverage_work_end@"$(date +%Y-%m-%d_%H:%M:%S)"@"$(hostname) >> log
	>>>

	runtime {
		docker: "oncowes/gatk4_beta1:latest"
		cpu: "1"
		memory:"10G"
		num_proc:"16"
		maxRetries:3
		mnt_db_dir:mnt_db_dir
		mnt_input_dir:mnt_input_dir
		mnt_out_dir:mnt_out_dir
		task_queue:task_queue
	}
	
	output{
		String tumor_Hets="${out_dir_path}/${task_ID}/2_bam2cnv/${tumor_name}_tumor.Hets.txt"
		String normal_Hets="${out_dir_path}/${task_ID}/2_bam2cnv/${normal_name}_normal.Hets.txt"
	}
}

task AllelicCNV_do {
	String task_ID
	String normal_name
	String normal_bam
	String tumor_name
	String tumor_bam
	String out_dir_path
	String tumor_tn_seg_calledseg
	String tumor_tn_tsv
	String tumor_Hets
	
	String mnt_db_dir
	String mnt_input_dir
	String mnt_out_dir
	String task_queue
	
	command <<<
		set -e
		
		if [ ! -d ${out_dir_path}/${task_ID}/2_bam2cnv ] ; then 
			mkdir -p ${out_dir_path}/${task_ID}/2_bam2cnv
		fi
		cd ${out_dir_path}/${task_ID}/2_bam2cnv
		
		echo -e "${task_ID} AllelicCNV_work_start@"$(date +%Y-%m-%d_%H:%M:%S)"@"$(hostname) >> log
		####step_7: AllelicCNV
		if [ ! -f ${tumor_name}.ACNV-sim-final.seg ] ; then
			/gatk/gatk-launch --javaOptions "-Xms5G -Xmx5G" AllelicCNV \
			--TMP_DIR ./tmp_dir \
			-S ${tumor_tn_seg_calledseg} \
			-TN ${tumor_tn_tsv} \
			-THET ${tumor_Hets} \
			-pre ${tumor_name}.ACNV \
			1>${tumor_name}.ACNV-sim-final.seg.stdout 2>${tumor_name}.ACNV-sim-final.seg.stderr
			if [[ $? -ne 0 && -f ${tumor_name}.ACNV-sim-final.seg ]] ; then rm ${tumor_name}.ACNV-sim-final.seg ; fi
		fi
		
		echo -e "${task_ID} AllelicCNV_work_end@"$(date +%Y-%m-%d_%H:%M:%S)"@"$(hostname) >> log
	>>>

	runtime {
		docker: "oncowes/gatk4_beta1:latest"
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
		String tumor_ACNV="${out_dir_path}/${task_ID}/2_bam2cnv/${tumor_name}.ACNV-sim-final.seg"
	}
}


task PlotACNVResults_do {
	String task_ID
	String normal_name
	String normal_bam
	String tumor_name
	String tumor_bam
	String out_dir_path
	String tumor_tn_tsv
	String tumor_Hets
	String tumor_ACNV
	String ref_dict
	
	String mnt_db_dir
	String mnt_input_dir
	String mnt_out_dir
	String task_queue
	
	command <<<
		set -e
		
		if [ ! -d ${out_dir_path}/${task_ID}/2_bam2cnv ] ; then 
			mkdir -p ${out_dir_path}/${task_ID}/2_bam2cnv
		fi
		cd ${out_dir_path}/${task_ID}/2_bam2cnv
		
		echo -e "${task_ID} PlotACNVResults_work_start@"$(date +%Y-%m-%d_%H:%M:%S)"@"$(hostname) >> log
		if [ ! -d cnvResult_PLOT ] ; then
				mkdir -p cnvResult_PLOT
		fi
		####step_8: PlotACNVResults
		if [ ! -f cnvResult_PLOT/${tumor_name}.PlotACNV_ACNV.png ] ; then
			/gatk/gatk-launch --javaOptions "-Xms5G -Xmx5G" PlotACNVResults \
			--TMP_DIR ./tmp_dir \
			--hets ${tumor_Hets} \
			--tangentNormalized ${tumor_tn_tsv} \
			--segments ${tumor_ACNV} \
			-SD ${ref_dict} \
			--output cnvResult_PLOT/ \
			--outputPrefix ${tumor_name}.PlotACNV \
			1>cnvResult_PLOT.stdout 2>cnvResult_PLOT.stderr
			if [[ $? -ne 0 && -f cnvResult_PLOT/${tumor_name}.PlotACNV_ACNV.png ]] ; then rm cnvResult_PLOT/${tumor_name}.PlotACNV_ACNV.png ; fi
		fi
		
		echo -e "${task_ID} PlotACNVResults_work_end@"$(date +%Y-%m-%d_%H:%M:%S)"@"$(hostname) >> log
	>>>

	runtime {
		docker: "oncowes/gatk4_beta1:latest"
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
		String PlotACNVResults_plot="${out_dir_path}/${task_ID}/2_bam2cnv/cnvResult_PLOT"
	}
}


task ConvertACNVResults_do {
	String task_ID
	String normal_name
	String normal_bam
	String tumor_name
	String tumor_bam
	String out_dir_path
	String tumor_tn_tsv
	String tumor_Hets
	String tumor_ACNV
	
	String mnt_db_dir
	String mnt_input_dir
	String mnt_out_dir
	String task_queue
	
	command <<<
		set -e
		
		if [ ! -d ${out_dir_path}/${task_ID}/2_bam2cnv ] ; then 
			mkdir -p ${out_dir_path}/${task_ID}/2_bam2cnv
		fi
		cd ${out_dir_path}/${task_ID}/2_bam2cnv
		
		echo -e "${task_ID} ConvertACNVResults_work_start@"$(date +%Y-%m-%d_%H:%M:%S)"@"$(hostname) >> log
		####step_9: ConvertACNVResults
		if [ ! -f ${tumor_name}.ACNV-sim-final.cnv.seg ] ; then
			/gatk/gatk-launch --javaOptions "-Xmx5G -Xmx5G" ConvertACNVResults \
			--tumorHets ${tumor_Hets} \
			--tangentNormalized ${tumor_tn_tsv} \
			--segments ${tumor_ACNV} \
			--outputDir ./ \
			1>${tumor_name}.ACNV-sim-final.cnv.seg.stdout 2>${tumor_name}.ACNV-sim-final.cnv.seg.stderr
			if [[ $? -ne 0 && -f ${tumor_name}.ACNV-sim-final.cnv.seg ]] ; then rm ${tumor_name}.ACNV-sim-final.cnv.seg ; fi
		fi
		echo -e "${task_ID} ConvertACNVResults_work_end@"$(date +%Y-%m-%d_%H:%M:%S)"@"$(hostname) >> log
	
		echo -e "${task_ID}\t${normal_name}\t${tumor_name}\t"$(readlink -f ${tumor_name}.ACNV-sim-final.cnv.seg) >log.out
		chmod 777 ./*
		
	>>>

	runtime {
		docker: "oncowes/gatk4_beta1:latest"
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
		String tumor_final_CNA="${out_dir_path}/${task_ID}/2_bam2cnv/${tumor_name}.ACNV-sim-final.cnv.seg"
		String log_file="${out_dir_path}/${task_ID}/2_bam2cnv/log"
	}
}


