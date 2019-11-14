####input
####bam list in format:
####            <taskID	normalName	normalBam	tumorName	tumorBam>


workflow bam2vcf{
	File bam_lists
	String out_dir_name
	
	String mnt_db_dir
	String mnt_input_dir
	String mnt_out_dir
	String task_queue

	String out_dir=mnt_out_dir+"/"+out_dir_name

	String? PON
        String chip_bed_path2=mnt_db_dir+"/chip_bed/qc.bed"
	String ref_path_2=mnt_db_dir+"/hs37d5_ref/hs37d5.fa"
	String gnomAD_path_2=mnt_db_dir+"/vcf_Mutect2_test/af-only-gnomad.raw.sites.b37.vcf.gz"
	String small_exac_common_path_2=mnt_db_dir+"/vcf_Mutect2_test/small_exac_common_3_b37.vcf.gz"
	Array [Array[String]] bam_lists_infos=read_tsv(bam_lists)
	
	scatter(bam_lists_info in bam_lists_infos){
		call getpileupsummaries_do{
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
				task_queue=task_queue,
				bed_path=chip_bed_path2,
				small_exac_common_path=small_exac_common_path_2
		}

		call calculatecontamination_do{
			input:
				task_ID=bam_lists_info[0],
				normal_name=bam_lists_info[1],
				normal_bam=bam_lists_info[2],
				tumor_name=bam_lists_info[3],
				tumor_bam=bam_lists_info[4],
				tumor_getpileupsummaries=getpileupsummaries_do.tumor_getpileupsummaries,
				normal_getpileupsummaries=getpileupsummaries_do.normal_getpileupsummaries,
				out_dir_path=out_dir,
				mnt_db_dir=mnt_db_dir,
				mnt_input_dir=mnt_input_dir,
				mnt_out_dir=mnt_out_dir,
				task_queue=task_queue
		}

	####
		call Mutect2_do{
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
				task_queue=task_queue,
				PON_path=PON,
				bed_path=chip_bed_path2,
				ref_path=ref_path_2,
				gnomAD_path=gnomAD_path_2
		}

		call LearnReadOrientationModel_FilterMutectCalls_do{
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
				task_queue=task_queue,
				ref_path=ref_path_2,
				f1r2_tar_gz=Mutect2_do.f1r2_tar_gz,
				somatic_raw_vcf_gz=Mutect2_do.somatic_raw_vcf_gz,
				tumol_calculatecontamination_table=calculatecontamination_do.tumor_calculatecontamination_table,
				somatic_raw_vcf_gz_stats=Mutect2_do.somatic_raw_vcf_gz_stats,
				somatic_raw_vcf_gz_tbi=Mutect2_do.somatic_raw_vcf_gz_tbi,
				segments_table=calculatecontamination_do.tumor_calculatecontamination_table
		}
		
		call filter_do{
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
				task_queue=task_queue,
				somatic_filtered_vcf_gz=LearnReadOrientationModel_FilterMutectCalls_do.somatic_filtered_vcf_gz
		}
	}
		
	output {
		####Array[Array[String]] vcf_out_infos=bam2vcf_do.vcf_out_info
		Array[String] filter_out_infos=filter_do.log_out
		Array[String] work_stat=filter_do.log_file
	}
}


task getpileupsummaries_do {
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

	String bed_path
	String small_exac_common_path
	
	command <<<
		set -e
		
		if [ ! -d ${out_dir_path}/${task_ID}/2_bam2vcf/tmp ] ; then 
			mkdir -p ${out_dir_path}/${task_ID}/2_bam2vcf/tmp
		fi
		
		cd ${out_dir_path}/${task_ID}/2_bam2vcf/

		echo -e "${task_ID} getpileupsummaries_work_start@"$(date +%Y-%m-%d_%H:%M:%S)"@"$(hostname) > log
		
		####step_1 tumor_getpileupsummaries
		if [ ! -f tumor_getpileupsummaries.table ] ; then 
			gatk GetPileupSummaries \
			-I ${tumor_bam} \
			-V ${small_exac_common_path} \
			-L ${bed_path} \
			-O ./tumor_getpileupsummaries.table \
			1>tumor_getpileupsummaries.table.stdout 2>tumor_getpileupsummaries.table.stderr
			 if [[ $? -ne 0 && -f tumor_getpileupsummaries.table ]] ; then rm tumor_getpileupsummaries.table ; fi
		fi
		
		####step_2 normal_getpileupsummaries
		if [ ! -f normal_getpileupsummaries.table ] ; then
			gatk GetPileupSummaries \
			-I ${normal_bam} \
			-V ${small_exac_common_path} \
			-L ${bed_path} \
			-O ./normal_getpileupsummaries.table \
			1>normal_getpileupsummaries.table.stdout 2>normal_getpileupsummaries.table.stderr
			if [[ $? -ne 0 && -f normal_getpileupsummaries.table ]] ; then rm normal_getpileupsummaries.table ; fi
		fi

		echo -e "${task_ID} getpileupsummaries_work_end@"$(date +%Y-%m-%d_%H:%M:%S)"@"$(hostname) >> log
	>>>
	runtime {
		docker:"oncowes/gatk4:latest"
		cpu:"1"
		memory:"5G"
		num_proc:"2"
		maxRetries:3
		mnt_db_dir:mnt_db_dir
		mnt_input_dir:mnt_input_dir
		mnt_out_dir:mnt_out_dir
		task_queue:task_queue
	}

	output{
		String tumor_getpileupsummaries="${out_dir_path}/${task_ID}/2_bam2vcf/tumor_getpileupsummaries.table"
		String normal_getpileupsummaries="${out_dir_path}/${task_ID}/2_bam2vcf/normal_getpileupsummaries.table"
	}
}

task calculatecontamination_do {
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
	String tumor_getpileupsummaries
	String normal_getpileupsummaries
	
	command <<<
		set -e

		if [ ! -d ${out_dir_path}/${task_ID}/2_bam2vcf/tmp ] ; then
			mkdir -p ${out_dir_path}/${task_ID}/2_bam2vcf/tmp
		fi

		cd ${out_dir_path}/${task_ID}/2_bam2vcf/

		echo -e "${task_ID} calculatecontamination_work_start@"$(date +%Y-%m-%d_%H:%M:%S)"@"$(hostname) >> log

		####step_3 tumorl_calculatecontamination
		if [[ ! -f tumor_calculatecontamination.table || ! -f segments.table ]] ; then
			gatk CalculateContamination \
			-I ${tumor_getpileupsummaries} \
			--matched ${normal_getpileupsummaries} \
			-O ./tumor_calculatecontamination.table \
			--tumor-segmentation segments.table \
			1>tumor_calculatecontamination.table.stdout 2>tumor_calculatecontamination.table.stderr
			if [[ $? -ne 0 && -f tumor_calculatecontamination.table ]] ; then rm tumor_calculatecontamination.table ; fi
		fi
		
		####step_4 normal_calculatecontamination
		if [[ ! -f normal_calculatecontamination.table || ! -f segments.table ]] ; then
			gatk CalculateContamination \
			-I ${normal_getpileupsummaries} \
			-O ./normal_calculatecontamination.table \
			1>normal_calculatecontamination.table.stdout 2>normal_calculatecontamination.table.stderr
			if [[ $? -ne 0 && -f normal_calculatecontamination.table ]] ; then rm normal_calculatecontamination.table ; fi
		fi
		echo -e "${task_ID} calculatecontamination_work_end@"$(date +%Y-%m-%d_%H:%M:%S)"@"$(hostname) >> log
	>>>
	
	runtime {
		docker:"oncowes/gatk4:latest"
		cpu: "1"
		memory:"5G"
		num_proc:"2"
		maxRetries:3
		mnt_db_dir:mnt_db_dir
		mnt_input_dir:mnt_input_dir
		mnt_out_dir:mnt_out_dir
		task_queue:task_queue
	}

	output{
		String tumor_calculatecontamination_table="${out_dir_path}/${task_ID}/2_bam2vcf/tumor_calculatecontamination.table"
		String normal_calculatecontamination_table="${out_dir_path}/${task_ID}/2_bam2vcf/normal_calculatecontamination.table"
		String segments_table="${out_dir_path}/${task_ID}/2_bam2vcf/segments.table"
	}
}

task Mutect2_do{
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
	String? PON_path
	String bed_path
	String ref_path
	String gnomAD_path

        command <<<
                set -e

                if [ ! -d ${out_dir_path}/${task_ID}/2_bam2vcf/tmp ] ; then
                        mkdir -p ${out_dir_path}/${task_ID}/2_bam2vcf/tmp
                fi

                cd ${out_dir_path}/${task_ID}/2_bam2vcf/

                echo -e "${task_ID} Mutect2_work_start@"$(date +%Y-%m-%d_%H:%M:%S)"@"$(hostname) >> log

		####step_5 Mutect2
		if [[ ! -f somatic_raw.vcf.gz || ! -f f1r2.tar.gz || ! -f tumor_normal_m2.bam ]] ; then
			gatk Mutect2 --java-options "-Djava.io.tmpdir=./tmp -Xms16g -Xmx16g" --native-pair-hmm-threads 8 \
			-R ${ref_path}  \
			-I ${tumor_bam} \
			-I ${normal_bam} \
			-O ./somatic_raw.vcf.gz \
			-L ${bed_path} \
			--germline-resource ${gnomAD_path} \
			-normal ${normal_name} \
			--af-of-alleles-not-in-resource 0.0000025 \
			-bamout tumor_normal_m2.bam \
			-A ClippingRankSumTest \
			--f1r2-tar-gz f1r2.tar.gz \
			${"-pon "+PON_path} \
			--tmp-dir ./tmp \
			1>somatic_raw.vcf.gz.stdout 2>somatic_raw.vcf.gz.stderr
			if [[ $? -ne 0 && -f somatic_raw.vcf.gz ]] ; then rm somatic_raw.vcf.gz ; fi
		fi
		echo -e "${task_ID} Mutect2_work_end@"$(date +%Y-%m-%d_%H:%M:%S)"@"$(hostname) >> log
        >>>

        runtime {
                docker:"oncowes/gatk4:latest"
                cpu: "1"
                memory:"16G"
                num_proc:"8"
                maxRetries:8
                mnt_db_dir:mnt_db_dir
                mnt_input_dir:mnt_input_dir
                mnt_out_dir:mnt_out_dir
                task_queue:task_queue
        }

        output{
                String somatic_raw_vcf_gz="${out_dir_path}/${task_ID}/2_bam2vcf/somatic_raw.vcf.gz"
		String somatic_raw_vcf_gz_tbi="${out_dir_path}/${task_ID}/2_bam2vcf/somatic_raw.vcf.gz.tbi"
		String f1r2_tar_gz="${out_dir_path}/${task_ID}/2_bam2vcf/f1r2.tar.gz"
		String somatic_raw_vcf_gz_stats="${out_dir_path}/${task_ID}/2_bam2vcf/somatic_raw.vcf.gz.stats"
		
        }
}

task LearnReadOrientationModel_FilterMutectCalls_do {
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

        String ref_path
	String f1r2_tar_gz
	String somatic_raw_vcf_gz
	String somatic_raw_vcf_gz_tbi
	String tumol_calculatecontamination_table
	String somatic_raw_vcf_gz_stats
	String segments_table
        command <<<
                set -e

                if [ ! -d ${out_dir_path}/${task_ID}/2_bam2vcf/tmp ] ; then
                        mkdir -p ${out_dir_path}/${task_ID}/2_bam2vcf/tmp
                fi

                cd ${out_dir_path}/${task_ID}/2_bam2vcf/

                echo -e "${task_ID} LearnReadOrientationModel_FilterMutectCalls_work_start@"$(date +%Y-%m-%d_%H:%M:%S)"@"$(hostname) >> log
		
		####step_6 LearnReadOrientationModel
		if [ ! -f artifact-priors.tar.gz ] ; then
			gatk LearnReadOrientationModel \
			-I ${f1r2_tar_gz} \
			-O ./artifact-priors.tar.gz \
			1>artifact-priors.tar.gz.stdout 2>artifact-priors.tar.gz.stderr
			if [[ $? -ne 0 && -f artifact-priors.tar.gz ]] ; then rm artifact-priors.tar.gz ; fi
		fi
		
		####step_7 FilterMutectCalls
		if [[ ! -f somatic_filtered.vcf.gz || ! -f artifact-priors.tar.gz ]] ; then
			gatk FilterMutectCalls \
			-V ${somatic_raw_vcf_gz} \
			-R ${ref_path}  \
			--contamination-table ${tumol_calculatecontamination_table} \
			--stats ${somatic_raw_vcf_gz_stats} \
			-O ./somatic_filtered.vcf.gz \
			--orientation-bias-artifact-priors artifact-priors.tar.gz \
			--tumor-segmentation ${segments_table} \
			--filtering-stats FMC.filtering.stats \
			1>somatic_filtered.vcf.gz.stdout 2>somatic_filtered.vcf.gz.stderr
			if [[ $? -ne 0 && -f somatic_filtered.vcf.gz ]] ; then rm somatic_filtered.vcf.gz ; fi
		fi
		echo -e "${task_ID} LearnReadOrientationModel_FilterMutectCalls_work_end@"$(date +%Y-%m-%d_%H:%M:%S)"@"$(hostname) >> log
        >>>

        runtime {
                docker:"oncowes/gatk4:latest"
                cpu: "1"
                memory:"5G"
                num_proc:"2"
                maxRetries:3
                mnt_db_dir:mnt_db_dir
                mnt_input_dir:mnt_input_dir
                mnt_out_dir:mnt_out_dir
                task_queue:task_queue
        }

        output{
               String somatic_filtered_vcf_gz="${out_dir_path}/${task_ID}/2_bam2vcf/somatic_filtered.vcf.gz"
        }
}

task filter_do {
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
        String somatic_filtered_vcf_gz

        command <<<
                set -e

                if [ ! -d ${out_dir_path}/${task_ID}/2_bam2vcf/tmp ] ; then
                        mkdir -p ${out_dir_path}/${task_ID}/2_bam2vcf/tmp
                fi

                cd ${out_dir_path}/${task_ID}/2_bam2vcf/

                echo -e "${task_ID} filter_work_start@"$(date +%Y-%m-%d_%H:%M:%S)"@"$(hostname) >> log
		
		####step_8 filter_PASS

		rm -rf ./tmp/
		
		if [ ! -f somatic_filtered_PASS.vcf ] ; then
			bcftools view -i 'FILTER="PASS"' ${somatic_filtered_vcf_gz} > somatic_filtered_PASS.vcf
			if [[ $? -ne 0 && -f somatic_filtered_PASS.vcf ]] ; then rm somatic_filtered_PASS.vcf ; fi
		fi
		
		echo -e "${task_ID} filter_work_end@"$(date +%Y-%m-%d_%H:%M:%S)"@"$(hostname) >> log		
		echo -e "${task_ID}\t${normal_name}\t${tumor_name}\t"$(readlink -f somatic_filtered_PASS.vcf) >log.out
		
	>>>

	runtime {
                docker:"oncowes/gatk4:latest"
                cpu: "1"
                memory:"2G"
                num_proc:"2"
                maxRetries:3
                mnt_db_dir:mnt_db_dir
                mnt_input_dir:mnt_input_dir
                mnt_out_dir:mnt_out_dir
                task_queue:task_queue
	}
	
	output{
		String log_out="${out_dir_path}/${task_ID}/2_bam2vcf/log.out"
		String log_file="${out_dir_path}/${task_ID}/2_bam2vcf/log"
	}
}
