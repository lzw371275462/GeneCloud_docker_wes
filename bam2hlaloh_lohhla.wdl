####envs:
####	1) docker:"oncowes/hlaloh:v3"
####	2) input file should be split by tab, not space


####input [bam list in format]:
####	<taskID	normalName	normalBam	tumorName	tumorBam	hlaOptiTypeTransResult>


workflow bam2hlaloh{
	File bam_lists
	String out_dir_name
	String? docker_user_name
	
	String mnt_db_dir                       ####e.g.:/mnt/GeneCloud/luozw/oncoWES_db/
	String mnt_input_dir
	String mnt_out_dir
	String task_queue
	String docker_user=select_first([docker_user_name,'$UID'])
	
	String out_dir=mnt_out_dir+"/"+out_dir_name
	
	####
	Array [Array[String]] bam_lists_infos=read_tsv(bam_lists)
	
	scatter(bam_lists_info in bam_lists_infos){
		call bam2hlaloh_do{
			input:
				task_ID=bam_lists_info[0],
				normal_name=bam_lists_info[1],
				normal_bam=bam_lists_info[2],			
				tumor_name=bam_lists_info[3],
				tumor_bam=bam_lists_info[4],
				hlaOptiType_result=bam_lists_info[5],
				out_dir_path=out_dir,
				mnt_db_dir=mnt_db_dir,
				mnt_input_dir=mnt_input_dir,
				mnt_out_dir=mnt_out_dir,
				task_queue=task_queue,
				docker_user=docker_user
		}
	}
	####
	output {
		Array[String] hlaloh_results=bam2hlaloh_do.hlaloh_result
		Array[String] hlaloh_log_out=bam2hlaloh_do.log_out
		Array[String] work_stat=bam2hlaloh_do.log_file
	}
}


task bam2hlaloh_do {
	String task_ID
	String normal_name
	String normal_bam
	String tumor_name
	String tumor_bam
	String hlaOptiType_result
	String out_dir_path
	
	String mnt_db_dir
	String mnt_input_dir
	String mnt_out_dir
	String task_queue
	String docker_user
	
	command <<<
		set -e
		
		if [ ! -d ${out_dir_path}/${task_ID}/2_bam2hlaLOH/ ] ; then 
			mkdir -p ${out_dir_path}/${task_ID}/2_bam2hlaLOH/
		fi
		cd ${out_dir_path}/${task_ID}/2_bam2hlaLOH/
		
		echo -e "${task_ID} bam2hlaloh_work_start@"$(date +%Y-%m-%d_%H:%M:%S)"@"$(hostname) > log

		outFile="${tumor_name}.10.DNA.HLAlossPrediction_CI.txt"

		export PATH=/opt/conda/bin/:$PATH
		####hla_loh
		if [ ! -f hlaloh.SUCCESS ] ; then
			/opt/conda/bin/Rscript /opt/lohhla/LOHHLAscript.R \
			--patientId ${tumor_name} \
			--CopyNumLoc FALSE \
			--outputDir ./ \
			--tumorBAMfile ${tumor_bam} \
			--hlaPath ${hlaOptiType_result} \
			--HLAfastaLoc /opt/lohhla/data/hlaFasta.fa \
			--mappingStep TRUE \
			--minCoverageFilter 10 \
			--fishingStep TRUE \
			--cleanUp TRUE  \
			--novoDir /opt/lohhla/bin/ \
			--HLAexonLoc  /opt/lohhla/data/hla.dat \
			--numMisMatch 1 \
			--plottingStep FALSE \
			--normalBAMfile ${normal_bam} \
			1>$outFile".stdout" 2>$outFile".stderr" \
			&& touch hlaloh.SUCCESS
		fi
		
		echo -e "${task_ID} bam2hlaloh_work_end@"$(date +%Y-%m-%d_%H:%M:%S)"@"$(hostname) >> log
		echo -e "${task_ID}\t${normal_name}\t${tumor_name}\t"$(readlink -f $outFile) > log.out
		
	>>>

	runtime {
		docker:"oncowes/hlaloh:v3"
		memory:"10G"
		num_proc:"8"
		maxRetries:3
		mnt_db_dir:mnt_db_dir
		mnt_input_dir:mnt_input_dir
		mnt_out_dir:mnt_out_dir
		task_queue:task_queue
		docker_user:docker_user
	}
	
	output{
		String hlaloh_result="${out_dir_path}/${task_ID}/2_bam2hlaLOH/${tumor_name}.10.DNA.HLAlossPrediction_CI.txt"
		String log_out="${out_dir_path}/${task_ID}/2_bam2hlaLOH/log.out"
		String log_file="${out_dir_path}/${task_ID}/2_bam2hlaLOH/log"
	}
}
