####input
####bam list in format:
####	<taskID		music2_bam2geneCov_list		 maf_file>

####music2_bam2geneCov_list in format:
####	<taskID normalName      normalBam       tumorName       tumorBam	music2_bam2geneCov_result>


workflow maf2MutSignatureGenes {
	File geneCov_list_file
	String out_dir_name="maf2MutSignatureGenes_outdir"
	String fdr_value=0.05
	
	String mnt_db_dir="/mnt/GeneCloud/luozw/oncoWES_db/"
	String mnt_input_dir
	String mnt_out_dir
	String task_queue="all.q"
	String docker_user='$UID'
	String out_dir=mnt_out_dir+"/"+out_dir_name

	String all_gene_tsv=mnt_db_dir+"/maf2signatureGene_db/all_gene.tsv"
	String ref_fa=mnt_db_dir+"/hs37d5_ref/hs37d5.fa"
	String oncoplot_script=mnt_db_dir+"/maf2signatureGene_db/wes_maftools_oncoplot.R"
	String smgsFormat_script=mnt_db_dir+"/maf2signatureGene_db/smgsFormat.pl"
	String Interact_script=mnt_db_dir+"/maf2signatureGene_db/maftools_somaticInteractions.R"
	
	
	Array [Array[String]] geneCov_infos=read_tsv(geneCov_list_file)
	scatter(geneCov_info in geneCov_infos){
		call get_geneCoverage_do {
			input:
				taskName=geneCov_info[0],
				cov_files=geneCov_info[1],
				all_gene_tsv=all_gene_tsv,
				ref_fa=ref_fa,
				out_dir=out_dir,
				mnt_db_dir=mnt_db_dir,
				mnt_input_dir=mnt_input_dir,
				mnt_out_dir=mnt_out_dir,
				task_queue=task_queue,
				docker_user=docker_user
		}
		
		call gene_mutation_rates_do {
			input:
				taskName=geneCov_info[0],
				maf_file=geneCov_info[2],
				all_gene_tsv=all_gene_tsv,
				ref_fa=ref_fa,
				out_dir=out_dir,
				bamlist=get_geneCoverage_do.bamlist,
				mnt_db_dir=mnt_db_dir,
				mnt_input_dir=mnt_input_dir,
				mnt_out_dir=mnt_out_dir,
				task_queue=task_queue,
				docker_user=docker_user
		}
		call test_by_fdr_do {
			input:
				taskName=geneCov_info[0],
				gene_mrs=gene_mutation_rates_do.gene_mrs,
				fdr_value=fdr_value,
				out_dir=out_dir,
				mnt_db_dir=mnt_db_dir,
				mnt_input_dir=mnt_input_dir,
				mnt_out_dir=mnt_out_dir,
				task_queue=task_queue,
				docker_user=docker_user
		}
		call plot_SigMutGenes_do {
			input:
				taskName=geneCov_info[0],
				maf_file=geneCov_info[2],
				smgs_simple=test_by_fdr_do.smgs_simple,
				oncoplot_script=oncoplot_script,
				smgsFormat_script=smgsFormat_script,
				Interact_script=Interact_script,
				out_dir=out_dir,
				mnt_db_dir=mnt_db_dir,
				mnt_input_dir=mnt_input_dir,
				mnt_out_dir=mnt_out_dir,
				task_queue=task_queue,
				docker_user=docker_user
		}
	}
	
	output{
		Array[String] SigMutGenes_results=plot_SigMutGenes_do.SigMutGenes_result
		Array[String] log=plot_SigMutGenes_do.log_file
	}
}

task get_geneCoverage_do {
	String taskName
	String cov_files
	String all_gene_tsv
	String ref_fa
	String out_dir
	String mnt_db_dir
	String mnt_input_dir
	String mnt_out_dir
	String task_queue
	String docker_user

	command <<<
		set -e
		if [ ! -d ${out_dir}/${taskName}/roi_covgs ] ; then
				mkdir -p  ${out_dir}/${taskName}/roi_covgs
		fi
		cd ${out_dir}/${taskName}/
		echo -e "${taskName} get_geneCoverage_do_work_start@"$(date +%Y-%m-%d_%H:%M:%S)"@"$(hostname) > log
		
		####step1: gene cov results summary
		cat ${cov_files}|cut -f 6|while read a ; do
			file_name=$(basename $a)
			if [ ! -f roi_covgs/$file_name ] ; then
				ln -s $a roi_covgs/
			fi
		done
		cat ${cov_files}|awk '{print $4"\t"$3"\t"$5}'  > bam.list
		
		export PATH=/usr/local/bin/miniconda3/bin/:$PATH
		####step1.Run bmr calc-covg again to get gene coverage
		if [ ! -f get_geneCoverage_do.SUCCESS ] ; then 
			music2 bmr calc-covg \
			--roi-file ${all_gene_tsv} \
			--reference-sequence ${ref_fa} \
			--bam-list ./bam.list \
			--output-dir ${out_dir} \
			1>get_geneCoverage_do.stdout 2>get_geneCoverage_do.stderr \
			&& touch get_geneCoverage_do.SUCCESS
		fi
		echo -e "${taskName} get_geneCoverage_do_work_end@"$(date +%Y-%m-%d_%H:%M:%S)"@"$(hostname) >> log
	>>>

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
		String bamlist ="${out_dir}/${taskName}/bam.list"
		String log_file="${out_dir}/${taskName}/log"
	}
}



task gene_mutation_rates_do {
	String taskName
	String maf_file
	String all_gene_tsv
	String ref_fa
	String bamlist
	
	String out_dir
	String mnt_db_dir
	String mnt_input_dir
	String mnt_out_dir
	String task_queue
	String docker_user

	command <<<
		set -e
		if [ ! -d ${out_dir}/${taskName}/ ] ; then
				mkdir -p  ${out_dir}/${taskName}/
		fi
		cd ${out_dir}/${taskName}/
		echo -e "${taskName} gene_mutation_rates_do_work_start@"$(date +%Y-%m-%d_%H:%M:%S)"@"$(hostname) >> log

		export PATH=/usr/local/bin/miniconda3/bin/:$PATH
		#step2.Run calc-bmr to measure overall and per-gene mutation rates. Give it extra memory, because it may need it
		if [ ! -f gene_mutation_rates_do.SUCCESS ] ; then
			music2 bmr calc-bmr \
			--roi-file ${all_gene_tsv} \
			--reference-sequence ${ref_fa} \
			--bam-list ${bamlist} \
			--maf-file ${maf_file} \
			--output-dir ${out_dir} \
			--show-skipped \
			1>gene_mutation_rates_do.stdout 2>gene_mutation_rates_do.stderr \
			&& touch gene_mutation_rates_do.SUCCESS
		fi
		echo -e "${taskName} gene_mutation_rates_do_work_end@"$(date +%Y-%m-%d_%H:%M:%S)"@"$(hostname) >> log
	>>>

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
		String gene_mrs ="${out_dir}/${taskName}/gene_mrs"
		String log_file="${out_dir}/${taskName}/log"
	}
}



task test_by_fdr_do {
	String taskName
	String gene_mrs
	String fdr_value
	String out_dir
	String mnt_db_dir
	String mnt_input_dir
	String mnt_out_dir
	String task_queue
	String docker_user

	command <<<
		set -e
		if [ ! -d ${out_dir}/${taskName}/ ] ; then
				mkdir -p  ${out_dir}/${taskName}/
		fi
		cd ${out_dir}/${taskName}/
		echo -e "${taskName} test_by_fdr_do_work_start@"$(date +%Y-%m-%d_%H:%M:%S)"@"$(hostname) >> log

		export PATH=/usr/local/bin/miniconda3/bin/:$PATH
		#step3.Run SMG test using an FDR threshold appropriate for these mutation rates
		if [ ! -f test_by_fdr_do.SUCCESS ] ; then
			music2 smg \
			--gene-mr-file ${gene_mrs} \
			--output-file smgs \
			--max-fdr ${fdr_value} \
			--processors 1 \
			1>test_by_fdr_do.stdout 2>test_by_fdr_do.stderr
			&& touch test_by_fdr_do.SUCCESS
		fi
		awk 'NR>1{print $1"\t"$12} BEGIN{print "gene\tq"}' smgs > smgs.simple
		
		echo -e "${taskName} test_by_fdr_do_work_end@"$(date +%Y-%m-%d_%H:%M:%S)"@"$(hostname) >> log
	>>>

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
		String smgs_simple = "${out_dir}/smgs.simple"
		String log_file="${out_dir}/log"
	}
}


task plot_SigMutGenes_do {
	String taskName
	String oncoplot_script
	String maf_file
	String smgs_simple
	String smgsFormat_script
	String Interact_script
	
	String out_dir
	String mnt_db_dir
	String mnt_input_dir
	String mnt_out_dir
	String task_queue
	String docker_user

	command <<<
		set -e
		if [ ! -d ${out_dir}/${taskName}/ ] ; then
				mkdir -p  ${out_dir}/${taskName}/
		fi
		cd ${out_dir}/${taskName}/
		echo -e "${taskName} plot_SigMutGenes_do_work_start@"$(date +%Y-%m-%d_%H:%M:%S)"@"$(hostname) >> log

		export PATH=/usr/local/bin/miniconda3/bin/:$PATH
		#step4.Stat and draw oncoplot
		/usr/local/bin/miniconda3/bin/Rscript \
		${oncoplot_script} \
		--maf ${maf_file} \
		--smg ${smgs_simple} \
		--outdir ./ \
		--outprefix SigMutGenes_Oncoplot 
		
		perl ${smgsFormat_script} \
		${maf_file} \
		${smgs_simple} \
		SigMutGenes.xls
		
		/mnt/X500/farmers/liuqf/libs/Rscript \
		${Interact_script} \
		--maf ${maf_file} \
		--smg ${smgs_simple} \
		--out_dir ./ \
		--outprefix SigMutGenes_Interactions
		
		echo -e "${taskName} plot_SigMutGenes_do_work_end@"$(date +%Y-%m-%d_%H:%M:%S)"@"$(hostname) >> log
	>>>

	runtime {
		docker:"oncowes/vcf2signature:v3"
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
		String SigMutGenes_result = "${out_dir}/SigMutGenes.xls"
		String log_file="${out_dir}/log"
	}
}
