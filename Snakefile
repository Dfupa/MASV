from datetime import datetime

date1 = str(datetime.now())
tmp = str.replace(date1," ",".") 
tmp2 = str.replace(tmp,":","")
date = str.replace(tmp2,"-","")

scripts_dir = config["Inputs"]["scripts_dir"]
logs_dir = config["Parameters"]["logs_dir"]
if not os.path.exists(logs_dir):
	os.makedirs(logs_dir)

# file wildcard 
ontfiles = config["Wildcards"]["ONT_fastqs"]

sample = config["General Parameters"]["sample_barcode"]
pipeline_version = "v" + str(config["General Parameters"]["version"])
workingdir = config["General Parameters"]["basedir"}

###############
#### RULES ####
###############

rule index_minimap2:
	input:
		ref = config["Inputs"]["reference_genome"] + {file}.*.fasta
	output:
		workingdir + "/index/minimap2.idx"
	threads: config["Minimap2 parameters"]["minimap2_cores"]
	conda: "MASV_pipeline.yml"

	shell:
			""minimap2 -t {threads} -ax map-ont -Y {input.ref} -d {output}"



rule mapping:
	input:
		reads = config["Inputs"]["ONT_reads_directory"] + "{file}.*.fastq.gz",
		aligner = config["Inputs"]["aligner_selection"]
		ref = config["Inputs"]["reference_genome"] + {file}.*.fasta
		indx = rules.index_minimap2.output
	output:
		trim1 = outdir + {file}.bam
		trim2 = config["Outputs"]["ILLUMINA_trim"] + "{file}.2_val_2.fq.gz",
	params:
		outdir = directory (config["Outputs"]["alignment_out"]),
		kmer_length_mm2 = config["Minimap2 parameters"]["minimap2_kmer_length"]
		match_score_mm2 = config["Minimap2 parameters"]["minimap2_match_score"],
		mismatch_score_mm2 = config["Minimap2 parameters"]["minimap2_mismatch_score"],
		gap_open_score_mm2 = config["Minimap2 parameters"]["minimap2_gap_open_score"],
		GT_AG_find_mm2 = config["Minimap2 parameters"]["GT_AG_find"],
		min_ident_ngmlr = config["Ngmlr parameters"]["ngmlr_min_ident"],
		match_score_ngmlr = config["Ngmlr parameters"]["ngmlr_match_score"],
		mismatch_score_ngmlr = config["Ngmlr parameters"]["ngmlr_mismatch_score"],
		gap_open_score_ngmlr = config["Ngmlr parameters"]["ngmlr_gap_open_score"],
		gap_extend_max_ngmlr = config["Ngmlr parameters"]["ngmlr_gap_extend_max"],
		gap_extend_min_ngmlr = config["Ngmlr parameters"]["ngmlr_gap_extend_min"],
		kmer_length_ngmlr = config["Ngmlr parameters"]["ngmlr_kmer_length"],
		kmer_skip_ngmlr = config["Ngmlr parameters"]["ngmlr_kmer_length"]

	log: 
		logs_dir + str(date) + ".{file}.alignment.log"

	threads: 
		minimap2_threads = config["Minimap2 parameters"]["minimap2_cores"]
		ngmlr_threads = config["Ngmlr parameters"]["ngmlr_cores"]

	#conda: "pipeline_env.yml"
  
	run:
		if aligner == "minimap2":												#If the selected aligner is minimap2
			shell("mkdir -p {params.outdir}; cd {params.outdir}; module purge; module load MINIMAP2/2.9; minimap2 -k {params.kmer_length_mm2} -A {params.match_score_mm2} -B {params.mismatch_score_mm2} -O {params.gap_open_score_mm2} -u {params.GT_AG_find_mm2} -ax map-ont {input.ref} {input.reads} > {output}")
			
  	elif aligner == "ngmlr":												#If the selected aligner is ngmlr
			shell("mkdir -p {params.outdir}; cd {params.outdir}; module purge; srun ngmlr -t {threads.ngmlr_threads} -x ont -i {params.min_ident_ngmlr} --match {params.match_score_ngmlr} --mismatch {params.mismatch_score_ngmlr} --gap-open {params.gap_open_score_ngmlr} --gap-extend-max {params.gap_extend_max_ngmlr} --gap-extend-min {params.gap_extend_min_ngmlr} --kmer-skip {params.kmer_skip_ngmlr} -k {params.kmer_length_ngmlr} -r {input.ref} -q {input.reads} -o {output}")
		else:
			shell("echo 'An error ocurred in the mapping step. Please, resubmit a valid aligner: minimap2, ngmlr.' > mapping.err; exit")


rule sort_bam:

shell.prefix("source ~/.bashrc; ")

#Author: Diego Fuentes
#Contact email: diegofupa@gmail.com
#Date:2019-10-16

import os
#import glob
#import gzip
from datetime import datetime
import sys
sys.path.append('/lib/src/')
from helper_find_file import find_files
from helper_find_file import which_tech

date1 = str(datetime.now())
tmp = str.replace(date1," ",".") 
tmp2 = str.replace(tmp,":","")
date = str.replace(tmp2,"-","")


##############
# PARAMETERS #
##############



#scripts_dir = config["Inputs"]["scripts_dir"]
logs_dir = config["Parameters"]["logs_dir"]
if not os.path.exists(logs_dir):
    os.makedirs(logs_dir)

# file wildcard 
ontfiles = config["Wildcards"]["ONT_fastqs"]

sample = config["General Parameters"]["sample_barcode"]
pipeline_version = "v" + str(config["General Parameters"]["version"])
workingdir = config["General Parameters"]["basedir"]

###############
#### RULES ####
###############

#RULES REGARDING SV CALLING

rule index_minimap2:
    input:
        ref = config["Inputs"]["reference_genome"] + "{file}.*.fasta"
    output: workingdir + str(date) + "/index/{sample}_minimap2.idx"
    params:
        kmer_length_mm2 = config["Minimap2"]["minimap2_kmer_length"]
    threads: config["Minimap2"]["minimap2_cores"]
    conda: "pipeline_env.yml"
    shell:
        "minimap2 -t {threads} -k {params.kmer_length_mm2} -ax map-ont -Y {input.ref} -d {output}"


rule mapping:
    input:
        reads = config["Inputs"]["ONT_reads_directory"] + "{file}.*.fastq.gz"
        aligner = config["Inputs"]["aligner_selection"],
        ref = config["Inputs"]["reference_genome"]
        #indx = rules.index_minimap2.output
    output:
        BAM = workingdir + str(date) + "/{params.outdir}/{sample}_{input.aligner}.bam"
    params:
        mapping_tech_minimap, mapping_tech_ngmlr = which_tech(parameter=config["General Parameters"]["seq_technology"])
        outdir = directory (config["Outputs"]["alignment_out"])
        match_score_mm2 = config["Minimap2"]["minimap2_match_score"]
        mismatch_score_mm2 = config["Minimap2"]["minimap2_mismatch_score"]
        gap_open_score_mm2 = config["Minimap2"]["minimap2_gap_open_score"]
        GT_AG_find_mm2 = config["Minimap2"]["GT_AG_find"]
        min_ident_ngmlr = config["Ngmlr"]["ngmlr_min_ident"]
        match_score_ngmlr = config["Ngmlr"]["ngmlr_match_score"]
        mismatch_score_ngmlr = config["Ngmlr"]["ngmlr_mismatch_score"]
        gap_open_score_ngmlr = config["Ngmlr"]["ngmlr_gap_open_score"]
        gap_extend_max_ngmlr = config["Ngmlr"]["ngmlr_gap_extend_max"]
        gap_extend_min_ngmlr = config["Ngmlr"]["ngmlr_gap_extend_min"]
        kmer_length_ngmlr = config["Ngmlr"]["ngmlr_kmer_length"]
        kmer_skip_ngmlr = config["Ngmlr"]["ngmlr_kmer_length"]

    log: logs_dir + str(date) + ".{file}.alignment.log"

    threads: 
        minimap2_threads = config["Minimap2"]["minimap2_cores"]
        ngmlr_threads = config["Ngmlr"]["ngmlr_cores"]

    conda: "pipeline_env.yml"
  
    run:
        if aligner == "minimap2":#If the selected aligner is minimap2
            shell("mkdir -p {params.outdir}; cd {params.outdir}; minimap2 -A {params.match_score_mm2} -B {params.mismatch_score_mm2} -O {params.gap_open_score_mm2} -u {params.GT_AG_find_mm2} -ax {params.mapping_tech_minimap} {input.ref} {input.reads} > {output}")

        elif aligner == "ngmlr":#If the selected aligner is ngmlr
            shell("mkdir -p {params.outdir}; cd {params.outdir}; srun ngmlr -t {threads.ngmlr_threads} -x {params.mapping_tech_ngmlr} -i {params.min_ident_ngmlr} --match {params.match_score_ngmlr} --mismatch {params.mismatch_score_ngmlr} --gap-open {params.gap_open_score_ngmlr} --gap-extend-max {params.gap_extend_max_ngmlr} --gap-extend-min {params.gap_extend_min_ngmlr} --kmer-skip {params.kmer_skip_ngmlr} -k {params.kmer_length_ngmlr} -r {input.ref} -q {input.reads} -o {output}")
        else:
            shell("echo 'An error ocurred in the mapping step. Please, resubmit a valid aligner: minimap2, ngmlr.' > mapping.err; exit")


rule sort_bam:
    input:
        BAM = rules.mapping.output.BAM
    output:
        outBAM = workingdir + str(date) + "/{rules.mapping.params.outdir}/{sample}_{rules.mapping.input.aligner}_sorted.bam",
        BAI = workingdir + str(date) + "/{params.outdir}/{sample}_{input.aligner}.bam.bai"
        
    log: logs_dir + str(date) + ".{file}.sort_bam.log"
        
    threads: config["Minimap2"]["minimap2_cores"]
    
    conda: "pipeline_env.yml"
        
    shell:
        "samtools sort -@ {threads} -O BAM -o {output.outBAM} {input.BAM} - && samtools index -@ {threads} {input.BAM} {output.BAI}"
        
"""rule bam_to_bedpe:
    input:
        BAM = rules.sort_bam.output.outBAM

    output: workingdir + str(date) + "/{rules.mapping.params.outdir}/{sample}_target.bedpe"
    
    log: logs_dir + str(date) + ".{file}.bam_to_bedpe.log"

    conda: "env.yml"
    shell:
        "bedtools bamtobed -bedpe -i {input.BAM}  > {output}"
"""       
rule alignment_stats:
    input:
        BAM = rules.mapping.sort_bam.output.outBAM
    output:
         "{rules.mapping.input.aligner}/alignment_stats/alignment_stats.txt"
    log:
        logs_dir + str(date) +"/alignment_stats/alignment_stats.log"
    shell:
        "python3 " + os.path.join(workflow.basedir, "lib/scr/alignment_stats.py") + \
            " -o {output} {input.BAM} 2> {log}

rule sv_calling:
    input:
        BAM = rules.sort_bam.output.outBAM
        svcaller = config["Inputs"]["svcaller_selection"]
        ref = config["Inputs"]["reference_genome"]
    output:
        BEDPE = workingdir + str(date) + "/{params.outdir}/{sample}_{input.svcaller}.bedpe"
        #VCF = workingdir + str(date) + "/{params.outdir}/{sample}_{input.svcaller}.VCF"
        outDIR = workingdir + str(date) + "/{params.outdir}/)
    params:
        outdir = directory (config["Outputs"]["sv_call_out"])
        min_sv_length_sniffles = config["Sniffles"]["sniffles_min_sv_length"]
        max_sv_length_sniffles = config["Sniffles"]["sniffles_max_sv_length"]
        min_read_length_sniffles = config["Sniffles"]["sniffles_min_read_length"]
        min_read_map_quality_sniffles =  config["Sniffles"]["sniffles_min_read_mapping_quality"]
        min_read_support_sniffles = config["Sniffles"]["sniffles_min_read_support"]
        num_reads_report_sniffles = config["Sniffles"]["sniffles_num_reads_report"]
        max_num_splits_sniffles = config["Sniffles"]["sniffles_num_reads_report"]
        genotype_sniffles = config["Sniffles"]["sniffles_genotype"]
        cluster_sniffles = config["Sniffles"]["sniffles_cluster"]
        min_homo_af_sniffles = config["Sniffles"]["sniffles_min_homo_af"]
        min_het_af_sniffles = config["Sniffles"]["sniffles_min_het_af"]
        min_sv_length_svim = config["Svim"]["svim_min_sv_length"]
        max_sv_lenght_svim = config["Svim"]["svim_max_sv_length"]
        min_read_length_svim = config["Svim"]["svim_min_read_length"]
        min_read_mapp_quality = config["Svim"]["svim_min_read_mapping_quality"]
        gap_tolerance_svim = config["Svim"]["svim_gap_tolerance"]
        overlap_tolerance_svim = config["Svim"]["svim_overlap_tolerance"]
        partition_max_distance_svim = config["Svim"]["svim_partition_max_distance"]
        sv_max_distance_svim = config["Svim"]["svim_sv_max_distance"]
        min_geno_score_svim = config["Svim"]["svim_min_geno_score"]
        homozygous_thresh_svim = config["Svim"]["svim_homozygous_thresh"]
        heterozygous_thresh_svim = config["Svim"]["svim_heterozygous_thresh"]
        min_depth_svim = config["Svim"]["svim_min_depth"]
        duplicat_as_insert = config["Svim"]["svim_duplicat_as_insert"]
        
    threads:
        sniffles_threads = config["Sniffles"]["sniffles_cores"]
        
    logs: logs_dir + str(date) + ".{file}.svcalling.log"
        
    run:
        if svcaller == "sniffles":#If the selected caller is sniffles
            shell("sniffles -m {input.BAM} -v {output} -s {params.min_read_support_sniffles} -r {params.min_read_length_sniffles} -l {params.min_sv_length_sniffles} -d {params.max_sv_length_sniffles} -q {params.min_read_map_quality_sniffles} -n {params.num_reads_report_sniffles} --genotype {params.genotype_sniffles} --min_homo_af {params.min_homo_af_sniffles} --min_het_af {params.min_het_af_sniffles} --cluster {params.cluster_sniffles} --report_read_strands -b {output.BEDPE}")
            
        if svcaller == "svim": #If the selected caller is svim
            shell("svim alignment {output.outDIR} {input.BAM} {input.ref} ")
                  
