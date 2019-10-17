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
	conda: "pipeline_env.yml"

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


