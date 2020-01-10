shell.prefix("source ~/.bashrc; ")

#Author: Diego Fuentes
#Contact email: diegofupa@gmail.com
#Date:2019-11-04

import os
#import glob
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

sample = config["General Parameters"]["sample_barcode"]
pipeline_version = "v" + str(config["General Parameters"]["version"])
workingdir = config["General Parameters"]["basedir"]

#Benchmark directory
benchmark_dir = workingdir + "/Benchmark/"
if not os.path.exists(benchmark_dir)
    os.makedirs(benchmark_dir)

#logs directory
logs_dir = config["Parameters"]["logs_dir"]
if not os.path.exists(logs_dir):
    os.makedirs(logs_dir)

# file wildcard 
ontfiles = config["Wildcards"]["ONT_fastqs"]


###################
# ALIGNMENT RULES #
###################

rule mapping:
    input:
        reads = expand(config["Inputs"]["ONT_reads_directory"] + "{ontfile}.*", ontfile=ontfiles.split(',')),
        aligner = config["Inputs"]["aligner_selection"],
        ref = config["Inputs"]["reference_genome"]
        
    output:
        protected(workingdir + str(date) + "/{params.outdir}/{sample}_{input.aligner}.bam")
        
    params:
        mapping_tech_minimap, mapping_tech_ngmlr = which_tech(parameter=config["General Parameters"]["seq_technology"]),
        outdir = directory (config["Outputs"]["alignment_out"]),
        match_score_mm2 = config["Minimap2"]["minimap2_match_score"],
        mismatch_score_mm2 = config["Minimap2"]["minimap2_mismatch_score"],
        gap_open_score_mm2 = config["Minimap2"]["minimap2_gap_open_score"],
        GT_AG_find_mm2 = config["Minimap2"]["GT_AG_find"],
        min_ident_ngmlr = config["Ngmlr"]["ngmlr_min_ident"],
        match_score_ngmlr = config["Ngmlr"]["ngmlr_match_score"],
        mismatch_score_ngmlr = config["Ngmlr"]["ngmlr_mismatch_score"],
        gap_open_score_ngmlr = config["Ngmlr"]["ngmlr_gap_open_score"],
        gap_extend_max_ngmlr = config["Ngmlr"]["ngmlr_gap_extend_max"],
        gap_extend_min_ngmlr = config["Ngmlr"]["ngmlr_gap_extend_min"],
        kmer_length_ngmlr = config["Ngmlr"]["ngmlr_kmer_length"],
        kmer_skip_ngmlr = config["Ngmlr"]["ngmlr_kmer_length"]
        
    logs:
        logs_dir + str(date) + ".{ontfile}.alignment.log"
        
    benchmark:
        benchmark_dir + str(date) + ".{ontfile}.alignment.benchmark.txt"

    threads: 
        minimap2_threads = config["Minimap2"]["minimap2_cores"],
        ngmlr_threads = config["Ngmlr"]["ngmlr_cores"]

    conda: "MASV_pipeline.yml"
  
    run:
        if aligner == "minimap2":   #If the selected aligner is minimap2
            shell("mkdir -p {params.outdir}; cd {params.outdir}; minimap2 -t {threads.minimap2_threads} \
            --MD -A {params.match_score_mm2} -B {params.mismatch_score_mm2} -O {params.gap_open_score_mm2} \
            -u {params.GT_AG_find_mm2} -ax {params.mapping_tech_minimap} {input.ref} {input.reads} | samtools sort \
            -@ {threads.minimap2_threads} -O BAM -o {output} 2> {logs}")

        elif aligner == "ngmlr":    #If the selected aligner is ngmlr
            shell("mkdir -p {params.outdir}; cd {params.outdir}; ngmlr -t {threads.ngmlr_threads} \
            -x {params.mapping_tech_ngmlr} -i {params.min_ident_ngmlr} --match {params.match_score_ngmlr} \
            --mismatch {params.mismatch_score_ngmlr} --gap-open {params.gap_open_score_ngmlr} \
            --gap-extend-max {params.gap_extend_max_ngmlr} --gap-extend-min {params.gap_extend_min_ngmlr} \
            --kmer-skip {params.kmer_skip_ngmlr} -k {params.kmer_length_ngmlr} -r {input.ref} \
            -q {input.reads} | samtools sort -@ {threads.minimap2_threads} -O BAM -o {output} 2> {logs}")
        else:
            shell("echo 'An error ocurred in the mapping step. Please, resubmit a valid aligner: minimap2, ngmlr.' > mapping.err; exit")


rule index_bam:
    input:
        BAM = rules.mapping.output
        
    output:
        BAI = workingdir + str(date) + "/{rules.mapping.params.outdir}/{sample}_{rules.mapping.input.aligner}_sorted.bam.bai"
        
    logs:
        logs_dir + str(date) + ".{ontfile}.index_bam.log"
        
    threads:
        config["Minimap2"]["minimap2_cores"]
    
    conda: "MASV_pipeline.yml"
        
    shell:
        "samtools index -@ {threads} {input.BAM} {output.BAI} 2> {logs}"
        

rule alignment_stats:
    input:
        BAM = rules.mapping.sort_bam.output.outBAM
    output:
         "{rules.mapping.input.aligner}/alignment_stats/alignment_stats.txt"
    logs:
        logs_dir + str(date) +".{ontfile}.alignment_stats.log"
        
    conda: "MASV_pipeline.yml"
        
    shell:
        "python3 " + os.path.join(workflow.basedir, "lib/scr/alignment_stats.py") + \
            " -o {output} {input.BAM} 2> {logs}
