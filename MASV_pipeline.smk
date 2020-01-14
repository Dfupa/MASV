shell.prefix("source ~/.bashrc; ")

#Author: Diego Fuentes
#Contact email: diegofupa@gmail.com
#Date:2020-01-05

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


include: "lib/rules/QC.smk"
include: "lib/rules/alignment.smk"
include: "lib/rules/bedtools_eval.smk"
include: "lib/rules/calling.smk"
#include: "lib/rules/truvari_eval.smk"


configfile: os.path.join(workflow.basedir, "lib/config/config.json") #As a default setting. The config will be named per user interest using MASV_get_config.py
    
    
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

##############
# MAIN RULES #
##############

rule all:
    input:
        expand(workingdir + str(date) + "/{outdir}/{sample}_{aligner}.bam",
              sample=config["General Parameters"]["sample_barcode"],
              outdir=config["Outputs"]["alignment_out"],
              aligner=config["Inputs"]["aligner_selection"]),
        expand(workingdir + str(date) + "/{outdir}/mosdepth/{sample}.mosdepth.global.dist.txt",
              sample=config["General Parameters"]["sample_barcode"],
              outdir=config["Outputs"]["alignment_out"])
        expand(workingdir + str(date) + "/{outdir}/nanoplot-qc",
              outdir=config["Outputs"]["alignment_out"])
        expand(str(date) + "_eval_stats_{svcaller}_feature.txt",
              svcaller=config["Inputs"]["svcaller_selection"])
    
    log:
        logs_dir + str(date) + ".rule_all.log" 

