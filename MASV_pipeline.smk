shell.prefix("source ~/.bashrc; ")

#Author: Diego Fuentes
#Contact email: diegofupa@gmail.com
#Date:2020-01-05

import os
#import glob
from datetime import datetime
import sys
sys.path.append('./lib/scr/')
from helper_functions import which_tech

date1 = str(datetime.now())
tmp = str.replace(date1," ",".") 
tmp2 = str.replace(tmp,":","")
date = str.replace(tmp2,"-","")

#First determine the config file
configfile: os.path.join(workflow.basedir, "lib/config/config.json") #As a default setting. The config will be named per user interest using MASV_get_config.py

##############
# PARAMETERS #
##############

#Global parameters
sample = config["Parameters"]["sample_barcode"]
pipeline_version = "v" + str(config["Parameters"]["version"])
workingdir = config["Parameters"]["basedir"]
outdir = config["Outputs"]["alignment_out"]
svout = config["Outputs"]["svcall_out"]
#aligner = config["Inputs"]["aligner_selection"]
svcaller = config["Inputs"]["svcaller_selection"]
aligner = "ngmlr"

#Benchmark directory
benchmark_dir = workingdir + "Benchmark/"
if not os.path.exists(benchmark_dir):
    os.makedirs(benchmark_dir)

#logs directory
logs_dir = config["Parameters"]["logs_dir"]
if not os.path.exists(logs_dir):
    os.makedirs(logs_dir)

# file wildcard
ontfiles = config["Wildcards"]["ONT_fastqs"]


################
# RULES STEPUP #
################

include: "lib/rules/alignment.smk"
include: "lib/rules/QC.smk"
include: "lib/rules/calling.smk"
#include: "lib/rules/bedtools_eval.smk"
#include: "lib/rules/truvari_eval.smk"

##############
# MAIN RULES #
##############


rule all:
    input:
        bam=expand(outdir + sample +"_"+aligner+".{ontfile}.bam",
        ontfile=ontfiles.split(',')),
        bai=expand(outdir + sample +"_"+ aligner+".{ontfile}.bam.bai",
        ontfile=ontfiles.split(',')),
        stats=expand(outdir + sample +"_"+ aligner + "/alignment_stats/{ontfile}.alignment_stats.txt",
        ontfile=ontfiles.split(',')),
        dist=expand(outdir + "mosdepth/"+ sample + "_"+ aligner+ ".{ontfile}.mosdepth.global.dist.txt",
        ontfile=ontfiles.split(',')),
        plot=expand(outdir + "mosdepth/" + sample  +"_"+aligner+ ".{ontfile}_global_plot.html",
        ontfile=ontfiles.split(',')),
        nanoplot=expand(outdir + str(date) + "_" + sample + ".{ontfile}_nanoplot-qc/",
        ontfile=ontfiles.split(',')),
        vcf = expand(svout+"/"+str(date)+"_"+sample+"_"+svcaller+".{ontfile}.vcf", ontfile=ontfiles.split(',')),
        svim = expand(svout+"/"+str(date)+"_"+sample+"_"+svcaller+".{ontfile}/{ontfile}_minscore_{minscore}.vcf",
        ontfile=ontfiles.split(','),
        minscore=config["Svim"]["svim_min_score"])

        #VCF = protected(outdir+sample+"_{input.svcaller}.{ontfile}.vcf"),
        #outDIR = directory(outdir+sample+"_{input.svcaller}.{ontfile}/")

        #expand(workingdir + str(date)  + "_" + sample + ".{ontfile}_eval_stats_{svcaller}_feature.txt",
        #ontfile=ontfiles.split(','),
        #svcaller=config["Inputs"]["svcaller_selection"])
    log:
        logs_dir + str(date) + "_" + sample +"_"+aligner+"_.all.rule.log"

if aligner == "minimap2":
    rule mapping_only:
        input:
            bam=expand(outdir + sample +"_"+aligner+".{ontfile}.bam",
            ontfile=ontfiles.split(','))
        #bai=expand(outdir + sample +"_"+ aligner+".{ontfile}.bam.bai",
        #ontfile=ontfiles.split(',')),
        #stats=expand(outdir + sample +"_"+aligner+"/alignment_stats/{ontfile}.alignment_stats.txt",
        #ontfile=ontfiles.split(',')),
        #dist=expand(outdir + "mosdepth/"+ sample + "_"+ aligner+ ".{ontfile}.mosdepth.global.dist.txt",
        #ontfile=ontfiles.split(',')),
        #plot=expand(outdir + "mosdepth/" + sample  +"_"+aligner+ ".{ontfile}_global_plot.html",
        #ontfile=ontfiles.split(',')),
        #nanoplot=expand(outdir + str(date) + "_" + sample + ".{ontfile}_nanoplot-qc/",
        #ontfile=ontfiles.split(','))
        log:
            logs_dir + str(date) + "_" + sample +"_"+aligner+"_just.mapping.rule.log"
