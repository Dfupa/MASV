shell.prefix("source ~/.bashrc; ")

#Author: Diego Fuentes
#Contact email: diegofupa@gmail.com
#Date:2020-01-05

import os
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
aligner = config["Inputs"]["aligner_selection"]



#Benchmark directory
benchmark_dir = workingdir + "Benchmark/"
if not os.path.exists(benchmark_dir):
    os.makedirs(benchmark_dir)

#logs directory
logs_dir = config["Parameters"]["logs_dir"]
if not os.path.exists(logs_dir):
    os.makedirs(logs_dir)

#mosdepth directory
mosdir= outdir+"mosdepth/"
if not os.path.exists(mosdir):
    os.makedirs(mosdir)


#svout dir
if not os.path.exists(svout+"/"):
    os.makedirs(svout+"/")

#sniffles dir
if not os.path.exists(svout+"/"+str(date)+"_Sniffles"):
    os.makedirs(svout+"/"+str(date)+"_Sniffles")

# file wildcard
ontfiles = config["Wildcards"]["ONT_reads_directory"] 


#svim filtering
minscore = config["Svim"]["svim_min_score"]

################
# RULES STEPUP #
################

include: "lib/rules/alignment.smk"
include: "lib/rules/QC.smk"
include: "lib/rules/calling.smk"
include: "lib/rules/bedtools_eval.smk"
#include: "lib/rules/truvari_eval.smk"

##############
# MAIN RULES #
##############

rule all:
    input:
        bam=expand(outdir + sample +"_"+ aligner + ".{ontfile}.bam",
        ontfile=ontfiles.split(',')),
        bai=expand(outdir + sample +"_"+ aligner + ".{ontfile}.bam.bai",
        ontfile=ontfiles.split(',')),
        stats=expand(outdir + sample +"_"+ aligner + "/alignment_stats/{ontfile}.alignment_stats.txt",
        ontfile=ontfiles.split(',')),
        dist=expand(outdir + "mosdepth/"+ sample + "_"+ aligner+ ".{ontfile}.mosdepth.global.dist.txt",
        ontfile=ontfiles.split(',')),
        plot=expand(outdir + "mosdepth/" + sample  +"_"+aligner+ ".{ontfile}_global_plot.html",
        ontfile=ontfiles.split(',')),
        nanoplot=expand(outdir + str(date) + "_" + sample + ".{ontfile}_nanoplot-qc/",
        ontfile=ontfiles.split(',')),
        snifflesvcf=expand(svout+"/"+str(date)+"_Sniffles/{ontfile}/"+str(date)+"_"+sample+"_sniffles.{ontfile}.vcf", ontfile=ontfiles.split(',')),
        svimvcf=expand(svout+"/"+str(date)+"_"+sample+"_svim.{ontfile}/final_results.vcf",
        ontfile=ontfiles.split(',')),
        svimfilt=expand(svout+"/"+str(date)+"_"+sample+"_svim.{ontfile}/final_results_minscore_{minscore}.vcf",
        ontfile=ontfiles.split(','),
        minscore=config["Svim"]["svim_min_score"]),
        reformatsvim=expand(svout+"/"+str(date)+"_"+sample+"_svim.{ontfile}/final_results_reformated.vcf", ontfile=ontfiles.split(',')),
        reformatsniffles=expand(svout+"/"+str(date)+"_Sniffles/{ontfile}/"+str(date)+"_"+sample+"_reformated_sniffles.{ontfile}.vcf", ontfile=ontfiles.split(',')),
        eval_sniffles=expand(svout+"/"+str(date)+"_Sniffles/{ontfile}/Eval_stats_sniffles_DEL.txt", ontfile=ontfiles.split(',')),
        eval_svim=expand(svout+"/"+str(date)+"_"+sample+"_svim.{ontfile}/Eval_stats_svim_DEL.txt", ontfile=ontfiles.split(','))


    log:
        logs_dir + str(date) + "_" + sample +"_"+aligner+".all.rule.log"

rule mapping_only:
    input:
        bam=expand(outdir + sample +"_"+ aligner + ".{ontfile}.bam",
        ontfile=ontfiles.split(',')),
        bai=expand(outdir + sample +"_"+ aligner + ".{ontfile}.bam.bai",
        ontfile=ontfiles.split(',')),
        stats=expand(outdir + sample +"_"+ aligner + "/alignment_stats/{ontfile}.alignment_stats.txt",
        ontfile=ontfiles.split(',')),
	dist=expand(outdir + "mosdepth/"+ sample + "_"+ aligner+ ".{ontfile}.mosdepth.global.dist.txt",
        ontfile=ontfiles.split(',')),
        plot=expand(outdir + "mosdepth/" + sample  +"_"+aligner+ ".{ontfile}_global_plot.html",
	ontfile=ontfiles.split(',')),
	nanoplot=expand(outdir + str(date) + "_" + sample + ".{ontfile}_nanoplot-qc/",
        ontfile=ontfiles.split(','))
    log:
        logs_dir + str(date) + "_" + sample +"_"+aligner+".mapping.only.rule.log"

rule sniffles:
    input:
        stats=expand(outdir + sample +"_"+ aligner + "/alignment_stats/{ontfile}.alignment_stats.txt",
        ontfile=ontfiles.split(',')),
        snifflesvcf=expand(svout+"/"+str(date)+"_Sniffles/{ontfile}/"+str(date)+"_"+sample+"_sniffles.{ontfile}.vcf", ontfile=ontfiles.split(',')))
    log:
        logs_dir + str(date) + "_" + sample +"_"+aligner+".sniffles.rule.log"

rule svim:
    input:
        stats=expand(outdir + sample +"_"+ aligner + "/alignment_stats/{ontfile}.alignment_stats.txt",
        ontfile=ontfiles.split(',')),
        svimvcf=expand(svout+"/"+str(date)+"_"+sample+"_svim.{ontfile}/final_results.vcf",
        ontfile=ontfiles.split(',')),
        svimfilt=expand(svout+"/"+str(date)+"_"+sample+"_svim.{ontfile}/final_results_minscore_{minscore}.vcf",
        ontfile=ontfiles.split(','),
	minscore=config["Svim"]["svim_min_score"])
    log:
        logs_dir + str(date) + "_" + sample +"_"+aligner+".svim.rule.log"

rule eval_sniffles:
    input:
        stats=expand(outdir + sample +"_"+ aligner + "/alignment_stats/{ontfile}.alignment_stats.txt",
        ontfile=ontfiles.split(',')),
        eval_sniffles=expand(svout+"/"+str(date)+"_Sniffles/{ontfile}/Eval_stats_sniffles_DEL.txt", ontfile=ontfiles.split(','))
    log:
        logs_dir + str(date) + "_" + sample +"_"+aligner+".eval_sniffles.rule.log"

rule eval_svim:
    input:
        stats=expand(outdir + sample +"_"+ aligner + "/alignment_stats/{ontfile}.alignment_stats.txt",
        ontfile=ontfiles.split(',')),
        eval_svim=expand(svout+"/"+str(date)+"_"+sample+"_svim.{ontfile}/Eval_stats_svim_DEL.txt", ontfile=ontfiles.split(','))
    log:
        logs_dir + str(date) + "_" + sample +"_"+aligner+".eval_svim.rule.log"

rule sanity_check:
    input:
        stats=expand(outdir + sample +"_"+ aligner + "/alignment_stats/{ontfile}.alignment_stats.txt",
	ontfile=ontfiles.split(',')),
        dist=expand(outdir + "mosdepth/"+ sample + "_"+ aligner+ ".{ontfile}.mosdepth.global.dist.txt",
        ontfile=ontfiles.split(',')),
        plot=expand(outdir + "mosdepth/" + sample  +"_"+aligner+ ".{ontfile}_global_plot.html",
        ontfile=ontfiles.split(',')),
        nanoplot=expand(outdir + str(date) + "_" + sample + ".{ontfile}_nanoplot-qc/",
        ontfile=ontfiles.split(','))
    log:
        logs_dir + str(date) + "_" + sample +"_"+aligner+".sanity.check.rule.log"
