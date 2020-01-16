import os

#We set the feature list to be used:

feature = ['<DEL>', '<INS>', '<INV>', '<DUP:TANDEM>'] #This can be edited manually to your needs

#####################
#BEDTOOLS EVALUATION#
#####################

rule eval_stats_sniffles:
    input:
        sniffles = lambda wildcards: expand(rules.sv_calling.output.VCF, ontfile=ontfiles.split(',')),
        truth = config["Inputs"]["hq_vcf"]
    output:
        protected(workingdir + str(date) + "_" + sample + "sniffles.{ontfile}_eval_stats_feature.txt")
    log:
        logs_dir + str(date) + "_" + sample + "sniffles.{ontfile}.eval_stats_bedtools.log"
    params:
        #feature = ['<DEL>', '<INS>', '<INV>', '<DUP:TANDEM>'] #This can be edited manually to your needs
        
    benchmark: 
        benchmark_dir + str(date) + "_" + sample + "sniffles.{ontfile}.eval.benchmark.txt"
        
    threads: 1
        
    conda: "MASV_pipeline.yml"
        
    run:
       for sv in feature:
           shell("python3 " + os.path.join(workflow.basedir, "lib/scr/eval_stats.py") + \
           "--truth {input.truth} --callset {input.sniffles} --sniffles True --svtype " + \
           str(sv) + " 2> {log}")
                      
                      
rule eval_stats_sniffles:
    input:
        svim = lambda wildcards: expand(rules.filter_svim.output, ontfile=ontfiles.split(',')),
        truth = config["Inputs"]["hq_vcf"]
    output:
        protected(workingdir + str(date) + "_" + sample + "svim.{ontfile}_eval_stats_feature.txt")
    log:
        logs_dir + str(date) + "_" + sample + "svim.{ontfile}.eval_stats_bedtools.log"
    params:
        iterator = 10 + int(config["Svim"]["svim_min_score"]),
        
    benchmark: 
        benchmark_dir + str(date) + "_" + sample + "svim.{ontfile}.eval.benchmark.txt"
        
    threads: 1
        
    conda: "MASV_pipeline.yml"
        
    run:
       for sv in params.feature:
           shell("python3 " + os.path.join(workflow.basedir, "lib/scr/eval_stats.py") + \
           "--truth {input.truth} --callset {input.SVIM} --plot True --svtype " + \
           str(sv) + " --iterator {params.iterator} 2> {log}")
                     

