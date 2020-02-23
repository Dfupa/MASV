################
# VCF REFORMAT #
################

rule vcf_reformat_sniffles:
    input:
        vcf = rules.sniffles_calling.output.VCF,
        ref = config["Inputs"]["reference_genome"] 
    output:
        protected(svout+"/"+str(date)+"_Sniffles/{ontfile}/"+str(date)+"_"+sample+"_reformated_sniffles.{ontfile}.vcf")
    log:
        logs_dir + str(date) + "_" + sample + "sniffles.{ontfile}.reformat.log"
    params:
        caller = 'sniffles'
           
    benchmark: 
        benchmark_dir + str(date) + "_" + sample + "sniffles.{ontfile}.reformat.txt"
        
    threads: 1
        
    conda: os.path.join(workflow.basedir, "MASV_pipeline.yml")
        
    shell:
       "python3 " + os.path.join(workflow.basedir, "lib/scr/insertion_fix.py") + \
       " -v {input.vcf} -g {input.ref} -o {output} --caller {params.caller} 2> {log}"
    
    
rule vcf_reformat_svim:
    input:
        vcf = rules.svim_calling.output.outVCF,
        ref = config["Inputs"]["reference_genome"] 
    output:
        protected(svout+"/"+str(date)+"_"+sample+"_svim.{ontfile}/final_results_reformated.vcf")
    log:
        logs_dir + str(date) + "_" + sample + "sniffles.{ontfile}.reformat.log"
    params:
        caller = 'svim'
           
    benchmark: 
        benchmark_dir + str(date) + "_" + sample + "svim.{ontfile}.reformat.txt"
        
    threads: 1
        
    conda: os.path.join(workflow.basedir, "MASV_pipeline.yml")
        
    shell:
       "python3 " + os.path.join(workflow.basedir, "lib/scr/insertion_fix.py") + \
       " -v {input.vcf} -g {input.ref} -o {output} --caller {params.caller} 2> {log}"
    
#####################
#BEDTOOLS EVALUATION#
#####################


rule eval_stats_sniffles:
    input:
        sniffles = rules.vcf_reformat_sniffles.output,
        truth = config["Inputs"]["hq_vcf"]
    output:
        protected(svout+"/"+str(date)+"_Sniffles/{ontfile}/Eval_stats_sniffles_DEL.txt")#Note that the DEL in the output is the feature in the params with a strip of < and >. It is necessary to find the correct name. It has to do with the original eval_stats.py
    log:
        logs_dir + str(date) + "_" + sample + "sniffles.{ontfile}.eval_stats_bedtools.log"
    params:
        feature = "'<DEL>'", #This can be edited manually to your needs
        iterator = 30,
        minsup = config["Sniffles"]["sniffles_min_read_support"]
        
    benchmark: 
        benchmark_dir + str(date) + "_" + sample + "sniffles.{ontfile}.eval.benchmark.txt"
        
    threads: 1
        
    conda: os.path.join(workflow.basedir, "MASV_pipeline.yml")
        
    shell:
        "python3 " + os.path.join(workflow.basedir, "lib/scr/eval_stats.py")+ " --truth {input.truth} --callset {input.sniffles} --plot --sniffles --svtype {params.feature} -i {params.iterator} -ms {params.minsup} 2> {log}"
                      
                      
rule eval_stats_svim:
    input:
        svim = rules.vcf_reformat_svim.output,
        truth = config["Inputs"]["hq_vcf"]
    output:
        protected(svout+"/"+str(date)+"_"+sample+"_svim.{ontfile}/Eval_stats_svim_DEL.txt") #same here
    log:
        logs_dir + str(date) + "_" + sample + "svim.{ontfile}.eval_stats_bedtools.log"
    params:
        feature = "'<DEL>'", #This can be edited manually to your needs
        iterator = 20
        
    benchmark: 
        benchmark_dir + str(date) + "_" + sample + "sniffles.{ontfile}.eval.benchmark.txt"
        
    threads: 1
        
    conda: os.path.join(workflow.basedir, "MASV_pipeline.yml")
        
    shell:
        "python3 " + os.path.join(workflow.basedir, "lib/scr/eval_stats.py") + " --truth {input.truth} --callset {input.svim} --plot --svtype {params.feature} -i {params.iterator} 2> {log}"
                     


