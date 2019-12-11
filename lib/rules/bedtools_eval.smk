import os

#####################
#BEDTOOLS EVALUATION#
#####################

rule eval_stats:
    input:
        SVIM = rules.filter_svim.output,
        SNIFFLES = rules.sv_calling.output.VCF,
        truth = config["Inputs"]["hq_vcf"]
    log:
        logs_dir + str(date) +".{ontfile}.eval_stats_bedtools.log"
    params:
        iterator = 10 + int(config["Svim"]["svim_mins_core"]),
        feature = ['<DEL>', '<INS>', '<INV>', '<DUP:TANDEM>'] #This can be edited manually to your needs
        
    benchmark: 
        benchmark_dir + str(date) + ".{ontfile}.eval.benchmark.txt"
        
    threads: 1
        
    conda: "pipeline_env.yml"
        
    run:
        if rules.sv_calling.input.svcaller == "svim":
            for sv in params.feature:
                shell("python3 " + os.path.join(workflow.basedir, "lib/scr/eval_stats.py") + \
                    "--truth {input.truth} --callset {input.SVIM} --plot True --svtype " + \
                     str(sv) + " --iterator {params.iterator} 2> {log}")
        if rules.sv_calling.input.svcaller == "sniffles":
            for sv in params.feature:
                shell("python3 " + os.path.join(workflow.basedir, "lib/scr/eval_stats.py") + \
                    "--truth {input.truth} --callset {input.SNIFFLES} --sniffles True--svtype " + \
                      str(sv) + " 2> {log}")
            

