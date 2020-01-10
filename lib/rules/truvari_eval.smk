########ALPHA VERSION##########
#Please, note that this code implementation is not yet 100% functional. It was scrapped for the main bulk of the project
#but it can be adapted to work under the right circustamces. Right now it is not included in the main snakemake pipeline.

##########The status of this rule.smk is WIP on hiatus ###############

####################
#TRUVARI EVALUATION#
####################

import gzip

rule truvari_conversion_svim:
    input:
        rules.filter_svim.output
        
    output:
        temp(rules.sv_calling.output.outDIR + "minscore_{rules.filter_svim.params.minscore,[0-9]+}.truvari.vcf")
    
    threads: 1
    
    conda: "MASV_pipeline.yml"
        
    logs:
        logs_dir + str(date) + "truvari_svim_conversion.log"
    
    run:
        if rules.sv_calling.input.svcaller == "svim":
            shell("cat {input} | sed 's/INS:NOVEL/INS/g' | sed 's/DUP:INT/INS/g' | sed 's/DUP:TANDEM/INS/g' | awk 'OFS=\"\\t\" {{ if($1 ~ /^#/) {{ print $0 }} else {{ if($5==\"<DEL>\" || $5==\"<INS>\") {{ print $1, $2, $3, $4, $5, $6, \"PASS\", $8, $9, $10 }} }} }}' > {output} 2> {logs}"
        else:
            shell("echo 'As svim wasn't selected this step is aborted. Continuing.' 2> {logs}")
            
rule sort_calls_sniffles:
    input:
        rules.sv_calling.output.VCF
    output:
        temp(rules.sv_calling.output.outDIR + "{sample}_{input.svcaller}.sorted.vcf")
        
    threads: 1
    
    conda: "MASV_pipeline.yml"
    
    logs:
        logs_dir + str(date) + ".{ontfile}.vcf_sorting.log"
    run:
        if rules.sv_calling.input.svcaller == "sniffles":
            shell("cat {input} | awk '$1 ~ /^#/ {print $0;next} {print $0 | \"sort -k1,1 -k2,2n\"}' > {input} && bcftools sort {input} > {output} 2> {logs}")
        else:
            shell("echo 'This step is for sniffles calls. Continuing.' 2> {logs}")

rule sort_calls_svim:
    input:
        rules.truvari_conversion_svim.output
    output:
        temp(rules.sv_calling.output.outDIR + "minscore_{rules.filter_svim.params.minscore,[0-9]+}.truvari.sorted.vcf")
    threads: 1
    
    conda: "MASV_pipeline.yml"
    
    logs:
        logs_dir + str(date) + "sorting_svim_minscore_{rules.filter_svim.params.minscore,[0-9]+}.truvari.log"
    
    shell:
        "bcftools sort {input} > {output} 2> {logs}"

rule bgzip_tabix_sniffles:
    input:
        rules.sort_calls_sniffles.output
    output:
        protected(rules.sv_calling.output.outDIR + "{sample}_{input.svcaller}.sorted.vcf.gz")
        
    conda: "MASV_pipeline.yml"
                  
    logs:
        logs_dir + str(date) + ".{ontfile}.bgzip_tabix_sniffles.log"
    
    shell:
        "bgzip {input} > {output} && tabix {output} 2> {logs}"


rule bgzip_tabix_svim:
    input:
        rules.sort_calls_svim.output
    output:
        protected(rules.sv_calling.output.outDIR + "minscore_{rules.filter_svim.params.minscore,[0-9]+}.truvari.sorted.vcf.gz")
    
    conda: "MASV_pipeline.yml"
    
    logs:
        logs_dir + str(date) + ".{ontfile}.bgzip_tabix_sniffles.log"
    shell:
        "bgzip {input} > {output} && tabix {output} 2> {logs}"
                  
rule truvari_eval:
    """Note that this implementation of truvari is limited without including a bed to delimit by regions or by
    using a p = 0.00 in order to exclude allele frequencies
    
    I could add this new parameters into the config file. Working in progress till implementation"""
                  
    input:       
        hq_vcf = config["Inputs"]["hq_vcf"],
        genome = rules.mapping.input.ref
        
    output:
        directory(workingdir + str(date) + params.outdir)
                  
    params:
        outdir = workingdir + str(date) + "/{rules.params.outdir}/truvari_eval
    
    conda: "MASV_pipeline.yml"
    
    logs: logs_dir + str(date) + ".{ontfile}.truvari_eval.log
                  
    benchmark:
        benchmark_dir + str(date) + ".{ontfile}.truvari_eval.benchmark.txt"
                  
    run:
        if rules.sv_calling.input.svcaller == "svim":
            shell("rm -rf {params.out_dir} && truvari -f {input.genome}\
                    -b {input.hq_vcf} -c rules.bgzip_tabix_svim.output -o {params.out_dir}\
                    --passonly -r 500 -p 0.00 2> {logs}")
                  
        if rules.sv_calling.input.svcaller == "sniffles":
            shell("rm -rf {params.out_dir} && truvari -f {input.genome}\
                    -b {input.hq_vcf} -c rules.bgzip_tabix_sniffles.output -o {params.out_dir}\
                    --passonly -r 500 -p 0.00 2> {logs}")
