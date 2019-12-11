
######MOSDEPTH STATS AND PLOT#######
#Main code is found in https://github.com/wdecoster/nano-snakemake/blob/master/rules/mosdepth.smk

rule mosdepth_get:
    input:
        BAM = rules.mapping.output,
        BAI = rules.index_bam.output.BAI
    threads: 4
    output:
        protected("{rules.mapping.input.aligner}/mosdepth/{sample}.mosdepth.global.dist.txt"),
        protected("{rules.mapping.input.aligner}/mosdepth/{sample}.regions.bed.gz"),
    params:
        windowsize = 500,
        prefix = "{sample}",
        aligner = "{rules.mapping.input.aligner}"
    log:
        logs_dir + str(date) +".{ontfile}.mosdepth_{sample}.log"
    benchmark:
        benchmark_dir + str(date) + ".{sample}.mosdepth.benchmark.txt"
        
    conda: "pipeline_env.yml"
    
    shell:
        "mkdir {rules.mapping.input.aligner}/mosdepth && \
                mosdepth --threads {threads} -n \
                  --by {params.windowsize} \
                  {params.aligner}/mosdepth/{params.prefix} {input.BAM} 2> {log}"


rule mosdepth_global_plot:
    input:
        "{rules.mapping.input.aligner}/mosdepth/{sample}.mosdepth.global.dist.txt"
    output:
        "{rules.mapping.input.aligner}/mosdepth_global_plot/global.html"
    log:
        "logs/{rules.mapping.input.aligner}/mosdepth/mosdepth_global_plot.log"
        
    conda: "pipeline_env.yml"
    
    shell:
        "python3 " + os.path.join(workflow.basedir, "lib/scr/mosdepth-plot-dist.py") + \
            " {input} -o {output} 2> {log}"
            
            
#####NANOPLOT####

rule nanoplot_qc:
    input:
        BAM = rules.mapping.output
    output:
        DIR = directory("{rules.mapping.input.aligner}/nanoplot-qc")
    params:
        sample = {sample}
        title = str(date) + "_{sample}"
        
    threads: config["Minimap2"]["minimap2_cores"]
    
    conda: "pipeline_env.yml"
    
    shell:
        "NanoPlot -t {threads} --bam {input.BAM} --raw -o {output.DIR} -p {params.sample}_ --N50 --title {params.title}"
