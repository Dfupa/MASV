#Main code is found in https://github.com/wdecoster/nano-snakemake/blob/master/rules/mosdepth.smk

######MOSDEPTH STATS AND PLOT#######

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
    shell:
        "python3 " + os.path.join(workflow.basedir, "lib/scr/mosdepth-plot-dist.py") + \
            " {input} -o {output} 2> {log}"
