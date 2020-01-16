
######MOSDEPTH STATS AND PLOT#######
#Main code is found in https://github.com/wdecoster/nano-snakemake/blob/master/rules/mosdepth.smk

rule mosdepth_get:
    input:
        BAM = lambda wildcards: expand(rules.mapping.output, ontfile=ontfiles.split(',')),
        BAI = lambda wildcards: expand(rules.index_bam.output.BAI, ontfile=ontfiles.split(','))
    threads: 4
    output:
        txt=protected(outdir+"mosdepth/"+ sample + "_" + aligner+ ".{ontfile}.mosdepth.global.dist.txt"),
        bed=protected(outdir+"mosdepth/"+ sample + "_" + aligner+ ".{ontfile}.regions.bed.gz")
    params:
        windowsize = 500,
        prefix = sample,
        dir = outdir
    log:
        logs_dir + str(date) + "_" + sample +".{ontfile}.mosdepth.log"
    benchmark:
        benchmark_dir + str(date) + "_" + sample +".{ontfile}.mosdepth.benchmark.txt"
        
    conda: "MASV_pipeline.yml"
    
    shell:
        "mkdir "+aligner+"/mosdepth && mosdepth --threads {threads} -n \
        --by {params.windowsize} {params.dir}mosdepth/{params.prefix} {input.BAM} 2> {log}"


rule mosdepth_global_plot:
    input:
        plot = lambda wildcards: expand(rules.mosdepth_get.output.txt, ontfile=ontfiles.split(','))
    output:
        protected(outdir+"mosdepth/{ontfile}_global_plot.html")
    log:
        logs_dir + str(date) + "_" + sample +".{ontfile}.mosdepth_global_plot.log"
        
    conda: "MASV_pipeline.yml"
    
    shell:
        "python3 " + os.path.join(workflow.basedir, "lib/scr/mosdepth-plot-dist.py") + \
        " {input.plot} -o {output} 2> {log}"
            
            
#####NANOPLOT####

rule nanoplot_qc:
    input:
        BAM = lambda wildcards: expand(rules.mapping.output, ontfile=ontfiles.split(','))
    output:
        DIR = directory(outdir+str(date)+"_"+sample+".{ontfile}_nanoplot-qc/")
    params:
        id = sample,
        title = str(date) + "_"+sample+".{ontfile}"
        
    threads: config["Minimap2"]["minimap2_cores "]
    
    conda: "MASV_pipeline.yml"
    
    shell:
        "NanoPlot -t {threads} --bam {input.BAM} --raw -o {output.DIR} -p {params.id}_ --N50 --title {params.title}"
