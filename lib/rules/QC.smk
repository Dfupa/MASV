
######MOSDEPTH STATS AND PLOT#######
#Main code is found in https://github.com/wdecoster/nano-snakemake/blob/master/rules/mosdepth.smk

rule mosdepth_get:
    input:
        BAM = outdir + sample +"_"+ aligner + ".{ontfile}.bam",
        BAI = rules.index_bam.output.BAI
    threads: 4
    output:
        txt=protected(outdir+"mosdepth/"+ sample + "_" + aligner+ ".{ontfile}.mosdepth.global.dist.txt"),
        bed=protected(outdir+"mosdepth/"+ sample + "_" + aligner+ ".{ontfile}.regions.bed.gz")
    params:
        windowsize = 500,
        prefix = outdir+"mosdepth/"+ sample + "_" + aligner+ ".{ontfile}"
    log:
        logs_dir + str(date) + "_" + sample +".{ontfile}.mosdepth.log"
    benchmark:
        benchmark_dir + str(date) + "_" + sample +".{ontfile}.mosdepth.benchmark.txt"
        
    conda: os.path.join(workflow.basedir, "MASV_pipeline.yml")
    
    shell:
        "mosdepth --threads {threads} -n \
         --by {params.windowsize} \
         {params.prefix} {input.BAM} 2> {log}"


rule mosdepth_global_plot:
    input:
        outdir+"mosdepth/{ontfile}.mosdepth.global.dist.txt"
    output:
        protected(outdir+"mosdepth/{ontfile}_global_plot.html")
    log:
        logs_dir + str(date) + "_" + sample +".{ontfile}.mosdepth_global_plot.log"
        
    conda: os.path.join(workflow.basedir, "MASV_pipeline.yml")
    
    shell:
        "python3 " + os.path.join(workflow.basedir, "lib/scr/mosdepth-plot-dist.py") + \
            " {input} -o {output} 2> {log}"
            
            
#####NANOPLOT####

rule nanoplot_qc:
    input:
        BAM = outdir + sample +"_"+aligner+".{ontfile}.bam",
        BAI = rules.index_bam.output.BAI
    output:
        DIR = directory(outdir+str(date)+"_"+sample+".{ontfile}_nanoplot-qc/")
    params:
        id = sample,
        title = str(date) + "_"+sample+".{ontfile}"
        
    threads: config["Minimap2"]["minimap2_cores "]
    
    conda: os.path.join(workflow.basedir, "MASV_pipeline.yml")
    
    shell:
        "NanoPlot -t {threads} --bam {input.BAM} --raw -o {output.DIR} -p {params.id}_ --N50 --title {params.title}"
