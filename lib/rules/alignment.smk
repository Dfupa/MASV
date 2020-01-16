
#Author: Diego Fuentes
#Contact email: diegofupa@gmail.com
#Date:2019-11-04


import os

###################
# ALIGNMENT RULES #
###################



rule mapping:
    input:
        reads = expand(config["Inputs"]["ONT_reads_directory"] + "{ontfile}.fastq.gz", ontfile=ontfiles.split(',')),
        ref = config["Inputs"]["reference_genome"]
    output:
        protected(outdir + sample +"_"+aligner+".{ontfile}.bam")
        
    params:
        mapping_tech_minimap = which_tech(parameter=config["Parameters"]["seq_technology"])[0], 
        mapping_tech_ngmlr = which_tech(parameter=config["Parameters"]["seq_technology"])[1],
        outdir = config["Outputs"]["alignment_out"],
        match_score_mm2 = config["Minimap2"]["minimap2_match_score"],
        mismatch_score_mm2 = config["Minimap2"]["minimap2_mismatch_score"],
        gap_open_score_mm2 = config["Minimap2"]["minimap2_gap_open_score"],
        GT_AG_find_mm2 = config["Minimap2"]["GT_AG_find"],
        min_ident_ngmlr = config["Ngmlr"]["ngmlr_min_ident"],
        match_score_ngmlr = config["Ngmlr"]["ngmlr_match_score"],
        mismatch_score_ngmlr = config["Ngmlr"]["ngmlr_mismatch_score"],
        gap_open_score_ngmlr = config["Ngmlr"]["ngmlr_gap_open_score"],
        gap_extend_max_ngmlr = config["Ngmlr"]["ngmlr_gap_extend_max"],
        gap_extend_min_ngmlr = config["Ngmlr"]["ngmlr_gap_extend_min"],
        kmer_length_ngmlr = config["Ngmlr"]["ngmlr_kmer_length"],
        kmer_skip_ngmlr = config["Ngmlr"]["ngmlr_kmer_length"],
        ngmlr_threads = config["Ngmlr"]["ngmlr_cores"]
    log:
        logs_dir + str(date) + "_" + sample +"_"+aligner+".{ontfile}.alignment.log"
        
    benchmark:
        benchmark_dir + str(date) + "_" + sample +"_"+aligner+".{ontfile}.alignment.benchmark.txt"

    threads: 
        config["Minimap2"]["minimap2_cores "],

    conda: "MASV_pipeline.yml"
  
    run:
        if aligner == "minimap2": #If the selected aligner is minimap2
            shell("mkdir -p {params.outdir}; cd {params.outdir}; minimap2 -t {threads} --MD -A {params.match_score_mm2} -B {params.mismatch_score_mm2} -u {params.GT_AG_find_mm2} -ax {params.mapping_tech_minimap} {input.ref} {input.reads} | samtools sort -@ {threads} -O BAM -o {output} 2> {log}")

        elif aligner == "ngmlr": #If the selected aligner is ngmlr
            shell("mkdir -p {params.outdir}; cd {params.outdir}; ngmlr -t {params.ngmlr_threads} -x {params.mapping_tech_ngmlr} -i {params.min_ident_ngmlr} \ 
            --match {params.match_score_ngmlr} --mismatch {params.mismatch_score_ngmlr} --gap-open {params.gap_open_score_ngmlr} --gap-extend-max {params.gap_extend_max_ngmlr} \ 
            --gap-extend-min {params.gap_extend_min_ngmlr} --kmer-skip {params.kmer_skip_ngmlr} -k {params.kmer_length_ngmlr} -r {input.ref} \ 
            -q {input.reads} | samtools sort-@ {threads} -O BAM -o {output} 2> {log}")
        else:
            shell("echo 'An error ocurred in the mapping step. Please, resubmit a valid aligner: minimap2, ngmlr.' > mapping.err; exit")


rule index_bam:
    input:
        BAM = lambda wildcards: expand(rules.mapping.output, ontfile=ontfiles.split(','))
        
    output:
        BAI = outdir + sample +"_"+ aligner+".{ontfile}.bam.bai"
        
    log:
        logs_dir + str(date) + "_" + aligner+".{ontfile}.index_bam.log"
        
    threads:
        config["Minimap2"]["minimap2_cores "]
    
    conda: "MASV_pipeline.yml"
        
    shell:
        "samtools index -@ {threads} {input.BAM} {output.BAI} 2> {log}"
        

rule alignment_stats:
    input:
        BAM = lambda wildcards: expand(rules.mapping.output, ontfile=ontfiles.split(','))
    output:
         outdir+sample+"_"+aligner+"/alignment_stats/{ontfile}.alignment_stats.txt"
    log:
        logs_dir + str(date) +".{ontfile}.alignment_stats.log"
        
    conda: "MASV_pipeline.yml"
    
    params: outdir + sample +"_"+aligner+"/alignment_stats/"
           
    shell:
        "mkdir -p {params}; python3 " + os.path.join(workflow.basedir, "lib/scr/alignment_stats.py") + " -o {output} {input.BAM} 2> {log}"
