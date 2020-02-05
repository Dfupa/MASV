
#Author: Diego Fuentes
#Contact email: diegofupa@gmail.com
#Date:2019-11-04


import os

###################
# ALIGNMENT RULES #
###################


if aligner == "minimap2":
    rule mapping_minimap2:
        input:
            reads = expand(config["Inputs"]["ONT_reads_directory"] + "{ontfile}.{format}", ontfile=ontfiles.split(','), format=["fastq", "fastq.gz"]),
            ref = config["Inputs"]["reference_genome"]
        output:
            protected(outdir + sample +"_"+aligner+".{ontfile}.bam")

        params:
            mapping_tech_minimap = which_tech(parameter=config["Parameters"]["seq_technology"])[0],
            outdir = config["Outputs"]["alignment_out"],
            match_score_mm2 = config["Minimap2"]["minimap2_match_score"],
            mismatch_score_mm2 = config["Minimap2"]["minimap2_mismatch_score"],
            gap_open_score_mm2 = config["Minimap2"]["minimap2_gap_open_score"],
            GT_AG_find_mm2 = config["Minimap2"]["GT_AG_find"],

        log:
            logs_dir + str(date) + "_" + sample +"_"+aligner+".{ontfile}.alignment_minimap2.log"

        benchmark:
            benchmark_dir + str(date) + "_" + sample +"_"+aligner+".{ontfile}.alignment_minimap2.benchmark.txt"

        threads: 
            config["Minimap2"]["minimap2_cores "],

        conda: os.path.join(workflow.basedir, "MASV_pipeline.yml")

        shell:
            "mkdir -p {params.outdir}; cd {params.outdir}; minimap2 -t {threads} --MD -A {params.match_score_mm2} -B {params.mismatch_score_mm2} -u {params.GT_AG_find_mm2} -ax {params.mapping_tech_minimap} {input.ref} {input.reads} | samtools sort -@ {threads} -O BAM -o {output} 2> {log}"


elif aligner == "ngmlr":
    rule mapping_ngmlr:
        input:
            reads = expand(config["Inputs"]["ONT_reads_directory"] + "{ontfile}.{format}", ontfile=ontfiles.split(','), format=["fastq", "fastq.gz"]),
            ref = config["Inputs"]["reference_genome"]
        output:
            protected(outdir + sample +"_"+ aligner + ".{ontfile}.bam")

        params:
            mapping_tech_ngmlr = which_tech(parameter=config["Parameters"]["seq_technology"])[1],
            min_ident_ngmlr = config["Ngmlr"]["ngmlr_min_ident"],
            match_score_ngmlr = config["Ngmlr"]["ngmlr_match_score"],
            mismatch_score_ngmlr = config["Ngmlr"]["ngmlr_mismatch_score"],
            gap_open_score_ngmlr = config["Ngmlr"]["ngmlr_gap_open_score"],
            gap_extend_max_ngmlr = config["Ngmlr"]["ngmlr_gap_extend_max"],
            gap_extend_min_ngmlr = config["Ngmlr"]["ngmlr_gap_extend_min"],
            kmer_length_ngmlr = config["Ngmlr"]["ngmlr_kmer_length"],
            kmer_skip_ngmlr = config["Ngmlr"]["ngmlr_kmer_length"],

        log:
            logs_dir + str(date) + "_" + sample +"_"+aligner+".{ontfile}.alignment_ngmlr.log"

        benchmark:
            benchmark_dir + str(date) + "_" + sample +"_"+aligner+".{ontfile}.alignment_ngmlr.benchmark.txt"

        threads: 
            config["Ngmlr"]["ngmlr_cores"]

        conda: os.path.join(workflow.basedir, "MASV_pipeline.yml")

        shell:
            "mkdir -p {params.outdir}; cd {params.outdir}; ngmlr -t {threads} -x {params.mapping_tech_ngmlr} -i {params.min_ident_ngmlr} --match {params.match_score_ngmlr} --mismatch {params.mismatch_score_ngmlr} --gap-open {params.gap_open_score_ngmlr} --gap-extend-max {params.gap_extend_max_ngmlr} --gap-extend-min {params.gap_extend_min_ngmlr} --kmer-skip {params.kmer_skip_ngmlr} -k {params.kmer_length_ngmlr} -r {input.ref} -q {input.reads} | samtools sort-@ {threads} -O BAM -o {output} 2> {log}"

rule index_bam:
    input:
        BAM = outdir + sample +"_"+ aligner + ".{ontfile}.bam"
        
    output:
        BAI = protected(outdir + sample +"_"+ aligner+".{ontfile}.bam.bai")
        
    log:
        logs_dir + str(date) + "_" + aligner+".{ontfile}.index_bam.log"
        
    threads:
        config["Minimap2"]["minimap2_cores "]
    
    conda: os.path.join(workflow.basedir, "MASV_pipeline.yml")
        
    shell:
        "samtools index -@ {threads} {input.BAM} {output.BAI} 2> {log}"
        

rule alignment_stats:
    input:
        BAM = outdir + sample +"_"+ aligner + ".{ontfile}.bam",
        BAI = outdir + sample +"_"+ aligner+".{ontfile}.bam.bai"
    output:
         outdir + sample + "_"+ aligner + "/alignment_stats/{ontfile}.alignment_stats.txt"
    log:
        logs_dir + str(date) +".{ontfile}.alignment_stats.log"
        
    conda: os.path.join(workflow.basedir, "MASV_pipeline.yml")
    
    params: outdir + sample +"_"+aligner+"/alignment_stats/"
           
    shell:
        "mkdir -p {params}; python3 " + os.path.join(workflow.basedir, "lib/scr/BAM_alignment_stats.py") + " -o {output} {input.BAM} 2> {log}"
