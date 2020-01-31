######################
## SV CALLING RULES ##
######################


rule sniffles_calling:
    input:
        BAM = outdir + sample +"_"+aligner+".{ontfile}.bam",
        BAI = outdir + sample +"_"+ aligner+".{ontfile}.bam.bai"
    output:
        VCF = protected(svout+"/"+str(date)+"_Sniffles/{ontfile}/"+str(date)+"_"+sample+"_sniffles.{ontfile}.vcf"),
    params:
        prefix = "{ontfile}",
        min_sv_length_sniffles = config["Sniffles"]["sniffles_min_sv_length"],
        max_sv_length_sniffles = config["Sniffles"]["sniffles_max_sv_length"],
        min_read_length_sniffles = config["Sniffles"]["sniffles_min_read_length"],
        min_read_map_quality_sniffles =  config["Sniffles"]["sniffles_min_read_mapping_quality"],
        num_reads_report_sniffles = config["Sniffles"]["sniffles_num_reads_report"],
        max_num_splits_sniffles = config["Sniffles"]["sniffles_num_reads_report"],
        genotype_sniffles = config["Sniffles"]["sniffles_genotype"],
        cluster_sniffles = config["Sniffles"]["sniffles_cluster"],#not currently added in the sell command line
        min_homo_af_sniffles = config["Sniffles"]["sniffles_min_homo_af"],
        min_het_af_sniffles = config["Sniffles"]["sniffles_min_het_af"],
        min_supp = config["Sniffles"]["sniffles_min_read_support"]
        
    threads:
        config["Sniffles"]["sniffles_cores"] 
        
    log: logs_dir + str(date) + "_" + sample + "_sniffles.{ontfile}.svcalling.log"
    
    benchmark: benchmark_dir + str(date) + "_" + sample + "_sniffles.{ontfile}.sv.caller.benchmark.txt"
        
    conda: os.path.join(workflow.basedir, "MASV_pipeline.yml")
        
    shell:
      "mkdir -p "+svout+"/"+str(date)+"_Sniffles/{params.prefix}; sniffles -m {input.BAM} -s {params.min_supp} -t {threads} -r {params.min_read_length_sniffles} \
       -l {params.min_sv_length_sniffles} -d {params.max_sv_length_sniffles} \
       -q {params.min_read_map_quality_sniffles} -n {params.num_reads_report_sniffles} \
       {params.genotype_sniffles} --min_homo_af {params.min_homo_af_sniffles} \
       --min_het_af {params.min_het_af_sniffles} --report_read_strands -v {output.VCF} 2> {log}"
            
rule svim_calling:
    input:
        BAM = outdir + sample +"_"+aligner+".{ontfile}.bam",
        BAI = outdir + sample +"_"+ aligner+".{ontfile}.bam.bai",
        ref = config["Inputs"]["reference_genome"]
    output:
        outVCF = protected(svout+"/"+str(date)+"_"+sample+"_svim.{ontfile}/final_results.vcf")
    params:        
        min_sv_length_svim = config["Svim"]["svim_min_sv_length"],
        max_sv_lenght_svim = config["Svim"]["svim_max_sv_length"],
        min_read_length_svim = config["Svim"]["svim_min_read_length"],
        min_read_mapp_quality = config["Svim"]["svim_min_read_mapping_quality"],
        gap_tolerance_svim = config["Svim"]["svim_gap_tolerance"],
        overlap_tolerance_svim = config["Svim"]["svim_overlap_tolerance"],
        partition_max_distance_svim = config["Svim"]["svim_partition_max_distance"],
        sv_max_distance_svim = config["Svim"]["svim_sv_max_distance"],
        min_geno_score_svim = config["Svim"]["svim_min_geno_score"],
        homozygous_thresh_svim = config["Svim"]["svim_homozygous_thresh"],
        heterozygous_thresh_svim = config["Svim"]["svim_heterozygous_thresh"],
        min_depth_svim = config["Svim"]["svim_min_depth"]

    threads: 1

    log: logs_dir + str(date) + "_" + sample + "_svim.{ontfile}.svcalling.log"

    benchmark: benchmark_dir + str(date) + "_" + sample + "_svim.{ontfile}.sv.caller.benchmark.txt"

    conda: os.path.join(workflow.basedir, "MASV_pipeline.yml")

    shell:   
        "svim alignment "+svout+"/"+str(date)+"_"+sample+"_svim.{wildcards.ontfile}/ {input.BAM} {input.ref} \
         --min_sv_size {params.min_sv_length_svim} --max_sv_size {params.max_sv_lenght_svim} \
         --min_mapq {params.min_read_mapp_quality} --distance_normalizer {params.min_read_length_svim} \
         --segment_gap_tolerance {params.gap_tolerance_svim} --segment_overlap_tolerance {params.overlap_tolerance_svim} \
         --trans_partition_max_distance {params.partition_max_distance_svim} --trans_sv_max_distance {params.sv_max_distance_svim} \
         --minimum_score {params.min_geno_score_svim} --homozygous_threshold {params.homozygous_thresh_svim} \
         --heterozygous_threshold {params.heterozygous_thresh_svim} --minimum_depth {params.min_depth_svim} \
         --duplications_as_insertions  2> {log}"

rule filter_svim:
    input:
        rules.svim_calling.output.outVCF        
    output:
        protected(svout+"/"+str(date)+"_"+sample+"_svim.{ontfile}/final_results_minscore_"+str(minscore)+".vcf")    
    threads: 1
    
    params:
        minscore = config["Svim"]["svim_min_score"]
   
    log:
        logs_dir + str(date) + "_" + sample + ".{ontfile}.svim_filtering.log"
            
    conda: os.path.join(workflow.basedir, "MASV_pipeline.yml")
        
    shell:
        "grep -v \"hom_ref\" {input} | awk '{{ if($1 ~ /^#/) {{ print $0 }} \
        else {{ if($6>={params.minscore}) {{ print $0 }} }} }}' > {output} 2> {log}"