######################
## SV CALLING RULES ##
######################

rule sv_calling:
    input:
        BAM = rules.sort_bam.output.outBAM,
        svcaller = config["Inputs"]["svcaller_selection"],
        ref = config["Inputs"]["reference_genome"]
    output:
        VCF = workingdir + str(date) + "/{params.outdir}/{sample}_{input.svcaller}.VCF",
        outDIR = workingdir + str(date) + "/{params.outdir}/)
    params:
        outdir = directory (config["Outputs"]["sv_call_out"]),
        min_sv_length_sniffles = config["Sniffles"]["sniffles_min_sv_length"],
        max_sv_length_sniffles = config["Sniffles"]["sniffles_max_sv_length"],
        min_read_length_sniffles = config["Sniffles"]["sniffles_min_read_length"],
        min_read_map_quality_sniffles =  config["Sniffles"]["sniffles_min_read_mapping_quality"],
        min_read_support_sniffles = config["Sniffles"]["sniffles_min_read_support"],
        num_reads_report_sniffles = config["Sniffles"]["sniffles_num_reads_report"],
        max_num_splits_sniffles = config["Sniffles"]["sniffles_num_reads_report"],
        genotype_sniffles = config["Sniffles"]["sniffles_genotype"],
        cluster_sniffles = config["Sniffles"]["sniffles_cluster"],
        min_homo_af_sniffles = config["Sniffles"]["sniffles_min_homo_af"],
        min_het_af_sniffles = config["Sniffles"]["sniffles_min_het_af"],
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
        heterozygous_thresh_svim = config["Svim"]["svim_heterozygous_thresh"]
        
    threads:
        sniffles_threads = config["Sniffles"]["sniffles_cores"] 
        
    logs: logs_dir + str(date) + ".{ontfile}.svcalling.log"
    
    benchmark: benchmark_dir + str(date) + ".{ontfile}.sv.caller.benchmark.txt"
        
    conda: "pipeline_env.yml"
        
    run:
        if svcaller == "sniffles":    #If the selected caller is sniffles
            shell("sniffles -m {input.BAM} -v {output} -s {params.min_read_support_sniffles} \ 
            -r {params.min_read_length_sniffles} -l {params.min_sv_length_sniffles} -d {params.max_sv_length_sniffles} \
            -q {params.min_read_map_quality_sniffles} -n {params.num_reads_report_sniffles} \
            --genotype {params.genotype_sniffles} --min_homo_af {params.min_homo_af_sniffles} \
            --min_het_af {params.min_het_af_sniffles} --cluster {params.cluster_sniffles} --report_read_strands \
            -v {output.VCF} 2> {log}")
            
        if svcaller == "svim":        #If the selected caller is svim
            shell("svim alignment --min_sv_size {params.min_sv_length_svim} --max_sv_size {params.max_sv_length_svim} \
            --min_mapq {params.min_read_mapp_quality} --distance_normalizer {params.min_read_length_svim} \
            --segment_gap_tolerance {params.gap_tolerance_svim} --segment_overlap_tolerance {params.overlap_tolerance_svim} \
            --trans_partition_max_distance{params.partition_max_distance_svim} --trans_sv_max_distance {params.sv_max_distance_svim} \
            --minimum_score {params.svim_min_geno_score} --homozygous_threshold {params.homozygous_thresh_svim} \
            --heterozygous_threshold {params.heterozygous_thresh_svim} --minimum_depth {params.min_depth_svim} \
            --duplications_as_insertions {output.outDIR} {input.BAM} {input.ref} 2> {log}")

rule filter_svim:
    input:
        rules.sv_calling.output.outDIR + "final_results.vcf"
        
    output:
        rules.sv_calling.output.outDIR + "minscore_{params.minscore,[0-9]+}.vcf" #Note that a numeric constraint is added
    
    threads: 1
    
    params:
        minscore = config["Svim"]["svim_mins_core"]
   
   log:
        logs: logs_dir + str(date) + ".{ontfile}.svim_filtering.log"
            
    conda: "pipeline_env.yml"
        
    run:
        if rules.sv_calling.input.svcaller == "svim"
            shell("grep -v \"hom_ref\" {input} | awk '{{ if($1 ~ /^#/) {{ print $0 }} else {{ if($6>={params.minscore}) {{ print $0 }} }} }}' > {output}")
        else:
            shell("echo 'The filter_svim step is aborted as it is not necessary for Sniffles calling. Continuing.'")
