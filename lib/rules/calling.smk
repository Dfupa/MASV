######################
## SV CALLING RULES ##
######################


rule sniffles_calling:
    input:
        BAM = lambda wildcards: expand(rules.mapping.output, ontfile=ontfiles.split(',')),
        ref = config["Inputs"]["reference_genome"]
    output:
        VCF = protected(svout+"/"+str(date)+"_"+sample+"_sniffles.{ontfile}.vcf"),
        #outDIR = directory(svout+"/"+str(date)+"_"+sample+"_{input.svcaller}.{ontfile}/")
    params:
        outdir = directory (config["Outputs"]["svcall_out"]),
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
        config["Sniffles"]["sniffles_cores"] 
        
    log: logs_dir + str(date) + "_" + sample + "_sniffles.{ontfile}.svcalling.log"
    
    benchmark: benchmark_dir + str(date) + "_" + sample + "_sniffles.{ontfile}.sv.caller.benchmark.txt"
        
    conda: "MASV_pipeline.yml"
        
    shell:
        #if svcaller == "sniffles":    #If the selected caller is sniffles
            "mkdir -p "+svout+"/; sniffles -m {input.BAM} -v {output} -s {params.min_read_support_sniffles} -t {threads} \ 
            -r {params.min_read_length_sniffles} -l {params.min_sv_length_sniffles} -d {params.max_sv_length_sniffles} \
            -q {params.min_read_map_quality_sniffles} -n {params.num_reads_report_sniffles} \
            --genotype {params.genotype_sniffles} --min_homo_af {params.min_homo_af_sniffles} \
            --min_het_af {params.min_het_af_sniffles} --cluster {params.cluster_sniffles} --report_read_strands \
            -v {output.VCF} 2> {log}"
            
        #if svcaller == "svim":        #If the selected caller is svim
            #shell("mkdir -p "+svout+"/; svim alignment --min_sv_size {params.min_sv_length_svim} --max_sv_size {params.max_sv_length_svim} \
            #--min_mapq {params.min_read_mapp_quality} --distance_normalizer {params.min_read_length_svim} \
            #--segment_gap_tolerance {params.gap_tolerance_svim} --segment_overlap_tolerance {params.overlap_tolerance_svim} \
            #--trans_partition_max_distance{params.partition_max_distance_svim} --trans_sv_max_distance {params.sv_max_distance_svim} \
            #--minimum_score {params.svim_min_geno_score} --homozygous_threshold {params.homozygous_thresh_svim} \
            #--heterozygous_threshold {params.heterozygous_thresh_svim} --minimum_depth {params.min_depth_svim} \
            #--duplications_as_insertions {output.outDIR} {input.BAM} {input.ref} 2> {log}")

rule svim_calling:
    input:
        BAM = lambda wildcards: expand(rules.mapping.output, ontfile=ontfiles.split(',')),
        ref = config["Inputs"]["reference_genome"]
    output:
        outDIR = directory(svout+"/"+str(date)+"_"+sample+"_svim.{ontfile}/")
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
        heterozygous_thresh_svim = config["Svim"]["svim_heterozygous_thresh"]

    threads: 1

    log: logs_dir + str(date) + "_" + sample + "_{input.svcaller}.{ontfile}.svcalling.log"

    benchmark: benchmark_dir + str(date) + "_" + sample + "_{input.svcaller}.{ontfile}.sv.caller.benchmark.txt"

    conda: "MASV_pipeline.yml"

    shell:   
        "mkdir -p "+svout+"/; svim alignment --min_sv_size {params.min_sv_length_svim} --max_sv_size {params.max_sv_length_svim} \
        --min_mapq {params.min_read_mapp_quality} --distance_normalizer {params.min_read_length_svim} \
        --segment_gap_tolerance {params.gap_tolerance_svim} --segment_overlap_tolerance {params.overlap_tolerance_svim} \
        --trans_partition_max_distance{params.partition_max_distance_svim} --trans_sv_max_distance {params.sv_max_distance_svim} \
        --minimum_score {params.svim_min_geno_score} --homozygous_threshold {params.homozygous_thresh_svim} \
        --heterozygous_threshold {params.heterozygous_thresh_svim} --minimum_depth {params.min_depth_svim} \
        --duplications_as_insertions {output.outDIR} {input.BAM} {input.ref} 2> {log}"

rule filter_svim:

    input:
        lambda wildcards: expand(rules.svim_calling.output.outDIR + "final_results.vcf", ontfile=ontfiles.split(','))
        
    output:
        protected(rules.svim_calling.output.outDIR + "{ontfile]_minscore_{params.minscore,[0-9]+}.vcf") #Note that a numeric constraint is added
    
    threads: 1
    
    params:
        minscore = config["Svim"]["svim_min_score"]
   
    log:
        logs_dir + str(date) + "_" + sample + ".{ontfile}.svim_filtering.log"
            
    conda: "MASV_pipeline.yml"
        
    run:
        if rules.sv_calling.input.svcaller == "svim":
            shell("grep -v \"hom_ref\" {input} | awk '{{ if($1 ~ /^#/) {{ print $0 }} \
            else {{ if($6>={params.minscore}) {{ print $0 }} }} }}' > {output} 2> {log}")
        else:
            shell("echo 'The filter_svim step is aborted as it is not necessary for Sniffles calling. Continuing.' 2> {log}")
