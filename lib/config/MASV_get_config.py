#!/usr/bin/env python3
import os
import json
import argparse
import sys
import re

#Author: Diego Fuentes
#Contact email: diegofupa@gmail.com
#Date:2019-10-01

ont_barcodes = []    #List to store the ont barcodes

#######################
###CONFIG FILE CLASS###
#######################
class CreateConfigurationFile(object):
    """Class which manages Configuration file Manager"""
      
    def __init__(self):
        """Class constructor of the pipeline"""
        #GENERAL PARAMETERS

        self.configFile = None                               #Name of the json configuration file to be created.
        self.version = 1                                     #Pipeline version
        self.logs_dir = "logs"                               #Directory to keep all the log files
        self.sample_barcode = None                           #Sample barcode 
        self.basedir = self.sample_barcode                   #Base directory for the pipeline run
        self.seq_technology = "nanopore"                     #Sequencing technology
        self.single = False                                  #Parameter that is going to be used for the helper function

        #INPUT PARAMETERS

        self.ONT_reads_directory = None          #Directory where the ont fastqs are stored
        self.reference_genome = None             #Reference genome provided in .fa or .fa.gz format
        self.hq_vcf = None                       #Provided high confidence .vcf.gz for truvari evaluation (Alpha)
        self.aligner_selection = "minimap2"      #Default aligner
        self.svcaller_selection = "svim"         #Default sv caller


        #OUTPUT PARAMETERS

        self.alignment_out = "Alignment"                     #Out directory of the alignment step
        self.sv_call_out =  "SV-calling"                     #Out directory of the sv calls
    
        #WILDCARD PARAMETER

        self.ONT_fastqs = None                               #List with basename of the ONT fastqs if any

        #MINIMAP2 PARAMETERS

        self.minimap2_cores = 4                  #Number of threads to run the minimap2 aligner
        self.minimap2_kmer_length = 15           #K-mer lenght in bases <10-28>
        self.minimap2_match_score = 2            #Match score
        self.minimap2_mismatch_score = 5         #Mismatch penalty
        self.minimap2_gap_open_score = 5         #Gap open penalty <4-24>
        self.GT_AG_find = "n"                    #how to find GT-AG. f:transcript strand, b:both strands, n:don't match GT-AG

        #NGMLR PARAMETERS

        self.ngmlr_cores = 8                     #Number of threads to run the ngmlr aligner
        self.ngmlr_min_ident = 0.75              #Alignments with an identity lower than this threshold will be discarded
        self.ngmlr_match_score = 2               #Match score
        self.ngmlr_mismatch_score = -5           #Mismatch score
        self.ngmlr_gap_open_score = -5           #Gap open score
        self.ngmlr_gap_extend_max = -5           #Gap extend max
        self.ngmlr_gap_extend_min = -1           #Gap extend min
        self.ngmlr_kmer_length = 13              #K-mer lenght in bases <10-15>
        self.ngmlr_kmer_skip = 2                 #Number of k-mers to skip when building the lookup table from the reference

        #SNIFFLES PARAMETERS

        self.sniffles_cores = 4                     #Number of threads to run the sniffles SV caller
        self.sniffles_min_sv_length = 40            #Minimum SV length
        self.sniffles_max_sv_length = 100000        #Maximum SV length
        self.sniffles_min_read_length = 1000        #Minimum read length
        self.sniffles_min_read_mapping_quality = 20 #Min mapping quality. Reads will lower mapping quality will be discarded
        self.sniffles_min_read_support = 10         #Minimum read support required to call a SV
        self.sniffles_num_reads_report = 0          #Report up to N reads that support the SV in the vcf file. -1: report all.
        self.sniffles_max_num_splits = 7            #Maximum number of splits per read to be still taken into account
        self.sniffles_genotype = False              #Enables Sniffles to compute the genotypes.
        self.sniffles_cluster = False               #Enables Sniffles to phase SVs that occur on the same reads
        self.sniffles_min_homo_af = 0.8             #Minimum variant allele frequency to be called as homozygous
        self.sniffles_min_het_af = 0.2              #Minimum variant allele frequency to be called as heterozygous

        #SVIM PARAMETERS

        self.svim_min_sv_length = 40                #Minimum SV length
        self.svim_max_sv_length = 100000            #Maximum SV length
        self.svim_min_read_length = 1000            #Minimum read length
        self.svim_min_read_mapping_quality = 20     #Min mapping quality. Reads will lower mapping quality will be discarded
        self.svim_gap_tolerance = 10                #Maximum tolerated gap between adjacent alignment segments
        self.svim_overlap_tolerance = 5             #Maximum tolerated overlap between adjacent alignment segments
        self.svim_partition_max_distance = 200      #Maximum distance in bp between translocation breakpoints in a partition
        self.svim_sv_max_distance = 500             #Maximum distance in bp between a translocation breakpoint and an SV signature to be combined
        self.svim_min_geno_score = 3                #Minimum score for genotyping
        self.svim_homozygous_thresh = 0.8           #Minimum variant allele frequency to be called as homozygous
        self.svim_heterozygous_thresh = 0.2         #Minimum variant allele frequency to be called as heterozygous
        self.svim_min_depth = 4                     #Minimum total read depth for genotyping
        self.svim_min_score = 10                    #Minimum quality score to filter SVIM
###
        #DICTIONARIES
        self.allParameters = {}
        self.generalParameters = {}
        self.inputParameters = {}
        self.outputParameters = {}
        self.wildcardParameters = {}
        self.minimap2Parameters = {}
        self.ngmlrParameters = {}
        self.snifflesParameters = {}
        self.svimParameters = {}
        
####

    def register_parameter(self, parser):
        """Register all parameters with the given
        argparse parser"""
        self.register_general(parser)
        self.register_input(parser)
        self.register_output(parser)
        self.register_wildcards(parser)
        self.register_minimap2(parser)
        self.register_ngmlr(parser)
        self.register_sniffles(parser)
        self.register_svim(parser)

    def register_general(self, parser):
        """Register all general parameters with the given
        argparse parser

        parser -- the argparse parser
        """
        general_group = parser.add_argument_group('General Parameters')
        general_group.add_argument('--configFile', dest="configFile", metavar="configFile", help='Configuration JSON to be generated. Default %s.' % self.configFile)
        general_group.add_argument('--version', type=int, dest="version", metavar="version", default=self.version, help='Pipeline run version. Default %s.' % self.version)
        general_group.add_argument('--logs-dir', dest="logs_dir", metavar="logs_dir", help='Directory to keep all the log files. Default sample_barcode id.')
        general_group.add_argument('--sample-barcode', dest="sample_barcode", metavar="sample_barcode", help='Sample barcode. Default %s.' % self.sample_barcode)
        general_group.add_argument('--basedir', dest="basedir", metavar="basedir", help='Base directory for the pipeline run. Default %s.' % self.basedir)
        general_group.add_argument('--single', dest="single", type=bool, default=self.single, help='Parameter used for the helper function find_files. Default %s.' % self.single)
        general_group.add_argument('--sequencing-technology', type=str, dest="seq_technology", metavar="seq_technology", default=self.seq_technology, help='Parameter used for determining the sequencing technology ("nanopore" or "pacbio"). Default %s.' % self.seq_technology)

    def register_input(self, parser):
        """Register all input parameters with the given
        argparse parser

        parser -- the argparse parser
        """
        input_group = parser.add_argument_group('Inputs')
        input_group.add_argument('--ont-reads-directory', dest="ONT_reads_directory", metavar="ONT_reads_directory", help='Directory where the ont fastqs are stored. Default %s.' % self.ONT_reads_directory)
        input_group.add_argument('--reference-genome', dest="reference_genome", metavar="reference_genome", help='Reference genome provided in .fa or .fa.gz format. Your path is  %s.' % self.reference_genome)
        input_group.add_argument('--hq-vcf', dest="hq_vcf", metavar="hq_vcf", help='Provided high confidence .vcf.gz for truvari evaluation (Alpha). Your path is  %s.' % self.hq_vcf) 
        input_group.add_argument('--aligner-selection', dest="aligner_selection", metavar="aligner_selection", default=self.aligner_selection, help='Selects the aligner to be used in the pipeline. Default "%s".' % self.aligner_selection)
        input_group.add_argument('--sv_caller-selection', dest="svcaller_selection",  metavar="svcaller_selection", default=self.svcaller_selection, help='Selects the SV caller to be used in the pipeline. Default "%s".' % self.svcaller_selection)


    
    def register_output(self, parser):
        """Register all output parameters with the given
        argparse parser

        parser -- the argparse parser
        """

        output_group = parser.add_argument_group('Outputs')
        output_group.add_argument('--alignment-out', dest="alignment_out", help='Out directory of the alignment step. Default "/%s"' % self.alignment_out)
        output_group.add_argument('--sv-call-out', dest="sv_call_out", help='Out directory of the sv calls. Default "/%s"' % self.sv_call_out)


    def register_wildcards(self, parser):
        """Register all wildcards parameters with the given
        argparse parser

        parser -- the argparse parser
        """
        wildcards_group = parser.add_argument_group('Wildcards')
        wildcards_group.add_argument('--ONT-fastqs', dest="ONT_fastqs", metavar="ONT_fastqs", help='List with basename of the ONT fastqs. Default %s' % self.ONT_fastqs)

    def register_minimap2(self, parser):
        """Register all minimap2 aligner parameters with the given
        argparse parser

        parser -- the argparse parser
        """
        minimap2_group = parser.add_argument_group('Minimap2 parameters')
        minimap2_group.add_argument('--minimap2-cores', type=int, dest="minimap2_cores", metavar="minimap2-cores", default=self.minimap2_cores, help='Number of threads to run the minimap2 aligner. Default %s.' % self.minimap2_cores)
        minimap2_group.add_argument('--minimap2-kmer-length', type=int, dest="minimap2_kmer_length", metavar="minimap2-kmer-length", default=self.minimap2_kmer_length, help='K-mer lenght in bases <10-28>. Default %s.' % self.minimap2_kmer_length)
        minimap2_group.add_argument('--minimap2-match-score', type = int, dest="minimap2_match_score", metavar="minimap2-match-score", default=self.minimap2_match_score, help='Match score. Default %s.' % self.minimap2_match_score)
        minimap2_group.add_argument('--minimap2-mismatch-penalty', type = int, dest="minimap2_mismatch_score", metavar="minimap2-mismatch-penalty", default=self.minimap2_mismatch_score, help='Mismatch penalty. Default %s.' % self.minimap2_mismatch_score)
        minimap2_group.add_argument('--minimap2-gap-open-penalty', type = int, dest="minimap2_gap_open_score", metavar="minimap2-gap-open-penalty", default=self.minimap2_gap_open_score, help='Gap open penalty <4-24>. Default %s.' % self.minimap2_gap_open_score)
        minimap2_group.add_argument('--GT-AG-find', type = str, dest="GT_AG_find", metavar="GT-AG-find", default=self.GT_AG_find, help='how to find GT-AG. f:transcript strand, b:both strands, n:dont match GT-AG. Default "%s".' % self.GT_AG_find)

    def register_ngmlr(self, parser):
        """Register all ngmlr aligner parameters with the given
        argparse parser

        parser -- the argparse parser
        """
        ngmlr_group = parser.add_argument_group('Ngmlr parameters')
        ngmlr_group.add_argument('--ngmlr-cores', type=int, dest="ngmlr_cores", metavar="ngmlr-cores", default=self.ngmlr_cores, help='Number of threads to run the ngmlr aligner. Default %s.' % self.ngmlr_cores)
        ngmlr_group.add_argument('--ngmlr-min-ident', type=float, dest="ngmlr_min_ident", metavar="ngmlr-min-ident", default=self.ngmlr_min_ident, help='Alignments with an identity lower than this threshold will be discarded. Default %s.' % self.ngmlr_min_ident)
        ngmlr_group.add_argument('--ngmlr-kmer-length', type=int, dest="ngmlr_kmer_length", metavar="ngmlr-kmer-length", default=self.ngmlr_kmer_length, help='K-mer lenght in bases <10-15>. Default %s.' % self.ngmlr_kmer_length)
        ngmlr_group.add_argument('--ngmlr-kmer-skip', type=int, dest="ngmlr_kmer_skip", metavar="ngmlr-kmer-skip", default=self.ngmlr_kmer_skip, help='Number of k-mers to skip when building the lookup table from the reference. Default %s.' % self.ngmlr_kmer_skip)
        ngmlr_group.add_argument('--ngmlr-match-score', type = int, dest="ngmlr_match_score", metavar="ngmlr-match-score", default=self.ngmlr_match_score, help='Match score. Default %s.' % self.ngmlr_mismatch_score)
        ngmlr_group.add_argument('--ngmlr-mismatch-score', type = int, dest="ngmlr_mismatch_score", metavar="ngmlr-mismatch-score", default=self.ngmlr_mismatch_score, help='Mismatch penalty. Default %s.' % self.ngmlr_mismatch_score)
        ngmlr_group.add_argument('--ngmlr-gap-open-score', type = int, dest="ngmlr_gap_open_score", metavar="ngmlr-gap-open-score", default=self.ngmlr_gap_open_score, help='Gap open penalty. Default %s.' % self.ngmlr_gap_open_score)
        ngmlr_group.add_argument('--ngmlr-gap-extend-max', type = int, dest="ngmlr_gap_extend_max", metavar="ngmlr-gap-extend-max", default=self.ngmlr_gap_extend_max, help='Gap extend max. Default %s.' % self.ngmlr_gap_extend_max)
        ngmlr_group.add_argument('--ngmlr-gap-extend-min', type = int, dest="ngmlr_gap_extend_min", metavar="ngmlr-gap-extend-min", default=self.ngmlr_gap_extend_min, help='Gap extend min. Default %s.' % self.ngmlr_gap_extend_min)

    def register_sniffles(self, parser):
        """Register all sniffles sv caller parameters with the given
        argparse parser

        parser -- the argparse parser
        """
        sniffles_group = parser.add_argument_group('Sniffles parameters')
        sniffles_group.add_argument('--sniffles-cores', type=int, dest="sniffles_cores", metavar="sniffles-cores", default=self.sniffles_cores, help='Number of threads to run the sniffles SV caller. Default %s.' % self.sniffles_cores)
        sniffles_group.add_argument('--sniffles-min-sv-length', type=int, dest="sniffles_min_sv_length", metavar="sniffles-min-sv-length", default=self.sniffles_min_sv_length, help='Minimum SV length. Default %s.' % self.sniffles_min_sv_length)
        sniffles_group.add_argument('--sniffles-max-sv-length', type=int, dest="sniffles_max_sv_length", metavar="sniffles-max-sv-length", default=self.sniffles_max_sv_length, help='Maximum SV length. Default %s.' % self.sniffles_max_sv_length)
        sniffles_group.add_argument('--sniffles-min-read-length', type=int, dest="sniffles_min_read_length", metavar="sniffles-min-read-length", default=self.sniffles_min_read_length, help='Minimum read length. Discard read if non of its segment is larger then this. Default %s.' % self.sniffles_min_read_length)
        sniffles_group.add_argument('--sniffles-min-read-mapping-quality', type=int, dest="sniffles_min_read_mapping_quality", metavar="sniffles-min-read-mapping-quality", default=self.sniffles_min_read_mapping_quality, help='Min mapping quality. Reads will lower mapping quality will be discarded. Default %s.' % self.sniffles_min_read_mapping_quality)
        sniffles_group.add_argument('--sniffles-min-read-support', type=int, dest="sniffles_min_read_support", metavar="sniffles-min-read-support", default=self.sniffles_min_read_support, help='Minimum read support required to call a SV. Default %s.' % self.sniffles_min_read_support)
        sniffles_group.add_argument('--sniffles-num-reads-report', type=int, dest="sniffles_num_reads_report", metavar="sniffles-num-reads-report", default=self.sniffles_num_reads_report, help='Report up to N reads that support the SV in the vcf file. -1: report all. Default %s.' % self.sniffles_num_reads_report)
        sniffles_group.add_argument('--sniffles-max-num-splits', type=int, dest="sniffles_max_num_splits", metavar="sniffles-max-num-splits", default=self.sniffles_max_num_splits, help='Maximum number of splits per read to be still taken into account. Default %s.' % self.sniffles_max_num_splits)
        sniffles_group.add_argument('--sniffles-genotype', type=bool, dest="sniffles_genotype", metavar="sniffles-genotype", default=self.sniffles_genotype, help='Enables Sniffles to compute the genotypes. Default "%s".' % self.sniffles_genotype)
        sniffles_group.add_argument('--sniffles-cluster', type=bool, dest="sniffles_cluster", metavar="sniffles-cluster", default=self.sniffles_cluster, help='Enables Sniffles to phase SVs that occur on the same reads. Default "%s".' % self.sniffles_cluster)
        sniffles_group.add_argument('--sniffles-min-homo-af', type=int, dest="sniffles_min_homo_af", metavar="sniffles-min-homo-af", default=self.sniffles_min_homo_af, help='Minimum variant allele frequency to be called as homozygous. Default %s.' % self.sniffles_min_homo_af)
        sniffles_group.add_argument('--sniffles-min-het-af', type=int, dest="sniffles_min_het_af", metavar="sniffles-min-het-af", default=self.sniffles_min_het_af, help='Minimum variant allele frequency to be called as heterozygous. Default %s.' % self.sniffles_min_het_af)


    def register_svim(self, parser):
        """Register all svim sv caller parameters with the given
        argparse parser

        parser -- the argparse parser
        """
        svim_group = parser.add_argument_group('Svim parameters')
        svim_group.add_argument('--svim-min-sv-length', type = int, dest="svim_min_sv_length", metavar="svim-min-sv-length", default=self.svim_min_sv_length, help='Minimum SV length to run the svim SV caller. Default %s.' % self.svim_min_sv_length)
        svim_group.add_argument('--svim-max-sv-length', type = int, dest="svim_max_sv_length", metavar="svim-max-sv-length", default=self.svim_max_sv_length, help='Maximum SV length to run the svim SV caller. Default %s.' % self.svim_max_sv_length)
        svim_group.add_argument('--svim-min-read-length', type = int, dest="svim_min_read_length", metavar="svim-min-read-length", default=self.svim_min_read_length, help='Minimum read length to run the svim SV caller. Default %s.' % self.svim_min_read_length)
        svim_group.add_argument('--svim-min-read-mapping-quality', type = int, dest="svim_min_read_mapping_quality", metavar="svim-min-read-mapping-quality", default=self.svim_min_read_mapping_quality, help='Minimum mapping quality. Reads will lower mapping quality will be discarded. Default %s.' % self.svim_min_read_mapping_quality)
        svim_group.add_argument('--svim-gap-tolerance', type = int, dest="svim_gap_tolerance", metavar="svim-gap-tolerance", default=self.svim_gap_tolerance, help='Maximum tolerated gap between adjacent alignment segments. Default %s' % self.svim_gap_tolerance)
        svim_group.add_argument('--svim-overlap-tolerance', type = int, dest="svim_overlap_tolerance", metavar="svim-overlap-tolerance", default=self.svim_overlap_tolerance, help='Maximum tolerated overlap between adjacent alignment segments. This parameter applies to overlaps on the reference and the read. Default %s.' % self.svim_overlap_tolerance)
        svim_group.add_argument('--svim-partition-max-distance', type = int, dest="svim_partition_max_distance", metavar="svim-partition-max-distance", default=self.svim_partition_max_distance, help='Maximum distance in bp between translocation breakpoints in a partition. Default %s.' % self.svim_partition_max_distance)
        svim_group.add_argument('--svim-sv-max-distance', type = int, dest="svim_sv_max_distance", metavar="svim-sv-max-distance", default=self.svim_sv_max_distance, help='Maximum distance in bp between a translocation breakpoint and an SV signature to be combined. Default %s.' % self.svim_sv_max_distance)
        svim_group.add_argument('--svim-min-geno-score', type = int, dest="svim_min_geno_score", metavar="svim-min-geno-score", default=self.svim_min_geno_score, help='Minimum score for genotyping.Default %s.' % self.svim_min_geno_score)
        svim_group.add_argument('--svim-homozygous-thresh', type = float, dest="svim_homozygous_thresh", metavar="svim-homozygous-thresh", default=self.svim_homozygous_thresh, help='Minimum variant allele frequency to be called as homozygous. Default %s.' % self.svim_homozygous_thresh)
        svim_group.add_argument('--svim-heterozygous-thresh', type = float, dest="svim_heterozygous_thresh", metavar="svim-heterozygous-thresh", default=self.svim_heterozygous_thresh, help='Minimum variant allele frequency to be called as heterozygous. Default %s.' % self.svim_heterozygous_thresh)
        svim_group.add_argument('--svim-min-depth', type = int, dest="svim_min_depth", metavar="svim-min-depth", default=self.svim_min_depth, help='Minimum total read depth for genotyping. Default %s.' % self.svim_min_depth)
        svim_group.add_argument('--svim-min-score', type = int, dest="svim_min_score", metavar="svim-min-score", default=self.svim_min_score, help='Minimum quality score used to filter SVIM results. Default "%s".' % self.svim_min_score)

####

    def check_parameters(self,args):
        """Check parameters consistency
            
        args -- set of parsed arguments"""

        if args.configFile==None:
            parser.print_help()
            sys.exit(-1)

        working_dir = os.getcwd() + "/"

        if args.sample_barcode == None:
            print("No sample_barcode specified. A barcode or identification is required")
            parser.print_help()
            sys.exit(-1)

        if args.basedir:
            args.basedir = os.path.abspath(args.basedir) + "/"
        else: 
            args.basedir = working_dir + "v" + str(args.version) + "/" + args.sample_barcode + "/"

        if args.logs_dir:
            args.logs_dir = os.path.abspath(args.logs_dir) + "/"
        else:
            args.logs_dir = args.basedir + self.logs_dir + "/"
            

        if args.ONT_reads_directory:
            args.ONT_reads_directory = os.path.abspath(args.ONT_reads_directory) + "/"
        else:
            args.ONT_reads_directory =  working_dir + "reads/ont/" + args.sample_barcode + "/"
        if not os.path.exists(args.ONT_reads_directory):
            print(args.ONT_reads_directory + " not found. The directory where the reads are located is required. Exiting now.")
            parser.print_help()
            sys.exit(-1)
            
        if args.reference_genome:
            args.reference_genome = os.path.abspath(args.reference_genome)+ "/"
        else:
            args.reference_genome = working_dir + "reference/genome.fa"

        if not os.path.exists(args.reference_genome):
            print("The reference genome has been not provided or it has not been found in "+args.reference_genome+". Exiting now")
            parser.print_help()
            sys.exit(-1)

        if args.hq_vcf:
            args.hq_vcf = os.path.abspath(args.hq_vcf)
        else:
            args.hq_vcf = working_dir + "truth_dataset/hq.sv.vcf.gz"
        if not os.path.exists(args.hq_vcf):
            print("The high confidence .vcf.gz has not been provided in the path "+args.hq_vcf+" . Note that the truvari evaluation will not be completed")

        if args.alignment_out:
            args.alignment_out = os.path.abspath(args.alignment_out) + "/"
        else:
            args.alignment_out = args.basedir + self.alignment_out + "/"

        if args.sv_call_out:
            args.sv_call_out = os.path.abspath(args.sv_call_out) + "/"
        else:
            args.sv_call_out = args.alignment_out  + self.sv_call_out + "/"

        ##Assign wildcards

        if args.ONT_fastqs == None:
            for r, d, f in os.walk(args.ONT_reads_directory):
                for file in f:
                    if re.search('.1.fastq.gz', file):
                        a = file.replace('.1.fastq.gz','')
                        ont_barcodes.append(a)
                        if args.ONT_fastqs == None:
                            args.ONT_fastqs = a
                        else:
                            args.ONT_fastqs += "," + a               

###

    def storeGeneralParameters(self,args):
        """Updates general parameters to the map of parameters to be store in a JSON file

        args -- set of parsed arguments
        """
        self.generalParameters["configFile"] = args.configFile
        self.generalParameters["version"] = args.version
        self.generalParameters["basedir"] = args.basedir
        self.generalParameters["logs_dir"] = args.logs_dir
        self.generalParameters["sample_barcode"] = args.sample_barcode
        self.generalParameters["single"] = args.single
        self.generalParameters["seq_technology"] = args.seq_technology
        self.allParameters["Parameters"] = self.generalParameters

    def storeInputParameters(self,args):
        """Updates input parameters to the map of parameters to be store in a JSON file

        args -- set of parsed arguments
        """

        self.inputParameters["aligner_selection"] = args.aligner_selection
        self.inputParameters["svcaller_selection"] = args.svcaller_selection
        self.inputParameters["ONT_reads_directory"] = args.ONT_reads_directory
        self.inputParameters["reference_genome"] = args.reference_genome
        self.inputParameters["hq_vcf"] = args.hq_vcf
        self.allParameters ["Inputs"] = self.inputParameters

    def storeOutputParameters(self,args):
        """Updates output parameters to the map of parameters to be store in a JSON file

        args -- set of parsed arguments
        """
        self.outputParameters["alignment_out"] = args.alignment_out
        self.outputParameters["svcall_out"] = args.sv_call_out
        self.allParameters ["Outputs"] = self.outputParameters

    def storeWildcardParameters(self,args):
        """Updates wildcard parameters to the map of parameters to be store in a JSON file

        args -- set of parsed arguments
        """
        self.wildcardParameters["ONT_fastqs"] = args.ONT_fastqs
        self.allParameters ["Wildcards"] = self.wildcardParameters

    def storeMinimap2Parameters(self,args):
        """Updates minimap2 aligner parameters to the map of parameters to be store in a JSON file

        args -- set of parsed arguments
        """
        self.minimap2Parameters["minimap2_cores "] = args.minimap2_cores
        self.minimap2Parameters["minimap2_kmer_length"] = args.minimap2_kmer_length
        self.minimap2Parameters["minimap2_match_score"] = args.minimap2_match_score
        self.minimap2Parameters["minimap2_mismatch_score"] = args.minimap2_mismatch_score
        self.minimap2Parameters["minimap2_gap_open_score"] = args.minimap2_gap_open_score
        self.minimap2Parameters["GT_AG_find"] = args.GT_AG_find
        self.allParameters ["Minimap2"] = self.minimap2Parameters



    def storeNgmlrParameters(self,args):
        """Updates Ngmlr aligner parameters to the map of parameters to be store in a JSON file

        args -- set of parsed arguments
        """
        self.ngmlrParameters["ngmlr_cores"] = args.ngmlr_cores
        self.ngmlrParameters["ngmlr_min_ident"] = args.ngmlr_min_ident
        self.ngmlrParameters["ngmlr_match_score"] = args.ngmlr_match_score
        self.ngmlrParameters["ngmlr_mismatch_score"] = args.ngmlr_mismatch_score
        self.ngmlrParameters["ngmlr_gap_open_score"] = args.ngmlr_gap_open_score
        self.ngmlrParameters["ngmlr_gap_extend_min"] = args.ngmlr_gap_extend_min
        self.ngmlrParameters["ngmlr_kmer_length"] = args.ngmlr_kmer_length
        self.ngmlrParameters["ngmlr_kmer_skip"] = args.ngmlr_kmer_skip
        self.allParameters ["Ngmlr"] = self.ngmlrParameters



    def storeSnifflesParameters(self,args):
        """Updates Sniffles SV caller parameters to the map of parameters to be store in a JSON file

        args -- set of parsed arguments
        """
        self.snifflesParameters["sniffles_cores"] = args.sniffles_cores
        self.snifflesParameters["sniffles_min_sv_length"] = args.sniffles_min_sv_length
        self.snifflesParameters["sniffles_min_sv_length"] = args.sniffles_min_sv_length
        self.snifflesParameters["sniffles_min_read_length"] = args.sniffles_min_read_length
        self.snifflesParameters["sniffles_min_read_mapping_quality"] = args.sniffles_min_read_mapping_quality
        self.snifflesParameters["sniffles_min_read_support"] = args.sniffles_min_read_support
        self.snifflesParameters["sniffles_num_reads_report"] = args.sniffles_num_reads_report
        self.snifflesParameters["sniffles_max_num_splits"] = args.sniffles_max_num_splits
        self.snifflesParameters["sniffles_genotype"] = args.sniffles_genotype
        self.snifflesParameters["sniffles_cluster"] = args.sniffles_cluster
        self.snifflesParameters["sniffles_min_homo_af"] = args.sniffles_min_homo_af
        self.snifflesParameters["sniffles_min_het_af"] = args.sniffles_min_het_af
        self.allParameters ["Sniffles"] = self.snifflesParameters


    def storeSvimParameters(self,args):
        """Updates Svim SV caller parameters to the map of parameters to be store in a JSON file

        args -- set of parsed arguments
        """
        self.svimParameters["svim_min_sv_length"] = args.svim_min_sv_length
        self.svimParameters["svim_max_sv_length"] = args.svim_max_sv_length
        self.svimParameters["svim_min_read_length"] = args.svim_min_read_length
        self.svimParameters["svim_min_read_mapping_quality"] = args.svim_min_read_mapping_quality
        self.svimParameters["svim_gap_tolerance"] = args.svim_gap_tolerance
        self.svimParameters["svim_overlap_tolerance"] = args.svim_overlap_tolerance
        self.svimParameters["svim_partition_max_distance"] = args.svim_partition_max_distance
        self.svimParameters["svim_sv_max_distance"] = args.svim_sv_max_distance
        self.svimParameters["svim_min_geno_score"] = args.svim_min_geno_score
        self.svimParameters["svim_homozygous_thresh"] = args.svim_homozygous_thresh
        self.svimParameters["svim_heterozygous_thresh"] = args.svim_heterozygous_thresh
        self.svimParameters["svim_min_depth"] = args.svim_min_depth
        self.svimParameters["svim_min_score"] = args.svim_min_score
        self.allParameters ["Svim"] = self.svimParameters


#####

#1.Create object class Configuration File
configManager = CreateConfigurationFile()

#2.Create object for argument parsinng
parser = argparse.ArgumentParser(prog="create_configuration_file",
                description="Create a configuration json file for the LRSV pipeline."
                )     

#2.1 Updates arguments and parsing
configManager.register_parameter(parser)

args = parser.parse_args()

#2.2 Check Parameters
configManager.check_parameters(args)

#3. store arguments to super map structure
configManager.storeGeneralParameters(args)
configManager.storeInputParameters(args)
configManager.storeOutputParameters(args)
configManager.storeWildcardParameters(args)
configManager.storeMinimap2Parameters(args)
configManager.storeNgmlrParameters(args)
configManager.storeSnifflesParameters(args)
configManager.storeSvimParameters(args)

        


###
#4. Store JSON file
with open(args.configFile, 'w') as of:
    json.dump(configManager.allParameters, of, indent=2)
