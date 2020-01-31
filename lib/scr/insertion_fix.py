#!/usr/bin/env python3

#Author: Diego Fuentes
#Contact email: diegofupa@gmail.com

"""
Description:
Script to fix the insertion calls end position based on sv length in order to be used by bedtools intersect. It takes into account the chromsize (contigsize) of the 
provided genome to make sure it does not exceed the contig max length.

Author: Diego Fuentes
Contact email: diegofupa@gmail.com
"""

import random
import re
import pybedtools
import logging
import os
import sys
import tempfile
from argparse import ArgumentParser, RawDescriptionHelpFormatter
from pysam import VariantFile
from pyfaidx import Faidx

######################
# Define the Objects #
######################

##First create the SV_Info object##
#SV_Info class is going to be used to process the VariantFile objects selecting the features necessary for the fix.

class SV_Info:

    def __init__(self, vcf_provided, caller):
        self.id = vcf_provided.id
        self.type = vcf_provided.info.get('SVTYPE', None)
        self.chr1 = vcf_provided.chrom
        self.pos1 = vcf_provided.pos
        self.pos2 = vcf_provided.stop
        if caller == "sniffles":
            self.length = abs(vcf_provided.info.get('SVLEN', None)) #abs because deletions are reported as negative length
            #self.read_support = vcf_provided.info.get('RE', None)
        else:
            self.length = vcf_provided.info.get('SVLEN', None)
            #self.read_support = vcf_provided.qual #Quality score of svim which is based on Read support

        self.vcf_record = vcf_provided

        if len(vcf_provided.samples.keys()) != 1:
            raise RuntimeError(
                "Currently only a single sample per file is supported")

    def __repr__(self):
        return "{} {} {} {} {} {} ".format(self.id, self.chr1,
                                                   self.pos1, self.pos2,
                                                   self.type, self.length)
##Then create the VCFile object##
#The VCF object will have two static methods to process the SVTYPE correctly, plus a method to fix the insertion coordinates where POS2 - POS1 <= 1

class VCFile:
    
    def __init__(self, vcf_path, genome, fix_param=False, caller="svim"):
        self._variants, self.header = self.read_vcf(vcf_path, genome, fixing=fix_param, caller="svim")

        
        if caller == "sniffles":

            SV_TYPES = ['DEL/INV', 'DUP', 'INV', 'TRA', 'BND', 'INVDUP', 'INS', 'DEL',
                        'DUP/INS', 'INV/INVDUP']
            SV_TYPE_NAMES = {'DEL': 'Deletion',
                             'INS': 'Insertion',
                             'INV': 'Inversion',
                             'TRA': 'Translocation',
                             'BND': 'Translocation',
                             'DUP': 'Duplication',
                             'DEL/INV': 'Deletion/Inversion',
                             'DUP/INS': 'Interspread Duplication',
                             'INVDUP': 'Inverted Duplication'}
            SV_BASE_TYPE =  {'DEL': 'DEL',
                             'INS': 'INS',
                             'INV': 'INV',
                             'TRA': 'TRA',
                             'BND': 'BND',
                             'DUP': 'DUP',
                             'DUP/INS': 'DUP',
                             'DEL/INV': 'INV',
                             'INVDUP': 'INV'}
        if caller == "svim":

            SV_TYPES = ['DEL', 'DUP:TANDEM', 'INV', 'DUP_INT', 'BND', 'INS']
            SV_TYPE_NAMES = {'DEL': 'Deletion',
                             'INS': 'Insertion',
                             'INV': 'Inversion',
                             'BND': 'Translocation',
                             'DUP:TANDEM': 'Tandem Duplication',
                             'DUP:INT': 'Interspread Duplication'}
            SV_BASE_TYPE =  {'DEL': 'DEL',
                             'INS': 'INS',
                             'INV': 'INV',
                             'BND': 'BND',
                             'DUP:TANDEM': 'DUP:TANDEM',
                             'DUP_INT': 'INS'}

    #Define two static methods

    @staticmethod
    def get_display_name(id):
        if id in VCFile.SV_TYPE_NAMES:
            return VCFile.SV_TYPE_NAMES[id]
        else:
            return id

    @staticmethod
    def get_base_type(sv_type):
        if sv_type in VCFile.SV_BASE_TYPE:
            return VCFile.SV_BASE_TYPE[sv_type]
        else:
            return sv_type

    
    def read_vcf(self, vcf_path, genome, fixing=False, caller="svim"):
        """
        This method changes INS, DUP_INT and DUP/INS stop ('END') coordinates if POS2 - POS1 <= 1.
        It takes into account the chromsize (contigsize) of the provided genome to make sure it does not
        exceed NEW_STOP > CONTIG_MAX_LENGTH.

        Params:
        - vcf_path: Str. VCF path provided to as an argument.
        - genome: Str. Genome path provided as an argument.
        - fixing: Boolean. When set to True takes into account the abovementioned SVs for fixing.
        - caller: Str. Selected caller provided as an argument. It is required by the SV_Info object.  
        """
        variants = []
        contig_lengths = {}
        vcf = VariantFile(vcf_path, "r")
        fa = Faidx(genome)
        for item in fa.index:
            contig_lengths[item] = int(fa.index[item].rlen)

        for item in vcf.fetch():
            sv_info = SV_Info(item, caller)

            if fixing and sv_info.type in ['INS', 'DUP_INT', 'DUP/INS']:
                if int(sv_info.pos2 - sv_info.pos1) <= 1:
                    logging.debug("Changing {}/{} to {}/{}".format(sv_info.pos1,
                                                                   sv_info.pos2,
                                                                   (sv_info.pos1 - round(int(sv_info.length)/2)),
                                                                   (sv_info.pos1 + round(int(sv_info.length)/2))))

                    for contig in contig_lengths:
                        if contig == sv_info.chr1:
                            max_len = contig_lengths[contig]
                            final = int(sv_info.pos1 + int(sv_info.length))
                            real_stop = sv_info.vcf_record.start + int(sv_info.length)

                            if final > max_len:
                                final = max_len
                            if real_stop > max_len:
                                real_stop = max_len


                            sv_info.pos2 = final
                            sv_info.vcf_record.stop = real_stop
                            #print("real start: "+str(sv_info.vcf_record.start)+" and real stop: "+str(real_stop)+" but the stop is "+str(sv_info.vcf_record.stop)+" and the length "+str(sv_info.vcf_record.rlen)+" and is a "+str(sv_info.vcf_record.info['SVTYPE']))
                else:
                    pass

            variants.append(sv_info)
        return variants, vcf.header            
                
    def write_vcf(self, vcfpath):
        vcf = VariantFile(vcfpath, 'w', header=self.header)
        for variant in self._variants:
            vcf.write(variant.vcf_record)
        vcf.close()


#######################
# Set the main script #
#######################

def main(argv=sys.argv[1:]):
    """
    Main script. Reformats a VCF file and fixes the INS, DUP_INT and DUP/INS stop ('END') coordinates.

    Params:
    - argv. Provided command line arguments
    """
    args = parse_args(argv=argv)

    #First two quick path.exists checks
    vcf_file = args.VCF
    if not os.path.exists(vcf_file):
        raise OSError("Could not find {}.".format(vcf_file))

    genome_file = args.GENOME
    if not os.path.exists(genome_file):
        raise OSError("Could not find {}.".format(genome_file))


    processed_file = VCFile(vcf_file, genome=genome_file, fix_param=True, caller=args.caller)

    processed_file.write_vcf(args.OUTPUT)
        
                    
def parse_args(argv):
    #usage = "Script to fix the insertion calls length in order to be used by bedtools intersect"
    parser = ArgumentParser(description=__doc__,
                            formatter_class=RawDescriptionHelpFormatter)
    parser.add_argument('-v', '--vcf',
                        dest='VCF',
                        type=str,
                        required=True,
                        help='Input VCF file')
    parser.add_argument('-g', '--genome',
                        dest='GENOME',
                        type=str,
                        required=True,
                        help='Input genome (.fa, .fasta) file')
    parser.add_argument('-o', '--output',
                        dest='OUTPUT',
                        type=str,
                        default="stdout",
                        help='Output VCF file')
    parser.add_argument("--caller",
                        dest="caller",
                        type=str,
                        default="svim",
                        help="Determine the type of caller in order to process the vcf file. Default \"svim\".")

    args = parser.parse_args(argv)
    return args


if __name__ == '__main__':
    main()
