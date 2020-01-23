import random
import re
import pybedtools
import logging
import os
import sys
import tempfile

from argparse import ArgumentParser, RawDescriptionHelpFormatter
A
from pysam import VariantFile
from pyfaidx import Faidx

#########First create the SV_Info class
class SV_Info:

    def __init__(self, vcf_provided, caller):
        self.id = vcf_provided.id
        self.type = vcf_provided.info.get('SVTYPE', None)
        self.chr1 = vcf_provided.chrom
        self.pos1 = vcf_provided.pos
        self.pos2 = vcf_provided.stop
        if caller == "sniffles":
            self.length = abs(vcf_provided.info.get('SVLEN', None))
            self.read_support = vcf_provided.info.get('RE', None)
            self.chr2 = vcf_provided.info.get('CHR2', None)
            #self.pos2 = vcf_provided.info.get('END', None)
        else:
            self.length = (vcf_provided.info.get('SVLEN', None))
            self.read_support = vcf_provided.qual #quality score of svim
            self.chr2 = None
            self.precise = None
            #self.pos2 = vcf_provided.info.get('END')


        if len(vcf_provided.samples.keys()) != 1:
            raise RuntimeError(
                "Currently only a single sample per file is supported")

    
        self.vcf_record = vcf_provided

    def __repr__(self):
        return "{} {} {} {} {} {} {} {} {}".format(self.chr1,
                                                   self.pos1, self.chr2,
                                                   self.pos2, self.type,
                                                   self.length, self.read_support,
                                                   self.precise)
class VCFile:
    
    def __init__(self, vcf_path, genome, fix_param=False, caller="svim"):
        self._variants, self.header = self.read_vcf(vcf_path, genome, fixing=fix_param, caller="svim")

        
        if caller == "sniffles":

            SV_TYPES = ['DEL/INV', 'DUP', 'INV', 'TRA', 'BND', 'INVDUP', 'INS', 'DEL',
                        'INV/INVDUP']
            SV_TYPE_NAMES = {'DEL': 'Deletion',
                             'INS': 'Insertion',
                             'INV': 'Inversion',
                             'TRA': 'Translocation',
                             'BND': 'Translocation',
                             'DUP': 'Duplication',
                             'DEL/INV': 'Deletion/Inversion',
                             'DUP/INS': 'Duplication',
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
    #def faidx_genome(self, genome):
    #    fa = Faidx(genome)
    #    dictio = {}
    #    for item in fa.index:
    #        dictio[item] = int(fa.index[item].rlen)
    #    return dictio
    
    def read_vcf(self, vcf_path, genome, fixing=False, caller="svim"):
        variants = []
        contig_lengths = {}
        vcf = VariantFile(vcf_path, "r")
        fa = Faidx(genome)
        for item in fa.index:
            contig_lengths[item] = int(fa.index[item].rlen)

        for item in vcf.fetch():
            sv_info = SV_Info(item, caller)

            if fixing and sv_info.type in ['INS', 'DUP_INT', 'DUP/INS']:
                logging.debug("Changing {}/{} to {}/{}".format(sv_info.pos1,
                                                                 sv_info.pos2,
                                                                 (sv_info.pos1 - round(int(sv_info.length)/2)),
                                                                 (sv_info.pos1 + round(int(sv_info.length)/2))))

                for contig in contig_lengths:
                    if contig == sv_info.chr1:
                        max_len = contig_lengths[contig]
                        initial = int(sv_info.pos1 - round(int(sv_info.length)/2))
                        #final = int(sv_info.pos1 + round(int(sv_info.length)/2))
                        final = int(sv_info.pos1 + int(sv_info.length))
                        real_start = sv_info.vcf_record.start - round(int(sv_info.length)/2)
                        #real_stop = sv_info.vcf_record.start + round(int(sv_info.length)/2)
                        real_stop = sv_info.vcf_record.start + int(sv_info.length)
                        if initial < 0:
                            initial = 0
                        if final > max_len:
                            final = max_len
                        if real_start < 0:
                            real_start = 0
                        if real_stop > max_len:
                            real_stop = max_len


                        #sv_info.pos1 = initial
                        sv_info.pos2 = final
                        sv_info.vcf_record.stop = real_stop
                        #sv_info.vcf_record.start = real_start
                        #print("real start: "+str(sv_info.vcf_record.start)+" and real stop: "+str(real_stop)+" but the stop is "+str(sv_info.vcf_record.stop)+" and the length "+str(sv_info.vcf_record.rlen)+" and is a "+str(sv_info.vcf_record.info['SVTYPE']))


                #sv_info.pos1 = int(sv_info.pos1 - round(int(sv_info.length)/2))
                #sv_info.pos2 = int(sv_info.pos1 + round(int(sv_info.length)/2))
                #sv_info.vcf_record.stop = sv_info.vcf_record.start + round(int(sv_info.length)/2)
                #sv_info.vcf_record.start = sv_info.vcf_record.start - round(int(sv_info.length)/2)
                #print(sv_info.vcf_record)

            variants.append(sv_info)
        return variants, vcf.header            
         

    #def check(self):
    #    for variant in self._variants:
    #        print(variant.pos2)
    #        if variant.pos1 >= variant.pos2 and variant.type != 'BND':
    #            logging.warning("POS1 >= POS2 for:")
    #            logging.warning("{}".format(variant))

                
    def write_vcf(self, vcfpath):
        vcf = VariantFile(vcfpath, 'w', header=self.header)
        for variant in self._variants:
            vcf.write(variant.vcf_record)
        vcf.close()
#########Now let's set the main functions of the script####



def main(argv=sys.argv[1:]):
    """
    Basic command line interface to cuecat.
    :param argv: Command line arguments
    :type argv: list
    :return: None
    :rtype: NoneType
    """
    args = parse_args(argv=argv)

    vcf_file = args.VCF
    if not os.path.exists(vcf_file):
        raise OSError("Could not find {}.".format(vcf_file))


    processed_file = VCFile(vcf_file, genome=args.GENOME, fix_param=True, caller=args.caller)

    # sniffles_file.fix_overlaping_ins_dup()

    #check_file(processed_file)

    processed_file.write_vcf(args.OUTPUT)
        
                

def check_file(vcf_file):
    vcf_file.check()
    
def parse_args(argv):
    usage = "Script to fix the insertion calls length in order to be used by bedtools intersect"
    parser = ArgumentParser(description=usage,
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
                        default="/dev/stdout",
                        help='Output VCF file')
    parser.add_argument("--caller",
                        dest="caller",
                        type=str,
                        help="Determine the type of caller in order to process the vcf file. Default \"svim\".")
    parser.add_argument("--check",
                        dest="CHECK",
                        action='store_true',
                        help="Check file for integrity.")

    args = parser.parse_args(argv)
    return args


if __name__ == '__main__':
    main()
