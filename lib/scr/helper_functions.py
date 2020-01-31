#!/usr/bin/env python3
import os
import glob
import sys


#Author: Diego Fuentes
#Contact email: diegofupa@gmail.com
#Date:2019-10-10

#######################
####HELPER FUNCTION####
#######################

#This function has been deprecated as right now. Leave it here for old documentation purposes

#def find_files(directory, single=False, pattern="*.fastq"):
    #"""
    # This function will find the files in a provided directory. Default pattern is *.fastq.
    
    #-directory. Input directory path
    #-single. Boolean parameter. Specifies just one file to look at.
    #-pattern. String parameter. It must include the file pattern with a wildcard
    #"""
    #if os.path.isfile(directory):
    #    return directory
    #else:
    #    print ("We were unable to find the provided path for the input directory. Exiting now.")
    #    sys.exit(-1)
    
    #files = []
    #for file in glob.glob(os.path.join(directory, pattern)):
    #    files.append(file)

    #if len(files) == 0:
    #    print("We were not able to find {} files in the provided directory {}. Trying to use the stored wildcard parameter".format(pattern, directory))
        

    #if single is True:
    #    if len(files) > 1:
    #        print("WARNING: Multiple {} files found in the {} directory. Returning the first one.".format(pattern, directory))
    #    return files[0]
    #else:
    #    return files

def which_tech(parameter="nanopore"):
    """
    This function selects the mapping parameters for either aligner based on the sequencing technology
    
    -parameter. String parameter. Default is "nanopore".
    """

    if parameter == "nanopore":
        minimap2_parameter = "map-ont"
        ngmlr_parameter = "ont"

    elif parameter == "pacbio":
        minimap2_parameter = "map-pb"
        ngmlr_parameter = "pacbio"

    else:
        print("Warning: the selected sequencing technology is not valid. ONT has been selected by default")
        minimap2_parameter = "map-ont"
        ngmlr_parameter = "ont"

    return minimap2_parameter, ngmlr_parameter
