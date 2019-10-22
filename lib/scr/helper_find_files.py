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

def find_files(directory, single=False, pattern="*.fastq"):
    """
    This function will find the files in a provided directory. Default pattern is *.fastq.
    
    -directory. Input directory path
    -single. Boolean parameter. Specifies just one file to look at.
    -pattern. String parameter. It must include the file pattern with a wildcard
    """
    if os.path.isfile(directory):
        return directory
    else:
        print ("We were unable to find the provided path for the input directory. Exiting now.")
        sys.exit(-1)
    
    files = []
    for file in glob.glob(os.path.join(directory, pattern)):
        files.append(file)

    if len(files) == 0:
        print("We were not able to find {} files in the provided directory {}. Trying to use the stored wildcard parameter".format(pattern, directory))
        

    if single is True:
        if len(files) > 1:
            print("WARNING: Multiple {} files found in the {} directory. Returning the first one.".format(pattern, directory))
        return files[0]
    else:
        return files
