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
    This function will find the files in a provided directory. Default is fastq.
    
    -directory. Input directory path
    -single. Boolean parameter. Specifies just one file to look at.
    -pattern. String parameter. It must include the file pattern with a wildcard
    """
    if os.path.isfile(folder):
        return folder
    else:
        print ()"We were unable to find the provided path for the input directory. Exiting now.")
        sys.exit(-1)
    
    files = []
    for file in glob.glob(os.path.join(folder, pattern)):
        files.append(file)

    if len(files) == 0:
        print("Could not find {} files in {}. Trying to use the stored wildcard parameter".format(pattern, folder))
        

    if single is True:
        if len(files) > 1:
            print("Warning: Multiple {} files found in {}. Returning the first one.".format(pattern, folder))
        return files[0]
    else:
        return files
