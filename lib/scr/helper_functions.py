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
