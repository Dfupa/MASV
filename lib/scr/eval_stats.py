#!/usr/bin/env python3
"""
Uses pybedtools to calculate the Recall and Precision of the benchmarking
"""
from __future__ import print_function

import io
import sys
import argparse
import os
import multiprocessing
import pandas as pd
import warnings
warnings.filterwarnings('ignore')

if __name__ == "__main__":

    ap = argparse.ArgumentParser(prog=os.path.basename(sys.argv[0]),
                                 usage=__doc__)
    ap.add_argument('--truth', required=True,
                    help='VCF/BED high confidence ("truth") dataset')
    ap.add_argument('--callset', required=True,
                    help='VCF/BED with the prediction callset')
    ap.add_argument('--svtype', default='<DEL>', type=str,
                    help='Define the feature type to filter. '
                         'Default is <DEL>')
    ap.add_argument('--processes', default=1, type=int,
                    help='Number of processes to use in parallel.')
    ap.add_argument('--sniffles', type=bool, default=False,
                   help='Parameter used to reformat the sniffles callset if provided')
    ap.add_argument('-v', '--verbose', action='store_true',
                    help='Verbose (goes to stderr)')
    args = ap.parse_args()

    hq = args.truth
    calls = args.callset
    feature = args.svtype
    sniffles = args.sniffles
    
    print("#############")
    print("#Step 1: cwd#")
    print("#############\n")
    owd = os.getcwd()
    nwd = os.path.dirname(calls)
    cwd = os.chdir(nwd)

    print("The original working dir was '"+str(owd)+"'. Now it is going to be changed to '"+str(nwd)+"'.\n")
    #Make sure the paths provided exist
    
    try:
        os.path.exists(hq)
    except IOerror:
        sys.exit(-1)
    
    try:
        os.path.exists(calls)
    except IOerror:
        sys.exit(-1)
            
    

    if args.processes > 3:
        print(
            "Only need 3 processes, resetting processes from {0} to 3".format(args.processes)
        )
        args.processes = 3
        
    def svim_reformat(callset, featuretype):
        """
        Reformat the svim vcf callset
        """
        print("###############################")
        print("#Step 2: Reformatting the file#")
        print("###############################\n")
        command_svim_1 = 'cat ' + os.path.join(callset) + ' | awk \'OFS="\t" {{ if($1 ~ /^#/) {{print $0}} else {{ if($5=="'+featuretype+'") {{ print $0 }} }}  }}\' > '+os.path.join(nwd, 'svim_reformated.vcf')
        print(os.popen(command_svim_1).read())
        print("Everything went smoothly. Bye!\n")      

    
    def sniffles_reformat(callset, featuretype):
        """
        Reformat the sniffles vcf callset
        """    
        if sniffles == True:
            print("###############################")
            print("#Step 2: Reformatting the file#")
            print("###############################\n")
            command_snif_1 = 'cat ' + os.path.join(callset) + ' | awk \'OFS="\t" {{ if($1 ~ /^#/) {{ print $0 }} }}\' > '+os.path.join(nwd, 'header.tmp.vcf')
            command_snif_2 = 'cat ' + os.path.join(callset) + ' | awk \'OFS="\t" {{ if($1 !~ /^#/) {{print $0}} }}\' | grep \"' + featuretype + '\" > '+os.path.join(nwd, 'subset.temp.vcf')
            command_snif_3 = 'cat '+os.path.join(nwd, 'header.tmp.vcf')+' '+os.path.join(nwd, 'subset.temp.vcf')+' > '+os.path.join(nwd, 'sniffles_reformated.vcf')
            print(os.popen(command_snif_1).read(), os.popen(command_snif_2).read(), os.popen(command_snif_3).read())
            print("Everything went smoothly. Bye!\n")
            print(os.popen('rm *.tmp.vcf').read())

    def recall_precision_stats(truth, callset):
        """
        Returns a table with the arithmetic 
        """
        print("############################################")
        print("#Step 3: Obtaining the precision and recall#")                                                                  
        print("############################################\n")
        #Total number of variants in the hq dataset
        command_hq_1 =  'cat ' + os.path.join(truth) + ' | awk \'OFS="\t" {{ if($1 !~ /^#/) {{print $0}} }}\' | wc -l'                           
        number_variants_hq = os.popen(command_hq_1).read()
                                                                           
        #Total number of variants called in the call dataset
        command_call_1 = 'cat ' + os.path.join(callset) + ' | awk \'OFS="\t" {{ if($1 !~ /^#/) {{print $0}} }}\' | wc -l'                                                                 
        number_variants_callset = os.popen(command_call_1).read()
                                                                           
        #Number of hits by the intersect -a truth -b callset
        command_out_1 = 'bedtools intersect -u -a '+os.path.join(truth)+' -b '+os.path.join(callset)+' | wc -l'
        truth_vs_call = os.popen(command_out_1).read()
                                                                           
        #Number of hits by the intersect -a callset -b truth
        command_out_1 = 'bedtools intersect -u -a '+os.path.join(callset)+' -b '+os.path.join(truth)+' | wc -l'
        call_vs_truth = os.popen(command_out_1).read()
                                                                           
        #Evaluation Metrics:                                                                   
        recall = int(truth_vs_call)/int(number_variants_hq) #TP/(TP+FN) considering that TP are any called sv that overlaps once with any sv of the truth dataset and being FN any called sv that do not overlap not even once
                                                                           
        precision = int(call_vs_truth)/int(number_variants_callset)
                                                                           
        df = pd.DataFrame([[recall],
                     [precision]], 
                      columns = ['Results in %'])
    
        df.insert(loc=0, column='Eval Metrics', value=['Sensitivity','Precision'])
        
        return df

    #Reformat the feature name for sniffles based on the default example
    if sniffles == True:
            feature = feature.replace("<", "")
            feature = feature.replace(">", "")
            reformat = sniffles_reformat(calls, feature)
            new_call = os.path.realpath('sniffles_reformated.vcf')                                                               
    #Reformat for svim                                                                       
    elif sniffles == False:
            reformat = svim_reformat(calls, feature)
            new_call = os.path.realpath('svim_reformated.vcf')
                                                                           
    results = recall_precision_stats(hq, new_call)
    
                                                                           
    file = open('eval_stats.txt', 'w')
    file.write("These are the results for the evaluation of "+os.path.basename(hq)+" truth dataset and "+os.path.basename(calls)+" callset dataset.\n")
    file.write("\n")
    file.write(results.to_string(index=False))
    file.write("\n\nAnd it is done! Bye.\n")
    file.close()
    print("Removing the temporary files, please wait ...\n")
    print(os.popen('rm *reformated.vcf').read())
    print("And it is done! Bye.\n")
