#!/usr/bin/env python3

#Author: Diego Fuentes
#Contact email: diegofupa@gmail.com
#Date:2019-12-01

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
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import numpy as np
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
    ap.add_argument('--plot', type=bool, default=False,
                   help='Parameter used to produce a plot of the eval metrics for different svim score filtering')
    ap.add_argument('-v', '--verbose', action='store_true',
                    help='Verbose (goes to stderr)')
    args = ap.parse_args()

    hq_path = args.truth
    calls_path = args.callset
    feature = args.svtype
    sniffles = args.sniffles
    plot = args.plot

    #Make sure the paths provided exist
    
    try:
        os.path.exists(hq_path)
    except IOerror:
        sys.exit(-1)
    
    try:
        os.path.exists(calls_path)
    except IOerror:
        sys.exit(-1)

    #We set the truth dataset path
    hq = os.path.abspath(hq_path)
    

    #The results will be stored where the predicted calls are originally stored
    print("#############")
    print("#Step 1: cwd#")
    print("#############\n")
    owd = os.getcwd()
    nwd = os.path.dirname(os.path.abspath(calls_path))
    print(nwd)
    if owd == nwd:
        pass
    else:
        cwd = os.chdir(nwd)
        nwd = os.getcwd()

    print("The original working dir was '"+str(owd)+"'. Results are going to be stored in '"+str(nwd)+"'.\n")

    #Provide only the file name for the callset
    calls = os.path.basename(os.path.abspath(calls_path))

            
    

    if args.processes > 3:
        print(
            "Only need 3 processes, resetting processes from {0} to 3".format(args.processes)
        )
        args.processes = 3
        
    def svim_reformat(callset, featuretype, tempfile = 'svim_reformated.vcf'):
        """
        Reformat the svim vcf callset
        """

        #command_svim_1 = 'cat ' + os.path.realpath(callset)+ ' | awk \'OFS="\\t" {{ if($1 ~ /^#/) {{print $0}}  }}\' > '+os.path.join(nwd, 'header.vcf')
        #command_svim_2 = 'cat ' + os.path.realpath(callset)+ ' | awk \'OFS="\\t" {{ if($5=="'+featuretype+'") {{ print a $1, $2, $3, $4, $5, $6, $7, $8, $9, $10 }} }}\' > '+ os.path.join(nwd, 'variant.vcf')
        #command_svim_3 = 'cat ' + os.path.join(nwd, 'header.vcf')+ ' ' + os.path.join(nwd, 'variant.vcf')+ ' > '+os.path.join(nwd, 'svim_reformated.vcf') 
        #print(os.popen(command_svim_1).read(), os.popen(command_svim_2).read(), os.popen(command_svim_3).read())
        #print("Everything went smoothly. Bye!\n")
        #print(os.popen('rm *.tmp.vcf').read())      
        command_svim_1 = 'cat ' + os.path.join(callset) + ' | awk \'OFS="\\t" {{ if($1 ~ /^#/) {{print $0}} else {{ if($5=="'+featuretype+'") {{ print $0 }} }}  }}\' > '+os.path.join(nwd, tempfile)
        print(os.popen(command_svim_1).read())
 
    
    def sniffles_reformat(callset, featuretype):
        """
        Reformat the sniffles vcf callset
        """    
        if sniffles == True:

            command_snif_1 = 'cat ' + os.path.abspath(callset) + ' | awk \'OFS="\\t" {{ if($1 ~ /^#/) {{ print $0 }} }}\' > '+os.path.join(nwd, 'header.tmp.vcf')
            command_snif_2 = 'cat ' + os.path.abspath(callset) + ' | awk \'OFS="\\t" {{ if($1 !~ /^#/) {{print $0}} }}\' | grep \"' + featuretype + '\" > '+os.path.join(nwd, 'subset.tmp.vcf')
            command_snif_3 = 'cat '+os.path.join(nwd, 'header.tmp.vcf')+' '+os.path.join(nwd, 'subset.temp.vcf')+' > '+os.path.join(nwd, 'sniffles_reformated.vcf')
            print(os.popen(command_snif_1).read(), os.popen(command_snif_2).read(), os.popen(command_snif_3).read())
            print(os.popen('rm *.tmp.vcf').read())

    def recall_precision_stats(truth, callset):
        """
        Returns a table with the arithmetic 
        """

        #Total number of variants in the hq dataset
        command_hq_1 =  'cat ' + os.path.abspath(truth) + ' | awk \'OFS="\\t" {{ if($1 !~ /^#/) {{print $0}} }}\' | wc -l'                           
        number_variants_hq = os.popen(command_hq_1).read()
                                                                  
        #Total number of variants called in the call dataset
        command_call_1 = 'cat ' + os.path.abspath(callset) + ' | awk \'OFS="\\t" {{ if($1 !~ /^#/) {{print $0}} }}\' | wc -l'                                                                 
        number_variants_callset = os.popen(command_call_1).read()
                                                                          
        #Number of hits by the intersect -a truth -b callset
        command_out_1 = 'bedtools intersect -u -a '+os.path.abspath(truth)+' -b '+os.path.abspath(callset)+' | wc -l'
        truth_vs_call = os.popen(command_out_1).read()
                                                                          
        #Number of hits by the intersect -a callset -b truth
        command_out_2 = 'bedtools intersect -u -a '+os.path.abspath(callset)+' -b '+os.path.abspath(truth)+' | wc -l'
        call_vs_truth = os.popen(command_out_2).read()
                                                                            
        #Evaluation Metrics:                                                                   
        recall = int(truth_vs_call)/int(number_variants_hq) #TP/(TP+FN) considering that TP are any called sv that overlaps once with any sv of the truth dataset and being FN any of the truth dataset SVs left to be called
                                                                           
        precision = int(call_vs_truth)/int(number_variants_callset)#TP/(TP+FP) considering that TP are any called sv that overlaps once with any sv of the truth dataste and being FP any of the other called svs that did not overlap not even once
        
        f1 = 2*((recall*precision)/(recall+precision))#harmonic mean of precision and recall
                                                                           
        df = pd.DataFrame([[recall],
                     [precision], [f1]], 
                      columns = ['Results in %'])
    
        df.insert(loc=0, column='Eval Metrics', value=['Sensitivity','Precision', 'F1'])
        
        return df
    
    def plot_svim_filtering(truth, callset, featuretype, iterations=21):
                                                                           
        
        performance_matrix = np.zeros((iterations,3))
        
        if sniffles == False and plot == True:
            
            print("############################################")
            print("#Optional Step: Creating the plot#")                                                                  
            print("############################################\n")
        
            for i in range(0, iterations, 1):
                subcommand = 'grep -v \"hom_ref\" '+os.path.realpath(callset)+' | awk  \'OFS="\\t" {{ if($1 ~ /^#/) {{ print $0 }} else {{ if($6>='+str(i)+') {{ print $0 }} }} }}\' > '+os.path.join(nwd, 'filter.temp.vcf')
                print(os.popen(subcommand).read())                                                                                                                                     
                new_calls = os.path.abspath('filter.temp.vcf')                                                           
                reformated = svim_reformat(new_calls, featuretype, tempfile = 'temp.vcf')
                reformated_call = os.path.abspath('temp.vcf')
                results = recall_precision_stats(truth, reformated_call)
                sensitivity = results.iloc[0]['Results in %']
                precision = results.iloc[1]['Results in %']    
                f1 = 2*((sensitivity*precision)/(sensitivity+precision))
                #Store it in the performance matrix
                performance_matrix[i, 0] = sensitivity
                performance_matrix[i, 1] = precision
                performance_matrix[i, 2] = f1
                #Remove temp files
                print(os.popen('rm filter.temp.vcf').read())
                print(os.popen('rm temp.vcf').read())
                print('*')
                
            
            #Let's start to plot
            print('\nNow lets plot it')
            X = performance_matrix[:, 0]#Sensitivity values
            Y = performance_matrix[:, 1]#Precision values
            F = performance_matrix[:, 2]#F1 values
            
            plt.figure(figsize=(20,14))
            #Subplot 1
            plt.subplot(2,2,1)
            plt.title("Positive feedback mRNA")
            line1, = plt.plot(Y, X, 'o-', color="b")
            plt.title("Evolution of the eval metrics based on SVIM score filtering (from >=0 to >=40, ascendant)")
            plt.xlabel('Precision')
            plt.xticks(np.array(Y))
            plt.ylabel('Sensitivity')
            plt.yticks(np.array(X))
            legend_handles = [ mlines.Line2D([], [], color='b', marker='o', \
                            markersize=10, label='Sensitivity-Precision trade-off')]
            plt.legend(handles=legend_handles, loc = 'best')

            #Subplot 2
            plt.subplot(2,2,2)
            plt.title("Precision based on the Svim Q score filtering")
            line1, = plt.plot(np.arange(0, iterations), Y, 'o-', color="r")
            plt.xlabel('Svim Q score selected')
            plt.xticks(np.arange(0, iterations))
            plt.ylabel('Precision')
            plt.yticks(np.array(Y))
            legend_handles = [ mlines.Line2D([], [], color='r', marker='o', \
                            markersize=10, label='Precision')]
            plt.legend(handles=legend_handles, loc = 'best')
            
            #Subplot 3
            plt.subplot(2,2,3)
            plt.title("Sensitivity based on the Svim Q score filtering")
            line1, = plt.plot(np.arange(0, iterations), X, 'o-', color="g")
            plt.xlabel('Svim Q score selected')
            plt.xticks(np.arange(0, iterations))
            plt.ylabel('Sensitivity')
            plt.yticks(np.array(X))
            legend_handles = [ mlines.Line2D([], [], color='g',marker='o', \
                            markersize=10, label='Sensitivity')]
            plt.legend(handles=legend_handles, loc = 'best')

            #Subplot 4
            plt.subplot(2,2,4)
            line1, = plt.plot(np.arange(0, iterations), F, 'o-', color="purple")
            plt.xlabel('Svim Q score selected')
            plt.xticks(np.arange(0, iterations))
            plt.ylabel('F1 score')
            plt.yticks(np.array(F))
            legend_handles = [ mlines.Line2D([], [], color='purple',marker='o', \
                            markersize=10, label='F1 score')]
            plt.legend(handles=legend_handles, loc = 'best')
                              
            #Save the plot
            plt.savefig(os.path.join(nwd, 'svim_filtering.png'))

            
            
        else:
            pass

    #Reformat the feature name for sniffles based on the default example
    if sniffles == True:
            feature = feature.replace("<", "")
            feature = feature.replace(">", "")
            
            print("###############################")
            print("#Step 2: Reformatting the file#")
            print("###############################\n")
            
            reformat = sniffles_reformat(calls, feature)
            new_call = os.path.realpath('sniffles_reformated.vcf')
            print("Everything went smoothly.\n")                                                               
    #Reformat for svim                                                                       
    elif sniffles == False:
        
            print("###############################")
            print("#Step 2: Reformatting the file#")
            print("###############################\n")
            
            reformat = svim_reformat(calls, feature)
            new_call = os.path.realpath('svim_reformated.vcf')
            print("Everything went smoothly.\n")
            
            
    print("############################################")
    print("#Step 3: Obtaining the precision and recall#")                                                                  
    print("############################################\n")                                                                       
    results = recall_precision_stats(hq, new_call)
    
    #This is an optional step, if the conditions are met, will produce a plot
    check_plot =plot_svim_filtering(hq, new_call, feature)    
          
    file = open('eval_stats.txt', 'w')
    file.write("These are the results for the evaluation of "+os.path.basename(hq)+" truth dataset and "+os.path.basename(calls)+" callset dataset.\n")
    file.write("\n")
    file.write(results.to_string(index=False))
    file.write("\n\nAnd it is done! Bye.\n")
    file.close()
    print("Removing the temporary files, please wait ...\n")
    print(os.popen('rm svim_reformated.vcf').read())
    print("And it is done! Bye.\n")
