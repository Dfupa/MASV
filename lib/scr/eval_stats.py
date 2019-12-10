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
from datetime import datetime
import warnings
warnings.filterwarnings('ignore')

date1 = str(datetime.now())
tmp = str.replace(date1," ",".") 
tmp2 = str.replace(tmp,":","")
date = str.replace(tmp2,"-","")

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
                   help='Parameter used to reformat the sniffles callset if provided. Default is \'False\'.')
    ap.add_argument('--plot', type=bool, default=False,
                   help='Parameter used to produce a plot of the eval metrics for different svim score filtering. Default is \'False\'.')
    ap.add_argument('--iterator', type=int, default=20,
                   help='Parameter used in the plotting step. It is the max svim filtering score range indicated from 0 to \'--iterator\'. Default is \'20\'.')
    ap.add_argument('-v', '--verbose', action='store_true',
                    help='Verbose (goes to stderr)')
    args = ap.parse_args()

    hq_path = args.truth
    calls_path = args.callset
    feature = args.svtype
    sniffles = args.sniffles
    plot = args.plot
    iterator = args.iterator

    #Make sure the paths provided exist
    
    try:
        os.path.exists(hq_path)
    except IOerror:
        sys.exit(-1)
    
    try:
        os.path.exists(calls_path)
    except IOerror:
        sys.exit(-1)

    if iterator >= 10:
        pass
    else:
        print("plot.py: error: argument --iterator: The iterator must surprass 10 to have an acceptable range")
        sys.exit(-1)

    #We set the truth dataset path
    hq = os.path.abspath(hq_path)
    

    #The results will be stored where the predicted calls are originally stored
    print("#############")
    print("#Step 1: cwd#")
    print("#############\n")
    owd = os.getcwd()
    nwd = os.path.dirname(os.path.abspath(calls_path))
    
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

     
        command_svim_1 = 'cat ' + os.path.join(callset) + ' | awk \'OFS="\\t" {{ if($1 ~ /^#/) {{print $0}} else {{ if($5=="'+featuretype+'") {{ print $0 }} }}  }}\' > '+os.path.join(nwd, tempfile)
        print(os.popen(command_svim_1).read())
 
    
    def sniffles_reformat(callset, featuretype, tempfile = 'sniffles_reformated.vcf'):
        """
        Reformat the sniffles vcf callset
        """   
        if sniffles == True:
            print(os.path.join(callset))
            command_snif_1 = 'cat ' + os.path.join(callset) + ' | awk \'OFS="\\t" {{ if($1 ~ /^#/) {{ print $0 }} }}\' > '+os.path.join('header.tmp.vcf')
            command_snif_2 = 'cat ' + os.path.join(callset) + ' | awk \'OFS="\\t" {{ if($1 !~ /^#/) {{print $0}} }}\' | grep \"' + featuretype + '\" > '+os.path.join(nwd, 'subset.tmp.vcf')
            command_snif_3 = 'cat '+os.path.join(nwd, 'header.tmp.vcf')+' '+os.path.join(nwd, 'subset.tmp.vcf')+' > '+os.path.join(nwd, tempfile)
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
    
    def plot_filtering(truth, callset, featuretype, iterations=20):
                                                                           
        
        performance_matrix = np.zeros((iterations,3))
        
        if sniffles == False and plot == True:
            dataset = "Svim"
            print("############################################")
            print("#Optional Step: Creating the plot#")                                                                  
            print("############################################\n")
        
            for i in range(0, iterations+1, 1):
                i = i - 1 #To avoid issues with Python's indexing
                subcommand = 'grep -v \"hom_ref\" '+os.path.realpath(callset)+' | awk  \'OFS="\\t" {{ if($1 ~ /^#/) {{ print $0 }} else {{ if($6>='+str(i)+') {{ print $0 }} }} }}\' > '+os.path.join(nwd, 'filter.temp.vcf')
                print(os.popen(subcommand).read())                                                                                                                                     
                new_calls = os.path.abspath('filter.temp.vcf')                                                           
                reformated = svim_reformat(new_calls, featuretype, tempfile = 'temp.vcf')
                reformated_call = os.path.abspath('temp.vcf')
                results = recall_precision_stats(truth, reformated_call)
                sensitivity = results.iloc[0]['Results in %']
                precision = results.iloc[1]['Results in %']    
                f1 = results.iloc[2]['Results in %']
                #Store it in the performance matrix
                performance_matrix[i, 0] = sensitivity
                performance_matrix[i, 1] = precision
                performance_matrix[i, 2] = f1
                #Remove temp files
                print(os.popen('rm filter.temp.vcf').read(), os.popen('rm temp.vcf').read(), 'Iteration: '+str(i+1)+'.')
                #print(os.popen('rm temp.vcf').read())
                
                
            
            #Let's start to plot
            print('\nNow lets plot it\n')
            X = performance_matrix[:, 0]#Sensitivity values
            Y = performance_matrix[:, 1]#Precision values
            F = performance_matrix[:, 2]#F1 values


            #Mean and standard deviation stored in a list. It will be added to the eval_stats.txt if plot == True
            result_list = []
            mean_recall = np.mean(X)
            result_list.append(mean_recall)
            mean_precision = np.mean(Y)
            result_list.append(mean_precision)
            mean_f1 = np.mean(F)
            result_list.append(mean_f1)
            std_recall = np.std(X)
            result_list.append(std_recall)
            std_precision = np.std(Y)
            result_list.append(std_precision)
            std_f1 = np.std(F)
            result_list.append(std_f1)
            

            #########PLOTS#########
            plt.figure(figsize=(20,14))
            plt.title("Results for "+str(dataset)+" based on a range of quality score filtering ("+str(iterations)+" iterations).")
            plt.grid()
            #Subplot 1
            plt.subplot(2,2,1)
            plt.title("Positive feedback mRNA")
            line1, = plt.plot(Y, X, 'o-', color="b")
            plt.xlabel('Precision')
            plt.xticks(np.arange(0, 0.4, 0.05))
            plt.ylabel('Sensitivity')
            plt.yticks(np.arange(0.6, 1.05, 0.05))
            legend_handles = [ mlines.Line2D([], [], color='b', marker='o', \
                            markersize=10, label='Sensitivity-Precision trade-off')]
            plt.legend(handles=legend_handles, loc = 'best')

            #Subplot 2
            plt.subplot(2,2,2)
            plt.title("Precision based on the Svim Q score filtering")
            line1, = plt.plot(np.arange(0, iterations), Y, 'o-', color="r")
            plt.fill_between(np.arange(0, iterations), Y - std_precision / np.sqrt(10), \
                         Y + std_precision / np.sqrt(10), alpha=0.1, color="green")
            plt.xlabel('Svim Q score selected')
            plt.xticks(np.arange(0, iterations))
            plt.ylabel('Precision')
            plt.yticks(np.arange(0, 0.45, 0.05))
            legend_handles = [ mlines.Line2D([], [], color='r', marker='o', \
                            markersize=10, label='Precision')]
            plt.legend(handles=legend_handles, loc = 'best')
            
            #Subplot 3
            plt.subplot(2,2,3)
            plt.title("Sensitivity based on the Svim Q score filtering")
            line1, = plt.plot(np.arange(0, iterations), X, 'o-', color="g")
            plt.fill_between(np.arange(0, iterations), X - std_recall / np.sqrt(10), \
                         X + std_recall / np.sqrt(10), alpha=0.1, color="green")
            plt.xlabel('Svim Q score selected')
            plt.xticks(np.arange(0, iterations))
            plt.ylabel('Sensitivity')
            plt.yticks(np.arange(0.6, 1.05, 0.05))
            legend_handles = [ mlines.Line2D([], [], color='g',marker='o', \
                            markersize=10, label='Sensitivity')]
            plt.legend(handles=legend_handles, loc = 'best')

            #Subplot 4
            plt.subplot(2,2,4)
            line1, = plt.plot(np.arange(0, iterations), F, 'o-', color="purple")
            plt.fill_between(np.arange(0, iterations), F - std_f1 / np.sqrt(10), \
                         F + std_f1 / np.sqrt(10), alpha=0.1, color="green")
            plt.xlabel('Svim Q score selected')
            plt.xticks(np.arange(0, iterations))
            plt.ylabel('F1 score')
            plt.yticks(np.arange(0, 0.65, 0.05))
            legend_handles = [ mlines.Line2D([], [], color='purple',marker='o', \
                            markersize=10, label='F1 score')]
            plt.legend(handles=legend_handles, loc = 'best')
                              
            #Save the plot
            plt.savefig(os.path.join(nwd, (str(dataset)+'_filtering'+str(featuretype)+'.png')))

            return result_list

        else:
            pass

    #Reformat the feature name for sniffles based on the default example
    if sniffles == True:
            dataset = "sniffles"
            if feature == '<DUP:TANDEM>':
                feature = '<DUP>'
            else:
                pass  
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
            dataset = "svim"
            print("###############################")
            print("#Step 2: Reformatting the file#")
            print("###############################\n")
            
            reformat = svim_reformat(calls, feature)
            new_call = os.path.realpath('svim_reformated.vcf')
            print("Everything went smoothly.\n")
            
            
    print("############################################")
    print("#Step 3: Obtaining the precision and recall#")                                                                  
    print("############################################\n")                                                                       
    try:
        results = recall_precision_stats(hq, new_call)
    except ValueError:
        sys.exit(-1)
    print("Everything went smoothly.\n")
    #This is an optional step, if the conditions are met, will produce a plot
    try:
        check_plot =plot_filtering(hq, new_call, feature, iterator)    
    except ValueError:
        sys.exit(-1)      
    file = open(str(date)+'_eval_stats_'+str(dataset)+'_'+str(feature)+'.txt', 'w')
    file.write("These are the results for the evaluation of "+os.path.basename(hq)+" truth dataset and "+os.path.basename(calls)+" callset dataset for the specific feature "+str(feature)+".\n")
    file.write("\n")
    if plot == True:
        file.write("Note that the added columns 'Mean values' and 'Std values' stand for the mean and standard deviation of the iterative loop using a range of 0 to "+str(iterator)+" Q scores.\n")
        results.insert(loc=2, column='Mean values', value=[check_plot[0],check_plot[1], check_plot[2]])
        results.insert(loc=3, column='Std values', value=[check_plot[3],check_plot[4], check_plot[5]])
        file.write("\n")
        file.write(results.to_string(index=False))
    else:
        file.write(results.to_string(index=False))
    file.write("\n\nAnd it is done! Bye.\n")
    file.close()
    print("\nRemoving the temporary files, please wait ...\n")
    if sniffles == False:
        print(os.popen('rm svim_reformated.vcf').read())
    if sniffles == True:
        print(os.popen('rm sniffles_reformated.vcf').read())
    print("And it is done! Bye.\n")
