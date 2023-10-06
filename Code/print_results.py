# -*- coding: utf-8 -*-
"""
Created on Wed Mar 21 15:46:19 2018

@author: kisha_000
"""

import matplotlib
#backend compatible with threading
matplotlib.use('Agg')
import os, fitz
from matplotlib import pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from scipy.cluster.hierarchy import dendrogram,cophenet
import traceback as tb
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

def print_latex_code(df_clusters,filename,cluster_number):
    latex_code = '\\textbf{Cluster '+str(cluster_number) + ' - '+ str(len(df_clusters)) + ' patients} \n' +\
                 ' \n\\vspace{3mm} \n \n' +\
                '\\begin{tabular}{cc} \n\hline \n'
    latex_code = latex_code + 'id\_patient & PE Temporal sequences \\\ \n\hline \n'
    for index,row in df_clusters.iterrows():
        latex_code = latex_code + str(row['id_patient']) +'\t & \t'+ str(row['aux_encode']) + '\t \\\ \n' 
    
    latex_code = latex_code + '\end{tabular} \n \n\\vspace{5mm} \n'

    return latex_code                            
    
#print the resulting clusters and the correspondent code for latex tables
def print_clusters(k,partition_found,df_encoded,filename):
        
    text_file = open(filename, "w")
    text_file.write(filename)
    text_file.write("\n")
    for c in range(0,k):
        text_file.write("\n")
        text_file.write("Cluster %s - %s elements" % (str(c+1) , len(partition_found[c])))
        text_file.write("\n")
        text_file.write("%s" % df_encoded.loc[partition_found[c]])
        text_file.write("\n \n")
        latex_code = print_latex_code(df_encoded.loc[partition_found[c]],filename,c+1)
        text_file.write("%s" % latex_code)
        text_file.write("\n \n")
        
    text_file.close()
    
def print_clusters_csv(k,partition_found,df_encoded,directory):
    #create directory
    try:
        if not os.path.exists('./'+directory):
            os.makedirs('./'+directory)
    except OSError:
        print ('Error: Creating directory. ' +  './'+directory)
        
    for c in range(0,k):
        cluster_name = 'Cluster '+ str(c+1) + ' - ' + str(len(partition_found[c])) + ' elements.csv'
        df_encoded.loc[partition_found[c]].to_csv(directory+cluster_name,encoding='utf-8', index=False)

def print_pdf_images(df_avgs, df_stds, gapPenalty, temporalPenalty, method, figurePartition):
    try:
        directory = './partial_analysis'
        if not os.path.exists(directory):
            os.makedirs(directory)
        pdfFileName = 'partial_analysis_' + str(round(gapPenalty, 2)) + '.pdf'      

        path = os.path.join(directory, pdfFileName)
        pdfFileOpen = PdfPages(path)
        pdfFileOpen.savefig(figurePartition)
        plt.close(figurePartition)
        
        fig1 = plt.figure(figsize=(10,5))
        ax = plt.gca()
        ax.xaxis.set_visible(False)
        ax.yaxis.set_visible(False)
        ax.axis('tight')
        ax.axis('off')
        colLabels=df_avgs.loc[:, df_avgs.columns != 'k_score_avg'].columns
        cell_text = []
        for row in range(len(df_avgs)):
            cell_text.append(df_avgs.iloc[row,0:-1].round(decimals=3))
        plt.title('Average values of five clustering indices \n gap: %.2f, Tp: %.2f, %s link' %(gapPenalty, temporalPenalty, method))
        plt.table(cellText=cell_text, colLabels=colLabels, loc='center',cellLoc='center',fontsize=20)
        pdfFileOpen.savefig(fig1)
        plt.close(fig1)

            #bar chart of standard deviation - standard deviation of all measures
            # Create a figure instance
        #    plt.figure(2)
        #    df_stds.loc[:,df_stds.columns != 'k'].plot.bar(figsize=(15,8))
        #    plt.title('Standard deviation of five measures versus number of clusters',fontsize=25)
        #    plt.xlabel('Number of clusters',labelpad=20,fontsize=20)    
        #    plt.ylabel('Standard deviation',labelpad=10,fontsize=20)    
        #    plt.xticks(size = 20)
        #    plt.yticks(size = 20)
        #    plt.show()
        
        fig2 = plt.figure(3)
        df_stds.loc[:,'Adjusted Rand'].plot.bar(figsize=(15,8),color='forestgreen')
        plt.title('Standard deviation of Adjusted Rand versus number of clusters \n gap: %.2f, Tp: %.2f, %s link' %(gapPenalty, temporalPenalty, method),fontsize=25)
        plt.xlabel('Number of clusters',labelpad=20,fontsize=15)    
        plt.ylabel('Standard deviation',labelpad=10,fontsize=15)    
        plt.xticks(size = 20)
        plt.yticks(size = 20)
        #plt.show()
        pdfFileOpen.savefig(fig2)
        plt.close(fig2)
        
        pdfFileOpen.close()
    except:
        print(tb.format_exc(), flush=True)
        #k=input("press close to exit")

def print_full_analysis(pdfFileName, gap_values):
    try:
        directory = '.\partial_analysis'
        fullAnalysisPdf = fitz.open()
        
        for gap in gap_values:
            partialFileName = 'partial_analysis_' + str(round(gap, 2)) + '.pdf'
            path = os.path.join(directory, partialFileName)
            if not os.path.isfile(path):
                continue
            with fitz.open(path) as partial_images:
                fullAnalysisPdf.insert_pdf(partial_images)
            os.remove(path)
        os.rmdir(directory)
        fullAnalysisPdf.save(pdfFileName)

    except:
        print(tb.format_exc(), flush=True)
        #k=input("press close to exit")
