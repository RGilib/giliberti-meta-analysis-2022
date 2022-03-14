#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 24 17:52:56 2021

@author: renatogiliberti
"""
import numpy as np
import pandas as pd
import statsmodels.stats.multitest
from scipy.stats import mannwhitneyu
from scipy.stats import fisher_exact
import argparse

parser = argparse.ArgumentParser(description='perform test to find if microbial species are biomarkers')
parser.add_argument('-d','--data',metavar='',required=True, help='the dataset on which perform the test')
parser.add_argument('-o','--output',metavar='',required=True, help='the path for the output list with pvalues')
parser.add_argument('-c1', '--condition1', metavar='', required=True, help = 'the condition to test')
parser.add_argument('-lvl', '--level', metavar ='', required = True, choices=['abundance','bool0','bool01','bool001','bool0001','bool00001'], help= 'the kind of data to test (relative abundance or boolean data)')
parser.add_argument('-c2', '--condition2', metavar='', required= False, default='control', help='the second condition to test')
parser.add_argument('-t', '--task', metavar='',required = True, default='study_condition', help='the classification task to test')
args=parser.parse_args()
if args.level=='abundance':
   
    data = pd.read_table(args.data ,header = None, index_col = 0)
    data = data.T 
    if 'ibd_huttenhower' in args.data:
        data = data.replace('CD',"UC")
    elif 'ibd_alm' in args.data:
        data = data.replace('CD',"UC")
    elif 'ibd_engstrand_maxee' in args.data:
        data = data.replace('CD',"UC")
    elif 'ra_littman'  in args.data:
        data = data.replace('RA',"PSA")
    elif 'mhe_zhang'  in args.data:
        data = data.replace('MHE',"CIRR")

    data_control = data[data[args.task] == args.condition2]
    data_disease = data[data[args.task]== args.condition1]
    keepcol=[col for col in data.columns if col.startswith('k__')]
    data_control = data_control[keepcol]
    data_disease = data_disease[keepcol]
    data = data[keepcol]
    columns = list(data) 
    p_val = []
    for column in columns :
        p_val.append(mannwhitneyu(pd.to_numeric(data_control[column],downcast='float'), pd.to_numeric(data_disease[column],downcast='float'), alternative='two-sided')[1])
        
    corr_pval= np.array(statsmodels.stats.multitest.multipletests(p_val, alpha=0.05, method = "fdr_bh")[1])
    bact = np.array(columns)
    
    corr_pval = np.vstack((bact, corr_pval)).T
    
    df = pd.DataFrame(corr_pval)
    for i in range(len(df)):
        df.iloc[i,1]=pd.to_numeric(df.iloc[i,1],downcast='float')
    df = df[df.iloc[:,1]<0.05]
    df.to_csv(args.output,index=False)

else:
    data = pd.read_table(args.data ,header = None, index_col = 0)
    data = data.T
    if 'ibd_huttenhower' in args.data:
        data = data.replace('CD',"UC")
    elif 'ibd_alm' in args.data:
        data = data.replace('CD',"UC")
    elif 'ibd_engstrand_maxee' in args.data:
        data = data.replace('CD',"UC")
    elif 'ra_littman'  in args.data:
        data = data.replace('RA',"PSA")
    elif 'mhe_zhang'  in args.data:
        data = data.replace('MHE',"CIRR")

    data_control = data[data[args.task] == args.condition2]
    data_disease = data[data[args.task] == args.condition1]
    keepcol=[col for col in data.columns if col.startswith('k__')]
    data_control = data_control[keepcol]
    data_disease = data_disease[keepcol]
    data = data[keepcol]
    columns = list(data) 
    p_val = []
    for column in columns:
        p_val.append(fisher_exact([[sum(data_control[column].values == "1.0") , sum(data_control[column].values == "0.0")], [sum(data_disease[column].values == "1.0"), sum(data_disease[column].values == "0.0")]])[1])               

    corr_pval= np.array(statsmodels.stats.multitest.multipletests(p_val, alpha=0.05, method = 'fdr_bh')[1])
     
    bact = np.array(columns)

    corr_pval = np.vstack((bact, corr_pval)).T

    df = pd.DataFrame(corr_pval)

    for i in range(len(df)):
        df.iloc[i,1]=pd.to_numeric(df.iloc[i,1],downcast='float')
    df = df[df.iloc[:,1]<0.05]
    df.to_csv(args.output,index=False)
