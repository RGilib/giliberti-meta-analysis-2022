#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 28 11:15:20 2022

@author: renatogiliberti
"""
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import scipy
from scipy.stats import spearmanr


import matplotlib as mpl
mpl.use('svg')
from matplotlib.patches import Patch
import colorsys
import seaborn as sns; sns.set_theme(color_codes=True)
import matplotlib.colors as mcolors
from distinctipy import distinctipy

def pvalue_s(df):
    pvals=[]
    for i in range(len(df.columns)):
        tmp=[]
        for j in range(len(df.columns)):
            t= df[[df.columns[i],df.columns[j]]].dropna()
            tmp.append(spearmanr(t.iloc[:,0],t.iloc[:,1])[1])
        pvals.append(tmp)
    return pd.DataFrame([i for i in pvals], index=df.columns,columns=df.columns)

dataset = ['JieZ_2017.txt','ChngKR_2016.txt','YeZ_2018.txt','RaymondF_2016.txt','QinN_2014.txt','FengQ_2015.txt','GuptaA_2019.txt','HanniganGD_2017.txt','ThomasAM_2018a.txt','ThomasAM_2018b.txt','VogtmannE_2016.txt','WirbelJ_2018.txt','YachidaS_2019.txt','YuJ_2015.txt','ZellerG_2014.txt','LiJ_2017.txt','IjazUZ_2017.txt','NielsenHB_2014.txt','GhensiP_2019_m.txt','GhensiP_2019.txt','Castro-NallarE_2015.txt','Heitz-BuschartA_2016.txt','KosticAD_2015.txt','KarlssonFH_2013.txt','QinJ_2012.txt','HMP_2012.txt']
task = ["study_condition","study_condition","study_condition","study_condition","study_condition","study_condition","study_condition","study_condition","study_condition","study_condition","study_condition","study_condition","study_condition","study_condition","study_condition","study_condition","study_condition","study_condition","study_condition","study_condition","study_condition","study_condition","study_condition","study_condition","study_condition",'body_site']
study_condition_controls = ["control","control","control","control","control","control","control","control","control","control","control","control","control","control","control","control","control","control","control","control","control","control","control","control","control",'stool']
study_condition_cases = ['ACVD','AD','BD','cephalosporins','cirrhosis','CRC','CRC','CRC','CRC','CRC','CRC','CRC','CRC','CRC','CRC','hypertension','IBD','IBD','mucositis','peri-implantitis','schizophrenia','T1D','T1D','T2D','T2D','oralcavity']
body_site = ['Gut','Skin','Gut','Gut','Gut','Gut','Gut','Gut','Gut','Gut','Gut','Gut','Gut','Gut','Gut','Gut','Gut','Gut','Oral','Oral','Oral','Gut','Gut','Gut','Gut',np.nan]
numbers=pd.DataFrame()
numbers_boolean = pd.DataFrame()

for i in range(len(dataset)):
    a=pd.read_csv('../data/shotgun/abundance/species/'+dataset[i],sep='\t',index_col='dataset_name').T
    a_cases =  a[a[task[i]]==study_condition_cases[i]]
    a_controls = a[a[task[i]]==study_condition_controls[i]]
    
    numbers.loc[dataset[i].strip('.txt')+' ('+study_condition_cases[i]+')','N_controls']=len(a_controls)
    numbers.loc[dataset[i].strip('.txt')+' ('+study_condition_cases[i]+')','N_cases']=len(a_cases)
    auc=pd.read_csv('../metrics/shotgun/Random_forest/abundance/species/auc.txt',sep='\t',header=None,index_col=0)
    numbers.loc[dataset[i].strip('.txt')+' ('+study_condition_cases[i]+')','AUC_abundance']= round(auc.loc[dataset[i].strip('.txt'),:][2],2)
    auc=pd.read_csv('../metrics/shotgun/Random_forest/bool0/species/auc.txt',sep='\t',header=None,index_col=0)
    numbers.loc[dataset[i].strip('.txt')+' ('+study_condition_cases[i]+')','AUC_presence/absence']= round(auc.loc[dataset[i].strip('.txt'),:][2],2)
    species = pd.read_csv('../results/biomarkers/shotgun/abundance/species/'+dataset[i].strip('.txt')+'.csv',sep=',',header=None)
    numbers.loc[dataset[i].strip('.txt')+' ('+study_condition_cases[i]+')','N_species_abundance']=len(species)
    species = pd.read_csv('../results/biomarkers/shotgun/bool0/species/'+dataset[i].strip('.txt')+'.csv',sep=',',header=None)
    numbers.loc[dataset[i].strip('.txt')+' ('+study_condition_cases[i]+')','N_species_presence/absence']=len(species)
    numbers.loc[dataset[i].strip('.txt')+' ('+study_condition_cases[i]+')','diseases']=study_condition_cases[i]
    auprc = pd.read_csv('../metrics/shotgun/Random_forest/abundance/species/auprc.txt',sep='\t',header=None,index_col=0)
    numbers.loc[dataset[i].strip('.txt')+' ('+study_condition_cases[i]+')','AUPRC_abundance']= round(auprc.loc[dataset[i].strip('.txt'),:][2],2)
    auprc = pd.read_csv('../metrics/shotgun/Random_forest/bool0/species/auprc.txt',sep='\t',header=None,index_col=0)
    numbers.loc[dataset[i].strip('.txt')+' ('+study_condition_cases[i]+')','AUPRC_presence/absence']= round(auprc.loc[dataset[i].strip('.txt'),:][2],2)
numbers=numbers.rename_axis('Dataset_name')
numbers_boolean.index = numbers.index
numbers_boolean.loc[:,'Abundance']=numbers.loc[:,'AUC_abundance']
numbers_boolean.loc[:,'Presence/Absence']=numbers.loc[:,'AUC_presence/absence']
lvl = ['bool00001','bool0001','bool001','bool01']

for i in range(len(lvl)):
    for j in range(len(dataset)):
        auc = pd.read_csv('../metrics/shotgun/Random_forest/'+lvl[i]+'/species/auc.txt',sep='\t',index_col=0,header=None)
        numbers_boolean.loc[dataset[j].strip('.txt')+' ('+study_condition_cases[j]+')',lvl[i]] = round(auc.loc[dataset[j].strip('.txt'),:][2],2)
for j in range(len(dataset)):
    numbers_boolean.loc[dataset[j].strip('.txt')+' ('+study_condition_cases[j]+')','disease']= study_condition_cases[j]
numbers_boolean.loc[:,'species_Abundance']=numbers['N_species_abundance']
numbers_boolean.loc[:,'species_Presence/absence'] = numbers['N_species_presence/absence']
for i in range(len(lvl)):
    for j in range(len(dataset)):
        species = pd.read_csv('../results/biomarkers/shotgun/'+lvl[i]+'/species/'+dataset[j].strip('.txt')+'.csv',sep='\t',header=None)
        numbers_boolean.loc[dataset[j].strip('.txt')+' ('+study_condition_cases[j]+')','species'+lvl[i]] = len(species)
newcols=['Threshold = 0.0001','Threshold = 0.001','Threshold = 0.01','Threshold = 0.1','Abundance','Presence/Absence','Threshold = 0.0001','Threshold = 0.001','Threshold = 0.01','Threshold = 0.1']
cols = list(numbers_boolean.columns)
unwanted_num = {'Abundance', 'Presence/Absence','disease'}
cols = [ele for ele in cols if ele not in unwanted_num]
rename_d= dict(zip(cols , newcols))
numbers_boolean=numbers_boolean.rename(columns=rename_d)
numbers_tass = pd.DataFrame()
taxonomical = ['genus','family','order']
numbers_tass.index = numbers_boolean.index
numbers_tass.loc[:,'Species Abundance'] = numbers_boolean.iloc[:,0]
numbers_tass.loc[:,'Species Presence/Absence'] = numbers_boolean.iloc[:,1]
for i in range(len(dataset)):
    for j in range(len(taxonomical)):
        auc_ab = pd.read_csv('../metrics/shotgun/Random_forest/abundance/'+taxonomical[j]+'/auc.txt',sep='\t',index_col=0,header=None)
        auc_pa = pd.read_csv('../metrics/shotgun/Random_forest/bool0/'+taxonomical[j]+'/auc.txt',sep='\t',index_col=0,header=None)
        numbers_tass.loc[dataset[i].strip('.txt')+' ('+study_condition_cases[i]+')',taxonomical[j]+' Abundance'] = round(auc_ab.loc[dataset[i].strip('.txt'),:][2] ,2)
        numbers_tass.loc[dataset[i].strip('.txt')+' ('+study_condition_cases[i]+')',taxonomical[j]+' Presence/absence'] = round(auc_pa.loc[dataset[i].strip('.txt'),:][2],2) 
numbers_tass.loc[:,'species_Species Abundance']=numbers_boolean.iloc[:,7]
numbers_tass.loc[:,'species_Species Presence/absence']=numbers_boolean.iloc[:,8]
for i in range(len(dataset)):
    for j in range(len(taxonomical)):
        species_ab = pd.read_csv('../results/biomarkers/shotgun/abundance/'+taxonomical[j]+'/'+dataset[i].strip('.txt')+'.csv',sep='\t',header=None)
        species_pa = pd.read_csv('../results/biomarkers/shotgun/bool0/'+taxonomical[j]+'/'+dataset[i].strip('.txt')+'.csv',sep='\t',header=None)
        numbers_tass.loc[dataset[i].strip('.txt')+' ('+study_condition_cases[i]+')',taxonomical[j]+'_species Abundance'] = len(species_ab)
        numbers_tass.loc[dataset[i].strip('.txt')+' ('+study_condition_cases[i]+')',taxonomical[j]+'_species Presence/absence'] = len(species_pa)
numbers_tass.loc[:,'diseases']=numbers.diseases

cols= numbers_tass.columns
new_cols=['Species Abundance', 'Species Presence/Absence', 'Genus Abundance',\
       'Genus Presence/Absence', 'Family Abundance', 'Family Presence/Absence',\
       'Order Abundance', 'Order Presence/Absence',\
       'Species Abundance', 'Species Presence/Absence', 'Genus Abundance',\
              'Genus Presence/Absence', 'Family Abundance', 'Family Presence/Absence',\
              'Order Abundance', 'Order Presence/Absence']
rename_di = dict(zip(cols,new_cols))
numbers_tass = numbers_tass.rename(columns=rename_di)
numbers_class=pd.DataFrame()
numbers_class.index=numbers_tass.index
numbers_class.loc[:,'ABUNDANCE RF']=numbers.loc[:,'AUC_abundance']
numbers_class.loc[:, 'PRESENCE/ABSENCE RF'] = numbers.loc[:,'AUC_presence/absence']

classifi=['lasso','enet','svm','lsvm']
classifier=['Lasso','ENet','LSVM','SVM']
lvl=['abundance','bool0']
for i in range(len(dataset)):
    for j in range(len(classifier)):
            auc_ab = pd.read_csv('../metrics/shotgun/'+classifi[j]+'/abundance/species/auc.txt',sep='\t',index_col=0,header=None)
            auc_pa = pd.read_csv('../metrics/shotgun/'+classifi[j]+'/bool0/species/auc.txt',sep='\t',index_col=0,header=None)
            numbers_class.loc[dataset[i].strip('.txt')+' ('+study_condition_cases[i]+')','ABUNDANCE '+classifier[j]]=round(auc_ab.loc[dataset[i].strip('.txt'),:][2],2)
            numbers_class.loc[dataset[i].strip('.txt')+' ('+study_condition_cases[i]+')','PRESENCE/ABSENCE '+classifier[j]]=round(auc_pa.loc[dataset[i].strip('.txt'),:][2],2)
numbers_class.loc[:,'diseases']=numbers.diseases

taxonomical = ['species','genus','family','order']
lvl = ['abundance','bool0','bool01','bool001','bool0001','bool00001']
dataset = ['ZellerG_2014.txt','YuJ_2015.txt','YachidaS_2019.txt','WirbelJ_2018.txt','YachidaS_2019.txt','VogtmannE_2016.txt','ThomasAM_2018b.txt','ThomasAM_2018a.txt','HanniganGD_2017.txt','GuptaA_2019.txt','FengQ_2015.txt']
numbers_lodo = pd.DataFrame()
for i in range(len(dataset)):
    for j in range(len(lvl)):
        for k in range(len(taxonomical)):
            auc=pd.read_csv('../metrics/shotgun/lodo/Random_forest/'+lvl[j]+'/'+taxonomical[k]+'/auc.txt', sep ='\t',index_col=0,header=None)
            numbers_lodo.loc[dataset[i].strip('.txt')+' (CRC)',taxonomical[k]+' '+lvl[j]]= round(auc.loc[dataset[i].strip('.txt'),:][2],2)
cols = numbers_lodo.columns
new_cols = ['Species Abundance','Genus Abundance',	'Family Abundance',	'Order Abundance'	,'Species Presence/Absence',	'Genus Presence/Absence',	'Family Presence/Absence',	'Order Presence/Absence',	'Species Thr = 0.1',	'Genus Thr = 0.1',\
            'Family Thr = 0.1',	'Order Thr = 0.1',	'Species Thr = 0.01',	'Genus Thr = 0.01'	,'Family Thr = 0.01'	,'Order Thr = 0.01',	'Species Thr = 0.001',	'Genus Thr = 0.001',	'Family Thr = 0.001',	'Order Thr = 0.001',	'Species Thr = 0.0001',	'Genus Thr = 0.0001',	'Family Thr = 0.0001',	'Order Thr = 0.0001']
rename_d = dict(zip(cols,new_cols))
numbers_lodo=numbers_lodo.rename(columns=rename_d)
lvl = ['Threshold = 0.0001','Threshold = 0.001','Threshold = 0.01','Threshold = 0.1']

number_16s= pd.DataFrame()
dataset = ['edd_singh','cdi_schubert','non_cdi_schubert','cdi_vincent_v3v5','cdi_youngster','ob_goodrich','ob_gordon_2008_v2','ob_zupancic','ob_ross','nash_ob_baker','crc_baxter','crc_zeller','crc_zhao','crc_xiang','ibd_gevers_2014','ibd_huttenhower','ibd_alm','ibd_engstrand_maxee','hiv_noguerajulian','hiv_dinh','hiv_lozupone','asd_son','autism_kb','t1d_alkanani','t1d_mejialeon','nash_chan','nash_ob_baker','ra_littman','mhe_zhang','par_scheperjans','nash_chan']
name = ['Singh_2015 (EDD)','Schubert_2014 (CDI)','Schubert_2014 (non-CDI)','Vincent_2013 (CDI)','Youngster_2014 (CDI)','Goodrich_2014 (OB)','Turnbaugh_2009 (OB)','Zupancic_2012 (OB)','Ross_2015 (OB)','Zhu_2013 (OB)','Baxter_2016 (CRC)','Zeller_2014 (CRC)','Wang_2012 (CRC)','Chen_2012 (CRC)','Gevers_2014 (IBD)','Morgan_2012 (IBD)','Papa_2012 (IBD)','Willing_2010 (IBD)','Noguera-Julian_2016 (HIV)','Dinh_2015 (HIV)','Lozupone_2013 (HIV)','Son_2015 (ASD)','Kang_2013 (ASD)','Alkanani_2015 (T1D)','Mejia-Leon_2014 (T1D)','Wang_2012 (CRC)','Zhu_2013 (NASH)','Scher_2013 (ART)','Zhang_2013 (LIV)','Scheperjans_2015 (PAR)','Wong_2013 (NASH)']
study_condition_controls = ["H","H","H","H","H","H","H","H","H","H","H","H","H","H","nonIBD","H","nonIBD","H","H","H","H","H","H","H","H","H","H","H","H","H",'H']
study_condition_cases = ["EDD", "CDI", "nonCDI", "CDI", "CDI", "OB", "OB", "OB", "OB", "nonNASH-OB", "CRC", "CRC", "CRC", "CRC", "CD", "UC", "UC", "UC", "HIV", "HIV", "HIV", "ASD", "ASD", "T1D", "T1D", "NASH","NASH","PSA","CIRR","PAR",'NASH']
export_file = ['edd_singh','cdi_schubert','non_cdi_schubert','cdi_vincent','cdi_youngster','ob_goodrich','ob_turnbaug','ob_zupancic','ob_ross','ob_zhu','crc_baxter','crc_zeller','crc_wang','crc_chen','ibd_gevers_2014','ibd_morgan','ibd_papa','ibd_willing','hiv_noguerajulian','hiv_dinh','hiv_lozupone','asd_son','asd_kang','t1d_alkanani','t1d_mejialeon','nash_wong','nash_zhu','art_scher','mhe_zhang','par_scheperjans','nash_wong']
subs1 = ["UC","UC","UC","PSA","CIRR"]
subs2 = ["CD","CD","CD","RA","MHE"]
for i in range(len(dataset)):
    data = (pd.read_table('../data/16s/abundance/'+ dataset[i]+".csv",header = None, index_col = 0))
    data = data.T         
    number_16s.loc[name[i],'N_controls']=len(data[data['study_condition']==study_condition_controls[i]])
    number_16s.loc[name[i],'N_cases']=len(data[data['study_condition']==study_condition_cases[i]])
    auc = pd.read_csv('../metrics/16s/Random_forest/abundance/auc.txt',sep='\t',header=None,index_col=0)
    number_16s.loc[name[i],'AUC_abundance']=auc.loc[export_file[i],:][2]
    auc = pd.read_csv('../metrics/16s/Random_forest/bool0/auc.txt',sep='\t',header=None,index_col=0)
    number_16s.loc[name[i],'AUC_presence/absence']=auc.loc[export_file[i],:][2]
    species = pd.read_csv('../results/biomarkers/16s/abundance/'+export_file[i]+'.csv',sep=',')
    number_16s.loc[name[i],'N_species_abundance'] = len(species)
    species = pd.read_csv('../results/biomarkers/16s/bool0/'+export_file[i]+'.csv',sep=',')
    number_16s.loc[name[i],'N_species_presence/absence'] = len(species)
    number_16s.loc[name[i],'diseases'] = study_condition_cases[i]
    auprc = pd.read_csv('../metrics/16s/Random_forest/abundance/auprc.txt',sep='\t',header=None,index_col=0)
    number_16s.loc[name[i],'AUPRC_abundance']=auprc.loc[export_file[i],:][2]
    auprc = pd.read_csv('../metrics/16s/Random_forest/bool0/auprc.txt',sep='\t',header=None,index_col=0)
    number_16s.loc[name[i],'AUPRC_presence/absence']=auprc.loc[export_file[i],:][2]
dataset = ['JieZ_2017.txt','ChngKR_2016.txt','YeZ_2018.txt','RaymondF_2016.txt','QinN_2014.txt','FengQ_2015.txt','GuptaA_2019.txt','HanniganGD_2017.txt','ThomasAM_2018a.txt','ThomasAM_2018b.txt','VogtmannE_2016.txt','WirbelJ_2018.txt','YachidaS_2019.txt','YuJ_2015.txt','ZellerG_2014.txt','LiJ_2017.txt','IjazUZ_2017.txt','NielsenHB_2014.txt','GhensiP_2019_m.txt','GhensiP_2019.txt','Castro-NallarE_2015.txt','Heitz-BuschartA_2016.txt','KosticAD_2015.txt','KarlssonFH_2013.txt','QinJ_2012.txt','HMP_2012.txt']
task = ["study_condition","study_condition","study_condition","study_condition","study_condition","study_condition","study_condition","study_condition","study_condition","study_condition","study_condition","study_condition","study_condition","study_condition","study_condition","study_condition","study_condition","study_condition","study_condition","study_condition","study_condition","study_condition","study_condition","study_condition","study_condition",'body_site']
study_condition_controls = ["control","control","control","control","control","control","control","control","control","control","control","control","control","control","control","control","control","control","control","control","control","control","control","control","control",'stool']
study_condition_cases = ['ACVD','AD','BD','cephalosporins','cirrhosis','CRC','CRC','CRC','CRC','CRC','CRC','CRC','CRC','CRC','CRC','hypertension','IBD','IBD','mucositis','peri-implantitis','schizophrenia','T1D','T1D','T2D','T2D','oralcavity']
body_site = ['Gut','Skin','Gut','Gut','Gut','Gut','Gut','Gut','Gut','Gut','Gut','Gut','Gut','Gut','Gut','Gut','Gut','Gut','Oral','Oral','Oral','Gut','Gut','Gut','Gut',np.nan]
#tab1
tab1 = pd.DataFrame()
for i in range(len(dataset)):
    tab1.loc[dataset[i].strip('.txt'),'body site'] = body_site[i]
    tab1.loc[dataset[i].strip('.txt'),'# controls'] = int(numbers.loc[dataset[i].strip('.txt')+' ('+study_condition_cases[i]+')','N_controls'])
    tab1.loc[dataset[i].strip('.txt'),'Cases']= study_condition_cases[i]
    tab1.loc[dataset[i].strip('.txt'),'# cases']= int(numbers.loc[dataset[i].strip('.txt')+' ('+study_condition_cases[i]+')','N_cases'])
tab1=tab1.rename(columns={0:'Dataset name'})

tab1.drop('HMP_2012').to_excel('../results/fig_and_tables/Table_1.xlsx')
#fig1
y= list(numbers.drop('HMP_2012 (oralcavity)').index)
x= list(numbers.drop('HMP_2012 (oralcavity)')['N_controls'])
x2=list(numbers.drop('HMP_2012 (oralcavity)')['N_cases'])
x3=list(numbers.drop('HMP_2012 (oralcavity)')['AUC_abundance'])
x4=list(numbers.drop('HMP_2012 (oralcavity)')['AUC_presence/absence'])
x5=list(numbers.drop('HMP_2012 (oralcavity)')['N_species_abundance'])
x6=list(numbers.drop('HMP_2012 (oralcavity)')['N_species_presence/absence'])
x7=list(numbers.drop('HMP_2012 (oralcavity)')['AUPRC_abundance'])
x8=list(numbers.drop('HMP_2012 (oralcavity)')['AUPRC_presence/absence'])
fig, (ax1, ax2, ax3, ax4)= plt.subplots(nrows=1, ncols=4,  sharey=True,  figsize=(30,30))
ax1.barh(y,x,label='Controls', color= 'darkcyan')
ax1.barh(y,x2,left =x, label='Cases', color='orange')
ax1.set_ylabel('Dataset',fontsize=20)
ax1.set_xlabel('Number of samples',fontsize=20)
ax1.legend(loc='upper right', prop={'size': 15}, labelspacing=1)
ax1.set_facecolor('white')
ax1.tick_params(axis='y', labelsize=20)
ax1.tick_params(axis='x', labelsize=20)
ax1.set_xlim(0, 600)
ax1.axvline(600, color='grey', alpha=0.5)
ax1.xaxis.grid(True, color='grey', linestyle='-', alpha=0.5)
ax1.set_title('A',fontsize=35, fontweight='bold')
ax1.patch.set_linewidth(5)



ax1.patch.set_edgecolor('grey')
ax2.hlines(y, 0, x3, linestyle='--')
ax2.hlines(y, x3, x4,  linestyle='--', color='r')
ax2.plot(x3, y, 'D',markersize=15, label= 'Relative abundance')
ax2.plot(x4, y, 'ro',markersize=11, label= 'Presence/Absence')
ax2.set_xlabel('AUC',fontsize=20)
ax2.legend(loc='best', prop={'size': 15}, labelspacing=1)
ax2.set_facecolor('white')
ax2.tick_params(axis='x', labelsize=20)
ax2.xaxis.grid(True, color='grey', linestyle='-', alpha=0.5)
ax2.set_xlim(0.5, 1)
ax2.axvline(1, color='grey', alpha=0.5)
ax2.set_title('B',fontsize=35, fontweight='bold')
ax2.patch.set_linewidth(5)
ax2.patch.set_edgecolor('grey')

ax3.hlines(y, 0, x7, linestyle='--')
ax3.hlines(y, x7, x8,  linestyle='--', color='r')
ax3.plot(x7, y, 'D',markersize=15, label= 'Relative abundance')
ax3.plot(x8, y, 'ro',markersize=11, label= 'Presence/Absence')
ax3.set_xlabel('AUPRC',fontsize=20)
ax3.legend(loc='best', prop={'size': 15}, labelspacing=1)
ax3.set_facecolor('white')
ax3.tick_params(axis='x', labelsize=20)
ax3.xaxis.grid(True, color='grey', linestyle='-', alpha=0.5)
ax3.set_xlim(0.5, 1)
ax3.axvline(1, color='grey', alpha=0.5)
ax3.set_title('C',fontsize=35, fontweight='bold')
ax3.patch.set_linewidth(5)
ax3.patch.set_edgecolor('grey')


ax4.hlines(y, -10, x5, linestyle='--')
ax4.hlines(y, x5, x6, linestyle='--', color='r')
ax3.hlines(y,x7,x8,linestyle='--', color='r')

ax4.plot(x5, y, 'D',markersize=15, label= 'Relative abundance')
ax4.plot(x6, y, 'ro',markersize=11, label= 'Presence/Absence')
ax4.set_xlabel('Species, q < 0.05',fontsize=20)
ax4.xaxis.grid(True, color='grey', linestyle='-', alpha=0.5)
ax4.legend(loc='best', prop={'size': 15}, labelspacing=1)
ax4.set_facecolor('white')
ax4.tick_params(axis='x', labelsize=20)
ax4.set_xlim(-10, 250)
ax4.axvline(250, color='grey', alpha=0.5)
ax4.set_title('D',fontsize=35, fontweight='bold')
ax4.patch.set_linewidth(5)
ax4.patch.set_edgecolor('grey')
fig.tight_layout()
plt.savefig('../results/fig_and_tables/Fig1.svg')
plt.savefig('../results/fig_and_tables/Fig1.png')

fig,(ax1)= plt.subplots(nrows=1, ncols=1,  figsize=(5,5))
ax1.scatter(x3,x4)
ax1.plot( [0,1],[0,1],"--r" )
ax1.title.set_text('AUC comparison')
ax1.set_xlabel('AUC relative abundance')
ax1.set_ylabel('AUC presence/absence')
ax1.set_xlim(0,1)
ax1.set_ylim(0,1)

plt.tight_layout()
plt.savefig('../results/fig_and_tables/FigS1.svg')
plt.savefig('../results/fig_and_tables/FigS1.png')

fig,(ax1)=plt.subplots(nrows=1,ncols=1,figsize=(10,10))
ax1.scatter(x3,x7,marker='D',label='Relative abundance')
ax1.scatter(x4,x8,marker='o',c='r',label='Presence/Absence')
corrauc_auprc=scipy.stats.spearmanr(x3,x7).correlation
corrauc_auprc2=scipy.stats.spearmanr(x4,x8).correlation
ax1.set_xlabel('AUC')
ax1.set_ylabel('AUPRC')
ax1.set_xlim(0.5,1)
ax1.set_ylim(0.5,1)
ax1.plot( [0,1],[0,1],"--r" )
ax1.title.set_text('AUC vs AUPRC')
ax1.legend(loc='best', prop={'size': 15}, labelspacing=1)
fig.tight_layout()

plt.savefig('../results/fig_and_tables/FigS2.svg')
plt.savefig('../results/fig_and_tables/FigS2.png')

#fig2
y= list(number_16s.index)
x= list(number_16s['N_controls'])
x2= list(number_16s['N_cases'])
x3= list(number_16s['AUC_abundance'])
x4= list(number_16s['AUC_presence/absence'])
x5= list(number_16s['N_species_abundance'])
x6= list(number_16s['N_species_presence/absence'])
x7=list(number_16s['AUPRC_abundance'])
x8=list(number_16s['AUPRC_presence/absence'])
fig, (ax1, ax2, ax3, ax4)= plt.subplots(nrows=1, ncols=4,  sharey=True,  figsize=(20,20))
ax1.barh(y,x,label='Controls', color= 'darkcyan')
ax1.barh(y,x2,left =x, label='Cases', color='orange')
ax1.set_ylabel('Dataset',fontsize=20)
ax1.set_xlabel('Number of samples',fontsize=20)
ax1.legend(loc='upper right', prop={'size': 15}, labelspacing=1)
ax1.set_facecolor('white')
ax1.tick_params(axis='y', labelsize=20)
ax1.tick_params(axis='x', labelsize=20)
ax1.set_xlim(0, 600)
ax1.axvline(600, color='grey', alpha=0.5)
ax1.xaxis.grid(True, color='grey', linestyle='-', alpha=0.5)
ax1.set_title('A',fontsize=35, fontweight='bold')
ax1.patch.set_linewidth(5)
ax1.patch.set_edgecolor('grey')


ax2.hlines(y, 0, x3, linestyle='--')
ax2.hlines(y, x3, x4,  linestyle='--', color='r')
ax2.plot(x3, y, 'D',markersize=15, label= 'Relative abundance')
ax2.plot(x4, y, 'ro',markersize=11, label= 'Presence/Absence')
ax2.set_xlabel('AUC',fontsize=20)
ax2.legend(loc='best', prop={'size': 15}, labelspacing=1)
ax2.set_facecolor('white')
ax2.tick_params(axis='x', labelsize=20)
ax2.xaxis.grid(True, color='grey', linestyle='-', alpha=0.5)
ax2.set_xlim(0.5, 1)
ax2.axvline(1, color='grey', alpha=0.5)
ax2.set_title('B',fontsize=35, fontweight='bold')
ax2.patch.set_linewidth(5)
ax2.patch.set_edgecolor('grey')

ax3.hlines(y, 0, x7, linestyle='--')
ax3.hlines(y, x7, x8,  linestyle='--', color='r')
ax3.plot(x7, y, 'D',markersize=15, label= 'Relative abundance')
ax3.plot(x8, y, 'ro',markersize=11, label= 'Presence/Absence')
ax3.set_xlabel('AUPRC',fontsize=20)
ax3.legend(loc='best', prop={'size': 15}, labelspacing=1)
ax3.set_facecolor('white')
ax3.tick_params(axis='x', labelsize=20)
ax3.xaxis.grid(True, color='grey', linestyle='-', alpha=0.5)
ax3.set_xlim(0.5, 1)
ax3.axvline(1, color='grey', alpha=0.5)
ax3.set_title('C',fontsize=35, fontweight='bold')
ax3.patch.set_linewidth(5)
ax3.patch.set_edgecolor('grey')

ax4.hlines(y, -10, x5, linestyle='--')
ax4.hlines(y, x5, x6, linestyle='--', color='r')
ax3.hlines(y,x7,x8,linestyle='--', color='r')

ax4.plot(x5, y, 'D',markersize=15, label= 'Relative abundance')
ax4.plot(x6, y, 'ro',markersize=11, label= 'Presence/Absence')
ax4.set_xlabel('Species, q < 0.05',fontsize=20)
ax4.xaxis.grid(True, color='grey', linestyle='-', alpha=0.5)
ax4.legend(loc='best', prop={'size': 15}, labelspacing=1)
ax4.set_facecolor('white')
ax4.tick_params(axis='x', labelsize=20)
ax4.set_xlim(-10, 110)
ax4.axvline(110, color='grey', alpha=0.5)
ax4.set_title('D',fontsize=35, fontweight='bold')
ax4.patch.set_linewidth(5)
ax4.patch.set_edgecolor('grey')
plt.tight_layout()
plt.savefig('../results/fig_and_tables/Fig2.svg')
plt.savefig('../results/fig_and_tables/Fig2.png')

fig,(ax1)= plt.subplots(nrows=1, ncols=1,  figsize=(5,5))
ax1.scatter(x3,x4)
ax1.plot( [0,1],[0,1],"--r" )
ax1.title.set_text('AUC comparison')
ax1.set_xlabel('AUC relative abundance')
ax1.set_ylabel('AUC presence/absence')
ax1.set_xlim(0,1)
ax1.set_ylim(0,1)

plt.tight_layout()
plt.savefig('../results/fig_and_tables/FigS4.svg')
plt.savefig('../results/fig_and_tables/FigS4.png')

#fig3
y= list(numbers_boolean.drop('HMP_2012 (oralcavity)').index)
x= list(numbers_boolean.drop('HMP_2012 (oralcavity)').iloc[:,:6]['Presence/Absence']-numbers_boolean.drop('HMP_2012 (oralcavity)').iloc[:,:6]['Abundance'])
x2=list(numbers_boolean.drop('HMP_2012 (oralcavity)').iloc[:,:6]['Threshold = 0.0001']-numbers_boolean.drop('HMP_2012 (oralcavity)').iloc[:,:6]['Abundance'])
x3=list(numbers_boolean.drop('HMP_2012 (oralcavity)').iloc[:,:6]['Threshold = 0.001']-numbers_boolean.drop('HMP_2012 (oralcavity)').iloc[:,:6]['Abundance'])
x4=list(numbers_boolean.drop('HMP_2012 (oralcavity)').iloc[:,:6]['Threshold = 0.01']-numbers_boolean.drop('HMP_2012 (oralcavity)').iloc[:,:6]['Abundance'])
x5=list(numbers_boolean.drop('HMP_2012 (oralcavity)').iloc[:,:6]['Threshold = 0.1']-numbers_boolean.drop('HMP_2012 (oralcavity)').iloc[:,:6]['Abundance'])

fig, (ax1, ax2)= plt.subplots(nrows=1, ncols=2,  sharey=True,  figsize=(20,20))
ax1.plot(x, y, 'D',markersize=20, label= 'Threshold = 0.0%')
ax1.plot(x2, y, 'ro',markersize=18, label= 'Threshold = 0.0001%')
ax1.plot(x3, y, 'yd',markersize=16, label= 'Threshold = 0.001%')
ax1.plot(x4, y, 'bs',markersize=14, label= 'Threshold = 0.01%')
ax1.plot(x5, y, 'kX',markersize=12, label= 'Threshold = 0.1%')
ax1.set_ylabel('Dataset',fontsize=20)
ax1.set_xlabel('AUC difference:\n Presence/absence - Rel. abundance profiles',fontsize=20)
ax1.legend(loc='upper right', prop={'size': 15}, labelspacing=1)
ax1.set_facecolor('white')
ax1.tick_params(axis='y', labelsize=20)
ax1.tick_params(axis='x', labelsize=20)
ax1.set_xlim(-0.5, 0.5)
ax1.xaxis.grid(True, color='grey', linestyle='-', alpha=0.5)
ax1.set_title('A',fontsize=35, fontweight='bold')
ax1.patch.set_linewidth(5)
ax1.patch.set_edgecolor('grey')
x= list(numbers_boolean.drop('HMP_2012 (oralcavity)').iloc[:,6:]['Presence/Absence']-numbers_boolean.drop('HMP_2012 (oralcavity)').iloc[:,6:]['Abundance'])
x2=list(numbers_boolean.drop('HMP_2012 (oralcavity)').iloc[:,6:]['Threshold = 0.0001']-numbers_boolean.drop('HMP_2012 (oralcavity)').iloc[:,6:]['Abundance'])
x3=list(numbers_boolean.drop('HMP_2012 (oralcavity)').iloc[:,6:]['Threshold = 0.001']-numbers_boolean.drop('HMP_2012 (oralcavity)').iloc[:,6:]['Abundance'])
x4=list(numbers_boolean.drop('HMP_2012 (oralcavity)').iloc[:,6:]['Threshold = 0.01']-numbers_boolean.drop('HMP_2012 (oralcavity)').iloc[:,6:]['Abundance'])
x5=list(numbers_boolean.drop('HMP_2012 (oralcavity)').iloc[:,6:]['Threshold = 0.1']-numbers_boolean.drop('HMP_2012 (oralcavity)').iloc[:,6:]['Abundance'])

ax2.plot(x, y, 'D',markersize=20, label= 'Threshold = 0.0%')
ax2.plot(x2, y, 'ro',markersize=18, label= 'Threshold = 0.0001%')
ax2.plot(x3, y, 'yd',markersize=16, label= 'Threshold = 0.001%')
ax2.plot(x4, y, 'bs',markersize=14, label= 'Threshold = 0.01%')
ax2.plot(x5, y, 'kX',markersize=12, label= 'Threshold = 0.1%')
ax2.set_xlabel('Significant taxa difference:\n Presence/absence - Rel. abundance profiles',fontsize=20)
ax2.legend(loc='best', prop={'size': 15}, labelspacing=1)
ax2.set_facecolor('white')
ax2.tick_params(axis='y', labelsize=20)
ax2.tick_params(axis='x', labelsize=20)
ax2.set_xlim(-170,170)
ax2.xaxis.grid(True, color='grey', linestyle='-', alpha=0.5)
ax2.set_title('B',fontsize=35, fontweight='bold')
ax2.patch.set_linewidth(5)
ax2.patch.set_edgecolor('grey')

plt.tight_layout()
plt.savefig('../results/fig_and_tables/Fig3.svg')
plt.savefig('../results/fig_and_tables/Fig3.png')

y= list(numbers_boolean.drop('HMP_2012 (oralcavity)').index)
x= list(numbers_boolean.iloc[:,:6].drop('HMP_2012 (oralcavity)')['Abundance'])
x1 = list(numbers_boolean.iloc[:,:6].drop('HMP_2012 (oralcavity)')['Presence/Absence'])
x2=list(numbers_boolean.iloc[:,:6].drop('HMP_2012 (oralcavity)')['Threshold = 0.0001'])
x3=list(numbers_boolean.iloc[:,:6].drop('HMP_2012 (oralcavity)')['Threshold = 0.001'])
x4=list(numbers_boolean.iloc[:,:6].drop('HMP_2012 (oralcavity)')['Threshold = 0.01'])
x5=list(numbers_boolean.iloc[:,:6].drop('HMP_2012 (oralcavity)')['Threshold = 0.1'])


#fig4

y= list(numbers_tass.drop('HMP_2012 (oralcavity)').index)
x= list(numbers_tass.iloc[:,:8].drop('HMP_2012 (oralcavity)')['Species Presence/Absence']-numbers_tass.iloc[:,:8].drop('HMP_2012 (oralcavity)')['Species Abundance'])
x2=list(numbers_tass.iloc[:,:8].drop('HMP_2012 (oralcavity)')['Genus Presence/Absence']-numbers_tass.iloc[:,:8].drop('HMP_2012 (oralcavity)')['Genus Abundance'])
x3=list(numbers_tass.iloc[:,:8].drop('HMP_2012 (oralcavity)')['Family Presence/Absence']-numbers_tass.iloc[:,:8].drop('HMP_2012 (oralcavity)')['Family Abundance'])
x4=list(numbers_tass.iloc[:,:8].drop('HMP_2012 (oralcavity)')['Order Presence/Absence']-numbers_tass.iloc[:,:8].drop('HMP_2012 (oralcavity)')['Order Abundance'])


fig, (ax1)= plt.subplots(nrows=1, ncols=1,  sharey=True,  figsize=(15,15))
ax1.plot(x, y, 'D',markersize=20, label= 'Species')
ax1.plot(x2, y, 'ro',markersize=18, label= 'Genus')
ax1.plot(x3, y, 'yd',markersize=16, label= 'Family')
ax1.plot(x4, y, 'bs',markersize=14, label= 'Order')
ax1.set_ylabel('Dataset',fontsize=20)
ax1.set_xlabel('AUC difference:\n Presence/absence - Rel. abundance profiles',fontsize=20)
ax1.legend(loc='best', prop={'size': 15}, labelspacing=1)
ax1.set_facecolor('white')
ax1.tick_params(axis='y', labelsize=20)
ax1.tick_params(axis='x', labelsize=20)
ax1.set_xlim(-0.5, 0.5)
ax1.xaxis.grid(True, color='grey', linestyle='-', alpha=0.5)
ax1.patch.set_linewidth(5)
ax1.patch.set_edgecolor('grey')

plt.tight_layout()
plt.savefig('../results/fig_and_tables/Fig4.svg')
plt.savefig('../results/fig_and_tables/Fig4.png')



#fig5



fig, (ax3)= plt.subplots(nrows=1, ncols=1,  sharey=True,  figsize=(20,20))
y= list(numbers_class.drop('HMP_2012 (oralcavity)').index)
x= list(numbers_class.drop('HMP_2012 (oralcavity)')['PRESENCE/ABSENCE Lasso']-numbers_class.drop('HMP_2012 (oralcavity)')['ABUNDANCE Lasso'])
x2=list(numbers_class.drop('HMP_2012 (oralcavity)')['PRESENCE/ABSENCE ENet']-numbers_class.drop('HMP_2012 (oralcavity)')['ABUNDANCE ENet'])
x3=list(numbers_class.drop('HMP_2012 (oralcavity)')['PRESENCE/ABSENCE LSVM']-numbers_class.drop('HMP_2012 (oralcavity)')['ABUNDANCE LSVM'])
x4=list(numbers_class.drop('HMP_2012 (oralcavity)')['PRESENCE/ABSENCE SVM']-numbers_class.drop('HMP_2012 (oralcavity)')['ABUNDANCE SVM'])
x5=list(numbers_class.drop('HMP_2012 (oralcavity)')['PRESENCE/ABSENCE RF']-numbers_class.drop('HMP_2012 (oralcavity)')['ABUNDANCE RF'])

ax3.plot(x, y, 'D',markersize=22, label= 'Lasso')
ax3.plot(x2, y, 'ro',markersize=20, label= 'ENet')
ax3.plot(x3, y, 'yd',markersize=18, label= 'LSVM')
ax3.plot(x4, y, 'bs',markersize=16, label= 'SVM')
ax3.plot(x5, y, 'kX',markersize=14, label= 'RFs')
ax3.set_xlabel('AUC difference: \n Presence/absence - Rel. abundance profiles',fontsize=20)
ax3.legend(loc='best', prop={'size': 15}, labelspacing=1)
ax3.set_facecolor('white')
ax3.tick_params(axis='y', labelsize=20)
ax3.tick_params(axis='x', labelsize=20)
ax3.set_xlim(-0.5, 0.5)
ax3.xaxis.grid(True, color='grey', linestyle='-', alpha=0.5)
ax3.patch.set_linewidth(5)
ax3.patch.set_edgecolor('grey')


plt.tight_layout()
plt.savefig('../results/fig_and_tables/Fig5.svg')
plt.savefig('../results/fig_and_tables/Fig5.png')

#tabS1
tabS1=pd.DataFrame()
for i in range(len(dataset )):
    tabS1.loc[name[i],'body_site']='Gut'
    tabS1.loc[name[i],'# controls']=number_16s.loc[name[i],'N_controls']
    tabS1.loc[name[i],'# cases']=number_16s.loc[name[i],'N_cases']
    tabS1.loc[name[i], 'Cases']= number_16s.loc[name[i],'diseases']
tabS1 =tabS1.rename(columns={0: 'Dataset name'})
tabS1.to_excel('../results/fig_and_tables/Table_S1.xlsx')

#tabS2
tabS2_auc=pd.DataFrame()
tabS2_auc.loc[:,'Dataset name'] = numbers.index
tabS2_auc= tabS2_auc.set_index('Dataset name')
tabS2_auc['Abundance']=numbers['AUC_abundance']
tabS2_auc.loc[:,'Presence/absence']=numbers.loc[:,'AUC_presence/absence']
colnames=['Thr = 0.0001','Thr = 0.001','Thr = 0.01','Thr = 0.1']
for i in range(len(colnames)):
    tabS2_auc.loc[:,colnames[i]]=numbers_boolean.iloc[:,i+2]
auc = pd.read_csv('../metrics/test/shotgun/Random_forest/bool0/species/auc.txt',sep='\t',header=None,index_col=0)
dataset = ['JieZ_2017.txt','ChngKR_2016.txt','YeZ_2018.txt','RaymondF_2016.txt','QinN_2014.txt','FengQ_2015.txt','GuptaA_2019.txt','HanniganGD_2017.txt','ThomasAM_2018a.txt','ThomasAM_2018b.txt','VogtmannE_2016.txt','WirbelJ_2018.txt','YachidaS_2019.txt','YuJ_2015.txt','ZellerG_2014.txt','LiJ_2017.txt','IjazUZ_2017.txt','NielsenHB_2014.txt','GhensiP_2019_m.txt','GhensiP_2019.txt','Castro-NallarE_2015.txt','Heitz-BuschartA_2016.txt','KosticAD_2015.txt','KarlssonFH_2013.txt','QinJ_2012.txt','HMP_2012.txt']

for i in range(len(dataset)):
    tabS2_auc.loc[dataset[i].strip('.txt')+ ' ('+study_condition_cases[i]+')','P-value Abundance vs Presence/absence']=round(auc.loc[dataset[i].strip('.txt')+'test.txt',:][2],2)
tabS2_f1=pd.DataFrame()
tabS2_precision=pd.DataFrame()
tabS2_recall=pd.DataFrame()
tabS2_auprc=pd.DataFrame()
colnames=['Abundance','Presence/absence','Thr = 0.0001','Thr = 0.001','Thr = 0.01','Thr = 0.1']
lvl=['abundance','bool0','bool00001','bool0001','bool001','bool01']
f1_t =pd.read_csv('../metrics/test/shotgun/Random_forest/bool0/species/f1.txt',sep='\t',header=None,index_col=0)
precision_t =pd.read_csv('../metrics/test/shotgun/Random_forest/bool0/species/precision.txt',sep='\t',header=None,index_col=0)
recall_t =pd.read_csv('../metrics/test/shotgun/Random_forest/bool0/species/recall.txt',sep='\t',header=None,index_col=0)

for i in range(len(colnames)):
    for j in range(len(dataset)):
        f1 = pd.read_csv('../metrics/shotgun/Random_forest/'+lvl[i]+'/species/f1.txt',sep='\t',index_col=0,header=None)
        tabS2_f1.loc[dataset[j].strip('.txt')+ ' ('+study_condition_cases[j]+')',colnames[i]] = f1.loc[dataset[j].strip('.txt'),:][2]
for j in range(len(dataset)):
    tabS2_f1.loc[dataset[j].strip('.txt')+ ' ('+study_condition_cases[j]+')','P-value Abundance vs Presence/absence']=f1_t.loc[dataset[j].strip('.txt')+'test.txt',:][2]
for i in range(len(colnames)):
    for j in range(len(dataset)):
        precision = pd.read_csv('../metrics/shotgun/Random_forest/'+lvl[i]+'/species/precision.txt',sep='\t',index_col=0,header=None)
        tabS2_precision.loc[dataset[j].strip('.txt')+ ' ('+study_condition_cases[j]+')',colnames[i]] = precision.loc[dataset[j].strip('.txt'),:][2]
for j in range(len(dataset)):
    tabS2_precision.loc[dataset[j].strip('.txt')+ ' ('+study_condition_cases[j]+')','P-value Abundance vs Presence/absence']=precision_t.loc[dataset[j].strip('.txt')+'test.txt',:][2]
for i in range(len(colnames)):
    for j in range(len(dataset)):
        recall = pd.read_csv('../metrics/shotgun/Random_forest/'+lvl[i]+'/species/recall.txt',sep='\t',index_col=0,header=None)
        tabS2_recall.loc[dataset[j].strip('.txt')+ ' ('+study_condition_cases[j]+')',colnames[i]] = recall.loc[dataset[j].strip('.txt'),:][2]
for j in range(len(dataset)):
    tabS2_recall.loc[dataset[j].strip('.txt')+ ' ('+study_condition_cases[j]+')','P-value Abundance vs Presence/absence']=recall_t.loc[dataset[j].strip('.txt')+'test.txt',:][2]
colnames=['Abundance','Presence/absence']
lvl=['abundance','bool0']
for i in range(len(colnames)):
    for j in range(len(dataset)):
        auprc = pd.read_csv('../metrics/shotgun/Random_forest/'+lvl[i]+'/species/auprc.txt',sep='\t',index_col=0,header=None)
        tabS2_auprc.loc[dataset[j].strip('.txt')+ ' ('+study_condition_cases[j]+')',colnames[i]] = auprc.loc[dataset[j].strip('.txt'),:][2]

tabS2_auc.drop('HMP_2012 (oralcavity)').to_excel('../results/fig_and_tables/Table_S2_auc.xlsx')
tabS2_f1.drop('HMP_2012 (oralcavity)').to_excel('../results/fig_and_tables/Table_S2_f1.xlsx')
tabS2_precision.drop('HMP_2012 (oralcavity)').to_excel('../results/fig_and_tables/Table_S2_precision.xlsx')
tabS2_recall.drop('HMP_2012 (oralcavity)').to_excel('../results/fig_and_tables/Table_S2_recall.xlsx')
tabS2_auprc.drop('HMP_2012 (oralcavity)').to_excel('../results/fig_and_tables/Table_S2_auprc.xlsx')
#tabS3 
learner=['Random_forest','lasso','enet','lsvm','svm']
colnames=['RF','LASSO','ENET','LSVM','SVM']
tabS3_ra=pd.DataFrame()
tabS3_pa=pd.DataFrame()
for i in range(len(learner)):
    auc_ab=pd.read_csv('../metrics/16s/'+learner[i]+'/abundance/auc.txt',sep='\t',header=None)
    auc_pa=pd.read_csv('../metrics/16s/'+learner[i]+'/bool0/auc.txt',sep='\t',header=None)

    for j in range((30)):
        tabS3_ra.loc[auc_ab[0].iloc[j],colnames[i]]=round(float(auc_ab.iloc[j,2]),2)
        tabS3_pa.loc[auc_pa[0].iloc[j],colnames[i]]=round(float(auc_pa.iloc[j,2]),2)
tabS3_ra.to_excel('../results/fig_and_tables/Table_S3_relative_abundances.xlsx')
tabS3_pa.to_excel('../results/fig_and_tables/Table_S3_presence-absence.xlsx')

#tabS4
tab_S4 = pd.DataFrame()
auc_p=[0.96,0.81,0.91,1,0.77,0.86,0.71,0.89,0.81,0.91,0.64,0.87,0.79,0.78,0.81,0.84,0.84]
dbs=['FengQ_2015.txt','GhensiP_2019_m.txt','GhensiP_2019.txt','GuptaA_2019.txt','HanniganGD_2017.txt','JieZ_2017.txt','KarlssonFH_2013.txt','LiJ_2017.txt','QinJ_2012.txt','QinN_2014.txt','ThomasAM_2018a.txt','ThomasAM_2018b.txt','WirbelJ_2018.txt','YachidaS_2019.txt','YeZ_2018.txt','YuJ_2015.txt','ZellerG_2014.txt']

auc = pd.read_csv('../metrics/shotgun/Random_forest/abundance/species/auc.txt',index_col=0,header=None,sep='\t')
for i in range(len(dbs)):
    tab_S4.loc[dbs[i].strip('.txt'),'AUC in original paper'] = auc_p[i]
    tab_S4.loc[dbs[i].strip('.txt'),'AUC in our_analysis'] = round(auc.loc[dbs[i].strip('.txt')][2],2)
tab_S4.to_excel('../results/fig_and_tables/Table_S4.xlsx')

#tabS5
db = ['Castro-NallarE_2015.csv',"ChngKR_2016.csv","FengQ_2015.csv","GhensiP_2019.csv","GhensiP_2019_m.csv","GuptaA_2019.csv","HanniganGD_2017.csv","Heitz-BuschartA_2016.csv","HMP_2012.csv","IjazUZ_2017.csv", "JieZ_2017.csv", "KarlssonFH_2013.csv", "KosticAD_2015.csv",  "LiJ_2017.csv", "NielsenHB_2014.csv", "QinJ_2012.csv", "QinN_2014.csv", "RaymondF_2016.csv", "ThomasAM_2018a.csv", "ThomasAM_2018b.csv", "VogtmannE_2016.csv", "WirbelJ_2018.csv", "YachidaS_2019.csv", "YeZ_2018.csv", "YuJ_2015.csv", "ZellerG_2014.csv"]
for i in range(len(db)):
    if i == 0:
        d1 = pd.read_table('../results/biomarkers/shotgun/abundance/species/'+db[i], sep=',', names=('species', 'pvalue'))
        d1=d1[1:]
        d1['dataset_name']=db[i]
    else:
        d = pd.read_table('../results/biomarkers/shotgun/abundance/species/'+db[i], sep=',', names=('species', 'pvalue'))
        d=d[1:]
        d['dataset_name']=db[i].strip('.csv')
        d1=d1.append(d)
d1=d1[d1['dataset_name']!='HMP_2012.txt']
tabS5 = pd.pivot_table(d1, values='pvalue', index='species', columns='dataset_name',fill_value=1)
tabS5.to_excel('../results/fig_and_tables/Table_S6.xlsx')



dataset = ['JieZ_2017.txt','ChngKR_2016.txt','YeZ_2018.txt','RaymondF_2016.txt','QinN_2014.txt','FengQ_2015.txt','GuptaA_2019.txt','HanniganGD_2017.txt','ThomasAM_2018a.txt','ThomasAM_2018b.txt','VogtmannE_2016.txt','WirbelJ_2018.txt','YachidaS_2019.txt','YuJ_2015.txt','ZellerG_2014.txt','LiJ_2017.txt','IjazUZ_2017.txt','NielsenHB_2014.txt','GhensiP_2019_m.txt','GhensiP_2019.txt','Castro-NallarE_2015.txt','Heitz-BuschartA_2016.txt','KosticAD_2015.txt','KarlssonFH_2013.txt','QinJ_2012.txt','HMP_2012.txt']
task = ["study_condition","study_condition","study_condition","study_condition","study_condition","study_condition","study_condition","study_condition","study_condition","study_condition","study_condition","study_condition","study_condition","study_condition","study_condition","study_condition","study_condition","study_condition","study_condition","study_condition","study_condition","study_condition","study_condition","study_condition","study_condition",'body_site']
study_condition_controls = ["control","control","control","control","control","control","control","control","control","control","control","control","control","control","control","control","control","control","control","control","control","control","control","control","control",'stool']
study_condition_cases = ['ACVD','AD','BD','cephalosporins','cirrhosis','CRC','CRC','CRC','CRC','CRC','CRC','CRC','CRC','CRC','CRC','hypertension','IBD','IBD','mucositis','peri-implantitis','schizophrenia','T1D','T1D','T2D','T2D','oralcavity']
body_site = ['Gut','Skin','Gut','Gut','Gut','Gut','Gut','Gut','Gut','Gut','Gut','Gut','Gut','Gut','Gut','Gut','Gut','Gut','Oral','Oral','Oral','Gut','Gut','Gut','Gut',np.nan]
level=['abundance','bool0','bool00001','bool0001','bool001','bool01']
colnames=['Abundance','Presence/absence','Thr = 0.0001','Thr = 0.001','Thr = 0.01','Thr = 0.1']
tabS6_thr=pd.DataFrame()
for i in range(len(dataset)):
    for j in range(len(colnames)):
        df = pd.read_csv('../results/biomarkers/shotgun/'+level[j]+'/species/'+dataset[i].strip('txt')+'csv',sep='\t')
        tabS6_thr.loc[dataset[i].strip('.txt')+ ' ('+study_condition_cases[i]+')',colnames[j]]=len(df)


tax=['species','genus','family','order']
tax_col=['SPECIES','GENUS','FAMILY','ORDER']
level=['abundance','bool0']
colnames=['Abundance','Presence/absence']

d1=pd.DataFrame()
for i in range(len(dataset)):
    for j in range(len(tax)):
        for k in range(len(level)):
            if i == 0 and j==0 and k==0 :
                df=pd.read_csv('../results/biomarkers/shotgun/'+level[k]+'/'+tax[j]+'/'+dataset[i].strip('txt')+'csv',sep=',')
                d1.loc[dataset[i].strip('.txt')+' ('+study_condition_cases[i]+')','tax_level']=tax_col[j]
                d1.loc[dataset[i].strip('.txt')+' ('+study_condition_cases[i]+')','level']=colnames[k]
                d1.loc[dataset[i].strip('.txt')+' ('+study_condition_cases[i]+')','N_species']=len(df)
            else:
                d=pd.DataFrame()
                df1=pd.read_csv('../results/biomarkers/shotgun/'+level[k]+'/'+tax[j]+'/'+dataset[i].strip('txt')+'csv',sep='\t')
                d.loc[dataset[i].strip('.txt')+' ('+study_condition_cases[i]+')','tax_level']=tax_col[j]
                d.loc[dataset[i].strip('.txt')+' ('+study_condition_cases[i]+')','level']=colnames[k]
                d.loc[dataset[i].strip('.txt')+' ('+study_condition_cases[i]+')','N_species']=len(df1)
                d1=d1.append(d)
tabS6=pd.pivot_table(d1,values='N_species',index=d1.index,columns=['tax_level','level'])
tabS6_thr.drop('HMP_2012 (oralcavity)').to_excel('../results/fig_and_tables/Table_S7_thr.xlsx')
tabS6.drop('HMP_2012 (oralcavity)').to_excel('../results/fig_and_tables/Table_S7.xlsx')


#tabS7
export_file = ['edd_singh','cdi_schubert','non_cdi_schubert','cdi_vincent','cdi_youngster','ob_goodrich','ob_turnbaug','ob_zupancic','ob_ross','ob_zhu','crc_baxter','crc_zeller','crc_wang','crc_chen','ibd_gevers_2014','ibd_morgan','ibd_papa','ibd_willing','hiv_noguerajulian','hiv_dinh','hiv_lozupone','asd_son','asd_kang','t1d_alkanani','t1d_mejialeon','nash_wong','nash_zhu','art_scher','mhe_zhang','par_scheperjans']
auc_p=[0.96,0.99,0.98,0.91,np.nan,0.67,0.84,0.44,0.49,0.86,0.77,0.82,0.90,0.78,0.71,0.81,0.84,0.66,0.67,0.22,0.92,0.39,0.76,0.71,0.77,0.68,0.93,0.62,.80,0.67]
tabS7_auc=pd.DataFrame()
tabS7_f1=pd.DataFrame()
tabS7_precision=pd.DataFrame()
tabS7_recall=pd.DataFrame()
auc_ab=pd.read_csv('../metrics/16s/Random_forest/abundance/AUC.txt',sep='\t',header=None,index_col=0)
auc_pa=pd.read_csv('../metrics/16s/Random_forest/bool0/AUC.txt',sep='\t',header=None,index_col=0)
f1_ab=pd.read_csv('../metrics/16s/Random_forest/abundance/f1.txt',sep='\t',header=None,index_col=0)
f1_pa=pd.read_csv('../metrics/16s/Random_forest/bool0/f1.txt',sep='\t',header=None,index_col=0)
precision_ab=pd.read_csv('../metrics/16s/Random_forest/abundance/precision.txt',sep='\t',header=None,index_col=0)
precision_pa=pd.read_csv('../metrics/16s/Random_forest/bool0/precision.txt',sep='\t',header=None,index_col=0)
recall_ab=pd.read_csv('../metrics/16s/Random_forest/abundance/recall.txt',sep='\t',header=None,index_col=0)
recall_pa=pd.read_csv('../metrics/16s/Random_forest/bool0/recall.txt',sep='\t',header=None,index_col=0)
auc_test = pd.read_csv('../metrics/test/16s/Random_forest/bool0/auc.txt',sep='\t',header=None,index_col=0)
f1_test = pd.read_csv('../metrics/test/16s/Random_forest/bool0/f1.txt',sep='\t',header=None,index_col=0)
precision_test = pd.read_csv('../metrics/test/16s/Random_forest/bool0/precision.txt',sep='\t',header=None,index_col=0)
recall_test = pd.read_csv('../metrics/test/16s/Random_forest/bool0/recall.txt',sep='\t',header=None,index_col=0)
tabS7_auprc=pd.DataFrame()
auprc_ab=pd.read_csv('../metrics/16s/Random_forest/abundance/auprc.txt',sep='\t',header=None,index_col=0)
auprc_pa=pd.read_csv('../metrics/16s/Random_forest/bool0/auprc.txt',sep='\t',header=None,index_col=0)

for i in range(len(export_file)):
    tabS7_auc.loc[export_file[i],'AUC reported in paper']=auc_p[i]
    tabS7_auc.loc[export_file[i],'AUC_abundance'] = round(auc_ab.loc[export_file[i],:][2],2)
    tabS7_auc.loc[export_file[i],'AUC_presence/absence'] = round(auc_pa.loc[export_file[i],:][2],2)
    tabS7_auc.loc[export_file[i],'p-value abundance vs presence/absence'] = round(auc_test.loc[export_file[i],:][2],2)
    tabS7_f1.loc[export_file[i],'f2_abundance'] = round(f1_ab.loc[export_file[i],:][2],2)
    tabS7_f1.loc[export_file[i],'f2_presence/absence'] = round(f1_pa.loc[export_file[i],:][2],2)
    tabS7_f1.loc[export_file[i],'p-value abundance vs presence/absence'] = round(f1_test.loc[export_file[i],:][2],2)
    tabS7_precision.loc[export_file[i],'precision_abundance'] = round(precision_ab.loc[export_file[i],:][2],2)
    tabS7_precision.loc[export_file[i],'precision_presence/absence'] = round(precision_pa.loc[export_file[i],:][2],2)
    tabS7_precision.loc[export_file[i],'p-value abundance vs presence/absence'] = round(precision_test.loc[export_file[i],:][2],2)
    tabS7_recall.loc[export_file[i],'recall_abundance'] = round(recall_ab.loc[export_file[i],:][2],2)
    tabS7_recall.loc[export_file[i],'recall_presence/absence'] = round(recall_pa.loc[export_file[i],:][2],2)
    tabS7_recall.loc[export_file[i],'p-value abundance vs presence/absence'] = round(recall_test.loc[export_file[i],:][2],2)
    tabS7_auprc.loc[export_file[i],'recall_abundance'] = round(auprc_ab.loc[export_file[i],:][2],2)
    tabS7_auprc.loc[export_file[i],'recall_presence/absence'] = round(auprc_pa.loc[export_file[i],:][2],2)
tabS7_auc.to_excel('../results/fig_and_tables/Table_S5_auc.xlsx')
tabS7_f1.to_excel('../results/fig_and_tables/Table_S5_f1.xlsx')
tabS7_precision.to_excel('../results/fig_and_tables/Table_S5_precision.xlsx')
tabS7_recall.to_excel('../results/fig_and_tables/Table_S5_recall.xlsx')
tabS7_auprc.to_excel('../results/fig_and_tables/Table_S5_auprc.xlsx')



#tabS8
export_file = ['edd_singh','cdi_schubert','non_cdi_schubert','cdi_vincent','cdi_youngster','ob_goodrich','ob_turnbaug','ob_zupancic','ob_ross','ob_zhu','crc_baxter','crc_zeller','crc_wang','crc_chen','ibd_gevers_2014','ibd_morgan','ibd_papa','ibd_willing','hiv_noguerajulian','hiv_dinh','hiv_lozupone','asd_son','asd_kang','t1d_alkanani','t1d_mejialeon','nash_wong','nash_zhu','art_scher','mhe_zhang','par_scheperjans','nash_wong']
name = ['Singh_2015 (EDD)','Schubert_2014 (CDI)','Schubert_2014 (non-CDI)','Vincent_2013 (CDI)','Youngster_2014 (CDI)','Goodrich_2014 (OB)','Turnbaugh_2009 (OB)','Zupancic_2012 (OB)','Ross_2015 (OB)','Zhu_2013 (OB)','Baxter_2016 (CRC)','Zeller_2014 (CRC)','Wang_2012 (CRC)','Chen_2012 (CRC)','Gevers_2014 (IBD)','Morgan_2012 (IBD)','Papa_2012 (IBD)','Willing_2010 (IBD)','Noguera-Julian_2016 (HIV)','Dinh_2015 (HIV)','Lozupone_2013 (HIV)','Son_2015 (ASD)','Kang_2013 (ASD)','Alkanani_2015 (T1D)','Mejia-Leon_2014 (T1D)','Wang_2012 (CRC)','Zhu_2013 (NASH)','Scher_2013 (ART)','Zhang_2013 (LIV)','Scheperjans_2015 (PAR)','Wong_2013 (NASH)']
tabS8=pd.DataFrame()
for i in range(len(name)):
    df = pd.read_csv('../results/biomarkers/16s/abundance/'+export_file[i]+'.csv',sep='\t')
    df2 = pd.read_csv('../results/biomarkers/16s/bool0/'+export_file[i]+'.csv',sep='\t')
    tabS8.loc[name[i],'abundance']=len(df)
    tabS8.loc[name[i],'presence/absence']=len(df2)
tabS8.to_excel('../results/fig_and_tables/Table_S8.xlsx')


tabS9 = pd.DataFrame()
auc_s =pd.read_csv('../metrics/shotgun/Random_forest/abundance/species/auc.txt',header=None,index_col=0,sep='\t')
auc_ge=pd.read_csv('../metrics/shotgun/Random_forest/abundance/genus/auc.txt',header=None,index_col=0,sep='\t')
auc_f=pd.read_csv('../metrics/shotgun/Random_forest/abundance/family/auc.txt',header=None,index_col=0,sep='\t')
auc_o=pd.read_csv('../metrics/shotgun/Random_forest/abundance/order/auc.txt',header=None,index_col=0,sep='\t')
dataset = ['JieZ_2017.txt','ChngKR_2016.txt','YeZ_2018.txt','RaymondF_2016.txt','QinN_2014.txt','FengQ_2015.txt','GuptaA_2019.txt','HanniganGD_2017.txt','ThomasAM_2018a.txt','ThomasAM_2018b.txt','VogtmannE_2016.txt','WirbelJ_2018.txt','YachidaS_2019.txt','YuJ_2015.txt','ZellerG_2014.txt','LiJ_2017.txt','IjazUZ_2017.txt','NielsenHB_2014.txt','GhensiP_2019_m.txt','GhensiP_2019.txt','Castro-NallarE_2015.txt','Heitz-BuschartA_2016.txt','KosticAD_2015.txt','KarlssonFH_2013.txt','QinJ_2012.txt','HMP_2012.txt']
study_condition_cases = ['ACVD','AD','BD','cephalosporins','cirrhosis','CRC','CRC','CRC','CRC','CRC','CRC','CRC','CRC','CRC','CRC','hypertension','IBD','IBD','mucositis','peri-implantitis','schizophrenia','T1D','T1D','T2D','T2D','oralcavity']
auc_svsg =pd.read_csv('../metrics/test/shotgun/Random_forest/genus/auc.txt',header=None,index_col=0,sep='\t')
auc_svsf =pd.read_csv('../metrics/test/shotgun/Random_forest/family/auc.txt',header=None,index_col=0,sep='\t')
auc_svso =pd.read_csv('../metrics/test/shotgun/Random_forest/order/auc.txt',header=None,index_col=0,sep='\t')
for i in range(len(dataset)):
    tabS9.loc[dataset[i].strip('.txt')+' ('+study_condition_cases[i]+')','AUC SPECIES']=round(auc_s.loc[dataset[i].strip('.txt')][2],2)
    tabS9.loc[dataset[i].strip('.txt')+' ('+study_condition_cases[i]+')','AUC GENUS']=round(auc_s.loc[dataset[i].strip('.txt')][2],2)
    tabS9.loc[dataset[i].strip('.txt')+' ('+study_condition_cases[i]+')','pvalue SPECIES vs GENUS']=round(auc_svsg.loc[dataset[i].strip('.txt')][2],2)
    tabS9.loc[dataset[i].strip('.txt')+' ('+study_condition_cases[i]+')','AUC FAMILY']=round(auc_f.loc[dataset[i].strip('.txt')][2],2)
    tabS9.loc[dataset[i].strip('.txt')+' ('+study_condition_cases[i]+')','pvalue SPECIES vs FAMILY']=round(auc_svsf.loc[dataset[i].strip('.txt')][2],2)
    tabS9.loc[dataset[i].strip('.txt')+' ('+study_condition_cases[i]+')','AUC ORDER']=round(auc_s.loc[dataset[i].strip('.txt')][2],2)
    tabS9.loc[dataset[i].strip('.txt')+' ('+study_condition_cases[i]+')','pvalue SPECIES vs ORDER']=round(auc_svso.loc[dataset[i].strip('.txt')][2],2)

tabS9.drop('HMP_2012 (oralcavity)').to_excel('../results/fig_and_tables/Table_S10.xlsx')
metrics = ['auc','f1','precision','recall']
lvl=['abundance','bool0']
colnames=['Abundance','Presence/absence']

for i in range(len(dataset)):
    for j in range(len(lvl)):
        for m in range(len(metrics)):
            if m == 0 and j == 0:
                d1 = pd.read_csv('../metrics/shotgun/Random_forest/'+lvl[j]+'/species/'+metrics[m]+'.txt',header = None,sep='\t')
                d1 = d1.iloc[:,:-1]
                d1['boolean']=colnames[j]
            else:
                d2 = pd.read_csv('../metrics/shotgun/Random_forest/'+lvl[j]+'/species/'+metrics[m]+'.txt',header = None,sep='\t')
                d2 = d2.iloc[:,:-1]

                d2['boolean']=colnames[j]
                d1 = d1.append(d2)
d1= d1.rename(columns={0:'Dataset name',1:'metrics',2:'values'})
d1['level']='species'
d1['values']=round(d1['values'],2)
tabS10_species=pd.pivot_table(d1,values='values',index='Dataset name',columns=['level','metrics','boolean'])

for i in range(len(dataset)):
    for j in range(len(lvl)):
        for m in range(len(metrics)):
            if m == 0 and j == 0:
                d1 = pd.read_csv('../metrics/shotgun/Random_forest/'+lvl[j]+'/genus/'+metrics[m]+'.txt',header = None,sep='\t')
                d1 = d1.iloc[:,:-1]
                d1['boolean']=colnames[j]
                d3 = pd.read_csv('../metrics/test/shotgun/Random_forest/bool0/genus/'+metrics[m]+'.txt',header = None,sep='\t')
                d3['boolean']='P-value Abundance vs Presence/absence'
                d1=d1.append(d3)
            else:
                 d2= pd.read_csv('../metrics/shotgun/Random_forest/'+lvl[j]+'/genus/'+metrics[m]+'.txt',header = None,sep='\t')
                 d2 = d2.iloc[:,:-1]

                 d2['boolean']=colnames[j]
                 d1 = d1.append(d2)
d1= d1.rename(columns={0:'Dataset name',1:'metrics',2:'values'})
d1['level']='genus'
d1['values']=round(d1['values'],2)
tabS10_genus=pd.pivot_table(d1,values='values',index='Dataset name',columns=['level','metrics','boolean'])

for i in range(len(dataset)):
    for j in range(len(lvl)):
        for m in range(len(metrics)):
            if m == 0 and j == 0:
                d1 = pd.read_csv('../metrics/shotgun/Random_forest/'+lvl[j]+'/family/'+metrics[m]+'.txt',header = None,sep='\t')
                d1 = d1.iloc[:,:-1]
                d1['boolean']=colnames[j]
                d3 = pd.read_csv('../metrics/test/shotgun/Random_forest/bool0/family/'+metrics[m]+'.txt',header = None,sep='\t')
                d3['boolean']='P-value Abundance vs Presence/absence'
                d1=d1.append(d3)
            else:
                 d2= pd.read_csv('../metrics/shotgun/Random_forest/'+lvl[j]+'/family/'+metrics[m]+'.txt',header = None,sep='\t')
                 d2 = d2.iloc[:,:-1]

                 d2['boolean']=colnames[j]
                 d1 = d1.append(d2)
d1= d1.rename(columns={0:'Dataset name',1:'metrics',2:'values'})
d1['level']='family'
d1['values']=round(d1['values'],2)
tabS10_family=pd.pivot_table(d1,values='values',index='Dataset name',columns=['level','metrics','boolean'])

for i in range(len(dataset)):
    for j in range(len(lvl)):
        for m in range(len(metrics)):
            if m == 0 and j == 0:
                d1 = pd.read_csv('../metrics/shotgun/Random_forest/'+lvl[j]+'/order/'+metrics[m]+'.txt',header = None,sep='\t')
                d1 = d1.iloc[:,:-1]
                d1['boolean']=colnames[j]
                d3 = pd.read_csv('../metrics/test/shotgun/Random_forest/bool0/order/'+metrics[m]+'.txt',header = None,sep='\t')
                d3['boolean']='P-value Abundance vs Presence/absence'
                d1=d1.append(d3)
            else:
                 d2= pd.read_csv('../metrics/shotgun/Random_forest/'+lvl[j]+'/order/'+metrics[m]+'.txt',header = None,sep='\t')
                 d2 = d2.iloc[:,:-1]

                 d2['boolean']=colnames[j]
                 d1 = d1.append(d2)
d1= d1.rename(columns={0:'Dataset name',1:'metrics',2:'values'})
d1['level']='order'
d1['values']=round(d1['values'],2)
tabS10_order=pd.pivot_table(d1,values='values',index='Dataset name',columns=['level','metrics','boolean'])

tabS10_species.to_excel('../results/fig_and_tables/Table_S11_species.xlsx')
tabS10_genus.to_excel('../results/fig_and_tables/Table_S11_genus.xlsx')
tabS10_family.to_excel('../results/fig_and_tables/Table_S11_family.xlsx')
tabS10_order.to_excel('../results/fig_and_tables/Table_S11_order.xlsx')






col_names=['Abundance',	'Presence/absence',	'Thr = 0.1',	'Thr = 0.01',	'Thr = 0.001',	'Thr = 0.0001']
filter_col = [col for col in numbers_lodo if col.startswith('Species')]
tabS11_species=numbers_lodo[filter_col]
rename_dict=dict(zip(tabS11_species.columns,col_names))
tabS11_species=tabS11_species.rename(columns=rename_dict)

filter_col = [col for col in numbers_lodo if col.startswith('Genus')]
tabS11_genus=numbers_lodo[filter_col]
rename_dict=dict(zip(tabS11_genus.columns,col_names))
tabS11_genus=tabS11_genus.rename(columns=rename_dict)

filter_col = [col for col in numbers_lodo if col.startswith('Family')]
tabS11_family=numbers_lodo[filter_col]
rename_dict=dict(zip(tabS11_family.columns,col_names))
tabS11_family=tabS11_species.rename(columns=rename_dict)

filter_col = [col for col in numbers_lodo if col.startswith('Order')]
tabS11_order=numbers_lodo[filter_col]
rename_dict=dict(zip(tabS11_species.columns,col_names))
tabS11_order=tabS11_species.rename(columns=rename_dict)
tabS11_species.to_excel('../results/fig_and_tables/Table_S12_species.xlsx')
tabS11_genus.to_excel('../results/fig_and_tables/Table_S12_genus.xlsx')
tabS11_family.to_excel('../results/fig_and_tables/Table_S12_family.xlsx')
tabS11_order.to_excel('../results/fig_and_tables/Table_S12_order.xlsx')

#tab_s12
dataset = ['JieZ_2017.txt','ChngKR_2016.txt','YeZ_2018.txt','RaymondF_2016.txt','QinN_2014.txt','FengQ_2015.txt','GuptaA_2019.txt','HanniganGD_2017.txt','ThomasAM_2018a.txt','ThomasAM_2018b.txt','VogtmannE_2016.txt','WirbelJ_2018.txt','YachidaS_2019.txt','YuJ_2015.txt','ZellerG_2014.txt','LiJ_2017.txt','IjazUZ_2017.txt','NielsenHB_2014.txt','GhensiP_2019_m.txt','GhensiP_2019.txt','Castro-NallarE_2015.txt','Heitz-BuschartA_2016.txt','KosticAD_2015.txt','KarlssonFH_2013.txt','QinJ_2012.txt','HMP_2012.txt']
task = ["study_condition","study_condition","study_condition","study_condition","study_condition","study_condition","study_condition","study_condition","study_condition","study_condition","study_condition","study_condition","study_condition","study_condition","study_condition","study_condition","study_condition","study_condition","study_condition","study_condition","study_condition","study_condition","study_condition","study_condition","study_condition",'body_site']
study_condition_controls = ["control","control","control","control","control","control","control","control","control","control","control","control","control","control","control","control","control","control","control","control","control","control","control","control","control",'stool']
study_condition_cases = ['ACVD','AD','BD','cephalosporins','cirrhosis','CRC','CRC','CRC','CRC','CRC','CRC','CRC','CRC','CRC','CRC','hypertension','IBD','IBD','mucositis','peri-implantitis','schizophrenia','T1D','T1D','T2D','T2D','oralcavity']
metrics=['auc','f1','precision','recall']
lvl=['abundance','bool0']
lvl2=['Abundance','Presence/absence']
learner = ['Random_forest','enet','lasso','lsvm','svm']
learner2 = ['Random forest','enet','lasso','lsvm','svm']
lista=[]
for i in range(len(dataset)):
    for m in range(len(lvl)):
        for j in range(len(metrics)):
            for k in range(len(learner)):
                dati = pd.read_csv('../metrics/shotgun/'+learner[k]+'/'+lvl[m]+'/species/'+metrics[j]+'.txt',sep='\t',index_col=0,header=None)
                lista.append((dataset[i].strip('.txt')+' ('+study_condition_cases[i]+')',metrics[j],lvl2[m],learner2[k],dati.loc[dataset[i].strip('.txt')][2]))
lista = pd.DataFrame(lista, columns = ['Dataset name','metrics','level','learner','value'])
tab_S12 = pd.pivot_table(lista,values='value',index='Dataset name',columns=['learner','metrics','level'])
tab_S12.to_excel('../results/fig_and_tables/Table_S13.xlsx')
#Fig6

y= list(numbers_lodo.index)
x= list(numbers_lodo['Species Abundance'])
x2= list(numbers_lodo['Species Presence/Absence'])
x3=list(numbers_lodo['Species Thr = 0.1'])
x4=list(numbers_lodo['Species Thr = 0.01'])
x5=list(numbers_lodo['Species Thr = 0.001'])
x6=list(numbers_lodo['Species Thr = 0.0001'])


fig, (ax1, ax2,ax3)= plt.subplots(nrows=1, ncols=3,  sharey=True,  figsize=(20,20))
ax1.plot(x2, y, 'D',markersize=22, label= 'Threshold = 0.0%')
ax1.plot(x6, y, 'ro',markersize=20, label= 'Threshold = 0.0001%')
ax1.plot(x5, y, 'yd',markersize=18, label= 'Threshold = 0.001%')
ax1.plot(x4, y, 'bs',markersize=16, label= 'Threshold = 0.01%')
ax1.plot(x3, y, 'kX',markersize=14, label= 'Threshold = 0.1%')
ax1.plot(x,y, '^', color='violet', markersize=12, label = 'Relative abundance')

ax1.set_ylabel('Dataset',fontsize=20)
ax1.set_xlabel('AUC',fontsize=20)
ax1.legend(loc='best', prop={'size': 15}, labelspacing=1)
ax1.set_facecolor('white')
ax1.tick_params(axis='y', labelsize=20)
ax1.tick_params(axis='x', labelsize=20)
ax1.set_xlim(0.40, 1)
ax1.xaxis.grid(True, color='grey', linestyle='-', alpha=0.5)
ax1.set_title('A',fontsize=35, fontweight='bold')
ax1.patch.set_linewidth(5)
ax1.patch.set_edgecolor('grey')


x= list(numbers_lodo['Genus Abundance']-numbers_lodo['Species Abundance'])
x2=list(numbers_lodo['Family Abundance']-numbers_lodo['Species Abundance'])
x3=list(numbers_lodo['Order Abundance']-numbers_lodo['Species Abundance'])


ax2.plot(x, y, 'D',markersize=22, label= 'Genus - Species')
ax2.plot(x2, y, 'ro',markersize=20, label= 'Family - Species')
ax2.plot(x3, y, 'yd',markersize=18, label= 'Order - Species')

ax2.set_xlabel('AUC difference:\n AUC Taxonomic level - AUC Species',fontsize=20)
ax2.legend(loc='best', prop={'size': 15}, labelspacing=1)
ax2.set_facecolor('white')
ax2.tick_params(axis='y', labelsize=20)
ax2.tick_params(axis='x', labelsize=20)
ax2.set_xlim(-0.5, 0.5)
ax2.xaxis.grid(True, color='grey', linestyle='-', alpha=0.5)
ax2.set_title('B',fontsize=35, fontweight='bold')
ax2.patch.set_linewidth(5)
ax2.patch.set_edgecolor('grey')

x= list(numbers_lodo['Species Presence/Absence']-numbers_lodo['Species Abundance'])
x2=list(numbers_lodo['Genus Presence/Absence']-numbers_lodo['Genus Abundance'])
x3=list(numbers_lodo['Family Presence/Absence']-numbers_lodo['Family Abundance'])
x4=list(numbers_lodo['Order Presence/Absence']-numbers_lodo['Order Abundance'])

ax3.plot(x, y, 'D',markersize=22, label= 'Species')
ax3.plot(x2, y, 'ro',markersize=20, label= 'Genus')
ax3.plot(x3, y, 'yd',markersize=18, label= 'Family')
ax3.plot(x4, y, 'bs',markersize=16, label= 'Order')
ax3.set_xlabel('AUC difference:\n Presence/absence - Rel. abundance profiles',fontsize=20)
ax3.legend(loc='best', prop={'size': 15}, labelspacing=1)
ax3.set_facecolor('white')
ax3.tick_params(axis='y', labelsize=20)
ax3.tick_params(axis='x', labelsize=20)
ax3.set_xlim(-0.5, 0.5)
ax3.xaxis.grid(True, color='grey', linestyle='-', alpha=0.5)
ax3.set_title('C',fontsize=35, fontweight='bold')
ax3.patch.set_linewidth(5)
ax3.patch.set_edgecolor('grey')


plt.tight_layout()
plt.savefig('../results/fig_and_tables/Fig6.svg')
plt.savefig('../results/fig_and_tables/Fig6.png')

#figS7

a = list(numbers_lodo.index)
aa = numbers.index.isin(a)
numbers=numbers[aa]
numbers_lodo.rename(columns={'Dataset Name':'Dataset_ID'}, inplace=True)
numbers_lodo3 = pd.merge(numbers_lodo,numbers,left_index=True, right_index=True, how='left')
y= list(numbers_lodo3.index)
x= list(numbers_lodo3['AUC_abundance'])
x2= list(numbers_lodo3['Species Abundance'])
x3=list(numbers_lodo3['AUC_presence/absence'])
x4=list(numbers_lodo3['Species Presence/Absence'])


fig, (ax1)= plt.subplots(nrows=1, ncols=1,  figsize=(15,20))
ax1.plot(x, y, 'D',markersize=22, label= 'CV-Abundance')
ax1.plot(x2, y, 'ro',markersize=20, label= 'Lodo-Abundance')
ax1.plot(x3, y, 'yd',markersize=18, label= 'CV-Presence/Absence')
ax1.plot(x4, y, 'bs',markersize=16, label= 'Lodo-Presence/Absence')

ax1.set_ylabel('Dataset',fontsize=20)
ax1.set_xlabel('AUC',fontsize=20)
ax1.legend(loc='best', prop={'size': 15},labelspacing=1)
ax1.set_facecolor('white')
ax1.tick_params(axis='y', labelsize=20)
ax1.tick_params(axis='x', labelsize=20)
ax1.set_xlim(0.40, 1)
ax1.xaxis.grid(True, color='grey', linestyle='-', alpha=0.5)
ax1.patch.set_linewidth(5)
ax1.patch.set_edgecolor('grey')




plt.tight_layout()
plt.savefig('../results/fig_and_tables/FigS9.svg')
plt.savefig('../results/fig_and_tables/FigS9.png')


#figS8

y= list(numbers_class.drop('HMP_2012 (oralcavity)').index)
x= list(numbers_class.drop('HMP_2012 (oralcavity)')['ABUNDANCE Lasso']-numbers_class.drop('HMP_2012 (oralcavity)')['ABUNDANCE RF'])
x2=list(numbers_class.drop('HMP_2012 (oralcavity)')['ABUNDANCE ENet']-numbers_class.drop('HMP_2012 (oralcavity)')['ABUNDANCE RF'])
x3=list(numbers_class.drop('HMP_2012 (oralcavity)')['ABUNDANCE LSVM']-numbers_class.drop('HMP_2012 (oralcavity)')['ABUNDANCE RF'])
x4=list(numbers_class.drop('HMP_2012 (oralcavity)')['ABUNDANCE SVM']-numbers_class.drop('HMP_2012 (oralcavity)')['ABUNDANCE RF'])


fig, (ax1, ax2)= plt.subplots(nrows=1, ncols=2,  sharey=True,  figsize=(20,20))
ax1.plot(x, y, 'D',markersize=20, label= 'Lasso')
ax1.plot(x2, y, 'ro',markersize=18, label= 'ENet')
ax1.plot(x3, y, 'yd',markersize=16, label= 'LSVM')
ax1.plot(x4, y, 'bs',markersize=14, label= 'SVM')
ax1.set_ylabel('Dataset',fontsize=20)
ax1.set_xlabel('AUC difference:\n Selected classifier - RF classifier',fontsize=20)
ax1.legend(loc='best', prop={'size': 15}, labelspacing=1)
ax1.set_facecolor('white')
ax1.tick_params(axis='y', labelsize=20)
ax1.tick_params(axis='x', labelsize=20)
ax1.set_xlim(-0.5, 0.5)
ax1.xaxis.grid(True, color='grey', linestyle='-', alpha=0.5)
ax1.set_title('A',fontsize=35, fontweight='bold')
ax1.patch.set_linewidth(5)
ax1.patch.set_edgecolor('grey')

x= list(numbers_class.drop('HMP_2012 (oralcavity)')['PRESENCE/ABSENCE Lasso']-numbers_class.drop('HMP_2012 (oralcavity)')['PRESENCE/ABSENCE RF'])
x2=list(numbers_class.drop('HMP_2012 (oralcavity)')['PRESENCE/ABSENCE ENet']-numbers_class.drop('HMP_2012 (oralcavity)')['PRESENCE/ABSENCE RF'])
x3=list(numbers_class.drop('HMP_2012 (oralcavity)')['PRESENCE/ABSENCE LSVM']-numbers_class.drop('HMP_2012 (oralcavity)')['PRESENCE/ABSENCE RF'])
x4=list(numbers_class.drop('HMP_2012 (oralcavity)')['PRESENCE/ABSENCE SVM']-numbers_class.drop('HMP_2012 (oralcavity)')['PRESENCE/ABSENCE RF'])

ax2.plot(x, y, 'D',markersize=20, label= 'Lasso')
ax2.plot(x2, y, 'ro',markersize=18, label= 'ENet')
ax2.plot(x3, y, 'yd',markersize=16, label= 'LSVM')
ax2.plot(x4, y, 'bs',markersize=14, label= 'SVM')
ax2.set_xlabel('AUC difference:\n Selected classifier - RF classifier',fontsize=20)
ax2.legend(loc='best', prop={'size': 15}, labelspacing=1)
ax2.set_facecolor('white')
ax2.tick_params(axis='y', labelsize=20)
ax2.tick_params(axis='x', labelsize=20)
ax2.set_xlim(-0.5, 0.5)
ax2.xaxis.grid(True, color='grey', linestyle='-', alpha=0.5)
ax2.set_title('B',fontsize=35, fontweight='bold')
ax2.patch.set_linewidth(5)
ax2.patch.set_edgecolor('grey')
plt.tight_layout()
plt.savefig('../results/fig_and_tables/FigS10.svg')
plt.savefig('../results/fig_and_tables/FigS10.png')

plt.close()

dbs = ["ChngKR_2016","FengQ_2015","GhensiP_2019","GhensiP_2019_m","GuptaA_2019","HanniganGD_2017","Heitz-BuschartA_2016","IjazUZ_2017", "JieZ_2017", "KarlssonFH_2013", "KosticAD_2015",  "LiJ_2017", "NielsenHB_2014", "QinJ_2012", "QinN_2014", "RaymondF_2016", "ThomasAM_2018a", "ThomasAM_2018b", "VogtmannE_2016", "WirbelJ_2018", "YachidaS_2019", "YeZ_2018", "YuJ_2015", "ZellerG_2014"]
study_condition_cases = ["AD", "CRC", "peri-implantitis", "mucositis", "CRC", "CRC", "T1D", "IBD", "ACVD", "T2D", "T1D", "hypertension", "IBD", "T2D", "cirrhosis", "cephalosporins", "CRC", "CRC", "CRC", "CRC", "CRC", "BD", "CRC", "CRC"]
level = ['abundance','bool0']
res = {}
for key in dbs:
    for value in study_condition_cases:
        res[key] = value
        study_condition_cases.remove(value)
        break  
study_condition_cases = ["AD", "CRC", "peri-implantitis", "mucositis", "CRC", "CRC", "T1D", "IBD", "ACVD", "T2D", "T1D", "hypertension", "IBD", "T2D", "cirrhosis", "cephalosporins", "CRC", "CRC", "CRC", "CRC", "CRC", "BD", "CRC", "CRC"]

disease= pd.DataFrame.from_dict(res,orient='index')
disease=disease.rename(columns={0:'disease'})

df = pd.read_csv('../results/biomarkers/shotgun/abundance/species/Castro-NallarE_2015.csv')
df2 = pd.read_csv('../results/biomarkers/shotgun/bool0/species/Castro-NallarE_2015.csv')
dm1 = pd.read_csv('../data/shotgun/abundance/species/Castro-NallarE_2015.txt', sep= '\t',index_col=0).T
dm2 = pd.read_csv('../data/shotgun/bool0/species/Castro-NallarE_2015.txt', sep= '\t',index_col=0).T
dm1 = dm1.loc[:, dm1.columns != 'study_condition'].astype(float).join(dm1['study_condition'])
dm2 = dm2.loc[:, dm2.columns != 'study_condition'].astype(float).join(dm2['study_condition'])
dm1 = dm1.groupby('study_condition').mean().T.reset_index()
dm2 = dm2.groupby('study_condition').mean().T.reset_index()
dm1.rename(columns={ dm1.columns[2]: "case" }, inplace = True)
dm2.rename(columns={ dm2.columns[2]: "case" }, inplace = True)
dm1=dm1.set_index('index')
dm2=dm2.set_index('index')
df=df.set_index('0')
df2=df2.set_index('0')

for i in list(df.index):
    if dm1.control.loc[i] < dm1.case.loc[i]:
        df['1'].loc[i] = df['1'].loc[i] * -1
for i in list(df2.index):
    if dm2.control.loc[i] < dm2.case.loc[i]:
        df2['1'].loc[i] = df2['1'].loc[i] * -1

df.rename(columns={'0':'species','1':'Castro-NallarE_2015_abundance'}, inplace=True)
df2.rename(columns={'0':'species','1':'Castro-NallarE_2015_presence/absence'}, inplace=True)

pval = pd.merge(df ,df2,how='outer',left_index=True,right_index=True)

for i in range(len(dbs)):
    for j in range(len(level)):
        df = pd.read_csv('../results/biomarkers/shotgun/'+level[j]+'/species/'+dbs[i]+'.csv')
        dm1 = pd.read_csv('../data/shotgun/'+lvl[j]+'/species/'+dbs[i]+'.txt', sep= '\t',index_col=0).T
        dm1 = dm1.loc[:, dm1.columns != 'study_condition'].astype(float).join(dm1['study_condition'])
        dm1 = dm1.groupby('study_condition').mean().T.reset_index()
        dm1.rename(columns={study_condition_cases[i]: "case" }, inplace = True)
        df=df.set_index('0')
        dm1=dm1.set_index('index')
        for a in list(df.index):
            if dm1.control.loc[a] < dm1.case.loc[a]:
                df['1'].loc[a] = df['1'].loc[a] * -1
        df.rename(columns={'1':dbs[i]+'_'+level[j]}, inplace=True)
        pval = pd.merge(pval, df , how='outer',left_index=True,right_index=True)


pval=pval.reset_index()
pval.rename(columns={'0':'species'}, inplace=True)
for column in pval.loc[:, pval.columns != 'species']:
    for j in range(len(pval)):
        if (pval[column].iloc[j] <= (-0.05) or pval[column].iloc[j] >= 0.05):
            pval[column].iloc[j] = None
pval = pval.set_index('species')

pval = pval.dropna(how='all')
pval = pval.fillna(1)
pval = pval.loc[:, ~pval.eq(1).all()]
pval.columns=[i.replace('_bool0','').replace('_abundance','') for i in pval.columns]
pval = pval.reset_index()
pval = pval.drop(columns=['Castro-NallarE_2015', "YeZ_2018"] )
meta1 = ['abundance','presence/absence']*18
meta2 = list(pval.loc[:, pval.columns != 'species'].columns)
meta3= list(pval.species)

pval["species"]=[i.split(".")[-1].replace('s__','').replace('_',' ') for i in pval.species]
pval.loc[pval.species == ' massiliensis', 'species'] = 'Collinsella massiliensis'
meta3 = pd.DataFrame(data=meta3)
meta3['class']=meta3[0]
species = list(meta3[0])
meta3["class"]=[i.split(".")[-5].replace('c__','').replace('_',' ') for i in meta3['class']]
meta3[0]=[i.split(".")[-1].replace('s__','').replace('_',' ') for i in meta3[0]]
meta3['class'].loc[51]=meta3['class'].loc[52]
meta3[0].iloc[51] = pval['species'].iloc[51]
meta = pd.DataFrame(meta1,meta2).reset_index().T
meta.columns=[pval.loc[:, pval.columns != 'species'].columns]
meta=meta.T
meta = meta.rename(columns={"index": "dataset",0:'presence/absence vs abundance'})
for i in meta.index:
    meta.loc[i,'disease'] = disease['disease'].loc[i]

N = 18
HSV_tuples = distinctipy.get_colors(N)
RGB_tuples = map(lambda x: colorsys.hsv_to_rgb(*x), HSV_tuples)
X = 26
HSV_tuples2 = distinctipy.get_colors(X)
RGB_tuples2 = map(lambda x: colorsys.hsv_to_rgb(*x), HSV_tuples2)
network_pal = sns.color_palette(HSV_tuples)
network_lut = dict(zip(meta.dataset.unique(), network_pal))
networks = meta.dataset
network_colors = pd.Series(networks).map(network_lut)
network_pal2 = sns.color_palette("Paired", len(meta['presence/absence vs abundance'].unique()))
network_lut2 = dict(zip(meta['presence/absence vs abundance'].unique(), network_pal2))
networks2 = meta['presence/absence vs abundance']
network_colors2 = pd.Series(networks2).map(network_lut2)
meta3 = meta3.sort_values(['class'])
network_pal3 = sns.color_palette(HSV_tuples2)
network_lut3 = dict(zip(meta3['class'].unique(), network_pal3))
networks3 = meta3['class']
network_colors3 = pd.Series(networks3).map(network_lut3)
network_colors3.index = meta3[0]
net_col ={'dataset':network_colors,'presence/absence vs abundance': network_colors2}
net_col = pd.DataFrame(data=net_col)
net_col = net_col.reset_index().drop(columns='level_0')

color1 = plt.cm.Reds(np.linspace(0., 1, 128))
color2 = plt.cm.Greens_r(np.linspace(0, 1, 128))
colors = np.vstack((color1, color2))
mymap = mcolors.LinearSegmentedColormap.from_list('my_colormap', colors)

pval2 = pval
pval2 = pval2.replace(1,np.nan)

pval2 = pval2.dropna(how='all').reset_index()
species2 = pval2.species
values=[]
pval2=pval2.set_index('species')
for i in list(pval2.index):
    values.append((i,sum(n < 0 for n in pval2.loc[i]),sum(n > 0 for n in pval2.loc[i])))
values= pd.DataFrame(data=values,columns=['species','N.Negativi','N.Positivi'])
percent=[]
for i in range(len(values)):
    percent.append((values.species.iloc[i],(values['N.Positivi'].iloc[i] / (values['N.Positivi'].iloc[i] + values['N.Negativi'].iloc[i]))*100,(values['N.Negativi'].iloc[i] / (values['N.Positivi'].iloc[i] + values['N.Negativi'].iloc[i]))*100,(values['N.Positivi'].iloc[i] / 32)*100,(values['N.Negativi'].iloc[i] / 32)*100))
    


pval3 = pval
pval3 = pval3.replace(1,np.nan)
pval3=pval3.set_index('species')
pval3[pval3<0] =0 
pval3[pval3>0] =100
pval3 = pval3.replace(np.nan,1)
species3 = list(pval3.index)
lista =pd.DataFrame(data=species3)
lista['diff']=0
ab=[0,2,4,6,8,10,12,14,16,18,20,22,24,26,28,30,32,34]
pa=[1,3,5,7,9,11,13,15,17,19,21,23,25,27,29,31,33,35]
new_header = pval3.columns[ab]
discord = pd.DataFrame(columns=new_header,index=pval3.index)


for i in range(len(pval3)):
    for j in range(18):
        if pval3.iloc[i,ab[j]]!=pval3.iloc[i,pa[j]]:
            lista['diff'].iloc[i]=lista['diff'].iloc[i]+1
        if (pval3.iloc[i,ab[j]] + pval3.iloc[i,pa[j]] == 100):
            discord.iloc[i,j] = 1
        else:
            discord.iloc[i,j] = 0
lista['diff']=(lista['diff']/18)*100
percent = pd.DataFrame(data=percent,columns=['species','percentuale positivi/significativi','percentuale negativi/significativi','percentuale positivi/totale','percentuale negativi/totale'])
lista.rename(columns={ 0: "species" }, inplace = True)
pval3 = pval3.T
pval3['disease'] = pval3.index.map(res)
counter_pos = pd.DataFrame(index=pval3.columns,columns=pval3['disease'].unique())
counter_pos = counter_pos.fillna(0)
for dis in pval3['disease'].unique():
    pcount = pval3[pval3['disease']==dis]
    pcount = pcount.T
    for j in pcount.index:
        if (pcount.loc[j,:]==100).any():
            counter_pos.loc[j,dis] = counter_pos.loc[j,dis] + 1
counter_neg = pd.DataFrame(index=pval3.columns,columns=pval3['disease'].unique())
counter_neg = counter_neg.fillna(0)
for dis in pval3['disease'].unique():
    pcount = pval3[pval3['disease']==dis]
    pcount = pcount.T
    for j in pcount.index:
        if (pcount.loc[j,:]==0).any():
            counter_neg.loc[j,dis] = counter_neg.loc[j,dis] + 1
pval3 = pval3.T
counter_pos.loc[:, 'somma'] = counter_pos.iloc[:, :].sum(axis=1)
counter_neg.loc[:, 'somma'] = counter_neg.iloc[:, :].sum(axis=1)

network_pal4 = sns.color_palette("Greens",11)
network_lut4 = dict(zip(sorted(percent['percentuale positivi/totale'].unique()), network_pal4))
networks4 = percent['percentuale positivi/totale']
networks4.index = percent.species

network_colors4 = pd.Series(networks4).map(network_lut4)
pval3= pval3.reset_index()
pval3 = pval3[pval3.index.isin(species2)]
network_colors3=network_colors3.to_frame()
net_col2=network_colors3.join(network_colors4)

network_pal5 = sns.color_palette("Reds",15)
network_lut5 = dict(zip(sorted(percent['percentuale negativi/totale'].unique()), network_pal5))
networks5 = percent['percentuale negativi/totale']
networks5.index = percent.species
network_colors5 = pd.Series(networks5).map(network_lut5)
net_col2=net_col2.join(network_colors5)
net_col2=net_col2.stack().unstack(fill_value=(1.0,1.0,1.0))
pval=pval.set_index('species')
net_col2.rename(columns={ 'percentuale positivi/totale': "Control enriched",'percentuale negativi/totale': 'Case enriched'  }, inplace = True)
new_rc_params = {'text.usetex': False,
"svg.fonttype": 'none'
}
mpl.rcParams.update(new_rc_params)


disease_percent=pd.DataFrame(index=counter_pos.index)
for i in disease_percent.index:
    disease_percent.loc[i,'disease_percent_pos'] = (counter_pos.somma.loc[i]/9)
    disease_percent.loc[i,'disease_percent_neg'] = (counter_neg.somma.loc[i]/9)
network_pal4 = sns.color_palette("Greens",7)
network_pal5 = sns.color_palette("Reds",7)

network_lut6 = dict(zip(sorted(disease_percent['disease_percent_pos'].unique()), network_pal4))
networks6= disease_percent['disease_percent_pos']
networks6.index = disease_percent.index
network_colors6 = pd.Series(networks6).map(network_lut6)
network_lut7 = dict(zip(sorted(disease_percent['disease_percent_neg'].unique()), network_pal5))
networks7= disease_percent['disease_percent_neg']
networks7.index = disease_percent.index
network_colors7 = pd.Series(networks7).map(network_lut7)
net_col3=net_col2['class'].to_frame()

net_col3=net_col3.join(network_colors6)
net_col3=net_col3.join(network_colors7)
net_col3=net_col3.rename(columns={'disease_percent_pos':'Percentage of diseases enriched in controls','disease_percent_neg':'Percentage of diseases enriched in cases'})
g =sns.clustermap(pval, col_colors=[net_col['dataset'],net_col['presence/absence vs abundance']] ,colors_ratio=(0.03,0.03),row_colors=net_col3 , cmap=mymap, figsize=(30 , 30), mask=(pval==1), vmin= -0.05, vmax=0.05 ,yticklabels=0,linewidth=0.01,linecolor="grey",square=True)
handles = [Patch(facecolor=network_lut3[name]) for name in network_lut3]
plt.legend(handles, network_lut3, title='class',
           bbox_to_anchor=(1, 1), bbox_transform=plt.gcf().transFigure, loc='upper right',fontsize=8)

plt.savefig('../results/fig_and_tables/FigS7.svg', format='svg',dpi=300 )
plt.savefig('../results/fig_and_tables/FigS7.pdf')

color1 = (1, 0, 0,1) 
color2 =(0, 1, 0,1) 
colors = np.vstack((color2, color1))
mymap = mcolors.LinearSegmentedColormap.from_list('my_colormap', colors,N=2)
discord = discord.astype(float)

discord=discord.loc[~(discord==0).all(axis=1)]
discord=discord.loc[:,~(discord==0).all(axis=0)]

g =sns.clustermap(discord ,colors_ratio=(0.03,0.03) , cmap=mymap, figsize=(30 , 30), vmin= 0, vmax=1 ,yticklabels=1,xticklabels=1,linewidth=0.05,linecolor="grey",square=True,col_cluster=False,row_cluster=False)
plt.savefig('../results/fig_and_tables/FigS8.svg', format='svg',dpi=300 )
plt.savefig('../results/fig_and_tables/FigS8.pdf')

number_of_points=pd.DataFrame()
db = ["FengQ_2015","GhensiP_2019","GuptaA_2019","IjazUZ_2017", "JieZ_2017", "NielsenHB_2014", "QinJ_2012", "QinN_2014", "RaymondF_2016", "WirbelJ_2018", "YachidaS_2019"]
lvl = ['abundance','bool0']
minimi = []
def pvalue_s(df):
    pvals=[]
    for i in range(len(df.columns)):
        tmp=[]
        for j in range(len(df.columns)):
            t= df[[df.columns[i],df.columns[j]]].dropna()
            tmp.append(spearmanr(t.iloc[:,0],t.iloc[:,1])[1])
        pvals.append(tmp)
    return pd.DataFrame([i for i in pvals], index=df.columns,columns=df.columns)
fig, axs = plt.subplots(4, 3,  figsize=(30,40))
for i in range(len(db)):
    abund = pd.read_csv('../results/biomarkers/shotgun/'+lvl[0]+'/species/'+db[i]+'.csv',index_col=0)
    prese = pd.read_csv('../results/biomarkers/shotgun/'+lvl[1]+'/species/'+db[i]+'.csv',index_col=0)
    tot = pd.merge(abund,prese,left_index=True,right_index=True)
    tot = tot.T
    tot = tot[tot.columns[(tot<=0.05).any()]]
    tot = tot.T
    corr = tot.corr(method = 'spearman')
    pvalue = pvalue_s(tot)

    axs = axs.ravel()
    axs[i].scatter(tot['1_x'],tot['1_y'])
    
    axs[i].plot( [0,1],[0,1],"--r" )
    axs[i].title.set_text('Dataset:'+db[i]+'\n (correlation = ' + str((round(corr.iloc[0,1],2))) +', p-value = '+"{:.2e}".format(pvalue.iloc[0,1])+', number of points = ' + str(len(tot))+')')
    axs[i].set_xlabel('P-value from Relative abundance', fontsize=13)
    axs[i].set_ylabel('P-value from Presence/Absence', fontsize=13)
    #axs[i].set_xlim(0,1)
    #axs[i].set_ylim(0,1)
    axs[i].set_xscale('log')
    axs[i].set_yscale('log')

    number_of_points.loc[db[i],'number_of_points']=len(tot)
    number_of_points.loc[db[i],'correlation']=corr.iloc[0,1]
plt.savefig('../results/fig_and_tables/FigS6.png')
plt.savefig('../results/fig_and_tables/FigS6.svg')
y= list(numbers_boolean.drop('HMP_2012 (oralcavity)').index)
x= list(numbers_boolean.iloc[:,:6].drop('HMP_2012 (oralcavity)')['Abundance'])
x1 = list(numbers_boolean.iloc[:,:6].drop('HMP_2012 (oralcavity)')['Presence/Absence'])
x2=list(numbers_boolean.iloc[:,:6].drop('HMP_2012 (oralcavity)')['Threshold = 0.0001'])
x3=list(numbers_boolean.iloc[:,:6].drop('HMP_2012 (oralcavity)')['Threshold = 0.001'])
x4=list(numbers_boolean.iloc[:,:6].drop('HMP_2012 (oralcavity)')['Threshold = 0.01'])
x5=list(numbers_boolean.iloc[:,:6].drop('HMP_2012 (oralcavity)')['Threshold = 0.1'])


fig = plt.figure(figsize=(50, 20),constrained_layout=True)
gs0 = fig.add_gridspec(3, 1)

spec1 = gs0[0].subgridspec(1, 5 )
spec2 = gs0[1].subgridspec(1, 4 )
spec3 = gs0[2].subgridspec(1, 5 )
ax1 = fig.add_subplot(spec1[0,0])
ax2 = fig.add_subplot(spec1[0,1])
ax3 = fig.add_subplot(spec1[0,2])
ax4 = fig.add_subplot(spec1[0,3])
ax5 = fig.add_subplot(spec1[0,4])
ax6 = fig.add_subplot(spec2[0,0])
ax7 = fig.add_subplot(spec2[0,1])
ax8 = fig.add_subplot(spec2[0,2])
ax9 = fig.add_subplot(spec2[0,3])
ax10 = fig.add_subplot(spec3[0,0])
ax11 = fig.add_subplot(spec3[0,1])
ax12 = fig.add_subplot(spec3[0,2])
ax13 = fig.add_subplot(spec3[0,3])
ax14 = fig.add_subplot(spec3[0,4])
                           
ax1.scatter(x,x1)
ax1.plot( [0,1],[0,1],"--r" )
ax1.title.set_text('AUC comparison Abundance vs Presence/Absence')
ax1.set_xlabel('AUC relative abundance')
ax1.set_ylabel('AUC presence/absence')
ax1.set_xlim(0,1)
ax1.set_ylim(0,1)
ax2.scatter(x,x2)
ax2.plot( [0,1],[0,1],"--r" )
ax2.title.set_text('AUC comparison Abundance vs Threshold = 0.0001')
ax2.set_xlabel('AUC relative abundance')
ax2.set_ylabel('AUC Threshold = 0.0001')
ax2.set_xlim(0,1)
ax2.set_ylim(0,1)
ax3.scatter(x,x3)
ax3.plot( [0,1],[0,1],"--r" )
ax3.title.set_text('AUC comparison Abundance vs Threshold = 0.001')
ax3.set_xlabel('AUC relative abundance')
ax3.set_ylabel('AUC Threshold = 0.001')
ax3.set_xlim(0,1)
ax3.set_ylim(0,1)
ax4.scatter(x,x4)
ax4.plot( [0,1],[0,1],"--r" )
ax4.title.set_text('AUC comparison Abundance vs Threshold = 0.01')
ax4.set_xlabel('AUC relative abundance')
ax4.set_ylabel('AUC Threshold = 0.01')
ax4.set_xlim(0,1)
ax4.set_ylim(0,1)
ax5.scatter(x,x5)
ax5.plot( [0,1],[0,1],"--r" )
ax5.title.set_text('AUC comparison Abundance vs Threshold = 0.1')
ax5.set_xlabel('AUC relative abundance')
ax5.set_ylabel('AUC Threshold = 0.1')
ax5.set_xlim(0,1)
ax5.set_ylim(0,1)
y= list(numbers_tass.drop('HMP_2012 (oralcavity)').index)
x= list(numbers_tass.drop('HMP_2012 (oralcavity)').iloc[:,:8]['Species Presence/Absence'])
x1= list(numbers_tass.drop('HMP_2012 (oralcavity)').iloc[:,:8]['Species Abundance'])
x2=list(numbers_tass.drop('HMP_2012 (oralcavity)').iloc[:,:8]['Genus Presence/Absence'])
x_2=list(numbers_tass.drop('HMP_2012 (oralcavity)').iloc[:,:8]['Genus Abundance'])
x3=list(numbers_tass.drop('HMP_2012 (oralcavity)').iloc[:,:8]['Family Presence/Absence'])
x_3=list(numbers_tass.drop('HMP_2012 (oralcavity)').iloc[:,:8]['Family Abundance'])
x4=list(numbers_tass.drop('HMP_2012 (oralcavity)').iloc[:,:8]['Order Presence/Absence'])
x_4=list(numbers_tass.drop('HMP_2012 (oralcavity)').iloc[:,:8]['Order Abundance'])
ax6.scatter(x1,x)
ax6.plot( [0,1],[0,1],"--r" )
ax6.title.set_text('AUC comparison Abundance vs Presence/Absence species level')
ax6.set_xlabel('AUC species level relative abundance')
ax6.set_ylabel('AUC species level presence/absence')
ax6.set_xlim(0,1)
ax6.set_ylim(0,1)
ax7.scatter(x_2,x2)
ax7.plot( [0,1],[0,1],"--r" )
ax7.title.set_text('AUC comparison Abundance vs Presence/Absence level genus')
ax7.set_xlabel('AUC genus level relative abundance')
ax7.set_ylabel('AUC genus level presence/absence')
ax7.set_xlim(0,1)
ax7.set_ylim(0,1)
ax8.scatter(x_3,x3)
ax8.plot( [0,1],[0,1],"--r" )
ax8.title.set_text('AUC comparison Abundance vs Presence/Absence level family')
ax8.set_xlabel('AUC family level relative abundance')
ax8.set_ylabel('AUC family level presence/absence')
ax8.set_xlim(0,1)
ax8.set_ylim(0,1)
ax9.scatter(x_4,x4)
ax9.plot( [0,1],[0,1],"--r" )
ax9.title.set_text('AUC comparison Abundance vs Presence/Absence level order')
ax9.set_xlabel('AUC order level relative abundance')
ax9.set_ylabel('AUC order level presence/absence')
ax9.set_xlim(0,1)
ax9.set_ylim(0,1)
x= list(numbers_class.drop('HMP_2012 (oralcavity)')['PRESENCE/ABSENCE Lasso'])
x_1= list(numbers_class.drop('HMP_2012 (oralcavity)')['ABUNDANCE Lasso'])

x2=list(numbers_class.drop('HMP_2012 (oralcavity)')['PRESENCE/ABSENCE ENet'])
x_2=list(numbers_class.drop('HMP_2012 (oralcavity)')['ABUNDANCE ENet'])

x3=list(numbers_class.drop('HMP_2012 (oralcavity)')['PRESENCE/ABSENCE LSVM'])
x_3=list(numbers_class.drop('HMP_2012 (oralcavity)')['ABUNDANCE LSVM'])

x4=list(numbers_class.drop('HMP_2012 (oralcavity)')['PRESENCE/ABSENCE SVM'])
x_4=list(numbers_class.drop('HMP_2012 (oralcavity)')['ABUNDANCE SVM'])

x5=list(numbers_class.drop('HMP_2012 (oralcavity)')['PRESENCE/ABSENCE RF'])
x_5=list(numbers_class.drop('HMP_2012 (oralcavity)')['ABUNDANCE RF'])

ax10.scatter(x_1,x)
ax10.plot( [0,1],[0,1],"--r" )
ax10.title.set_text('AUC comparison Abundance vs Presence/Absence Lasso classifier')
ax10.set_xlabel('AUC relative abundance')
ax10.set_ylabel('AUC presence/absence')
ax10.set_xlim(0,1)
ax10.set_ylim(0,1)
ax11.scatter(x_2,x2)
ax11.plot( [0,1],[0,1],"--r" )
ax11.title.set_text('AUC comparison Abundance vs Presence/Absence ENet classifier')
ax11.set_xlabel('AUC relative abundance')
ax11.set_ylabel('AUC presence/absence')
ax11.set_xlim(0,1)
ax11.set_ylim(0,1)
ax12.scatter(x_3,x3)
ax12.plot( [0,1],[0,1],"--r" )
ax12.title.set_text('AUC comparison Abundance vs Presence/Absence LSVM classifier')
ax12.set_xlabel('AUC relative abundance')
ax12.set_ylabel('AUC presence/absence')
ax12.set_xlim(0,1)
ax12.set_ylim(0,1)
ax13.scatter(x_4,x4)
ax13.plot( [0,1],[0,1],"--r" )
ax13.title.set_text('AUC comparison Abundance vs Presence/Absence SVM classifier')
ax13.set_xlabel('AUC relative abundance')
ax13.set_ylabel('AUC presence/absence')
ax13.set_xlim(0,1)
ax13.set_ylim(0,1)
ax14.scatter(x_5,x5)
ax14.plot( [0,1],[0,1],"--r" )
ax14.title.set_text('AUC comparison Abundance vs Presence/Absence RF classifier')
ax14.set_xlabel('AUC relative abundance')
ax14.set_ylabel('AUC presence/absence')
ax14.set_xlim(0,1)
ax14.set_ylim(0,1)
plt.savefig('../results/fig_and_tables/FigS3.svg')
plt.savefig('../results/fig_and_tables/FigS3.png')

mean_reads=pd.read_csv('../results/fig_and_tables/meanreads.csv',index_col=1).iloc[:,-2:]
corr = mean_reads.corr(method='spearman')
pvalue = pvalue_s(mean_reads)
fig, (ax1)= plt.subplots(nrows=1, ncols=1, figsize=(10,10))
ax1.scatter(mean_reads['1'],mean_reads['2'])
ax1.set_xlim(mean_reads['1'].min()-1000000,mean_reads['1'].max()+1000000)
ax1.set_ylim(mean_reads['2'].min()-10,mean_reads['2'].max()+50)
ax1.title.set_text('Average number of reads vs number of significant species \n correlation = '+"{:.2e}".format(corr.iloc[0,1]))
ax1.set_xlabel('Average number of reads')
ax1.set_ylabel('Number of differentially abundant species')
plt.tight_layout()
plt.savefig('../results/fig_and_tables/FigS5.svg')
plt.savefig('../results/fig_and_tables/FigS5.png')
