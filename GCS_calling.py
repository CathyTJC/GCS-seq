### GCS Calling
# Author: Cathy/Juechun Tang

import os
import copy
import random
import matplotlib.pyplot as plt
import numpy as np
from collections import Counter
import pandas as pd

#######
# Parses txt file with coverage at each genomic position
# Input: format pos\t start\t end
# Extract 1) the 3' end counts into a list 2) total number of counts
#######
def start_end_parsing(file):
    print('Now is processing: ' + str(file))
    NE=[]
    Total_cov=0
    filein=open(file, 'r')
    for line in filein:
        line=line.rstrip().split('\t')
        if line[1] not in ['Pos']:  
            NE.append(int(line[2]))
            Total_cov+=int(line[2])
    filein.close()
    return NE, Total_cov

def read_file(file1, file2, file3,file4):
    Flag = start_end_parsing(file1)
    FlagIC = start_end_parsing(file2)
    Flagless = start_end_parsing(file3)
    FlaglessIC= start_end_parsing(file4)
    return Flag, FlagIC, Flagless, FlaglessIC

#######
# Normalization sequencing depth of input control, input list of [counts, total counts]
# return a list of normalized values
####### 
def depth_normIC(NE, eps = 1e-11):
    tmp = np.array(NE[0])
    print('Normalize IC by sequencing depth: ', NE[1])
    # print('No. of 0s: ', np.sum(np.argwhere(NE[0]==0).astype(int)))
    return np.maximum(tmp,eps)/NE[1], NE[1]
   
#######
# Input: FLAG, FLAG IC, FLAGless, FLAGless IC files
# 1) Normalization to Input control
#######  
def normalization(Flag, FlagIC, Flagless, FlaglessIC, eps = 0.00000000001):
    print('Normalize by sequencing depth...')
    Flag_norm = depth_normIC(Flag, eps)[0]
    Flagless_norm = depth_normIC(Flagless, eps)[0]
    FlagIC_norm = depth_normIC(FlagIC, eps)[0]
    FlaglessIC_norm = depth_normIC(FlaglessIC, eps)[0]

    Flag_normed = Flag_norm/FlagIC_norm
    Flagless_normed = Flagless_norm/FlaglessIC_norm

    return Flag_normed,Flagless_normed

#######
# Get gyrase cleavage site 
####### 
def call_cleavage(seq, alpha = 0.1, beta = 0.0):
    GCS = {}
    input = np.array(seq)
    length = len(input)-5
    high_mean = (input[:length]+input[5:length+5])/2
    rough_mean = (input[1:length+1]+input[2:length+2]+ input[3:length+3]+ input[4:length+4])/4
    high_min = np.minimum(input[:length],input[5:length+5])
    rough_max = input[1:length+1] 
    for i in range(2,5):
        rough_max = np.maximum(input[i:length+i],rough_max)
    positions  = np.argwhere((high_mean - rough_mean >=alpha) & (high_min- rough_max > beta)).squeeze()
    for pos in positions:
        GCS[pos] = (input[pos],input[pos+5], high_mean[pos] - rough_mean[pos], high_min[pos]- rough_max[pos])
    return GCS

#######
# Obtain GCS that occrued in the 3 replicates
####### 
def get_HF_GCS(drug, rep1,rep2,rep3,threshold =3):
    mg_list = list(rep1.keys())+list(rep2.keys())+list(rep3.keys())
    ct_list = dict(Counter(mg_list))
    GCS_final = dict([x for x in ct_list.items()if x[1]>=threshold ]) # select GCS frequency >2
    print(str(drug) + '\t' + str(len(GCS_final))) 
    GCS_dict ={}
    for key, freq in GCS_final.items():
        GCS_dict[key]=[]
        for rep in [rep1,rep2,rep3]:
            if key in rep.keys():
                GCS_dict[key].append([x for x in rep[key]])
            else:
                    GCS_dict[key].append(['NA','NA','NA'])
        GCS_dict[key].append([freq])  
    GCS_dict = dict(sorted(GCS_dict.items()))
    return GCS_dict

#######
# get GCS dictionary for each drug
# Needs normalized read counts Flag_normed and Flagless_normed prepared first
def cleavage_wrap(folder, alpha = 0.1, beta = 0.0, thereshold = 3):
    GCS_control = {}
    GCS = {}
    for drug in folder:
        for i in range(1,4):
            GCS_control[(drug,i)] = call_cleavage(Flagless_normed[(drug,i)],alpha, beta)
            GCS[(drug,i)] = call_cleavage(Flag_normed[(drug,i)],alpha, beta)
    for drug in folder:
        for i in range(1,4):
            print(drug+ str(i)+ ':', len(GCS_control[(drug,i)]), len(GCS[(drug,i)]))        
    HF_GCS = {}
    HF_GCS_control = {}
    for drug in folder:
        HF_GCS[drug] = get_HF_GCS(drug, GCS[(drug,1)], GCS[(drug,2)], GCS[(drug,3)])
        HF_GCS_control[drug] = get_HF_GCS(drug, GCS_control[(drug,1)], GCS_control[(drug,2)], GCS_control[(drug,3)])
    
    return HF_GCS, HF_GCS_control


def write_HF_GCS(dict,path_out):
        for drug, sub_dict in dict.items():
                file_name = drug +'.txt'
                completeName = os.path.join(path_out, file_name)
                GCSs_out=open(completeName, 'w')
                GCSs_out.write('pos' +'\t'+'count_i_rep1'+'\t'+'count_i+5_rep1'+'\t'+'dif_mean_rep1'+'\t'+ 
                'dif_min_max_rep1'+'\t'+'count_i_rep2'+'\t'+'count_i+5_rep2'+'\t'+'dif_mean_rep2'+'\t'+'dif_min_max_rep2'+'\t'+
                'count_i_rep3'+'\t'+'count_i+5_rep3'+'\t'+'dif_mean_rep3'+'\t'+'dif_min_max_rep3'+'\t' +'frequency'+'\n')
                for pos, value in sub_dict.items():
                        GCSs_out.write(str(pos+1) +'\t')
                        for rep in value:
                                for num in rep:
                                        GCSs_out.write(str(num) + '\t')
                        GCSs_out.write('\n')
                GCSs_out.close() 


####################### 
# Initialization (the input and output directory need to be changed according to the location of the coverage file and desired output directory)
pwd="/Volumes/TJC/GCS_calling_files/coverage_drug/"
path_out = '/Volumes/TJC/GCS_calling_files/HF_GCS_0408/'
if not os.path.exists(path_out):
        os.makedirs(path_out) 
folder = ['MOXI','LEVO','NOR','CIP','GEMI'] #drug used
#folder = ['4Mu-LEVO']
idx = [1,2,3] #replicates
eps =[0.00000000001]
 
Flag_normed = {}
Flagless_normed = {}
Flag = {}
FlagIC = {}
Flagless = {}
FlaglessIC = {}

# read in files 
for drug in folder:
    if drug in Flag_normed.keys():
        continue
    for i in idx:
        file1 = pwd + drug +"/Flag" + str(i) + ".txt"
        file2 = pwd+ drug+"/FlagIC"+ str(i)+".txt"
        file3 = pwd+ drug+"/Flagless"+str(i)+".txt"
        file4 = pwd+ drug+"/FlaglessIC"+str(i)+".txt"

        if ((drug,i) not in Flag.keys()) or ((drug,i) not in FlagIC.keys()) or ((drug,i) not in Flagless.keys()) or ((drug,i) not in FlaglessIC.keys()):
            (Flag[(drug,i)], FlagIC[(drug,i)], Flagless[(drug,i)], FlaglessIC[(drug,i)]) = read_file(file1, file2, file3,file4)

# normalization
for drug in folder:
    for i in idx:
        Flag_normed[(drug,i)],Flagless_normed[(drug,i)] = normalization(Flag[(drug,i)], FlagIC[(drug,i)], Flagless[(drug,i)], FlaglessIC[(drug,i)],eps)

#call cleavge
HF_GCS, HF_GCS_control =cleavage_wrap(folder, alpha=0.04,beta = 0.0)

#save file
write_HF_GCS(HF_GCS,path_out)


##############################
## Add in strength info and save the dataframe again
pwd =  '/Volumes/TJC/GCS_calling_files/HF_GCS_0408/' #dir to the GCS folders
path_out = '/Volumes/TJC/GCS_calling_files/HF_GCS_strength_0408/'
if not os.path.exists(path_out):
        os.makedirs(path_out)

drug_list = ['MOXI', 'LEVO', 'CIP', 'GEMI', 'NOR']
#drug_list = ['4Mu-LEVO','LEVO']
GCS_HF_count = {drug : pd.read_csv(pwd + drug + '.txt',sep = '\t',index_col = False) for drug in drug_list}
rep = ['rep1','rep2','rep3']

for drug, df in GCS_HF_count.items():
   for replicate in rep:
      df['avg_count_'+replicate] = (df['count_i_'+replicate]+ df['count_i+5_'+replicate])/2
      df['cleavage_percentage_'+replicate] =df['dif_mean_'+replicate]/df['avg_count_'+replicate]
      df['strength_'+replicate] = df['avg_count_'+replicate]*df['cleavage_percentage_'+replicate]
   df['avg_strength'] = (df['strength_rep1']+df['strength_rep2']+df['strength_rep3'])/3
  
   #save the df again
   df.to_csv(path_out+ drug+ '.txt', index = False)
