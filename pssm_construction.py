# author: Juechun Tang
# 02-24-21 updated
## 10.2022 updated
# construct motifs for MOXI, NOR, CIP, GEMI, and LEVO treated samples

import os
import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq
from Bio import motifs
import pandas as pd

#######
#FASTA sequences parsing.
#######
def obtain_seq(seq_path):
    seq_oi=open(seq_path, 'r')
    for record in SeqIO.parse(seq_oi, "fasta"):
        sequence=str(record.seq)
    seq_oi.close()      
    return sequence, record.id

#######
# GCSs data parsing, returns a dictionary in the format of of {GCS position:strength}
# take dicionary as input, user define name and input_path
#######
def GCSs_parsing(input_dict,col_num):
    GCSs_sets_dict={}
    for k, v in input_dict.items():
        GCSs_dict={}
        filein=open(v, 'r')
        for line in filein:
            line=line.rstrip().split(',')
            if line[0] not in ['pos']:
                GCSs_dict[int(line[0])]=float(line[col_num])
            else:
                continue
        filein.close()
        GCSs_sets_dict[k]=GCSs_dict
        print('Number of trusted GCSs for ' + str(k) + ' : ' + str(len(GCSs_dict)))
    return GCSs_sets_dict

# strength normalization only used if we use a wegihted motif construction
def strength_normalization (GCS_dict):
    GCS_norm_set_dict ={}
    for key, sub_dict in GCS_dict.items():
        GCS_norm_dict ={}
        least_val = min(sub_dict.values())
        print(least_val)
        for pos, strength in sub_dict.items():
            if strength < 1e+10:
                GCS_norm_dict[pos]=round(strength/least_val)
            else:
                continue
        GCS_norm_set_dict[key]=GCS_norm_dict
    return GCS_norm_set_dict
    
    
#######
# Motif construction
#######

# if weight = 1, the motif considers the relative abundance of each GCS
def motif_construct(GCSs_dict, genome, weight):
    win_width_l=5
    win_width_r=13
    seqs=[]
    for k, v in GCSs_dict.items():
        seq = genome[(int(k)-1-win_width_l):(int(k)-1+win_width_r)]
        if weight == 0:
            seqs.append(seq) 
        if weight == 1:    #Returns (weighted if weight == 1) sequences under the GCSs.
            [seqs.append(seq) for i in range(v)]       
    print('Number of sequences for motif construction: ' + str(len(seqs)))
    print('Len of sequences for motif construction: ' + str(len(seqs[0])))
    #Create PWM and PSSM.
    background={'A': 0.2461995, 'C': 0.2542337, 'G': 0.2536619, 'T': 0.2459049} #for MU-ori strain
    instances=[]
    for seq in seqs:
        instances.append(Seq(seq.upper()))
    m = motifs.create(instances)
    pwm=m.counts.normalize(pseudocounts=0.5)
    pssm=pwm.log_odds(background)  
    print('Consensus sequence: ' )
    print(m.consensus)
    return seqs,pwm, m


def wrap_motif_construct(Source_genome_path,GCSs_files_paths, weight,path_out,col_num):
    Source_sequence=obtain_seq(Source_genome_path)[0]
    GCSs_sets_dict = GCSs_parsing(GCSs_files_paths,col_num)
    if weight ==1:
        GCSs_sets_dict = strength_normalization (GCSs_sets_dict)
    for k, v in GCSs_sets_dict.items():
        GCSs_dict = v
        print('Processing samples from: ' + k )
        seqs,pwm, m = motif_construct(GCSs_dict, Source_sequence, weight)
        # make weblogo
        m.weblogo(path_out+ k+"_motif.pdf",format="pdf",show_fineprint = False,logo_title = k)

        # save the seqs for  seq logo if use onlihne web logo
        with open(path_out+  k + "_motif_seq.txt", "w") as output:
                for item in seqs:
                    output.write("%s\n" % item)
        # save the pwm
        df = pd.DataFrame(pwm)
        df.to_csv(path_out+k +'_pwm.csv')
    print('Finish motif construction!')
    return seqs,pwm, m


#Input: GCS input
pwd = '/Volumes/TJC/ChIP_project/GCS_calling_files/HF_GCS_strength_0408/TOP_strength_GCSs/'
GCSs_input={'MOXI': pwd+"TopGCS_MOXI_sort.txt",'LEVO': pwd+"TopGCS_LEVO_sort.txt", 'NOR': pwd+ "TopGCS_NOR_sort.txt",'CIP': pwd+"TopGCS_CIP_sort.txt",
           'GEMI': pwd+"TopGCS_GEMI_sort.txt"}

#Input: path to the mapping genome FASTA.
fasta_path_input = "/Volumes/TJC/ChIP_project/GCS_calling_files/Mu_ori_mu_insert_MG1655.fa"

#Output path
Output_dir="/Volumes/TJC/ChIP_project/GCS_calling_files/motif_topStrength_sort/"
if not os.path.exists(Output_dir):
    os.makedirs(Output_dir)     
seqs,pwm, m = wrap_motif_construct(fasta_path_input,GCSs_input, 1,Output_dir,1)


Output_dir="/Volumes/TJC/ChIP_project/GCS_calling_files/motif_topStrength_sort/noweight/"
if not os.path.exists(Output_dir):
    os.makedirs(Output_dir)     
seqs,pwm, m = wrap_motif_construct(fasta_path_input,GCSs_input, 0,Output_dir,1)
