#motif scoring function
# 2019, Author:Cathy Tang
# 2022 updated
# The script is used to scan log-odds score of a given seqence or acrross the genome to predict cleavage probability 

import os
import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq
from Bio import motifs
import pandas as pd
from random import shuffle

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
# High Fidelity GCSs data parsing, returns a dictionary of {GCS position:strength}
#######
def GCSs_parsing(input_dict):
    GCSs_sets_dict={}
    for k, v in input_dict.items():
        GCSs_dict={}
        filein=open(v, 'r')
        for line in filein:
            line=line.rstrip().split(',')
            if line[0] not in ['pos']:
                GCSs_dict[int(line[0])]=float(line[23])
            else:
                continue
        filein.close()
        GCSs_sets_dict[k]=GCSs_dict
        print('Number of trusted GCSs for ' + str(k) + ' : ' + str(len(GCSs_dict)))
    return GCSs_sets_dict


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
    win_width_l=63
    win_width_r=67
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
    return seqs,pwm, m,pssm

def wrap_motif_construct(Source_genome_path,GCSs_files_paths, weight):
    Source_sequence=obtain_seq(Source_genome_path)[0]
    GCSs_sets_dict = GCSs_parsing(GCSs_files_paths)
    if weight ==1:
        GCSs_sets_dict = strength_normalization (GCSs_sets_dict)
    for k, v in GCSs_sets_dict.items():
        GCSs_dict = v
        print('Processing samples from: ' + k )
        seqs,pwm, m,pssm = motif_construct(GCSs_dict, Source_sequence, weight)
    return seqs,pwm, m,pssm

def scan_score(pssm,sequence,Output_path,target_name,win_width_l=63, win_width_r=67,genome = True):
    if genome == True:
       sequence=obtain_seq(fasta_path_input)[0]
    #Scans forward sequence
    test_seq=Seq(str(sequence.upper()))
    whole_genome_scores=pssm.calculate(test_seq) 
    outfile=open(Output_path + target_name + '_scan_forward_with_combined_motif.txt', 'w')
    
    dict_f={}
    #If test_seq len is equal to pssm, pssm.calculate() returns float32 but not the list.
    if len(test_seq)==win_width_l+win_width_r:
        one_pos_score=[]
        one_pos_score.append(whole_genome_scores)
        whole_genome_scores=one_pos_score
        dict_f[63+2]=[whole_genome_scores, test_seq[63+1]]
    for i in range(len(whole_genome_scores)):
        outfile.write(str(i+63+2) + '\t' + str(whole_genome_scores[i]) + '\t'+ str(test_seq[i+63+1]) + '\n')
        dict_f[i+63+2]=[whole_genome_scores[i], test_seq[i+63+1]]
    outfile.close()

    #Scans reverse complement sequence
    test_seq_rc=Seq(str(sequence.upper())).reverse_complement()
    whole_genome_scores_rc=pssm.calculate(test_seq_rc)   
    outfile_rc=open(Output_path + target_name + '_scan_rc_with_combined_motif.txt', 'w')
    dict_rc={}
    
    #If test_seq_rc len is equal to pssm, pssm.calculate() returns float32 but not the list.
    if len(test_seq_rc)==win_width_l+win_width_r:
        one_pos_score_rc=[]
        one_pos_score_rc.append(whole_genome_scores_rc)
        whole_genome_scores_rc=one_pos_score_rc 
        dict_rc[67-1]=[whole_genome_scores_rc, test_seq_rc[67-1]]
    for i in range(len(whole_genome_scores_rc)):
        outfile_rc.write(str(i+67-1) + '\t' + str(whole_genome_scores_rc[-i-1]) + '\t'+ str(test_seq_rc[-(i+67-1)]) + '\n')
        dict_rc[i+67-1]=[whole_genome_scores_rc[-i-1], test_seq_rc[-(i+67-1)]]
    outfile_rc.close()     
    return dict_f, dict_rc,pssm


#Input: GCS input (user defined)
pwd = '/GCS_extended_version/'
GCSs_input={'LEVO': pwd+"LEVO.txt"}
#Input: path to the mapping genome FASTA.
fasta_path_input = "/Mu_ori_mu_insert_MG1655.fa"

#Output path (user defined)
Output_dir="/motif_construct/"
if not os.path.exists(Output_dir):
    os.makedirs(Output_dir)     
seqs,pwm, m,pssm = wrap_motif_construct(fasta_path_input,GCSs_input, 1)

#Output path for score (user defined)
Output_path="/Volumes/TJC/GCS_calling_files/motif_construct/score/"
if not os.path.exists(Output_path):
    os.makedirs(Output_path)   

sequence = 'ttgagccaggaatacattgaagacaaagaagtcacattgacaaagttaagtagcggccgccgccttctggaggcgttgctgatccttattgtcctgtttgccgtctggttgatggctgccttactaagctttaacccttcggaccccagctggtcgcaaacggcctggcatgaacctatccataatttaggtgggatgcccggtgcgtggttggcagatacgctgttctttatttttggcgtgatggcttacaccattcccgtcattattgtcggcggttgttggtttgcctggcgtcatcagtccagcgacgaatacattgattattttgccgtttcgctacgcatcattggcgttttggcgctcatccttacctcctgtggtctggcggcaatcaacgctgacgatatctggtattttgcctccggtggcgtcattggcagcttactaagcactacgctacaaccactgctacacagtagcgggggaactattgcgctgctctgcgtttgggcagcgggcctgacgttgttcaccggttggtcatgggtgaccattgctgaaaaactcggcggctggattttaaacattctcaccttcgccagtaatcgtacccgtcgcgatgatacctgggtcgatgaagatgagtatgaagacgacgaagagtatgaagatgaaaatcacggcaaacagcatgaatcacgccgtgcccgtattcttcgcggcgcgctagcgcgtcgtaaacggttggcggaaaaattcattaatccgatggggcggcaaacagacgctgcgttgttctccggtaagcggatggatgatgacgaagagattacctacactgcacgcggtgtggctgccgacccggacgacgtcctattttcgggcaatcgtgcaacgcagccagaatatgacgaatacgatccattattaaacggtgcgccaattaccgaacctgtcgctgtagcagctgctgctaccacggcgacacaaagctgggctgcgccggttgaacctgtgactcagacgccgcctgttgcctctgttgatgttccacctgcgcaacctacagtagcctggcagcctgtaccgggtccacaaacgggagagccggttattgctcctgcaccggaaggttacccacagcagtcacaatatgcgcagcctgcagtgcaatataatgagccgctgcaacaaccagtacagccgcagcagccgtattatgcacctgcagctgaacaacctgcgcaacagccgtattatgcccctgcgccagaacaaccggtggcaggtaacgcctggcaagccgaagagcagcaatccacttttgctccacagtctacataccagactgagcaaacttatcagcagccagccgctcaggagccgttgtaccaacagccgcaacccgttgaacagcagcctgttgtggagcctgaacccgttgtagaagagacaaaacccgcgcgtccgccgctttactactttgaagaagttgaagagaagcgagcccgtgaacgtgaacaacttgcggcctggtatcaaccgattccagaaccggttaaagaaccagaaccgatcaaatcttcgctgaaagcaccttctgttgcagcagtacctccagtagaagccgctgccgctgtttccccgctggcatctggcgtgaaaaaagcgacactggcgacgggggctgccgcaaccgttgccgcgccagtcttcagtctggcaaatagcggtggaccgcgtcctcaggtcaaagaggggattggtccgcagttgccacgaccgaaacgtatccgcgtgccaactcgtcgtgaactggcgtcttacggtattaagctgccctcacagcgtgcggcggaagaaaaagcccgtgaagcccagcgcaatcagtacgattctggcgatcagtacaacgatgatgaaatcgatgcgatgcagcaggatgaactggcacgtcagttcgcccagacacagcagcaacgctatggcgaacagtatcaacatgatgtgcccgtaaacgcagaagatgcagatgctgcggcagaggctgaactggctcgtcagtttgcgcaaactcaacaacaacgttattccggcgaacaaccggctggggcgaatccgttctcgctggatgattttgaattttcgccaatgaaagcgttgctggatgatggtccacacgaaccgttgtttacgccaattgttgaacctgtacagcagccgcaacaaccggttgcaccgcagcagcaatatcagcagccgcaacaaccagttccgccgcagccgcagtatcagcagccacaacagccggttgcgccgcagccacaatatcagcagccgcaacaaccggttgcgccacagcagcaatatcagcagccgcaacaaccggttgcgccgcagcagcagtatcagcagccacaacagccagttgcgccacaaccgcaggataccctgcttcatccgctgttgatgcgtaatggcgacagccgtccgttgcataaaccgacgacgccgctgccttctctggatttgctgacaccgccgccgagcgaagtggagccggtagatacctttgcgcttgaacaaatggctcgcctggtggaagcgcgtctggctgatttccgtattaaagccgatgtcgtcaattactctccggggccggttatcactcgctttgaattgaacctggcaccgggcgtaaaagcggcgcgcatttctaacttgtcacgggaccttgcccgttcactttcgacggtggcggtgcgtgtcgttgaagttattcctggcaaaccctatgtaggtctggagttaccgaataaaaaacgacaaaccgtttatctgcgcgaagttttggataacgccaaattccgcgataatccgtcgccattaaccgtggtgctgggtaaagatatcgccggtgagccggtggttgccgatctggcgaaaatgccgcacttgttggttgcggggactaccggttccggtaaatctgtcggtgtgaacgcgatgatcctgagcatgctttataaagcacagccagaagatgtgcgtttcatcatgatcgacccgaaaatgctggagctttcggtttatgaaggcattccgcatctgttaacggaagtcgttactgatatgaaagatgccgccaacgcgctgcgctggtgtgttaacgagatggagcgtcggtataaactgatgtctgcgctgggtgtgcgtaatctggcgggttataacgaaaaaattgctgaagccgatcgcatgatgcgtccgattccagacccgtactggaagccgggtgacagtatggatgcccagcatccggtgctgaaaaaagaaccatacattgtggtgttggttgacgaatttgccgacctgatgatgacggtaggtaaaaaagtggaagagctgatagcacgtctggcgcaaaaagcccgtgccgcgggtatccacctcgtactggcaactcagcgtccatcggttgatgttattactggtctgattaaagcgaatattccgacccgtatcgcctttaccgtatccagtaagattgactcacgtaccattcttgatcaggctggcgcggaatcactgctgggtatgggggatatgctctactctgggccgaactccacgttgccggtacgtgtccatggtgcttttgttcgcgatcaggaagttcatgccgtggtgcaggactggaaagcgcgtggtAgAccacagtatgttgatggcatcacctccgacagcgaaagcgaaggtggtgcgggtggtttcgatggcgctgaagaactggatccgttgttcgatcaggcggtgcagtttgtcactgaaaaacgcaaagcgtcaatttctggcgtacagcgtcagttccgcattggttataaccgtgcagcgcgtattatcgaacagatggaagcgcaggggattgtcagcgaacaggggcacaacggtaatcgtgaagtgctggccccaccgccgtttgactaa'
dict_f, dict_rc,pssm = scan_score(pssm,sequence=sequence,Output_path=Output_path,target_name='whole_genome',
win_width_l=63, win_width_r=67,genome = False)

## Mu Scrambled scanning
def shuffle_word(word):
    word = list(word)
    shuffle(word)
    return ''.join(word)

#creat shuffled Mu sites
dict = {} # Mu: score
Mu_scramble = []
Mu_SGS="GATTCATACACCGTTAAATACCGGTTTAAAAATCCCGTGGCGCGTTTTAAAAAATCTGTGCGGGTGATTTTATGCCTGATTCTGTTTATTGCCTCAGAGCGGCGCTGACGCGTTTTCTGATGGCATCAAAAAT"
Mu_l_60 = "AAGGAAGATAAAACGGGATTCATACACCGTTAAATACCGGTTTAAAAATCCCGTGGCGCG"
Mu_r_60 = "TTTCTGATGGCATCAAAAATTTCCTGTTCCCCGGTCTTATCCAGCCCCATATAAGGACGC"
Mu_Mid_69="TTTTAAAAAATCTGTGCGGGTGATTTTATGCCTGATTCTGTTTATTGCCTCAGAGCGGCGCTGACGCGT"

# for i in range(200):
#     shuffled = shuffle_word(Mu_Mid_69)
#     Mu_scram = Mu_l_60 + shuffled + Mu_r_60
#     Mu_scramble.append(Mu_scram)




