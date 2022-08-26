import os
import scipy
from random import shuffle
import numpy as np


# read in the hf-gcs with strength
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

def track_transform(GCS_dict):
    GCS_coverage_dict ={}
    for key, sub_dict in GCS_dict.items():
        coverage = [0]*4642893
        for pos, strength in sub_dict.items():
            coverage[pos-1] = strength
        GCS_coverage_dict[key]=coverage
    return GCS_coverage_dict

def write_track_file(path_out, GCS_dict):
    for drug, val_list in GCS_dict.items():
        file_name = drug +'_track.txt'
        completeName = os.path.join(path_out, file_name)
        track_out=open(completeName, 'w')
        for val in val_list:
            track_out.write(str(val)+'\n')
        track_out.close() 

pwd = '/GCS/GCS_extended_version/'  #need to change to proper directory
GCSs_input={'MOXI': pwd+"MOXI.txt",'LEVO': pwd+"LEVO.txt", 'NOR': pwd+ "NOR.txt",'CIP': pwd+"CIP.txt",
           'GEMI': pwd+"GEMI.txt"}

path_out = '/DNAPlotter_track/' #user defined
if not os.path.exists(path_out):
        os.makedirs(path_out)

GCSs_sets_dict = GCSs_parsing(GCSs_input)
GCSs_sets_dict = strength_normalization (GCSs_sets_dict)
GCS_coverage_dict = track_transform(GCSs_sets_dict)
write_track_file(path_out, GCS_coverage_dict)

