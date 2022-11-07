#03-28
#read in coverage file and save as sepearte files

import os
import pandas as pd
pwd="/Volumes/TJC/GCS_calling_files/coverage/"

path_out = '/Volumes/TJC/GCS_calling_files/coverage_drug/'
if not os.path.exists(path_out):
    os.makedirs(path_out)

path = {}
path['LEVO'] = pwd + 'LEVO.tabular'
path['MOXI'] = pwd + 'MOXI-NOR.tabular'
path['NOR'] = pwd + 'MOXI-NOR.tabular'
path['CIP'] = pwd + 'CIP-GEMI.tabular'
path['GEMI'] = pwd + 'CIP-GEMI.tabular'
path['4Mu-LEVO'] = pwd + '4Mu_LEVO.tabular'

column={}
HEAD = ['Seqname', 'Pos']
BODY = ['Flagless1','Flagless2','Flagless3','FlaglessIC1','FlaglessIC2','FlaglessIC3', 
'Flag1','Flag2','Flag3','FlagIC1','FlagIC2','FlagIC3']
BODY_COPY = [ name+'_copy' for name in BODY]
BODY2 = ['Flag1','Flag2','Flag3','FlagIC1','FlagIC2','FlagIC3', 
'Flagless1','Flagless2','Flagless3','FlaglessIC1','FlaglessIC2','FlaglessIC3']
BODY2_COPY = [ name+'_copy' for name in BODY2]
column['LEVO'] = HEAD + BODY2
column['NOR'] = HEAD + BODY + BODY_COPY
column['MOXI'] = HEAD + BODY + BODY_COPY
column['CIP'] = HEAD + BODY2 + BODY2_COPY
column['GEMI'] = HEAD + BODY2 + BODY2_COPY
column['4Mu-LEVO'] = HEAD + BODY2 + BODY2_COPY

# for drug in ['MOXI','NOR', 'LEVO', 'CIP', 'GEMI']:
for drug in ['4Mu-LEVO']:
    print('Now is processing: ' + str(path[drug]))
    table = pd.read_table(path[drug], sep = '\t',header = None)
    table.columns= column[drug]
    for s in ['Flag', 'Flagless']:
        for surfix in ['', 'IC']:
            for i in range(1,4):
                name = s + surfix + str(i)
                column_name = name+ '_copy' if drug in ['GEMI', 'NOR','4Mu-LEVO'] else name
                if not os.path.exists(path_out+drug):
                    os.makedirs(path_out+drug)
                table[['Pos',column_name]].to_csv(path_out+ drug + '/' +name+ '.txt', sep = '\t')

