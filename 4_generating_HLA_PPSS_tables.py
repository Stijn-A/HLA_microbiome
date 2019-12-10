#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov  5 14:30:49 2019

@author: stijn
"""
from itertools import combinations
import json

alleles_per_pseudo_seq = {}
with open('/home/stijn/stijn2/groningen/HLA_ABC_pseudo.dat', 'r') as f:
    for line in f:
        line = line.strip('\n')
        words = line.split(' ')
        if words[1] not in alleles_per_pseudo_seq.keys():
            alleles_per_pseudo_seq[words[1]] = [words[0]]
        else:
            alleles_per_pseudo_seq[words[1]].append(words[0])

pseudo_seq2HLA = {}
with open('/home/stijn/stijn2/groningen/uniq_pseudo_seq_alleles.txt') as f:
    for line in f:
        line = line.split(' ')
        HLA_allele = line[0]
        pseudo_seq = line[1]
        pseudo_seq = pseudo_seq.strip('\n')
        pseudo_seq2HLA[pseudo_seq] = HLA_allele

binder_2_A = {}
binder_2_B = {}
binder_2_C = {}
for k, v in pseudo_seq2HLA.items():
    
    if v.startswith('HLA-A'):
        binder_2_A[k] = []
        with open('/home/stijn/stijn2/groningen/peptides_HLA_A/' + v) as f:
            for line in f:
                if line.startswith('    1'):        #this are the lines with a predicted peptide affinity
                    line = line.strip('\n').strip(' <= SB').strip(' <= WB') #remove extra's so last item is %Rank
                    line_split = line.split(' ')
                    last_item = float((line_split[-1::])[0]) #the %Rank binding affinities
                    if last_item <= 2:  # %Rank below 2 is predicted to bind the HLA allele (weak + strong binders)
                        binder_2_A[k].append(line_split[9]) #the peptide
                        
    if v.startswith('HLA-B'):
        binder_2_B[k] = []
        with open('/home/stijn/stijn2/groningen/peptides_HLA_B/' + v) as f:
            for line in f:
                if line.startswith('    1'):        #this are the lines with a predicted peptide affinity
                    line = line.strip('\n').strip(' <= SB').strip(' <= WB') #remove extra's so last item is %Rank
                    line_split = line.split(' ')
                    last_item = float((line_split[-1::])[0]) #the %Rank binding affinities
                    if last_item <= 2:  # %Rank below 2 is predicted to bind the HLA allele (weak + strong binders)
                        binder_2_B[k].append(line_split[9]) #the peptide
                        
    if v.startswith('HLA-C'):
        binder_2_C[k] = []
        with open('/home/stijn/stijn2/groningen/peptides_HLA_C/' + v) as f:
            for line in f:
                if line.startswith('    1'):        #this are the lines with a predicted peptide affinity
                    line = line.strip('\n').strip(' <= SB').strip(' <= WB') #remove extra's so last item is %Rank
                    line_split = line.split(' ')
                    last_item = float((line_split[-1::])[0]) #the %Rank binding affinities
                    if last_item <= 2:  # %Rank below 2 is predicted to bind the HLA allele (weak + strong binders)
                        binder_2_C[k].append(line_split[9]) #the peptide
'''
with open('/home/stijn/stijn2/9_python_files/binder_2_A.txt', 'w') as f:
    json.dump(binder_2_A, f)
    
with open('/home/stijn/stijn2/9_python_files/binder_2_B.txt', 'w') as f:
    json.dump(binder_2_B, f)

with open('/home/stijn/stijn2/9_python_files/binder_2_C.txt', 'w') as f:
    json.dump(binder_2_C, f)
'''
#calculate jaccard Predicted Peptidome Similarity between MHC pseudo sequences

with open('/home/stijn/stijn2/9_python_files/binder_2_A.txt') as f:
    binder_2_A = json.load(f)
    
with open('/home/stijn/stijn2/9_python_files/binder_2_B.txt') as f:
    binder_2_B = json.load(f)

with open('/home/stijn/stijn2/9_python_files/binder_2_C.txt') as f:
    binder_2_C = json.load(f)


jaccard_2_A = {}
for i in combinations(binder_2_A.keys(), 2):
    same_items = list(set(binder_2_A[i[0]]).intersection(binder_2_A[i[1]]))
    nr_peptides_shared = len(same_items)
    total_peptides = len(binder_2_A[i[0]]) + len(binder_2_A[i[1]])
    jaccard_index = nr_peptides_shared/(total_peptides - nr_peptides_shared)
    key = i[0] + '-' + i[1]
    jaccard_2_A[key] = jaccard_index

jaccard_2_B = {}
for i in combinations(binder_2_B.keys(), 2):
    same_items = list(set(binder_2_B[i[0]]).intersection(binder_2_B[i[1]]))
    nr_peptides_shared = len(same_items)
    total_peptides = len(binder_2_B[i[0]]) + len(binder_2_B[i[1]])
    jaccard_index = nr_peptides_shared/(total_peptides - nr_peptides_shared)
    key = i[0] + '-' + i[1]
    jaccard_2_B[key] = jaccard_index

jaccard_2_C = {}
for i in combinations(binder_2_C.keys(), 2):
    same_items = list(set(binder_2_C[i[0]]).intersection(binder_2_C[i[1]]))
    nr_peptides_shared = len(same_items)
    total_peptides = len(binder_2_C[i[0]]) + len(binder_2_C[i[1]])
    jaccard_index = nr_peptides_shared/(total_peptides - nr_peptides_shared)
    key = i[0] + '-' + i[1]
    jaccard_2_C[key] = jaccard_index

'''
with open('/home/stijn/stijn2/9_python_files/jaccard_2_A_groningen.txt', 'w') as f:
    json.dump(jaccard_2_A, f)
    
with open('/home/stijn/stijn2/9_python_files/jaccard_2_B_groningen.txt', 'w') as f:
    json.dump(jaccard_2_B, f)

with open('/home/stijn/stijn2/9_python_files/jaccard_2_C_groningen.txt', 'w') as f:
    json.dump(jaccard_2_C, f)
'''
#for HLA gene A
with open('/home/stijn/stijn2/9_python_files/jaccard_2_A_groningen.txt') as f:
    jaccard_2_A = json.load(f)

allele2pseudo_A = {} #make dictionary with all alleles translated to pseudo seq.
for k, v in alleles_per_pseudo_seq.items():
    for i in v:
        if i.startswith('HLA-A'):
            allele2pseudo_A[i] = k
sorted_allele2pseudo_A = {}
for i in sorted(allele2pseudo_A.keys()):
    sorted_allele2pseudo_A[i] = allele2pseudo_A[i]

HLA_A_jaccard = {}
for i1 in sorted_allele2pseudo_A.items():
    for i2 in sorted_allele2pseudo_A.items():
        string_A = i1[0] + '-' + i2[0]
        if i1[1] == i2[1]:
            HLA_A_jaccard[string_A] = 1
        else:
            string_pseudo_1 = i1[1] + '-' + i2[1]
            #print(string_pseudo_1)
            string_pseudo_2 = i2[1] + '-' + i1[1]
            if string_pseudo_1 in jaccard_2_A.keys():
                HLA_A_jaccard[string_A] = jaccard_2_A[string_pseudo_1]
            else:
                HLA_A_jaccard[string_A] = jaccard_2_A[string_pseudo_2]
        
A_keys = sorted_allele2pseudo_A.keys()
M_A = []
for i, k1 in enumerate(A_keys):
    M_A.append([])
    for k2 in A_keys:
        string_1 = k1 + '-' + k2
        #print(string_1)
        M_A[i].append(HLA_A_jaccard[string_1])
        
with open ('/home/stijn/stijn2/groningen/table_HLA_A.txt', 'w') as f: 
    f.write('HLA-alleles' + '\t')
    for allele in A_keys:
        f.write(allele + '\t')
    f.write('\n')    
    for i, allele in enumerate(A_keys):
        f.write(allele + '\t' + '\t'.join([str(score) for score in M_A[i]]) + '\n' )

# For HLA gene B
with open('/home/stijn/stijn2/9_python_files/jaccard_2_B_groningen.txt') as f:
    jaccard_2_B = json.load(f)

allele2pseudo_B = {} #make dictionary with all alleles translated to pseudo seq.
for k, v in alleles_per_pseudo_seq.items():
    for i in v:
        if i.startswith('HLA-B'):
            allele2pseudo_B[i] = k
sorted_allele2pseudo_B = {}
for i in sorted(allele2pseudo_B.keys()):
    sorted_allele2pseudo_B[i] = allele2pseudo_B[i]

HLA_B_jaccard = {}
for i1 in sorted_allele2pseudo_B.items():
    for i2 in sorted_allele2pseudo_B.items():
        string_B = i1[0] + '-' + i2[0]
        if i1[1] == i2[1]:
            HLA_B_jaccard[string_B] = 1
        else:
            string_pseudo_1 = i1[1] + '-' + i2[1]
            #print(string_pseudo_1)
            string_pseudo_2 = i2[1] + '-' + i1[1]
            if string_pseudo_1 in jaccard_2_B.keys():
                HLA_B_jaccard[string_B] = jaccard_2_B[string_pseudo_1]
            else:
                HLA_B_jaccard[string_B] = jaccard_2_B[string_pseudo_2]
        
B_keys = sorted_allele2pseudo_B.keys()
M_B = []
for i, k1 in enumerate(B_keys):
    M_B.append([])
    for k2 in B_keys:
        string_1 = k1 + '-' + k2
        #print(string_1)
        M_B[i].append(HLA_B_jaccard[string_1])
        
with open ('/home/stijn/stijn2/groningen/table_HLA_B.txt', 'w') as f: 
    f.write('HLA-alleles' + '\t')
    for allele in B_keys:
        f.write(allele + '\t')
    f.write('\n')    
    for i, allele in enumerate(B_keys):
        f.write(allele + '\t' + '\t'.join([str(score) for score in M_B[i]]) + '\n' )
        
# For HLA gene C
with open('/home/stijn/stijn2/9_python_files/jaccard_2_C_groningen.txt') as f:
    jaccard_2_C = json.load(f)

allele2pseudo_C = {} #make dictionary with all alleles translated to pseudo seq.
for k, v in alleles_per_pseudo_seq.items():
    for i in v:
        if i.startswith('HLA-C'):
            allele2pseudo_C[i] = k
sorted_allele2pseudo_C = {}
for i in sorted(allele2pseudo_C.keys()):
    sorted_allele2pseudo_C[i] = allele2pseudo_C[i]

HLA_C_jaccard = {}
for i1 in sorted_allele2pseudo_C.items():
    for i2 in sorted_allele2pseudo_C.items():
        string_C = i1[0] + '-' + i2[0]
        if i1[1] == i2[1]:
            HLA_C_jaccard[string_C] = 1
        else:
            string_pseudo_1 = i1[1] + '-' + i2[1]
            string_pseudo_2 = i2[1] + '-' + i1[1]
            if string_pseudo_1 in jaccard_2_C.keys():
                HLA_C_jaccard[string_C] = jaccard_2_C[string_pseudo_1]
            else:
                HLA_C_jaccard[string_C] = jaccard_2_C[string_pseudo_2]
        
C_keys = sorted_allele2pseudo_C.keys()
M_C = []
for i, k1 in enumerate(C_keys):
    M_C.append([])
    for k2 in C_keys:
        string_1 = k1 + '-' + k2
        #print(string_1)
        M_C[i].append(HLA_C_jaccard[string_1])
        
with open ('/home/stijn/stijn2/groningen/table_HLA_C.txt', 'w') as f: 
    f.write('HLA-alleles' + '\t')
    for allele in C_keys:
        f.write(allele + '\t')
    f.write('\n')    
    for i, allele in enumerate(C_keys):
        f.write(allele + '\t' + '\t'.join([str(score) for score in M_C[i]]) + '\n' )