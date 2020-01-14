#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec  9 12:08:53 2019

@author: stijn
Stijn P. Andeweg
s.p.andeweg@students.uu.nl


This script contains six parts:
(1) reading all the BAM files into python and saving them using json (line 60)
        - is already done and because it takes a long time now we use the saved 
        dictionaries with json
(2) select the datasets on the HLA class I or II profile  (line 191)
(3) extract bacterial reads and filter datasets (line 1258)
(4) calculate the alpha diversity (line 1395)
(5) calculate the beta diversity (line 1956)
(6) PPSS (2361)

Results:                                Line
    Figure 1                            1112(A) 1021(B)
    Figure 2                            1541(A) 1587(B) 1780(C) 1826 (D)
    Figure 3                            1671
    Figure 4                            2511(A) 2550(B) 2572(C) 2602(D)
    Supplementary Figure 1              2627 (A) 2678(B) 2728(C)
    Supplementary Figure 2              2777
    Supplementary Figure 3              1244
    Supplementary Figure 4              568 (A-C), 748 (D-E), 995(F)
    Supplementary Figure 5              1541(A) 1616 (B) 1780(C) 1826 (D)
    Supplementary Figure 6              1904(A-B) 1671(C)
    Supplementary Figure 7              2119(A) 2157(B) 2181(C) 2218(D)
    Supplementary Figure 8              2308(A) 2339(B)
    Individual statistics main text     2807
    
    
For loading and saving of files my local directories are used.
Change /home/stijn/otherdirectories/ to own directories for reproducability.

"""

import json
import subprocess
from itertools import chain
from operator import itemgetter
import operator
import numpy as np
import sys
sys.path.append('/home/stijn/.local/lib/python3.5/site-packages')
#import dendropy
sys.path.append('/home/stijn/Downloads/EMDUnifrac-master/src/')
import EMDUnifrac as EMDU
from scipy import stats
from skbio.diversity import alpha
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.markers import TICKDOWN
import matplotlib
from itertools import product
from itertools import combinations

'''
part (1) reading all the BAM files into python and saving them using json



# functions make a list of all the bam files in given directory 'path'
def bam_loc_list(path):
# makes list of bam file names
    result = subprocess.run(['ls',  path], stdout=subprocess.PIPE)
    output_string1 = str(result.stdout)
    output_string = output_string1.strip('\"b')[1:-1]
    output_list = output_string.split("\\n")
    output_list = list(filter(None, output_list))
    for i in output_list:
        if i.endswith('bai'):
            output_list.remove(i)
    return output_list


# run samtools view to make readable text file from bam file
def samtools_view(path, bam_name):        
    bash_cmd1 = 'samtools view -q 20 '
    right_arrow = ' > '
    txt ='_txt'
    a = [bash_cmd1 + path + bam_name + right_arrow + path + bam_name + txt] 
    subprocess.call(a, shell=True)


# open the samtools_view text file and load into python
def open_bam_txt(path, fname):
    with open(path + fname, 'r') as textFile:
        lines = [line.split() for line in textFile]
        bam_lines = [item[:11] for item in lines]
        return bam_lines


# function to remove the samtools_view created text file after loading it into python
def rm_bam_txt(path, bam_name):
    bash_cmd2 = 'rm -rf '
    txt ='_txt'
    b = [bash_cmd2 + path + bam_name + txt]
    subprocess.call(b, shell=True)


# combine samtools_view, open_bam_txt and open_bam_txt to load file into
# python and remove unnecessary text file    
def bam_to_bam_txt(path, bam_name):
    samtools_view(path, bam_name)
    txt ='_txt'
    fname = bam_name + txt
    bam_txt = open_bam_txt(path, fname)
    rm_bam_txt(path, bam_name)
    return bam_txt


# create dictionary with as key the reference genome location and as value the 
# read numbers mapping
def d_Qname_RNAME(bam_txt):
    d1 = {}
    for Qname, flag, RNAME, POS, MAPQ, CIGAR, MRNM, MPOS, ISIZE, SEQQuery, QUAL\
    in bam_txt:
        d1.setdefault(Qname, []).append(RNAME)
    return d1


# Check how many primary reads map to a location, and devide the multimapping
# reads according to the nr of primary mapping reads to each location. 
# Returns dictionary with number of reads per location
def ConCarneSortingScheme(d1):
    location_one_map = {}
    for key, value in d1.items():
        number_mapping_locations = len(value)
        if number_mapping_locations == 1:
            for item in value:
                if item not in location_one_map.keys():
                    location_one_map[item] = 0
                    location_one_map[item] += 1
                else:
                    location_one_map[item] += 1
        if number_mapping_locations > 1:
            for item in value:
                if item not in location_one_map.keys():
                    location_one_map[item] = 0  
    reads_more_locations = {}
    for key, value in d1.items():
        number_mapping_locations = len(value)
        if number_mapping_locations > 1:
            total_per_value = 0
            for item in value:
                total_per_value += location_one_map[item]
            for item in value:
                if total_per_value > 0:    
                    reads_more_locations.setdefault(key, [])\
                    .append((item, location_one_map[item]/total_per_value))
                else:
                    reads_more_locations.setdefault(key, [])\
                    .append((item, 1/len(value)))
    second_locations = {}
    for key, value in reads_more_locations.items():
        for items in value:
            if items[0] not in second_locations.keys():
                second_locations[items[0]] = 0
                second_locations[items[0]] += items[1]
            else:
                second_locations[items[0]] += items[1]
    dict3 = {}
    for k, v in chain(location_one_map.items(), second_locations.items()):
        if k not in dict3.keys():
            dict3[k] = 0
            dict3[k] += v
        else:
            dict3[k] += v
    sorted_dict3 = sorted(dict3.items(), key=operator.itemgetter(1))
    return sorted_dict3


path= '/home/location/BAM_files/'   #change to location where mapped sequence
# files to SILVA and IMGT/HLA are (sorted BAM files)
input_list = bam_loc_list(path)   #makes a list of all present bam files
sep = '.'
BAM_files_dict = {}
for item in input_list:
    i1 = item.split(sep, 1)[0]
    if i1 not in BAM_files_dict.keys():
        bam_txt = bam_to_bam_txt(path, item)    #creates a bam text file
        d1 = d_Qname_RNAME(bam_txt)             #creates a dictionary
        dict4 = ConCarneSortingScheme(d1)       #sorts the multi-mapping reads
        BAM_files_dict[i1] = dict4              #creates nested dictionary 
                                        #dict[SRA_acc][reference_genome_location]
                                        #  = number of mapped reads
        

(2) select the datasets on the HLA profile / perform HLA typing
'''


#This are the BAM files loaded into python of the current study. (Part 1 done)
with open('/home/stijn/stijn2/9_python_files/all_sample_dict.txt') as f:
    all_sample_dict = json.load(f)
    
    
#remove multiple occurances of the same person from data
person_a = ['SRR8204380', 'SRR8204362', 'SRR8204372', 'SRR8204374']
person_b = ['SRR8204377', 'SRR8204369', 'SRR8204378']
person_c = ['SRR8204365', 'SRR8204347', 'SRR8204348']
person_d = ['SRR8204345']
remove = []
for k, v in all_sample_dict.items():
    if k in person_a:
        remove.append(k) 
    if k in person_b:
        remove.append(k) 
    if k in person_c:
        remove.append(k) 
    if k in person_d:
        remove.append(k) 
for k in remove: del all_sample_dict[k]


# substract HLA reads from dataset and create dictionary with only HLA reads
def only_HLA_dict(dir_sample):
    hla_reads_dict = {}
    for k, v in dir_sample.items():
        hla_reads_dict[k] = {}
        for a, b in v.items():
            if a.startswith('HLA'):
                hla_reads_dict[k][a] = b
    return hla_reads_dict
hla_dict = only_HLA_dict(all_sample_dict)


#creates HLA dict for 2-digit classification from the data containing 2,4,6,8 digits.
def HLA_collapse_digits(dir_bam_HLA, n): # for 2-digits n=1, 4-digits n =2, for 6 digits n=3
    HLA_4digit_dict = {}
    for key, value in dir_bam_HLA.items():
        HLA_4digit_dict[key] = {}
        for HLA, reads in value.items():
            splitHLA = HLA.split(":")
            HLA_digit = ":".join(splitHLA[:n])
            if HLA_digit not in HLA_4digit_dict[key].keys():
                HLA_4digit_dict[key][HLA_digit] = reads
            else:
                HLA_4digit_dict[key][HLA_digit] += reads
    return(HLA_4digit_dict)
hla_2digit_dict = HLA_collapse_digits(hla_dict, 1) #create dictionary with HLA two-digits display 


# seperate the HLA reads in a nested dictionary per HLA gene
def HLA_per_gene_dict(dict_HLA_4digit):
    dict_HLA_4digit_perHLA = {}
    for acc, HLA in dict_HLA_4digit.items():
        dict_HLA_4digit_perHLA[acc] = {}
        for HLA, reads in HLA.items():
            HLAgene = HLA.split('*')[0]
            if HLAgene not in dict_HLA_4digit_perHLA[acc].keys():
               dict_HLA_4digit_perHLA[acc][HLAgene] = {}
               dict_HLA_4digit_perHLA[acc][HLAgene][HLA] = reads
            else:
                dict_HLA_4digit_perHLA[acc][HLAgene][HLA] = reads
    return dict_HLA_4digit_perHLA
dict_HLA_2digit_perHLA = HLA_per_gene_dict(hla_2digit_dict)


# calculate for each HLA gene the R score
def calc_new_R_dict(dict_HLA_4digit_perHLA):
    HLA_allele_R_dict = {}
    for acc, value in dict_HLA_4digit_perHLA.items():
        HLA_allele_R_dict[acc] = {}
        for HLAgene, HLAalleles in value.items():
           alleles_sorted = sorted(HLAalleles.items(), key=lambda kv: kv[1])
           top2HLA = 0
           restHLA = 1                 #pseudo count
           for i in alleles_sorted[-2:]:
               top2HLA += i[1]
           for i in alleles_sorted[:-2]:
               restHLA += i[1]
           HLA_R = top2HLA/restHLA
           HLA_allele_R_dict[acc][HLAgene] = HLA_R
    return HLA_allele_R_dict
HLA_allele_R_dict_2d = calc_new_R_dict(dict_HLA_2digit_perHLA)    


# make a dict. with only the HLA Class I or Class II genes.
def select_classIorII_HLAR_perHLA(dict_HLA_4digit_perHLA, class_I):
    dict_HLA_4digit_perHLAcI = {}
    for k, v in dict_HLA_4digit_perHLA.items():
        dict_HLA_4digit_perHLAcI[k] = {}
        for k2, v2 in v.items():
            if k2 in class_I:
                dict_HLA_4digit_perHLAcI[k][k2] = v2
    return dict_HLA_4digit_perHLAcI
class_I = ['HLA-A', 'HLA-B', 'HLA-C']
dict_HLA_2digit_perHLAcI = select_classIorII_HLAR_perHLA(dict_HLA_2digit_perHLA,\
                                                         class_I)

class_II = ['HLA-DQB1', 'HLA-DRB1', 'HLA-DPB1']
dict_HLA_2digit_perHLAcII = select_classIorII_HLAR_perHLA(dict_HLA_2digit_perHLA,\
                                                          class_II)


# select only the R scores for the class I HLA genes
def make_ClassI_genesR_dict_f(HLA_allele_R_dict):
    class_I = ['HLA-A', 'HLA-B', 'HLA-C']
    ClassI_genes_dict = {}
    for k, v in HLA_allele_R_dict.items():
        ClassI_genes_dict[k] = {}
        for i in class_I:
            if i in v.keys():
                R_score = HLA_allele_R_dict[k][i]
                ClassI_genes_dict[k][i] = R_score
    
    ClassI_genes_dict_f = {}
    for k, v in ClassI_genes_dict.items():
        if not v == {}:
            ClassI_genes_dict_f[k] = v
    return ClassI_genes_dict_f
ClassI_genes_dict_f_2d = make_ClassI_genesR_dict_f(HLA_allele_R_dict_2d)


# select only the R scores for the class II HLA genes
def make_ClassII_genesR_dict_f(HLA_allele_R_dict, class_II):
    ClassII_genes_dict = {}
    for k, v in HLA_allele_R_dict.items():
        ClassII_genes_dict[k] = {}
        for i in class_II:
            if i in v.keys():
                R_score = HLA_allele_R_dict[k][i]
                ClassII_genes_dict[k][i] = R_score
    ClassII_genes_dict_f = {}
    for k, v in ClassII_genes_dict.items():
        if not v == {}:
            ClassII_genes_dict_f[k] = v
    return ClassII_genes_dict_f

ClassII_genes_dict_f_2d = make_ClassII_genesR_dict_f(HLA_allele_R_dict_2d, class_II)


# Selection for R >= 2 for all class I MHC genes.
selected_dict_CI_2d = {}
selected_HLARdict_CI_2d = {}
class_I = ['HLA-A', 'HLA-B', 'HLA-C']
for k, v in ClassI_genes_dict_f_2d.items():
    if (all(i in  v.keys() for i in class_I)):
        if (all( v2 >= 2 for v2 in v.values())):    #minimal R score
            dicti = all_sample_dict[k]
            dictHLA = dict_HLA_2digit_perHLAcI[k]
            selected_dict_CI_2d[k] = dicti
            selected_HLARdict_CI_2d[k] = dictHLA


# Selection for R >= 2 for all class II MHC genes.
selected_dict_CII_2d = {}
selected_HLARdict_CII_2d = {}
for k, v in ClassII_genes_dict_f_2d.items():   #change dict 2d or 4d
    count = 0
    for v2 in v.values():
        if v2 >= 2:     #minimal R score
            count += 1
    if count >= 3:      #number of genes present
        dicti = all_sample_dict[k]
        dictiHLA = dict_HLA_2digit_perHLAcII[k]
        selected_dict_CII_2d[k] = dicti
        selected_HLARdict_CII_2d[k] = dictiHLA


# take the selected samples on the R score and pick their top2 HLA hits to 
# create a HLA profile class I
HLA_profile_c1 = {}
for k, v in selected_HLARdict_CI_2d.items():    #change 2d or 4d
    HLA_profile_c1[k] = {}
    for k2, v2 in v.items():
        HLA_reads = []
        for k3, v3 in v2.items():
            HLA_reads.append((k3, v3))
        highest_HLA = max(HLA_reads, key=itemgetter(1))
        HLA_reads.remove(highest_HLA)
        if len(HLA_reads) > 0:
            cutoff = v2[highest_HLA[0]] * 0.2
            if max(HLA_reads, key=itemgetter(1))[1] >= cutoff:      
                #second hit of heterozygote need to be 20% of the first hit
                second_highest_HLA = max(HLA_reads, key=itemgetter(1))
            elif v2[highest_HLA[0]] >= 10:  #homozygotes need to have 10 reads
                                            #at least
                second_highest_HLA = highest_HLA
            else:
                second_highest_HLA = 'No allele'
        elif v2[highest_HLA[0]] >= 10:
            second_highest_HLA = highest_HLA
        else:
                second_highest_HLA = 'No allele'
        HLA_profile_c1[k][k2] = highest_HLA, second_highest_HLA
        
        
#take the selected samples on the R score and pick their top2 HLA hits to 
# create a HLA profile class II
HLA_profile_c2 = {}
for k, v in selected_HLARdict_CII_2d.items():    
    HLA_profile_c2[k] = {}
    for k2, v2 in v.items():
        HLA_reads = []
        for k3, v3 in v2.items():
            HLA_reads.append((k3, v3))
        highest_HLA = max(HLA_reads, key=itemgetter(1))
        HLA_reads.remove(highest_HLA)
        if len(HLA_reads) > 0:
            cutoff = v2[highest_HLA[0]] * 0.2
            if max(HLA_reads, key=itemgetter(1))[1] >= cutoff:      
                #second hit of heterozygote need to be 20% of the first hit
                second_highest_HLA = max(HLA_reads, key=itemgetter(1))
            elif v2[highest_HLA[0]] >= 10:  #homozygotes need to have 10 reads 
                                            #at least
                second_highest_HLA = highest_HLA
            else:
                second_highest_HLA = 'No allele'
        elif v2[highest_HLA[0]] >= 10:
            second_highest_HLA = highest_HLA
        else:
                second_highest_HLA = 'No allele'
        HLA_profile_c2[k][k2] = highest_HLA, second_highest_HLA


'''
Figures on HLA only Figure 1, Supplement figure 4
Preparation for NMDP comparison HLA class I
'''

# calculate HLA allele frequencies
HLA_counts_c1 = {}
for k, v in HLA_profile_c1.items():
    for k2, v2 in v.items():
        for k in v2:
            if k[0] not in HLA_counts_c1.keys():
                HLA_counts_c1[k[0]] = 1
            else:
                HLA_counts_c1[k[0]] += 1

#caclulate relative allele frequencies of current study
counts = 0
for k, v in HLA_counts_c1.items():
    counts += v
total_alleles_p_gene = counts/3              
HLA_frequency_c1 = {}
for k, v in HLA_counts_c1.items():
    frequency = v / total_alleles_p_gene
    HLA_frequency_c1[k] = frequency

def getHLAgenedict(HLA_frequency_dict, HLAgene):
    HLA_A = {}
    for k,v in HLA_frequency_dict.items():
        if k.startswith(HLAgene):
            HLA_A[k] = v
    return HLA_A
HLA_A_frequency = getHLAgenedict(HLA_frequency_c1, 'HLA-A')
HLA_B_frequency = getHLAgenedict(HLA_frequency_c1, 'HLA-B')
HLA_C_frequency = getHLAgenedict(HLA_frequency_c1, 'HLA-C')


# open NMDP files and read the European cohort frequencies
def openAlleleFrequencyNMDP(file_location, HLAgene, EUR_freq):
    dfA_xls = pd.read_excel(file_location)
    df1 = dfA_xls[[HLAgene, EUR_freq, 'EUR_rank']]
    A_dict_numbers = df1.set_index(HLAgene).transpose().to_dict(orient='list')
    
    A_dict_freq_4d = {}
    A_dict_rank = {}
    for k, v in A_dict_numbers.items():
        k = k.strip('g')
        k = k.strip('N')
        HLA = 'HLA-' + HLAgene + '*' + k[:2] + ':' + k[2:]
        A_dict_freq_4d[HLA] = v[0]
        A_dict_rank[HLA] = v[1]
    return A_dict_freq_4d

file_location = "/home/stijn/Documents/7_HLA_allelefrequencies/A.xls"
A_dict_freq_4d = openAlleleFrequencyNMDP("/home/stijn/Documents/7_HLA_allelefrequencies/A.xls", 'A','EUR_freq')
B_dict_freq_4d = openAlleleFrequencyNMDP("/home/stijn/Documents/7_HLA_allelefrequencies/B.xls", 'B','EUR_freq')
C_dict_freq_4d = openAlleleFrequencyNMDP("/home/stijn/Documents/7_HLA_allelefrequencies/C.xls", 'C','EUR_freq')


# translate four-digit frequencies to two-digit frequencies.
def d2tod4FreqDict(A_dict_freq_4d):
    A_dict_freq_2d = {}
    for k, v in A_dict_freq_4d.items():
        k = k.split(':')
        k = k[0]
        if k not in A_dict_freq_2d.keys():
            A_dict_freq_2d[k] = v
        else:
            A_dict_freq_2d[k] += v
    return A_dict_freq_2d
A_dict_freq_2d = d2tod4FreqDict(A_dict_freq_4d)
B_dict_freq_2d = d2tod4FreqDict(B_dict_freq_4d)
C_dict_freq_2d = d2tod4FreqDict(C_dict_freq_4d)


# make a complite list of all the digits present in both populations
def make_allele_list(database_dict, Mydata_dict):
    allele_list = []
    for k in database_dict.keys():
        k = k.split('*')
        k = k[1]
        allele_list.append(k)
    for i in Mydata_dict.keys():
        digits = i.split('*')
        
        digits = digits[1]
        digits = digits.strip('Q')
        digits = digits.strip('N')
        if digits not in allele_list:
            allele_list.append(digits)
    allele_list.sort()
    return allele_list

allele_list_A = make_allele_list(A_dict_freq_2d, HLA_A_frequency)
allele_list_B = make_allele_list(B_dict_freq_2d, HLA_B_frequency)
allele_list_C = make_allele_list(C_dict_freq_2d, HLA_C_frequency)


#add the frequencies from both populations  for HLA-A, B and C.
#data HLA-A       
labels_A = []
frequencies_mydata_A = []
frequencies_database_A = []
for i in allele_list_A:
    HLA = 'HLA-A*' + i
    labels_A.append(HLA)
    if HLA in HLA_A_frequency.keys():
        frequencies_mydata_A.append(HLA_A_frequency[HLA])
    else:
        frequencies_mydata_A.append(0)
    if HLA in A_dict_freq_2d.keys():
        frequencies_database_A.append(A_dict_freq_2d[HLA])
    else:
        frequencies_database_A.append(0)

#data HLA-B
labels_B = []
frequencies_mydata_B = []
frequencies_database_B = []
for i in allele_list_B:
    HLA = 'HLA-B*' + i
    labels_B.append(HLA)
    if HLA in HLA_B_frequency.keys():
        frequencies_mydata_B.append(HLA_B_frequency[HLA])
    else:
        frequencies_mydata_B.append(0)
    if HLA in B_dict_freq_2d.keys():
        frequencies_database_B.append(B_dict_freq_2d[HLA])
    else:
        frequencies_database_B.append(0)

#data HLA-C    C_dict_freq_2d    HLA_C_frequency
labels_C = []
frequencies_mydata_C = []
frequencies_database_C = []
for i in allele_list_C:
    HLA = 'HLA-C*' + i
    labels_C.append(HLA)
    if HLA in HLA_C_frequency.keys():
        frequencies_mydata_C.append(HLA_C_frequency[HLA])
    else:
        frequencies_mydata_C.append(0)
    if HLA in C_dict_freq_2d.keys():
        frequencies_database_C.append(C_dict_freq_2d[HLA])
    else:
        frequencies_database_C.append(0)


'''
Supplementary Figure 4A-C 'HLA frequencies compared between current study and 
National Marrow Donor Program (NMDP).'
'''

#plotting HLA-A
plt.figure(figsize=(10,8), dpi=100)
x = np.arange(len(labels_A))
width = 0.35

fig, ax = plt.subplots()
rects1 = ax.bar(x - width/2, frequencies_mydata_A, width, label='Current study')
rects2 = ax.bar(x + width/2, frequencies_database_A, width, label='NMDP')

ax.set_ylabel('Population frequency', fontsize=16)
ax.set_title('HLA-A', fontsize=16)
ax.set_xticks(x)
ax.set_xticklabels(labels_A, rotation=90)
ax.legend()

fig.tight_layout()
plt.savefig('/home/stijn/figures/figures_03_12_19_svg/HLA_A_freq_NMDP.svg')
plt.show()

#plotting HLA-B
plt.figure(figsize=(10,8), dpi=100)
x = np.arange(len(labels_B))

fig, ax = plt.subplots()
rects1 = ax.bar(x - width/2, frequencies_mydata_B, width, label='Current study')
rects2 = ax.bar(x + width/2, frequencies_database_B, width, label='NMDP')

ax.set_ylabel('Population frequency', fontsize=16)
ax.set_title('HLA-B', fontsize=16)
ax.set_xticks(x)
ax.set_xticklabels(labels_B, rotation=90)
ax.legend()

fig.tight_layout()
plt.savefig('/home/stijn/figures/figures_03_12_19_svg/HLA_B_freq_NMDP.svg')
plt.show()

#plotting HLA-C
plt.figure(figsize=(10,8), dpi=100)
x = np.arange(len(labels_C))

fig, ax = plt.subplots()
rects1 = ax.bar(x - width/2, frequencies_mydata_C, width, label='Current study')
rects2 = ax.bar(x + width/2, frequencies_database_C, width, label='NMDP')

ax.set_ylabel('Population frequency', fontsize=16)
ax.set_title('HLA-C', fontsize=16)
ax.set_xticks(x)
ax.set_xticklabels(labels_C, rotation=90)
ax.legend()

fig.tight_layout()
plt.savefig('/home/stijn/figures/figures_03_12_19_svg/HLA_C_freq_NMDP.svg')
plt.show()


'''
Preparation for NMDP comparison HLA class II Supplementary Figure 4D-E
'''

#calculating allele counts current study
HLA_counts_c2 = {}
for k, v in HLA_profile_c2.items():
    for k2, v2 in v.items():
        for k in v2:
            if k[0] not in HLA_counts_c2.keys():
                HLA_counts_c2[k[0]] = 1
            else:
                HLA_counts_c2[k[0]] += 1

#calculate total alleles
countsDQB1 = 0
countsDPB1 = 0
countsDRB1 = 0
for k, v in HLA_counts_c2.items():
    if k.startswith('HLA-DQB1'):
        countsDQB1 += v
    if k.startswith('HLA-DPB1'):
        countsDPB1 += v
    if k.startswith('HLA-DRB1'):
        countsDRB1 += v

#calculate relative allele frequencies
HLA_frequency_c2 = {}
for k, v in HLA_counts_c2.items():
    if k.startswith('HLA-DQB1'):
        frequency = v /  countsDQB1
    if k.startswith('HLA-DPB1'):
        frequency = v / countsDPB1
    if k.startswith('HLA-DRB1'):
        frequency = v / countsDRB1
    HLA_frequency_c2[k] = frequency

    
def HLAgeneC2dict(HLA_frequency_c2, HLAgene):
    HLAdict = {}
    for k,v in HLA_frequency_c2.items():
        if k.startswith(HLAgene):
            HLAdict[k] = v
    return HLAdict

HLA_DQB1_frequency_dict = HLAgeneC2dict(HLA_frequency_c2, 'HLA-DQB1')
HLA_DRB1_frequency_dict = HLAgeneC2dict(HLA_frequency_c2, 'HLA-DRB1')       

#open data from the NMDP             
def openAlleleFrequencyNMDPc2(file_location, HLAgene):
    dfA_xls = pd.read_excel(file_location)
    df1 = dfA_xls[[HLAgene, 'EUR_freq', 'EUR_rank']]
    A_dict_numbers = df1.set_index(HLAgene).transpose().to_dict(orient='list')
    
    A_dict_freq_4d = {}
    A_dict_rank = {}
    for k, v in A_dict_numbers.items():
        k = k.strip('g')
        k = k.strip('N')
        HLA = 'HLA-' + HLAgene + '*' + k[:2] + ':' + k[2:]
        A_dict_freq_4d[HLA] = v[0]
        A_dict_rank[HLA] = v[1]
    return A_dict_freq_4d
DQB1_dict_freq_4d = openAlleleFrequencyNMDPc2("/home/stijn/Documents/7_HLA_allelefrequencies/DQB1.xls", 'DQB1')
DRB1_dict_freq_4d = openAlleleFrequencyNMDPc2("/home/stijn/Documents/7_HLA_allelefrequencies/DRB1.xls", 'DRB1')

#transform four-digits to two-digits
def d2tod4FreqDict(A_dict_freq_4d):
    A_dict_freq_2d = {}
    for k, v in A_dict_freq_4d.items():
        k = k.split(':')
        k = k[0]
        if k not in A_dict_freq_2d.keys():
            A_dict_freq_2d[k] = v
        else:
            A_dict_freq_2d[k] += v
    return A_dict_freq_2d
DQB1_dict_freq_2d = d2tod4FreqDict(DQB1_dict_freq_4d)
DRB1_dict_freq_2d = d2tod4FreqDict(DRB1_dict_freq_4d)

#data HLA-DQB1   2D

allele_list_DQB1 = make_allele_list(DQB1_dict_freq_2d, DQB1_dict_freq_2d)

labels_DQ = []
frequencies_mydata_DQ = []
frequencies_database_DQ = []
for i in allele_list_DQB1:
    HLA = 'HLA-DQB1*' + i
    labels_DQ.append(HLA)
    if HLA in HLA_DQB1_frequency_dict.keys():
        frequencies_mydata_DQ.append(HLA_DQB1_frequency_dict[HLA])
    else:
        frequencies_mydata_DQ.append(0)
    if HLA in DQB1_dict_freq_2d.keys():
        frequencies_database_DQ.append(DQB1_dict_freq_2d[HLA])
    else:
        frequencies_database_DQ.append(0)


#data HLA-DRB1      2D

allele_list_DRB1 = make_allele_list(DRB1_dict_freq_2d, DRB1_dict_freq_2d)
labels_DR = []
frequencies_mydata_DR = []
frequencies_database_DR = []
for i in allele_list_DRB1:
    HLA = 'HLA-DRB1*' + i
    labels_DR.append(HLA)
    if HLA in HLA_DRB1_frequency_dict.keys():
        frequencies_mydata_DR.append(HLA_DRB1_frequency_dict[HLA])
    else:
        frequencies_mydata_DR.append(0)
    if HLA in DRB1_dict_freq_2d.keys():
        frequencies_database_DR.append(DRB1_dict_freq_2d[HLA])
    else:
        frequencies_database_DR.append(0)


'''
Supplementary Figure 4D-E 'HLA frequencies compared between current study and 
National Marrow Donor Program (NMDP).'
'''


#plotting HLA-DQB1
plt.figure(figsize=(10,8), dpi=100)
x = np.arange(len(labels_DQ))
width = 0.35

fig, ax = plt.subplots()
rects1 = ax.bar(x - width/2, frequencies_mydata_DQ, width, label='Current study')
rects2 = ax.bar(x + width/2, frequencies_database_DQ, width, label='NMDP')

ax.set_ylabel('Population frequency', fontsize=16)
ax.set_title('HLA-DQB1', fontsize=16)
ax.set_xticks(x)
ax.set_xticklabels(labels_DQ, rotation=90)
ax.legend()

fig.tight_layout()
plt.savefig('/home/stijn/figures/figures_03_12_19_svg/HLA_DQB1_freq_NMDP.svg')
plt.show()   


#plotting HLA-DRB1
plt.figure(figsize=(10,8), dpi=100)
x = np.arange(len(labels_DR))
width = 0.35

fig, ax = plt.subplots()
rects1 = ax.bar(x - width/2, frequencies_mydata_DR, width, label='Current study')
rects2 = ax.bar(x + width/2, frequencies_database_DR, width, label='NMDP')

ax.set_ylabel('Population frequency', fontsize=16)
ax.set_title('HLA-DRB1', fontsize=16)
ax.set_xticks(x)
ax.set_xticklabels(labels_DR, rotation=90)
ax.legend()

fig.tight_layout()
plt.savefig('/home/stijn/figures/figures_03_12_19_svg/HLA_DRB1_freq_NMDP.svg')
plt.show()


'''
Data preparation Figure 1B and Supplementary Figure 4F
'''

#classify the genes as homozygote, hetorozygote or one-allele known
zygote_mhc_profile_c1 = {}
for k, v in HLA_profile_c1.items():   
    zygote_mhc_profile_c1[k]  = {}
    for k2, v2 in v.items():
        if v2[0][0] == v2[1][0]:
            zygote_mhc_profile_c1[k][k2] = 'homozygote'
        elif v2[1] == 'No allele':
            zygote_mhc_profile_c1[k][k2] = 'one_allele'
        else:
            zygote_mhc_profile_c1[k][k2] = 'heterozygote'


counts_zygote_A = {}
counts_zygote_B = {}
counts_zygote_C = {}
for k, v in zygote_mhc_profile_c1.items():
    for k2, v2 in v.items():
        if k2 == 'HLA-A':
            if v2 not in counts_zygote_A.keys():
                counts_zygote_A[v2] = 1
            else:
                counts_zygote_A[v2] += 1
        if k2 == 'HLA-B':
            if v2 not in counts_zygote_B.keys():
                counts_zygote_B[v2] = 1
            else:
                counts_zygote_B[v2] += 1
        if k2 == 'HLA-C':
            if v2 not in counts_zygote_C.keys():
                counts_zygote_C[v2] = 1
            else:
                counts_zygote_C[v2] += 1

#calculate relative zygosity frequency
#A
relative_freq_zygosity_A = {}
totalcounts = counts_zygote_A['homozygote'] +counts_zygote_A['heterozygote']
relative_freq_zygosity_A['homozygote'] = counts_zygote_A['homozygote']/totalcounts
relative_freq_zygosity_A['heterozygote'] = counts_zygote_A['heterozygote']/totalcounts
#B
relative_freq_zygosity_B = {}
totalcounts = counts_zygote_B['homozygote'] + counts_zygote_B['heterozygote']
relative_freq_zygosity_B['homozygote'] = counts_zygote_B['homozygote']/totalcounts
relative_freq_zygosity_B['heterozygote'] = counts_zygote_B['heterozygote']/totalcounts
#C
relative_freq_zygosity_C = {}
totalcounts = counts_zygote_C['homozygote'] + counts_zygote_C['heterozygote']
relative_freq_zygosity_C['homozygote'] = counts_zygote_C['homozygote']/totalcounts
relative_freq_zygosity_C['heterozygote'] = counts_zygote_C['heterozygote']/totalcounts


df_A = pd.DataFrame(list(relative_freq_zygosity_A.items()), columns = ['zygosity', 'relative_frequency'])
category_A = ['A', 'A']
df_A['HLA_gene'] = category_A 
df_B = pd.DataFrame(list(relative_freq_zygosity_B.items()), columns = ['zygosity', 'relative_frequency'])
category_B = ['B', 'B']
df_B['HLA_gene'] = category_B 
df_C = pd.DataFrame(list(relative_freq_zygosity_C.items()), columns = ['zygosity', 'relative_frequency'])
category_C = ['C', 'C']
df_C['HLA_gene'] = category_C 


#data from NMDP
obs_freq_A = {'homozygote':0.1521, 'heterozygote': 1-0.1521}
df_A_NMDP = pd.DataFrame(list(obs_freq_A.items()), columns = ['zygosity', 'relative_frequency'])
category_A = ['A\nNMDP', 'A\nNMDP']
df_A_NMDP['HLA_gene'] = category_A

obs_freq_B = {'homozygote':0.1038, 'heterozygote': 1-0.1038}
df_B_NMDP = pd.DataFrame(list(obs_freq_B.items()), columns = ['zygosity', 'relative_frequency'])
category_B = ['B\nNMDP', 'B\nNMDP']
df_B_NMDP['HLA_gene'] = category_B

obs_freq_C = {'homozygote':0.0702, 'heterozygote': 1-0.0702}
df_C_NMDP = pd.DataFrame(list(obs_freq_C.items()), columns = ['zygosity', 'relative_frequency'])
category_C = ['C\nNMDP', 'C\nNMDP']
df_C_NMDP['HLA_gene'] = category_C


#background 'total'
df_A["total"] = df_A['relative_frequency'][0] + df_A['relative_frequency'][1]
df_A.drop(0, inplace=True)
df_B["total"] = df_B['relative_frequency'][0] + df_B['relative_frequency'][1]
df_B.drop(0, inplace=True)
df_C["total"] = df_C['relative_frequency'][0] + df_C['relative_frequency'][1]
df_C.drop(0, inplace=True)
df_A_NMDP["total"] = df_A_NMDP['relative_frequency'][0] + df_A_NMDP['relative_frequency'][1]
df_A_NMDP.drop(0, inplace=True)
df_B_NMDP["total"] = df_B_NMDP['relative_frequency'][0] + df_B_NMDP['relative_frequency'][1]
df_B_NMDP.drop(0, inplace=True)
df_C_NMDP["total"] = df_C_NMDP['relative_frequency'][0] + df_C_NMDP['relative_frequency'][1]
df_C_NMDP.drop(0, inplace=True)

df_c1 = df_A.append(df_A_NMDP, ignore_index=True)\
.append(df_B, ignore_index=True).append(df_B_NMDP, ignore_index=True)\
.append(df_C, ignore_index=True).append(df_C_NMDP, ignore_index=True)


#zygosity analysis class II
zygote_mhc_profile_c2 = {}
for k, v in HLA_profile_c2.items():   
    zygote_mhc_profile_c2[k]  = {}
    for k2, v2 in v.items():
        if v2[0][0] == v2[1][0]:
            zygote_mhc_profile_c2[k][k2] = 'homozygote'
        elif v2[1] == 'No allele':
            zygote_mhc_profile_c2[k][k2] = 'one_allele'
        else:
            zygote_mhc_profile_c2[k][k2] = 'heterozygote'


counts_zygote_DRB1 = {}
counts_zygote_DQB1 = {}
counts_zygote_DPB1 = {}
for k, v in zygote_mhc_profile_c2.items():
    for k2, v2 in v.items():
        if k2 == 'HLA-DRB1':
            if v2 not in counts_zygote_DRB1.keys():
                counts_zygote_DRB1[v2] = 1
            else:
                counts_zygote_DRB1[v2] += 1
        if k2 == 'HLA-DQB1':
            if v2 not in counts_zygote_DQB1.keys():
                counts_zygote_DQB1[v2] = 1
            else:
                counts_zygote_DQB1[v2] += 1
        if k2 == 'HLA-DPB1':
            if v2 not in counts_zygote_DPB1.keys():
                counts_zygote_DPB1[v2] = 1
            else:
                counts_zygote_DPB1[v2] += 1

#DRB1
relative_freq_zygosity_DRB1 = {}
totalcounts = counts_zygote_DRB1['homozygote'] +counts_zygote_DRB1['heterozygote']
relative_freq_zygosity_DRB1['homozygote'] = counts_zygote_DRB1['homozygote']/totalcounts
relative_freq_zygosity_DRB1['heterozygote'] = counts_zygote_DRB1['heterozygote']/totalcounts
#DQB1
relative_freq_zygosity_DQB1 = {}
totalcounts = counts_zygote_DQB1['homozygote'] + counts_zygote_DQB1['heterozygote']
relative_freq_zygosity_DQB1['homozygote'] = counts_zygote_DQB1['homozygote']/totalcounts
relative_freq_zygosity_DQB1['heterozygote'] = counts_zygote_DQB1['heterozygote']/totalcounts
#DPB1
relative_freq_zygosity_DPB1 = {}
totalcounts = counts_zygote_DPB1['homozygote'] + counts_zygote_DPB1['heterozygote']
relative_freq_zygosity_DPB1['homozygote'] = counts_zygote_DPB1['homozygote']/totalcounts
relative_freq_zygosity_DPB1['heterozygote'] = counts_zygote_DPB1['heterozygote']/totalcounts


df_DRB1 = pd.DataFrame(list(relative_freq_zygosity_DRB1.items()), columns = ['zygosity', 'relative_frequency'])
category_A = ['DRB1', 'DRB1']
df_DRB1['HLA_gene'] = category_A 
df_DQB1 = pd.DataFrame(list(relative_freq_zygosity_DQB1.items()), columns = ['zygosity', 'relative_frequency'])
category_B = ['DQB1', 'DQB1']
df_DQB1['HLA_gene'] = category_B 
df_DPB1 = pd.DataFrame(list(relative_freq_zygosity_DPB1.items()), columns = ['zygosity', 'relative_frequency'])
category_C = ['DPB1', 'DPB1']
df_DPB1['HLA_gene'] = category_C 

#data NMDP HLA class II
#source https://bioinformatics.bethematchclinical.org/workarea/downloadasset.aspx?id=6406       03-09-2019
obs_freq_DRB1 = {'homozygote':0.0914 , 'heterozygote': 1-0.0914}
df_DRB1_NMDP = pd.DataFrame(list(obs_freq_DRB1.items()), columns = ['zygosity', 'relative_frequency'])
category_A = ['DRB1\nNMDP', 'DRB1\nNMDP']
df_DRB1_NMDP['HLA_gene'] = category_A

obs_freq_DQB1 = {'homozygote':0.1405, 'heterozygote': 1-0.1405}
df_DQB1_NMDP = pd.DataFrame(list(obs_freq_DQB1.items()), columns = ['zygosity', 'relative_frequency'])
category_B = ['DQB1\nNMDP', 'DQB1\nNMDP']
df_DQB1_NMDP['HLA_gene'] = category_B

#background 'total'
df_DRB1["total"] = df_DRB1['relative_frequency'][0] + df_DRB1['relative_frequency'][1]
df_DRB1.drop(0, inplace=True)
df_DQB1["total"] = df_DQB1['relative_frequency'][0] + df_DQB1['relative_frequency'][1]
df_DQB1.drop(0, inplace=True)
df_DPB1["total"] = df_DPB1['relative_frequency'][0] + df_DPB1['relative_frequency'][1]
df_DPB1.drop(0, inplace=True)
df_DRB1_NMDP["total"] = df_DRB1_NMDP['relative_frequency'][0] + df_DRB1_NMDP['relative_frequency'][1]
df_DRB1_NMDP.drop(0, inplace=True)
df_DQB1_NMDP["total"] = df_DQB1_NMDP['relative_frequency'][0] + df_DQB1_NMDP['relative_frequency'][1]
df_DQB1_NMDP.drop(0, inplace=True)

#supplement plot
df_c2 = df_DRB1.append(df_DRB1_NMDP, ignore_index=True).append(df_DQB1, ignore_index=True)\
.append(df_DQB1_NMDP, ignore_index=True).append(df_DPB1, ignore_index=True)

df_c1_c2 = df_c1.append(df_c2, ignore_index = True)

df_ABC_DR_DQ_DP = df_A.append(df_B, ignore_index=True).append(df_C, ignore_index=True)\
    .append(df_DRB1, ignore_index=True).append(df_DQB1, ignore_index=True)\
    .append(df_DPB1, ignore_index=True)
df_ABC_DR_DQ_DP_NMDP = df_A_NMDP.append(df_B_NMDP, ignore_index=True)\
    .append(df_C_NMDP, ignore_index=True).append(df_DRB1_NMDP, ignore_index=True)\
    .append(df_DQB1_NMDP, ignore_index=True)

'''
Supplementary Figure 4F 'HLA frequencies compared between current study and 
National Marrow Donor Program (NMDP).'
'''

plt.figure(figsize=(10,8), dpi=100)
sns.set_context({"figure.figsize": (6, 8)})
x_pos = [0,3,6,9,12,15] 
plt.bar(x_pos, df_ABC_DR_DQ_DP.total, label='homozygote', color = 'skyblue')
bottom_plot = plt.bar(x_pos, df_ABC_DR_DQ_DP.relative_frequency, label = 'heterozygote', \
                      color = 'darkblue')
x_pos = [1,4,7,10,13] 
plt.bar(x_pos, df_ABC_DR_DQ_DP_NMDP.total, label='homozygote', color = 'navajowhite')
bottom_plot = plt.bar(x_pos, df_ABC_DR_DQ_DP_NMDP.relative_frequency, label = 'heterozygote', \
                      color = 'darkorange')
topbar = plt.Rectangle((0,0),1,1,fc="orange", edgecolor = 'none')
bottombar = plt.Rectangle((0,0),1,1,fc='blue',  edgecolor = 'none')

plt.legend(bbox_to_anchor=(1, 1.05), frameon=False)
x_pos = [0,1,3,4,6,7,9,10,12,13,15] 
plt.xticks(x_pos, df_c1_c2.HLA_gene) 
plt.xlabel("HLA genes")
plt.ylabel("Population frequency", fontsize=14)
plt.savefig('/home/stijn/figures/figures_03_12_19_svg/HLA_hzhz_freq_NMDP.svg')
plt.show()

'''
Figure 1B 'Overview of the dataset content per group selected on having a 
complete profile for the HLA class I or II genes.'
'''
df_c2 = df_A.append(df_B, ignore_index=True).append(df_C, ignore_index=True)\
.append(df_DRB1, ignore_index=True).append(df_DQB1, ignore_index=True)\
.append(df_DPB1, ignore_index=True)
sns.set_style("white")
sns.set_context({"figure.figsize": (6, 8)})

sns.barplot(x = df_c2.HLA_gene, y=df_c2.total , color = "salmon")
bottom_plot = sns.barplot(x = df_c2.HLA_gene, y=df_c2.relative_frequency, color = "teal")
topbar = plt.Rectangle((0,0),1,1,fc="salmon", edgecolor = 'none')
bottombar = plt.Rectangle((0,0),1,1,fc='teal',  edgecolor = 'none')
l = plt.legend([bottombar, topbar], ['heterozygote', 'homozygote'], loc='upper right', bbox_to_anchor=(1, 1.02), ncol = 2, prop={'size':14})
l.draw_frame(False)
sns.despine(left=True)
#bottom_plot.set_xlabel("Classical HLA genes")

for item in ([bottom_plot.xaxis.label, bottom_plot.yaxis.label] +
             bottom_plot.get_xticklabels() + bottom_plot.get_yticklabels()):
    item.set_fontsize(16)

plt.ylabel('Population frequency', fontsize=16)
plt.savefig('/home/stijn/figures/figures_03_12_19_svg/HLA_c1_c2_htrz_mhzgt_freq.svg')
plt.show()

'''
Data preperation for Figure 1A 'Overview of the dataset content per group 
selected on having a complete profile for the HLA class I or II genes.'
'''

reads_HLA_A = []
reads_HLA_B = []
reads_HLA_C = []
nr_A_genes = 0
nr_B_genes = 0
nr_C_genes = 0
for k, v in HLA_profile_c1.items():
    for k2, v2 in v.items():
        if k2 == 'HLA-A':
            reads_A = 0
            for i in v2:
                if len(i) == 2:
                    reads_A += i[1]
                    nr_A_genes += 1
            reads_HLA_A.append(reads_A)
        if k2 == 'HLA-B':
            reads_B = 0
            for i in v2:
                if len(i) == 2:
                    reads_B += i[1]
                    nr_B_genes += 1
            reads_HLA_B.append(reads_B)
        if k2 == 'HLA-C':
            reads_C = 0
            for i in v2:
                if len(i) == 2:
                    reads_C += i[1]
                    nr_C_genes += 1
            reads_HLA_C.append(reads_C)

reads_HLA_DPB1 = []
reads_HLA_DRB1 = []
reads_HLA_DQB1 = []
nr_DPB1_genes = 0
nr_DRB1_genes = 0
nr_DQB1_genes = 0
for k, v in HLA_profile_c2.items():
    for k2, v2 in v.items():
        if k2 == 'HLA-DPB1':
            reads_DPB1 = 0
            for i in v2:
                if len(i) == 2:
                    reads_DPB1 += i[1]
                    nr_DPB1_genes += 1
            reads_HLA_DPB1.append(reads_DPB1)
        if k2 == 'HLA-DRB1':
            reads_DRB1 = 0
            for i in v2:
                if len(i) == 2:
                    reads_DRB1 += i[1]
                    nr_DRB1_genes += 1
            reads_HLA_DRB1.append(reads_DRB1)
        if k2 == 'HLA-DQB1':
            reads_DQB1 = 0
            for i in v2:
                if len(i) == 2:
                    reads_DQB1 += i[1]
                    nr_DQB1_genes += 1
            reads_HLA_DQB1.append(reads_DQB1)
'''
Figure 1A 'Overview of the dataset content per group selected on having a 
complete profile for the HLA class I or II genes.'
'''
plt.figure(figsize=(10,8))
ax1 = plt.subplot(2, 2, 1)
plt.title('HLA class I \n(N=' + str(len(HLA_profile_c1.keys())) + ' datasets)', fontsize=16)
plt.boxplot(reads_HLA_A, positions=np.array([0.5]), showfliers=True)
plt.boxplot(reads_HLA_B, positions=np.array([1]), showfliers=True)
plt.boxplot(reads_HLA_C, positions=np.array([1.5]), showfliers=True)
plt.xlim(0,2)
plt.ylim(100,5000)
ax1.spines['bottom'].set_visible(False)
ax1.xaxis.tick_top()
ax1.tick_params(axis = 'x', labeltop='off', top=False, labelcolor = 'white') 
ax2 = plt.subplot(2, 2, 3)
plt.boxplot(reads_HLA_A, positions=np.array([0.5]), showfliers=True)
plt.boxplot(reads_HLA_B, positions=np.array([1]), showfliers=True)
plt.boxplot(reads_HLA_C, positions=np.array([1.5]), showfliers=True)
plt.xlim(0,2)
plt.ylim(0,100)
ax2.spines['top'].set_visible(False)
ax2.xaxis.tick_bottom()
plt.ylabel('                  reads per dataset', fontsize=16)
# plot the dashed lines for discontinious axes
d = .015  # how big to make the diagonal lines in axes coordinates
# arguments to pass to plot, just so we don't keep repeating them
kwargs = dict(transform=ax1.transAxes, color='k', clip_on=False)
ax1.plot((-d, +d), (-d, +d), **kwargs)        # top-left diagonal
ax1.plot((1 - d, 1 + d), (-d, +d), **kwargs)  # top-right diagonal

kwargs.update(transform=ax2.transAxes)  # switch to the bottom axes
ax2.plot((-d, +d), (1 - d, 1 + d), **kwargs)  # bottom-left diagonal
ax2.plot((1 - d, 1 + d), (1 - d, 1 + d), **kwargs)  # bottom-right diagonal
ticks = ['HLA-A','HLA-B','HLA-C']
plt.xticks([0.5, 1, 1.5], ticks, fontsize=13)
# FOR HLA CLASS II
ax1 = plt.subplot(2, 2, 2)
plt.title('HLA class II \n(N=' + str(len(HLA_profile_c2.keys())) + ' datasets)', fontsize=16)
plt.boxplot(reads_HLA_DPB1, positions=np.array([0.5]), showfliers=True)
plt.boxplot(reads_HLA_DRB1, positions=np.array([1]), showfliers=True)
plt.boxplot(reads_HLA_DQB1, positions=np.array([1.5]), showfliers=True)
plt.xlim(0,2)
plt.ylim(100,7000)

ax1.spines['bottom'].set_visible(False)
ax1.xaxis.tick_top()
ax1.tick_params(axis = 'x', labeltop='off', top=False, labelcolor = 'white') 
ax2 = plt.subplot(2, 2, 4)
plt.boxplot(reads_HLA_DPB1, positions=np.array([0.5]), showfliers=True)
plt.boxplot(reads_HLA_DRB1, positions=np.array([1]), showfliers=True)
plt.boxplot(reads_HLA_DQB1, positions=np.array([1.5]), showfliers=True)
plt.xlim(0,2)
plt.ylim(0,100)
ax2.spines['top'].set_visible(False)
ax2.xaxis.tick_bottom()
# plot the dashed lines for discontinious axes
d = .015  # how big to make the diagonal lines in axes coordinates
# arguments to pass to plot, just so we don't keep repeating them
kwargs = dict(transform=ax1.transAxes, color='k', clip_on=False)
ax1.plot((-d, +d), (-d, +d), **kwargs)        # top-left diagonal
ax1.plot((1 - d, 1 + d), (-d, +d), **kwargs)  # top-right diagonal

kwargs.update(transform=ax2.transAxes)  # switch to the bottom axes
ax2.plot((-d, +d), (1 - d, 1 + d), **kwargs)  # bottom-left diagonal
ax2.plot((1 - d, 1 + d), (1 - d, 1 + d), **kwargs)  # bottom-right diagonal
plt.tight_layout()
ticks = ['HLA-DPB1','HLA-DRB1','HLA-DQB1']
plt.xticks([0.5, 1, 1.5], ticks, fontsize=13)
plt.savefig('/home/stijn/figures/figures_03_12_19_svg/mapping_HLA_reads.svg')
plt.show()


'''
Data preparation for Supplementary Figure 3
'''


def total_reads_mapping(sample_dict):
    dict_total_HLA_p_sample = {}
    for key, value in sample_dict.items():
        dict_total_HLA_p_sample[key] = 0
        for locations, reads in value.items():
            dict_total_HLA_p_sample[key] += reads    
    return dict_total_HLA_p_sample

t_reads_mapping_hla_dict = total_reads_mapping(hla_dict)


samples_more_0 = 0
samples_more_1 = 0
samples_more_50 = 0
samples_more_100 = 0
samples_more_150 = 0
samples_more_200 = 0
samples_more_250 = 0
samples_more_300 = 0
samples_more_350 = 0
samples_more_400 = 0
samples_more_450 = 0
for k, v in t_reads_mapping_hla_dict.items():
    if v >= 0:
        samples_more_0 += 1
    if v >= 1:
        samples_more_1 += 1
    if v >= 50:
        samples_more_50 += 1
    if v >= 100:
        samples_more_100 += 1
    if v >= 150:
        samples_more_150 += 1
    if v >= 200:
        samples_more_200 += 1
    if v >= 250:
        samples_more_250 += 1
    if v >= 300:
        samples_more_300 += 1
    if v >= 350:
        samples_more_350 += 1
    if v >= 400:
        samples_more_400 += 1
    if v >= 450:
        samples_more_450 += 1

objects = ('all datasets','HLA reads \n>= 1', 'HLA reads \n>= 50', 'HLA reads \n>= 100', 'HLA reads \n>= 150', 'HLA reads \n>= 200',\
           'HLA reads \n>= 250', 'HLA reads \n>= 300', 'HLA reads \n>= 350', 'HLA reads \n>= 400',\
           'HLA reads \n>= 450')
y_pos = np.arange(len(objects))
number_of_samples = [samples_more_0, samples_more_1, samples_more_50, samples_more_100, \
                     samples_more_150, samples_more_200, samples_more_250, \
                     samples_more_300, samples_more_350, samples_more_400, \
                     samples_more_450]
'''
Supplementary Figure 3 'Reads mapping to HLA genes in all samples'
'''


plt.figure(figsize=(10,8))
plt.barh(y_pos, number_of_samples, align='center', alpha=0.5)
plt.yticks(y_pos, objects)
plt.xlabel('Number of datasets left')
plt.title('datasets left over on strictness of HLA read count')
plt.savefig('/home/stijn/figures/figures_03_12_19_svg/samples_left_HLA_reads.svg')
plt.show()


'''
(3) Microbiome preparation; extract bacterial reads and filter datasets
''' 


#translating SILVA accession numbers to SILVA taxonomy ID's    
def importSilvaAcc2taxid(file_location):
    silvaAcc2taxid = {}
    with open(file_location, "r") as taxonomy_file:
        for line in taxonomy_file:
            line = line.split('\t')
            acc = line[0]     #remove .split('.')[0] with real data
            taxid = line[1].strip('\n')
            silvaAcc2taxid[acc] = taxid
    return silvaAcc2taxid
file_location = "/home/stijn/stijn2/9_python_files/silva_taxonomy/tax_slv_ssu_132.acc_taxid"
silvaAcc2taxid = importSilvaAcc2taxid(file_location)


#take only the reads mapping to the SILVA SSU database
def only_SILVA_dict(sample_dict):
    SILVA_all_sample_dict = {}
    for k, v in sample_dict.items():
        SILVA_all_sample_dict[k] = {}
        for location, reads in v.items():
            if not location.startswith('HLA'):
                SILVA_all_sample_dict[k][location] = reads
    return SILVA_all_sample_dict
SILVA_all_sample_dd = only_SILVA_dict(all_sample_dict)


def importPrimaryAcc2Domain(file_location):
    silvaPrimaryAcc2Domain = {}
    with open(file_location, "r") as taxonomy_file:
        for line in taxonomy_file:
            line = line.split('\t')
            path = line[3]     #remove .split('.')[0] with real data
            path = path.split(';')
            path = list(filter(None, path))
            #length_path = len(path)
            PrimaryAcc = line[0]
            domain = path[:1]
            silvaPrimaryAcc2Domain[PrimaryAcc] = domain
    return silvaPrimaryAcc2Domain

#extrac bacterial SILVA reads
file_location = "/home/stijn/stijn2/9_python_files/silva_taxonomy/taxmap_slv_ssu_ref_132-corrected.txt"
silvaPrimaryAcc2Domain = importPrimaryAcc2Domain(file_location)


#take only the reads mapping to bacterial SSU in the SILVA database
def only_bacteria_dict(SILVA_all_sample_dd, silvaPrimaryAcc2Domain):
    bact_only_dd = {}
    for k, v in SILVA_all_sample_dd.items():
        bact_only_dd[k] = {}
        for k2, v2 in v.items():
            PrimaryAcc = k2.split('.')
            PrimaryAcc = PrimaryAcc[0]
            if silvaPrimaryAcc2Domain[PrimaryAcc] == ['Bacteria']:
                bact_only_dd[k][k2] = v2
    return bact_only_dd
bact_only_dd = only_bacteria_dict(SILVA_all_sample_dd, silvaPrimaryAcc2Domain)


# selecting samples on bacterial reads and fraction per species. 
# Taking mapping locations with relative frequency of 0.0001 and sample with >1000 reads.
# min_frac = 0.0001
# minimal_reads_sample = 1000
def sample_selection(bact_all_sample_dd, min_frac, minimal_reads_sample):  
    mapped_reads  = total_reads_mapping(bact_all_sample_dd)
    bact_restricted_dict = {}
    for k, v in bact_all_sample_dd.items():
        total_reads = 0
        bact_restricted_dict[k] = {}
        for location, reads in v.items():
            if reads > mapped_reads[k]*min_frac:    #or use 0.001
                total_reads += reads
                bact_restricted_dict[k][location] = reads
        if total_reads < minimal_reads_sample:
            del bact_restricted_dict[k]
    return (bact_restricted_dict)
#each bacteria at least 0.0005 fraction of the reads and each dataset at least 1000 bacterial reads
bact_only_dd_selected = sample_selection(bact_only_dd, 0.0001, 1000) 


#Calculate alpha diversity on genus level. 
#create dictionary for SILVA accession to genus name translation.
def importPrimaryAcc2Genus(file_location):
    silvaPrimaryAcc2Genus = {}
    with open(file_location, "r") as taxonomy_file:
        for line in taxonomy_file:
            line = line.split('\t')
            path = line[3]     #remove .split('.')[0] with real data
            path = path.split(';')
            path = list(filter(None, path))
            #length_path = len(path)
            PrimaryAcc = line[0]
            genus = path[-1:]
            for word in genus:
                if not word[0].isupper() & word[1].islower():
                    genus = path[-2:]
                    genus = '_'.join(genus)
                    genus = genus.replace(' ', '_')
                elif word.startswith('Subgroup'):
                    genus = path[-2:]
                    genus = '_'.join(genus)
                    genus = genus.replace(' ', '_')
                else:
                    genus = '_'.join(genus)
            if genus[-1].isdigit() and genus[-2] == ' ':
                genus = genus[:-2]
            genus = genus.replace(' ', '_')
            #print(PrimaryAcc, length_path, genus)
            silvaPrimaryAcc2Genus[PrimaryAcc] = genus
    return silvaPrimaryAcc2Genus
file_location = "/home/stijn/stijn2/9_python_files/silva_taxonomy/taxmap_slv_ssu_ref_132-corrected.txt"
silvaPrimaryAcc2Genus = importPrimaryAcc2Genus(file_location)


#translates the SILVA accession number to the genus name
def Acc_translated_to_genus_dict(bact_all_sample_dd, silvaPrimaryAcc2Genus):
    bact_genus_all_sample_dd = {}
    for k, v in bact_all_sample_dd.items():
        bact_genus_all_sample_dd[k] = {}
        for k2, v2 in v.items():
            PrimaryAcc = k2.split('.')
            PrimaryAcc = PrimaryAcc[0]
            genus = silvaPrimaryAcc2Genus[PrimaryAcc]
            #genus = ' '.join(genus[:])
            if genus not in bact_genus_all_sample_dd[k].keys():
                bact_genus_all_sample_dd[k][genus] = v2
            else:
                bact_genus_all_sample_dd[k][genus] += v2
    return(bact_genus_all_sample_dd)
bact_only_dd_selected_genus = Acc_translated_to_genus_dict(bact_only_dd_selected, silvaPrimaryAcc2Genus)


'''
part (4) Alpha/within Sample diversity compared to zygosity analysis 
'''

#zygote_mhc_profile_c1 was produced earlier
#only take complete HLA haplotypes
complete_hla_profile = {}
for k, v in zygote_mhc_profile_c1.items():
    if 'one_allele' not in v.values():
        complete_hla_profile[k] = v


#from collections import Counter
htrzyg_c1 = {}
for k, v in complete_hla_profile.items():
    htrzygt_count = list(v.values()).count('heterozygote')
    if htrzygt_count == 0:
        htrzyg_c1[k] = 0
    if htrzygt_count == 1:
        htrzyg_c1[k] = 1
    if htrzygt_count == 2:
        htrzyg_c1[k] = 2
    if htrzygt_count == 3:
        htrzyg_c1[k] = 3


data = {}
data['Richness'] = []
data['Evenness'] = []
heterogenes = []
for k, v in htrzyg_c1.items():
    richness = alpha.margalef(list(bact_only_dd_selected_genus[k].values()))
    data['Richness'].append(richness)
    evenness = alpha.heip_e(list(bact_only_dd_selected_genus[k].values()))
    data['Evenness'].append(evenness)
    string = 'hetero- \nzygote\n gene(s) ' + str(v)
    heterogenes.append(string)
df = pd.DataFrame.from_dict(data)

df['heterogenes'] = pd.Series(heterogenes)
df0 = df[df['heterogenes'] == 'hetero- \nzygote\n gene(s) 0']
df1 = df[df['heterogenes'] == 'hetero- \nzygote\n gene(s) 1']
df2 = df[df['heterogenes'] == 'hetero- \nzygote\n gene(s) 2']
df3 = df[df['heterogenes'] == 'hetero- \nzygote\n gene(s) 3']


# nice significance bar for figures.
#from https://stackoverflow.com/questions/33873176/how-to-add-significance-levels-on-bar-graph-using-pythons-matplotlib
def significance_bar(start,end,height,displaystring,linewidth = 1.2,markersize = 8,boxpad  =0.3,fontsize = 8,color = 'k'):
    # draw a line with downticks at the ends
    plt.plot([start,end],[height]*2,'-',color = color,lw=linewidth,marker = TICKDOWN,markeredgewidth=linewidth,markersize = markersize)
    # draw the text with a bounding box covering up the line
    plt.text(0.5*(start+end),height,displaystring,ha = 'center',va='center',bbox=dict(facecolor='1.', edgecolor='none',boxstyle='Square,pad='+str(boxpad)),size = fontsize)


def check_pval(p):
    if p>=0.05:
        displaystring = r'n.s.'
    elif p<0.0001:
        displaystring = r'***'
    elif p<0.001:
        displaystring = r'**'
    else:
        displaystring = r'*'
    return displaystring

# nice heatmap for statistics
# source of heatmap: https://matplotlib.org/3.1.1/gallery/images_contours_and_fields/image_annotated_heatmap.html
def heatmap(data, row_labels, col_labels, ax=None,
            cbar_kw={}, cbarlabel="", **kwargs):
    """
    Create a heatmap from a numpy array and two lists of labels.
    """

    if not ax:
        ax = plt.gca()

    # Plot the heatmap
    im = ax.imshow(data, **kwargs)

    # Create colorbar
    cbar = ax.figure.colorbar(im, ax=ax, **cbar_kw)
    cbar.ax.set_ylabel(cbarlabel, rotation=-90, va="bottom")

    # We want to show all ticks...
    ax.set_xticks(np.arange(data.shape[1]))
    ax.set_yticks(np.arange(data.shape[0]))
    # ... and label them with the respective list entries.
    ax.set_xticklabels(col_labels)
    ax.set_yticklabels(row_labels)

    # Let the horizontal axes labeling appear on top.
    ax.tick_params(top=True, bottom=False,
                   labeltop=True, labelbottom=False)

    # Rotate the tick labels and set their alignment.
    plt.setp(ax.get_xticklabels(), rotation=-30,  ha="right",
             rotation_mode="anchor")

    # Turn spines off and create white grid.
    for edge, spine in ax.spines.items():
        spine.set_visible(False)

    ax.set_xticks(np.arange(data.shape[1]+1)-.5, minor=True)
    ax.set_yticks(np.arange(data.shape[0]+1)-.5, minor=True)
    ax.grid(which="minor", color="w", linestyle='-', linewidth=3)
    ax.tick_params(which="minor", bottom=False, left=False)

    return im, cbar

def annotate_heatmap(im, data=None, valfmt="{x:.2f}",
                     textcolors=["black", "white"],
                     threshold=None, **textkw):
    """
    A function to annotate a heatmap.
    """

    if not isinstance(data, (list, np.ndarray)):
        data = im.get_array()

    # Normalize the threshold to the images color range.
    if threshold is not None:
        threshold = im.norm(threshold)
    else:
        threshold = im.norm(data.max())/2.

    # Set default alignment to center, but allow it to be
    # overwritten by textkw.
    kw = dict(horizontalalignment="center",
              verticalalignment="center")
    kw.update(textkw)

    # Get the formatter in case a string is supplied
    if isinstance(valfmt, str):
        valfmt = matplotlib.ticker.StrMethodFormatter(valfmt)

    # Loop over the data and create a `Text` for each "pixel".
    # Change the text's color depending on the data.
    texts = []
    for i in range(data.shape[0]):
        for j in range(data.shape[1]):
            kw.update(color=textcolors[int(im.norm(data[i, j]) > threshold)])
            text = im.axes.text(j, i, valfmt(data[i, j], None), **kw)
            texts.append(text)

    return texts
'''
Figure 2A 'The influence of HLA zygosity on gut microbiome alpha diversity.'
Supplementary Figure 5A
'''


plt.figure(figsize=(10,8), dpi=100)
ticks = ['homozygous','hetero- \nzygote\n gene(s) 1','hetero- \nzygote\n gene(s) 2','hetero- \nzygote\n gene(s) 3']
plt.subplot(211)
plt.boxplot(df0['Richness'], positions=np.array([0.5]))
plt.boxplot(df1['Richness'], positions=np.array([1]))
plt.boxplot(df2['Richness'], positions=np.array([1.5]))
plt.boxplot(df3['Richness'], positions=np.array([2]))
plt.xlim(0,2.4)
plt.ylabel('Richness', fontsize=16)
plt.xticks([])
#significance bars
significance_bar(0.5,1,85, check_pval(stats.ks_2samp(df0['Richness'],df1['Richness'], alternative='less')[1]))
plt.text(2.1, 81, 'P='+ str(stats.ks_2samp(df0['Richness'],df1['Richness'], alternative='less')[1])[0:5])

significance_bar(0.5,1.5,90, check_pval(stats.ks_2samp(df0['Richness'],df2['Richness'], alternative='less')[1]))
plt.text(2.1, 86, 'P='+ str(stats.ks_2samp(df0['Richness'],df2['Richness'], alternative='less')[1])[0:5])

significance_bar(0.5,2,95, check_pval(stats.ks_2samp(df0['Richness'],df3['Richness'], alternative='less')[1]))
plt.text(2.1, 91, 'P='+ str(stats.ks_2samp(df0['Richness'],df3['Richness'], alternative='less')[1])[0:5])

#significance_bar(1,1.5,100, check_pval(stats.ks_2samp(df1['Richness'],df2['Richness'], alternative='less')[1]))
#plt.text(2.1, 96, 'P='+ str(stats.ks_2samp(df1['Richness'],df2['Richness'], alternative='less')[1])[0:5])

#significance_bar(1,2,105, check_pval(stats.ks_2samp(df1['Richness'],df3['Richness'], alternative='less')[1]))
#plt.text(2.1, 101, 'P='+ str(stats.ks_2samp(df1['Richness'],df3['Richness'], alternative='less')[1])[0:5])

plt.subplot(212)
plt.boxplot(df0['Evenness'], positions=np.array([0.5]))
plt.boxplot(df1['Evenness'], positions=np.array([1]))
plt.boxplot(df2['Evenness'], positions=np.array([1.5]))
plt.boxplot(df3['Evenness'], positions=np.array([2]))
plt.ylabel('Evenness', fontsize=16)
plt.xlim(0,2.4)

plt.xticks([0.5, 1, 1.5, 2], ticks, fontsize=16)
plt.text(0.5, 0.59, 'N='+ str(len(df0)))
plt.text(1, 0.59, 'N='+ str(len(df1)))
plt.text(1.5, 0.59, 'N='+ str(len(df2)))
plt.text(2, 0.59, 'N='+ str(len(df3)))

plt.subplot(211)
plt.violinplot(df0['Richness'], positions=np.array([0.5]), showmeans=True)
plt.violinplot(df1['Richness'], positions=np.array([1]), showmeans=True)
plt.violinplot(df2['Richness'], positions=np.array([1.5]), showmeans=True)
plt.violinplot(df3['Richness'], positions=np.array([2]), showmeans=True)
plt.xlim(0,2.4)
plt.xticks([])
plt.subplot(212)
plt.violinplot(df0['Evenness'], positions=np.array([0.5]), showmeans=True)
plt.violinplot(df1['Evenness'], positions=np.array([1]), showmeans=True)
plt.violinplot(df2['Evenness'], positions=np.array([1.5]), showmeans=True)
plt.violinplot(df3['Evenness'], positions=np.array([2]), showmeans=True)
plt.xlim(0,2.4)
plt.savefig('/home/stijn/figures/figures_03_12_19_svg/alpha_htzgt_genes_c1.svg')
plt.show()

'''
Figure 2B 'The influence of HLA zygosity on gut microbiome alpha diversity.'
'''


list_df_A = [df0, df1, df2, df3]
Richness_P = []
Evenness_P = []
for i in range(len(list_df_A)):
    P_values_Richness = []
    P_values_Evenness = []
    for j in range(len(list_df_A)):
        P_values_Richness.append(stats.ks_2samp(list_df_A[i]['Richness'],\
                                       list_df_A[j]['Richness'], alternative='greater')[1])
        P_values_Evenness.append(stats.ks_2samp(list_df_A[i]['Evenness'],\
                                       list_df_A[j]['Evenness'], alternative='greater')[1])
    Richness_P.append(P_values_Richness)
    Evenness_P.append(P_values_Evenness)
fig, ax = plt.subplots()
im, cbar = heatmap(np.array(Richness_P), ticks, ticks, ax=ax,
                   cmap="Reds", cbarlabel="P-value")
texts = annotate_heatmap(im, valfmt="{x:.6f}")

fig.tight_layout()
plt.title('Richness two sample Kolmogorov-Smirnov\n item x-axis < item y-axis')
plt.savefig('/home/stijn/figures/figures_03_12_19_svg/heatmapKS_c1_richness.svg')
plt.show()


'''
Supplementary Figure 5B 'The influence of HLA zygosity on gut microbiome alpha diversity.'
'''


fig, ax = plt.subplots()
im, cbar = heatmap(np.array(Evenness_P), ticks, ticks, ax=ax,
                   cmap="Reds", cbarlabel="P-value")
texts = annotate_heatmap(im, valfmt="{x:.6f}")

fig.tight_layout()
plt.title('Evenness two sample Kolmogorov-Smirnov\n item x-axis < item y-axis')
plt.savefig('/home/stijn/figures/figures_03_12_19_svg/heatmapKS_c1_evenness.svg')
plt.show()


'''
Data processing Figure 3 and Supplementary Figure 6C
'''

#check differance per homozygotic gene
homozygote_genes_dict = {}
for k, v in zygote_mhc_profile_c1.items():              #complete_hla_profile used before
    homozygote_genes_dict[k] = []
    for k2, v2 in v.items():
        if k2 == 'HLA-A' and v2 == 'homozygote':
            homozygote_genes_dict[k].append('HLA-A')
        if k2 == 'HLA-B' and v2 == 'homozygote':
            homozygote_genes_dict[k].append('HLA-B')
        if k2 == 'HLA-C' and v2 == 'homozygote':
            homozygote_genes_dict[k].append('HLA-C')


data = {}
data['Richness'] = []
data['Evenness'] = []
HLA_gene = []
for k, v in homozygote_genes_dict.items():
    richness = alpha.margalef(list(bact_only_dd_selected_genus[k].values()))
    data['Richness'].append(richness)
    evenness = alpha.heip_e(list(bact_only_dd_selected_genus[k].values()))
    data['Evenness'].append(evenness)
    HLA_gene.append(v)
df = pd.DataFrame.from_dict(data)

df['HLA_gene'] = pd.Series(HLA_gene)


A = df.HLA_gene.apply(lambda x: 'HLA-A' in x)
df1 = df[A]
B = df.HLA_gene.apply(lambda x: 'HLA-B' in x)
df2 = df[B]
C = df.HLA_gene.apply(lambda x: 'HLA-C' in x)
df3 = df[C]

'''
Figure 3 and Supplementary Figure 6C
'''
plt.figure(figsize=(10,8), dpi=100)
plt.subplot(211)
plt.boxplot(df1['Richness'], positions=np.array([0.5]))
plt.boxplot(df2['Richness'], positions=np.array([1]))
plt.boxplot(df3['Richness'], positions=np.array([1.5]))
significance_bar(0.5,1,85, check_pval(stats.ks_2samp(df2['Richness'],df1['Richness'], alternative='less')[1]))
significance_bar(1,1.5,95, check_pval(stats.ks_2samp(df2['Richness'],df3['Richness'], alternative='less')[1]))
plt.text(1.6, 81, 'P-value='+ str(stats.ks_2samp(df2['Richness'],df1['Richness'], alternative='less')[1])[0:6])
plt.text(1.6, 93, 'P-value='+ str(stats.ks_2samp(df2['Richness'],df3['Richness'], alternative='less')[1])[0:6])
plt.xlim(0,2)
plt.ylabel('Richness', fontsize=16)

plt.xticks([])
plt.subplot(212)
plt.boxplot(df1['Evenness'], positions=np.array([0.5]))
plt.boxplot(df2['Evenness'], positions=np.array([1]))
plt.boxplot(df3['Evenness'], positions=np.array([1.5]))
significance_bar(0.5,1,0.57, check_pval(stats.ks_2samp(df2['Evenness'],df1['Evenness'], alternative='less')[1]))

significance_bar(1,1.5,0.67, check_pval(stats.ks_2samp(df2['Evenness'],df3['Evenness'], alternative='less')[1]))
plt.text(1.6, 0.56, 'P-value='+ str(stats.ks_2samp(df2['Evenness'],df1['Evenness'], alternative='less')[1])[0:6])
#plt.text(1.6, 0.60, 'P-value='+ str(stats.ks_2samp(df1['Evenness'],df3['Evenness'])[1])[0:6])
plt.text(1.6, 0.64, 'P-value='+ str(stats.ks_2samp(df2['Evenness'],df3['Evenness'], alternative='less')[1])[0:6])
plt.ylabel('Evenness', fontsize=16)
plt.xlim(0,2)
ticks = ['HLA-A','HLA-B','HLA-C']
plt.xticks([0.5, 1, 1.5], ticks, fontsize=16)
plt.text(0.4, 0.75, 'N='+ str(len(df1)))
plt.text(0.9, 0.75, 'N='+ str(len(df2)))
plt.text(1.4, 0.75, 'N='+ str(len(df3)))

plt.subplot(2, 1, 1)
plt.violinplot(df1['Richness'], positions=np.array([0.5]), showmeans=True)
plt.violinplot(df2['Richness'], positions=np.array([1]), showmeans=True)
plt.violinplot(df3['Richness'], positions=np.array([1.5]), showmeans=True)
plt.xlim(0,2)
plt.xticks([])
plt.subplot(2, 1, 2)
plt.violinplot(df1['Evenness'], positions=np.array([0.5]), showmeans=True)
plt.violinplot(df2['Evenness'], positions=np.array([1]), showmeans=True)
plt.violinplot(df3['Evenness'], positions=np.array([1.5]), showmeans=True)
plt.xlim(0,2)
plt.savefig('/home/stijn/figures/figures_03_12_19_svg/alpha_hmzgt_genes_c1.svg')
plt.show()

#text in figure difference between HLA-A and C
stats.ks_2samp(df1['Evenness'],df3['Evenness'])
stats.ks_2samp(df1['Richness'],df3['Richness'])

'''
Data processing Supplementary Figure 5 and Figure 2C
'''

    
zygote_mhc_profile_c2 = {}
for k, v in HLA_profile_c2.items():   #Separate the genes!
    zygote_mhc_profile_c2[k]  = {}
    for k2, v2 in v.items():
        if v2[0][0] == v2[1][0]:
            zygote_mhc_profile_c2[k][k2] = 'homozygote'
        elif v2[1] == 'No allele':
            zygote_mhc_profile_c2[k][k2] = 'one_allele'
        else:
            zygote_mhc_profile_c2[k][k2] = 'heterozygote'

complete_hla_profile = {}
for k, v in zygote_mhc_profile_c2.items():
    if 'one_allele' not in v.values():
        complete_hla_profile[k] = v
        

htrzyg_c2 = {}
for k, v in complete_hla_profile.items():
    htrzygt_count = list(v.values()).count('heterozygote')
    if htrzygt_count == 0:
        htrzyg_c2[k] = 0
    if htrzygt_count == 1:
        htrzyg_c2[k] = 1
    if htrzygt_count == 2:
        htrzyg_c2[k] = 2
    if htrzygt_count == 3:
        htrzyg_c2[k] = 3
        
        
data = {}
data['Richness'] = []
data['Evenness'] = []
heterogenes = []
for k, v in htrzyg_c2.items():
    richness = alpha.margalef(list(bact_only_dd_selected_genus[k].values()))
    data['Richness'].append(richness)
    evenness = alpha.heip_e(list(bact_only_dd_selected_genus[k].values()))
    data['Evenness'].append(evenness)
    string = 'hetero- \nzygote\n gene(s) ' + str(v)
    heterogenes.append(string)
df = pd.DataFrame.from_dict(data)


df['heterogenes'] = pd.Series(heterogenes)
df0 = df[df['heterogenes'] == 'hetero- \nzygote\n gene(s) 0']
df1 = df[df['heterogenes'] == 'hetero- \nzygote\n gene(s) 1']
df2 = df[df['heterogenes'] == 'hetero- \nzygote\n gene(s) 2']
df3 = df[df['heterogenes'] == 'hetero- \nzygote\n gene(s) 3']
#dfs1 = pd.melt(df, id_vars = 'HLA_gene')

'''
Supplementary Figure 5C 'Alpha diversity over hetero- & homozygotic HLA class II genes.'
Figure 2C
'''

plt.figure(figsize=(10,8), dpi=100)
plt.subplot(211)
plt.boxplot(df0['Richness'], positions=np.array([0.5]))
plt.boxplot(df1['Richness'], positions=np.array([1]))
plt.boxplot(df2['Richness'], positions=np.array([1.5]))
plt.boxplot(df3['Richness'], positions=np.array([2]))

plt.xlim(0,2.4)
plt.ylabel('Richness', fontsize=16)
plt.title('Alpha diversity over amount of HLA class II heterozygous genes')
plt.xticks([])
plt.subplot(212)
plt.boxplot(df0['Evenness'], positions=np.array([0.5]))
plt.boxplot(df1['Evenness'], positions=np.array([1]))
plt.boxplot(df2['Evenness'], positions=np.array([1.5]))
plt.boxplot(df3['Evenness'], positions=np.array([2]))

plt.ylabel('Evenness', fontsize=16)
plt.xlim(0,2.4)
ticks = ['homozygous','hetero- \nzygote\n gene(s) 1','hetero- \nzygote\n gene(s) 2','hetero- \nzygote\n gene(s) 3']
plt.xticks([0.5, 1, 1.5, 2], ticks, fontsize=16)
plt.text(0.5, 0.75, 'N='+ str(len(df0)))
plt.text(1, 0.75, 'N='+ str(len(df1)))
plt.text(1.5, 0.75, 'N='+ str(len(df2)))
plt.text(2, 0.75, 'N='+ str(len(df3)))

plt.subplot(211)
plt.violinplot(df0['Richness'], positions=np.array([0.5]), showmeans=True)
plt.violinplot(df1['Richness'], positions=np.array([1]), showmeans=True)
plt.violinplot(df2['Richness'], positions=np.array([1.5]), showmeans=True)
plt.violinplot(df3['Richness'], positions=np.array([2]), showmeans=True)
plt.xlim(0,2.4)
plt.xticks([])
plt.subplot(212)
plt.violinplot(df0['Evenness'], positions=np.array([0.5]), showmeans=True)
plt.violinplot(df1['Evenness'], positions=np.array([1]), showmeans=True)
plt.violinplot(df2['Evenness'], positions=np.array([1.5]), showmeans=True)
plt.violinplot(df3['Evenness'], positions=np.array([2]), showmeans=True)
plt.xlim(0,2.4)
plt.savefig('/home/stijn/figures/figures_03_12_19_svg/alpha_htzgt_genes_c2.svg')
plt.show()

'''
Figure 2D and Supplementary Figure 5D
'''
list_df_B = [df0, df1, df2, df3]
Richness_P = []
Evenness_P = []
for i in range(len(list_df_B)):
    P_values_Richness = []
    P_values_Evenness = []
    for j in range(len(list_df_B)):
        P_values_Richness.append(stats.ks_2samp(list_df_B[i]['Richness'],\
                                       list_df_B[j]['Richness'], alternative='greater')[1])
        P_values_Evenness.append(stats.ks_2samp(list_df_B[i]['Evenness'],\
                                       list_df_B[j]['Evenness'], alternative='greater')[1])
    Richness_P.append(P_values_Richness)
    Evenness_P.append(P_values_Evenness)
plt.figure(figsize=(10,8))
fig, ax = plt.subplots()
im, cbar = heatmap(np.array(Richness_P), ticks, ticks, ax=ax,
                   cmap="Reds", cbarlabel="P-value")
texts = annotate_heatmap(im, valfmt="{x:.6f}")

fig.tight_layout()
plt.title('Richness two sample Kolmogorov-Smirnov\n item x-axis < item y-axis')
plt.savefig('/home/stijn/figures/figures_03_12_19_svg/heatmapKS_c2_richness.svg')
plt.show()

plt.figure(figsize=(10,8))
fig, ax = plt.subplots()
im, cbar = heatmap(np.array(Evenness_P), ticks, ticks, ax=ax,
                   cmap="Reds", cbarlabel="P-value")
texts = annotate_heatmap(im, valfmt="{x:.6f}")

fig.tight_layout()
plt.title('Evenness two sample Kolmogorov-Smirnov\n item x-axis < item y-axis')
plt.savefig('/home/stijn/figures/figures_03_12_19_svg/heatmapKS_c2_evenness.svg')
plt.show()

'''
Supplementary Figure 6A-B
'''

#check differance per homozygotic gene class II
homozygote_genes_dictc2 = {}
for k, v in zygote_mhc_profile_c2.items():           #zygote_mhc_profile_c2   or complete_hla_profile
    homozygote_genes_dictc2[k] = []
    for k2, v2 in v.items():
        if k2 == 'HLA-DPB1' and v2 == 'homozygote':
            homozygote_genes_dictc2[k].append('HLA-DPB1')
        if k2 == 'HLA-DRB1' and v2 == 'homozygote':
            homozygote_genes_dictc2[k].append('HLA-DRB1')
        if k2 == 'HLA-DQB1' and v2 == 'homozygote':
            homozygote_genes_dictc2[k].append('HLA-DQB1')
            

data = {}
data['Richness'] = []
data['Evenness'] = []
HLA_gene = []
for k, v in homozygote_genes_dictc2.items():
    richness = alpha.margalef(list(bact_only_dd_selected_genus[k].values()))
    data['Richness'].append(richness)
    evenness = alpha.heip_e(list(bact_only_dd_selected_genus[k].values()))
    data['Evenness'].append(evenness)
    HLA_gene.append(v)
df = pd.DataFrame.from_dict(data)

df['HLA_gene'] = pd.Series(HLA_gene)
DPB1 = df.HLA_gene.apply(lambda x: 'HLA-DPB1' in x)
df1 = df[DPB1]
DRB1 = df.HLA_gene.apply(lambda x: 'HLA-DRB1' in x)
df2 = df[DRB1]
DQB1 = df.HLA_gene.apply(lambda x: 'HLA-DQB1' in x)
df3 = df[DQB1]



'''
Supplementary Figure 6A-B 'Effect of homozygous HLA genes on the alpha diversity of the microbiome'
'''

plt.figure(figsize=(10,8), dpi=100)
plt.subplot(2, 1, 1)
plt.boxplot(df1['Richness'], positions=np.array([0.5]))
plt.boxplot(df2['Richness'], positions=np.array([1]))
plt.boxplot(df3['Richness'], positions=np.array([1.5]))
significance_bar(0.5,1,85, check_pval(stats.ks_2samp(df1['Richness'],df2['Richness'])[1]))
significance_bar(0.5,1.5,90, check_pval(stats.ks_2samp(df1['Richness'],df3['Richness'])[1]))
significance_bar(1,1.5,95, check_pval(stats.ks_2samp(df2['Richness'],df3['Richness'])[1]))
plt.text(1.7, 83, 'P-value='+ str(stats.ks_2samp(df1['Richness'],df2['Richness'])[1])[0:6])
plt.text(1.7, 88, 'P-value='+ str(stats.ks_2samp(df1['Richness'],df3['Richness'])[1])[0:6])
plt.text(1.7, 93, 'P-value='+ str(stats.ks_2samp(df2['Richness'],df3['Richness'])[1])[0:6])
plt.xlim(0,2)
plt.ylabel('Richness', fontsize=16)
#plt.title('Alpha diversity over HLA class II homozygote genes')
plt.xticks([])
plt.subplot(2, 1, 2)
plt.boxplot(df1['Evenness'], positions=np.array([0.5]))
plt.boxplot(df2['Evenness'], positions=np.array([1]))
plt.boxplot(df3['Evenness'], positions=np.array([1.5]))
significance_bar(0.5,1,0.71, check_pval(stats.ks_2samp(df1['Evenness'],df2['Evenness'])[1]))
significance_bar(0.5,1.5,0.77, check_pval(stats.ks_2samp(df1['Evenness'],df3['Evenness'])[1]))
significance_bar(1,1.5,0.83, check_pval(stats.ks_2samp(df2['Evenness'],df3['Evenness'])[1]))
plt.text(1.7, 0.68, 'P-value='+ str(stats.ks_2samp(df1['Evenness'],df2['Evenness'])[1])[0:6])
plt.text(1.7, 0.74, 'P-value='+ str(stats.ks_2samp(df1['Evenness'],df3['Evenness'])[1])[0:6])
plt.text(1.7, 0.80, 'P-value='+ str(stats.ks_2samp(df2['Evenness'],df3['Evenness'])[1])[0:6])
plt.ylabel('Evenness', fontsize=16)
plt.xlim(0,2)
ticks = ['HLA-DQB1','HLA-DRB1','HLA-DPB1']
plt.xticks([0.5, 1, 1.5], ticks, fontsize=16)
plt.text(0.4, 0.95, 'N='+ str(len(df1)))
plt.text(0.9, 0.95, 'N='+ str(len(df2)))
plt.text(1.4, 0.95, 'N='+ str(len(df3)))

plt.subplot(2, 1, 1)
plt.violinplot(df1['Richness'], positions=np.array([0.5]), showmeans=True)
plt.violinplot(df2['Richness'], positions=np.array([1]), showmeans=True)
plt.violinplot(df3['Richness'], positions=np.array([1.5]), showmeans=True)
plt.xlim(0,2)
plt.xticks([])
plt.subplot(2, 1, 2)
plt.violinplot(df1['Evenness'], positions=np.array([0.5]), showmeans=True)
plt.violinplot(df2['Evenness'], positions=np.array([1]), showmeans=True)
plt.violinplot(df3['Evenness'], positions=np.array([1.5]), showmeans=True)
plt.xlim(0,2)
plt.savefig('/home/stijn/figures/figures_03_12_19_svg/alpha_hmzgt_genes_c2.svg')
plt.show()


'''
part (5) Beta diversity/between sample diversity 
Functions for beta diversity
'''


#load SILVA tree for weighted UniFrac
with open("/home/stijn/stijn2/9_python_files/silva_taxonomy/tax_slv_ssu_132.tre", 'r') as newick_silva_tree:
    a = newick_silva_tree.read().replace('\n', '')
    #tree editing branch length is 1.0
    b = a.split(')')
    c = list(':1.0)'.join(b))
    d = ''.join(c)
        
    e = d.split(',')
    f = list(':1.0,'.join(e))
    g = ''.join(f)
    
    tree_str = g.replace(';', ':1.0;')
(Tint, lint, nodes_in_order) = EMDU.parse_tree(tree_str)  


#tree is imported from the SILVA database
with open('/home/stijn/stijn2/9_python_files/tree_str.txt', 'w') as f:
    json.dump(tree_str, f)


#calculate weighted UniFrac distance from dict with SILVA TaxIDs
def calculate_UniFrac(sample1, sample2, sample_dictionary):
    envs = {}
    for k, v in sample_dictionary[sample1].items():
        #taxid = silvaAcc2taxid[k]
        if k not in envs.keys():
            envs[k] = {sample1:v}
        else: 
            envs[k].update({sample1:v})
        
    for k, v in sample_dictionary[sample2].items():
        if k not in envs.keys():
            envs[k] = {sample2:v}
        else: 
            envs[k].update({sample2:v})

    for k, v in envs.items():
        if sample1 not in v:
                envs[k].update({sample1:0})
        if sample2 not in v:
                envs[k].update({sample2:0})       


    (envs_prob_dict, samples) = EMDU.parse_envs(envs, nodes_in_order)  # Parse the environments.
    (Z , F, diffab) = EMDU.EMDUnifrac_weighted_flow(Tint, lint, nodes_in_order, envs_prob_dict[sample1], envs_prob_dict[sample2])  

    return Z


#because the silva taxonomy tree uses taxonomy ID's we need to translate the SILVA accession nr's to taxIDs
def Acc2Taxid(bact_restricted_dict, silvaAcc2taxid):    #is also species level to genus level
    #file_location = "/home/stijn/stijn2/9_python_files/silva_taxonomy/tax_slv_ssu_132.acc_taxid"
    bact_restricted_taxid_dict = {}
    for k, v in bact_restricted_dict.items():       
        bact_restricted_taxid_dict[k] = {}          
        for k2, v2 in v.items():
            taxid = silvaAcc2taxid[k2]
            if taxid not in bact_restricted_taxid_dict[k].keys():
                bact_restricted_taxid_dict[k][taxid] = v2
            else:
                bact_restricted_taxid_dict[k][taxid] += v2
    return bact_restricted_taxid_dict
bact_restricted_taxid_dict = Acc2Taxid(bact_only_dd_selected, silvaAcc2taxid)


#Calculate spearman rank correlation 
def spearmanr(d1, d2):
    d1_z = {key1: d1[key1] for key1 in d1}
    d2_z = {key2: d2[key2] for key2 in d2}
    
    for key2 in d2:
        if key2 not in d1_z:
            d1_z[key2] = 0
    
    for key1 in d1:
        if key1 not in d2_z:
            d2_z[key1] = 0
    list_d1_z = []
    list_d2_z = []
    for k, v in d1_z.items():
        list_d1_z.append(v)
        list_d2_z.append(d2_z[k])
    spearman = stats.spearmanr(list_d1_z, list_d2_z)[0]
    return spearman


#make a list of sample (SRA acc number) pairs
def calculate_list_of_pairs(sample_dict):
    list_of_acc = []
    for i in sample_dict.keys():
        list_of_acc.append(i)
    
    list_of_pairs = []
    for p1 in range(len(list_of_acc)):
        for p2 in range(p1+1, len(list_of_acc)):
            list_of_pairs.append([list_of_acc[p1],list_of_acc[p2]])
    return list_of_pairs
HLA_profile_c1_pairs = calculate_list_of_pairs(HLA_profile_c1)


#make list of alleles per individual
HLA_profile_list = {}
for k, v in HLA_profile_c1.items():
    HLA_profile_list[k] = []
    for k2, v2 in v.items():
        HLA_profile_list[k].append(v2[0][0])
        if not v2[1][0] == 'N':
            HLA_profile_list[k].append(v2[1][0])


#compare lists of alleles between indivduals and calculate the overlap.         
matching_alleles = {}   
for k in HLA_profile_c1_pairs:
    same_items = list(set(HLA_profile_list[k[0]]) & set(HLA_profile_list[k[1]]))
    nmr_alleles = len(same_items)
    k2 = k[0] + '-' + k[1]
    matching_alleles[k2] = nmr_alleles

#make a dict with weighted UniFrac Distance between each sample pair    
UniFrac_c1 = {}
for i in HLA_profile_c1_pairs:
    UniFrac_distance = calculate_UniFrac(i[0], i[1], bact_restricted_taxid_dict)
    key = i[0] + '-' + i[1]
    UniFrac_c1[key] = UniFrac_distance

#make a dict with Microbiome Spearman Rank Correlation between each sample pair
Spearman_c1 = {}
for i in HLA_profile_c1_pairs:
    Spearman_cor = spearmanr(bact_restricted_taxid_dict[i[0]],\
                             bact_restricted_taxid_dict[i[1]])
    key = i[0] + '-' + i[1]
    Spearman_c1[key] = Spearman_cor


data = {}
data['weighted_Unifrac'] = []
data['SpearmanR'] = []
data['Acc_pair'] = []
data['matching_alleles'] = []

for k, v in matching_alleles.items():
    data['weighted_Unifrac'].append(UniFrac_c1[k])
    data['SpearmanR'].append(Spearman_c1[k])
    data['Acc_pair'].append(k)
    data['matching_alleles'].append(v)
df = pd.DataFrame.from_dict(data)


df0 = df[df['matching_alleles'] == 0]
df1 = df[df['matching_alleles'] == 1]
df2 = df[df['matching_alleles'] == 2]
df3 = df[df['matching_alleles'] == 3]
df4 = df[df['matching_alleles'] == 4]
df5 = df[df['matching_alleles'] == 5]
df6 = df[df['matching_alleles'] == 6]

'''
Supplementary Figure 7A 'Microbiome beta diversity for samples with shared HLA class I alleles.'
'''
#SpearmanR figure

plt.figure(figsize=(10,8), dpi=100)
plt.boxplot(df0['SpearmanR'], positions=np.array([0.5]))
plt.boxplot(df1['SpearmanR'], positions=np.array([1]))
plt.boxplot(df2['SpearmanR'], positions=np.array([1.5]))
plt.boxplot(df3['SpearmanR'], positions=np.array([2]))
plt.boxplot(df4['SpearmanR'], positions=np.array([2.5]))
plt.boxplot(df5['SpearmanR'], positions=np.array([3]))
plt.boxplot(df6['SpearmanR'], positions=np.array([3.5]))
plt.xlim(0,3.9)

plt.ylabel('SpearmanR', fontsize=16)
ticks = ['0 alleles \nshared', '1 allele \nshared', '2 alleles \nshared', '3 alleles \nshared',\
         '4 alleles \nshared','5 alleles \nshared', '6 alleles \nshared']
plt.xticks([0.5, 1, 1.5, 2, 2.5, 3, 3.5], ticks, fontsize=16)
plt.text(0.5, 1.05, 'N='+ str(len(df0)))
plt.text(1, 1.05, 'N='+ str(len(df1)))
plt.text(1.5, 1.05, 'N='+ str(len(df2)))
plt.text(2, 1.05, 'N='+ str(len(df3)))
plt.text(2.5, 1.05, 'N='+ str(len(df4)))
plt.text(3, 1.05, 'N='+ str(len(df5)))
plt.text(3.5, 1.05, 'N='+ str(len(df6)))

plt.violinplot(df0['SpearmanR'], positions=np.array([0.5]), showmeans=True)
plt.violinplot(df1['SpearmanR'], positions=np.array([1]), showmeans=True)
plt.violinplot(df2['SpearmanR'], positions=np.array([1.5]), showmeans=True)
plt.violinplot(df3['SpearmanR'], positions=np.array([2]), showmeans=True)
plt.violinplot(df4['SpearmanR'], positions=np.array([2.5]), showmeans=True)
plt.violinplot(df5['SpearmanR'], positions=np.array([3]), showmeans=True)
#plt.violinplot(df6['SpearmanR'], positions=np.array([3.5]), showmeans=True)   plt.violinplot gives error with empty list
plt.xlim(0,3.9)
plt.savefig('/home/stijn/figures/figures_03_12_19_svg/alleles_shared_beta_div_spearman.svg')
plt.show()

'''
Supplementary Figure 7B 'Microbiome beta diversity for samples with shared HLA class I alleles.'

'''
list_df_C = [df0, df1, df2, df3, df4, df5]
Spearman_P = []
for i in range(len(list_df_C)):
    P_values_Spearman = []
    for j in range(len(list_df_C)):
        P_values_Spearman.append(stats.ks_2samp(list_df_C[i]['SpearmanR'],\
                                       list_df_C[j]['SpearmanR'], alternative='greater')[1])
    Spearman_P.append(P_values_Spearman)

plt.figure(figsize=(10,8))
fig, ax = plt.subplots()
im, cbar = heatmap(np.array(Spearman_P), ticks, ticks, ax=ax,
                   cmap="Reds", cbarlabel="P-value")
texts = annotate_heatmap(im, valfmt="{x:.4f}")

fig.tight_layout()
plt.title('Spearman two sample Kolmogorov-Smirnov\n item x-axis < item y-axis')
plt.savefig('/home/stijn/figures/figures_03_12_19_svg/heatmapKS_c1_sharedalles_spearmanR.svg')
plt.show()

'''
Supplementary Figure 7C 'Microbiome beta diversity for samples with shared HLA class I alleles.'
'''

plt.figure(figsize=(10,8), dpi=100)
plt.boxplot(df0['weighted_Unifrac'], positions=np.array([0.5]))
plt.boxplot(df1['weighted_Unifrac'], positions=np.array([1]))
plt.boxplot(df2['weighted_Unifrac'], positions=np.array([1.5]))
plt.boxplot(df3['weighted_Unifrac'], positions=np.array([2]))
plt.boxplot(df4['weighted_Unifrac'], positions=np.array([2.5]))
plt.boxplot(df5['weighted_Unifrac'], positions=np.array([3]))
plt.boxplot(df6['weighted_Unifrac'], positions=np.array([3.5]))
plt.xlim(0,3.9)

plt.ylabel('weighted UniFrac', fontsize=16)
ticks = ['0 alleles \nshared', '1 allele \nshared', '2 alleles \nshared', '3 alleles \nshared',\
         '4 alleles \nshared','5 alleles \nshared', '6 alleles \nshared']
plt.xticks([0.5, 1, 1.5, 2, 2.5, 3, 3.5], ticks, fontsize=16)
plt.text(0.5, 9.5, 'N='+ str(len(df0)))
plt.text(1, 9.5, 'N='+ str(len(df1)))
plt.text(1.5, 9.5, 'N='+ str(len(df2)))
plt.text(2, 9.5, 'N='+ str(len(df3)))
plt.text(2.5, 9.5, 'N='+ str(len(df4)))
plt.text(3, 9.5, 'N='+ str(len(df5)))
plt.text(3.5, 9.5, 'N='+ str(len(df6)))

plt.violinplot(df0['weighted_Unifrac'], positions=np.array([0.5]), showmeans=True)
plt.violinplot(df1['weighted_Unifrac'], positions=np.array([1]), showmeans=True)
plt.violinplot(df2['weighted_Unifrac'], positions=np.array([1.5]), showmeans=True)
plt.violinplot(df3['weighted_Unifrac'], positions=np.array([2]), showmeans=True)
plt.violinplot(df4['weighted_Unifrac'], positions=np.array([2.5]), showmeans=True)
plt.violinplot(df5['weighted_Unifrac'], positions=np.array([3]), showmeans=True)
#plt.violinplot(df6['weighted_Unifrac'], positions=np.array([3.5]), showmeans=True)  plt.violinplot gives error with empty list
plt.xlim(0,3.9)
plt.savefig('/home/stijn/figures/figures_03_12_19_svg/alleles_shared_beta_div_weighted_Unifrac.svg')
plt.show()

'''
Supplementary Figure 7D 'Microbiome beta diversity for samples with shared HLA class I alleles.'
'''

list_df_C = [df0, df1, df2, df3, df4, df5]
Spearman_P = []
for i in range(len(list_df_C)):
    P_values_Spearman = []
    for j in range(len(list_df_C)):
        P_values_Spearman.append(stats.ks_2samp(list_df_C[i]['weighted_Unifrac'],\
                                       list_df_C[j]['weighted_Unifrac'], alternative='greater')[1])
    Spearman_P.append(P_values_Spearman)

plt.figure(figsize=(10,8))
fig, ax = plt.subplots()
im, cbar = heatmap(np.array(Spearman_P), ticks, ticks, ax=ax,
                   cmap="Reds", cbarlabel="P-value")
texts = annotate_heatmap(im, valfmt="{x:.4f}")

fig.tight_layout()
plt.title('weighted UniFrac two sample Kolmogorov-Smirnov\n item x-axis < item y-axis')
plt.savefig('/home/stijn/figures/figures_03_12_19_svg/heatmapKS_c1_sharedalles_weighted_Unifrac.svg')
plt.show()

'''
Data processing for Supplementary Figure 8 Supertypes.
'''

#allele to supertypes translation
supertypes_dict = {}    #origin: https://bmcimmunol.biomedcentral.com/articles/10.1186/1471-2172-9-1#additional-information
with open('/home/stijn/Documents/7_HLA_allelefrequencies/classIsupertypes.txt') as f:
    for line in f:
        line = line.split('\t')
        if not line[1] == 'Unclassified':
            supertypes_dict[line[0]] = line[1]
        
#remove unclassified
supertypes_dict_2d = {}
for k,v in supertypes_dict.items():
    k = k[0:4]
    supertypes_dict_2d[k] = v


supertypes_profile_list = {} #fix HLA-C
for k, v in HLA_profile_c1.items():
    supertypes_profile_list[k] = []
    for k2, v2 in v.items():
        #print(k2)
        if k2.startswith(('HLA-A', 'HLA-B')):
            hla = v2[0][0]
            hla = hla[4::]
            supertype = supertypes_dict_2d[hla]
            supertypes_profile_list[k].append(supertype)
            if not v2[1][0] == 'N':
                hla = v2[1][0]
                hla = hla[4::]
                supertype = supertypes_dict_2d[hla]
                supertypes_profile_list[k].append(supertype)

matching_supertype = {}   
for k in HLA_profile_c1_pairs:
    same_items = list(set(supertypes_profile_list[k[0]]) & \
                      set(supertypes_profile_list[k[1]]))
    nmr_alleles = len(same_items)
    k2 = k[0] + '-' + k[1]
    matching_supertype[k2] = nmr_alleles


data = {}
#data['Aitchison'] = []
data['weighted_Unifrac'] = []
data['SpearmanR'] = []
data['Acc_pair'] = []
data['matching_supertypes'] = []

for k, v in matching_supertype.items():
#    data['Aitchison'].append(classI_aitch_distance_TaxID_dict[k][0])
    data['weighted_Unifrac'].append(UniFrac_c1[k])
    data['SpearmanR'].append(Spearman_c1[k])
    data['Acc_pair'].append(k)
    data['matching_supertypes'].append(v)

df = pd.DataFrame.from_dict(data)

df0 = df[df['matching_supertypes'] == 0]
df1 = df[df['matching_supertypes'] == 1]
df2 = df[df['matching_supertypes'] == 2]
df3 = df[df['matching_supertypes'] == 3]
df4 = df[df['matching_supertypes'] == 4]

'''
Supplementary Figure 8A
'''
#spearmanR
plt.figure(figsize=(10,8), dpi=100)
plt.boxplot(df0['SpearmanR'], positions=np.array([0.5]))
plt.boxplot(df1['SpearmanR'], positions=np.array([1]))
plt.boxplot(df2['SpearmanR'], positions=np.array([1.5]))
plt.boxplot(df3['SpearmanR'], positions=np.array([2]))
plt.boxplot(df4['SpearmanR'], positions=np.array([2.5]))
#plt.title('HLA-A/B supertypes shared')
plt.xlim(0,2.9)
plt.ylabel('SpearmanR')
ticks = ['0 supertypes \nshared', '1 supertypes \nshared', '2 supertypes \nshared', '3 supertypes \nshared',\
         '4 supertypes \nshared']
plt.xticks([0.5, 1, 1.5, 2, 2.5, 3, 3.5], ticks)
plt.text(0.5, 1.05, 'N='+ str(len(df0)))
plt.text(1, 1.05, 'N='+ str(len(df1)))
plt.text(1.5, 1.05, 'N='+ str(len(df2)))
plt.text(2, 1.05, 'N='+ str(len(df3)))
plt.text(2.5, 1.05, 'N='+ str(len(df4)))

plt.violinplot(df0['SpearmanR'], positions=np.array([0.5]), showmeans=True)
plt.violinplot(df1['SpearmanR'], positions=np.array([1]), showmeans=True)
plt.violinplot(df2['SpearmanR'], positions=np.array([1.5]), showmeans=True)
plt.violinplot(df3['SpearmanR'], positions=np.array([2]), showmeans=True)
plt.violinplot(df4['SpearmanR'], positions=np.array([2.5]), showmeans=True)
plt.xlim(0,2.9)
plt.savefig('/home/stijn/figures/figures_03_12_19_svg/supertypes_spearman_c1.svg')
plt.show()

'''
Supplementary Figure 8B
'''
list_df_D = [df0, df1, df2, df3, df4]
Spearman_P = []
for i in range(len(list_df_D)):
    P_values_Spearman = []
    for j in range(len(list_df_D)):
        P_values_Spearman.append(stats.ks_2samp(list_df_D[i]['SpearmanR'],\
                                       list_df_D[j]['SpearmanR'], alternative='greater')[1])
    Spearman_P.append(P_values_Spearman)

fig, ax = plt.subplots()
im, cbar = heatmap(np.array(Spearman_P), ticks, ticks, ax=ax,
                   cmap="Reds", cbarlabel="P-value")
texts = annotate_heatmap(im, valfmt="{x:.4f}")

fig.tight_layout()
plt.title('Spearman two sample Kolmogorov-Smirnov\n item x-axis < item y-axis')
plt.savefig('/home/stijn/figures/figures_03_12_19_svg/heatmapKS_c1_supertypes_spearmanR.svg')
plt.show()

'''
part (6) PPSS
Data processing for PPSS
'''

#lists of alleles found in the dataset and used to predict peptide affinity to.
list_alleles_netMHCpan = []
with open('/home/stijn/stijn2/11_netMHCpan/list_A') as f:
    for line in f:
        line = line.strip('\n')
        list_alleles_netMHCpan.append(line)
with open('/home/stijn/stijn2/11_netMHCpan/list_B') as f:
    for line in f:
        line = line.strip('\n')
        list_alleles_netMHCpan.append(line)
with open('/home/stijn/stijn2/11_netMHCpan/list_C') as f:
    for line in f:
        line = line.strip('\n')
        list_alleles_netMHCpan.append(line)

# create a list of binding peptides for each allele from NetMHCpan4.0 output
Binder_2 = {}
for i in list_alleles_netMHCpan:
    Binder_2[i] = []
    with open('/home/stijn/stijn2/11_netMHCpan/' + i) as f:
        for line in f:
            if line.startswith('    1'):
                line = line.strip('\n').strip(' <= SB').strip(' <= WB')
                line_split = line.split(' ')
                last_item = float((line_split[-1::])[0])
                if last_item <= 2:
                    Binder_2[i].append(line_split[9])

Binder_2_A = {}
Binder_2_B = {}
Binder_2_C = {}
for k, v in Binder_2.items():
    if k.startswith('HLA-A'):
        Binder_2_A[k] = v
    if k.startswith('HLA-B'):
        Binder_2_B[k] = v
    if k.startswith('HLA-C'):
        Binder_2_C[k] = v
    

#Jaccard index between alleles
jaccard_2_A = {}
for i in combinations(Binder_2_A.keys(), 2):
    same_items = list(set(Binder_2_A[i[0]]).intersection(Binder_2_A[i[1]]))
    nr_peptides_shared = len(same_items)
    total_peptides = len(Binder_2_A[i[0]]) + len(Binder_2_A[i[1]])
    jaccard_index = nr_peptides_shared/(total_peptides - nr_peptides_shared)
    jaccard_2_A[(i[0][:7],i[1][:7])] = jaccard_index
    jaccard_2_A[(i[0][:7],i[0][:7])] = 1
    jaccard_2_A[(i[1][:7],i[1][:7])] = 1
jaccard_2_B = {}
for i in combinations(Binder_2_B.keys(), 2):
    same_items = list(set(Binder_2_B[i[0]]).intersection(Binder_2_B[i[1]]))
    nr_peptides_shared = len(same_items)
    total_peptides = len(Binder_2_B[i[0]]) + len(Binder_2_B[i[1]])
    jaccard_index = nr_peptides_shared/(total_peptides - nr_peptides_shared)
    jaccard_2_B[(i[0][:7], i[1][:7])] = jaccard_index
    jaccard_2_B[(i[0][:7],i[0][:7])] = 1
    jaccard_2_B[(i[1][:7],i[1][:7])] = 1
jaccard_2_C = {}
for i in combinations(Binder_2_C.keys(), 2):
    same_items = list(set(Binder_2_C[i[0]]).intersection(Binder_2_C[i[1]]))
    nr_peptides_shared = len(same_items)
    total_peptides = len(Binder_2_C[i[0]]) + len(Binder_2_C[i[1]])
    jaccard_index = nr_peptides_shared/(total_peptides - nr_peptides_shared)
    jaccard_2_C[(i[0][:7], i[1][:7])] = jaccard_index
    jaccard_2_C[(i[0][:7],i[0][:7])] = 1
    jaccard_2_C[(i[1][:7],i[1][:7])] = 1


#PPSS score per gene between individuals.
jaccard_score_A = {} 
jaccard_score_B = {}  
jaccard_score_C = {}    
for k2 in HLA_profile_c1_pairs:
    k = k2[0] + '-' + k2[1]
    jaccard_A = 0
    HLA_A1 = []
    HLA_A2 = []
    for i in HLA_profile_list[k2[0]]:
        if i.startswith('HLA-A'):
            i = i.replace('*', '')
            HLA_A1.append(i)
    for i in HLA_profile_list[k2[1]]:
        if i.startswith('HLA-A'):
            i = i.replace('*', '')
            HLA_A2.append(i)
    for b in list(product(HLA_A1, HLA_A2)):
        sorted_b = (tuple(sorted(b)))
        jaccard_A += jaccard_2_A[sorted_b]
    jaccard_score_A[k] = jaccard_A
    #for B
    jaccard_B = 0
    HLA_B1 = []
    HLA_B2 = []
    for i in HLA_profile_list[k2[0]]:
        if i.startswith('HLA-B'):
            i = i.replace('*', '')
            HLA_B1.append(i)
    for i in HLA_profile_list[k2[1]]:
        if i.startswith('HLA-B'):
            i = i.replace('*', '')
            HLA_B2.append(i)
    for b in list(product(HLA_B1, HLA_B2)):
        sorted_b = (tuple(sorted(b)))
        jaccard_B += jaccard_2_B[sorted_b]
    jaccard_score_B[k] = jaccard_B
    #for C
    jaccard_C = 0
    HLA_C1 = []
    HLA_C2 = []
    for i in HLA_profile_list[k2[0]]:
        if i.startswith('HLA-C'):
            i = i.replace('*', '')
            HLA_C1.append(i)
    for i in HLA_profile_list[k2[1]]:
        if i.startswith('HLA-C'):
            i = i.replace('*', '')
            HLA_C2.append(i)
    for b in list(product(HLA_C1, HLA_C2)):
        sorted_b = (tuple(sorted(b)))
        jaccard_C += jaccard_2_C[sorted_b]
    jaccard_score_C[k] = jaccard_C
  
#PPSS score between individuals.      
jaccard_total_dict = {}
for k in HLA_profile_c1_pairs:
    k2 = k[0] + '-' + k[1]
    jaccard_total = jaccard_score_A[k2] + jaccard_score_B[k2] + jaccard_score_C[k2]
    jaccard_total_dict[k2] = jaccard_total
 

data_jaccard = {}
data_jaccard['weighted_Unifrac'] = []
data_jaccard['SpearmanR'] = []
data_jaccard['Acc_pair'] = []
data_jaccard['jaccard'] = []
for k, v in jaccard_total_dict.items():
    data_jaccard['weighted_Unifrac'].append(UniFrac_c1[k])
    data_jaccard['SpearmanR'].append(Spearman_c1[k])
    data_jaccard['Acc_pair'].append(k)
    data_jaccard['jaccard'].append(v)
df = pd.DataFrame.from_dict(data_jaccard)


'''
Figure 4A 'Individuals presenting similar microbial human gut peptides on their 
HLA have a more similar microbiome'
'''

df1 = df[df['jaccard'] < 4]
df2_select =  df[df['jaccard'] < 6]
df2 = df2_select[df2_select['jaccard'] >= 4]
df3_select =  df[df['jaccard'] < 8]
df3 = df3_select[df3_select['jaccard'] >= 6]
df4 = df[df['jaccard'] >= 8]

#SpearmanR
plt.figure(figsize=(10,8), dpi=100)
#plt.boxplot(df0['SpearmanR'], positions=np.array([0.5]))
plt.boxplot(df1['SpearmanR'], positions=np.array([0.5]))
plt.boxplot(df2['SpearmanR'], positions=np.array([1]))
plt.boxplot(df3['SpearmanR'], positions=np.array([1.5]))
plt.boxplot(df4['SpearmanR'], positions=np.array([2]))
plt.xlim(0,2.4)
plt.ylabel('SpearmanR', fontsize=16)
#plt.violinplot(df0['SpearmanR'], positions=np.array([0.5]), showmeans=True)
plt.violinplot(df1['SpearmanR'], positions=np.array([0.5]), showmeans=True)
plt.violinplot(df2['SpearmanR'], positions=np.array([1]), showmeans=True)
plt.violinplot(df3['SpearmanR'], positions=np.array([1.5]), showmeans=True)
plt.violinplot(df4['SpearmanR'], positions=np.array([2]), showmeans=True)
ticks = ['PPSS < 4', 'PPSS 4-6', 'PPSS 6-8',\
         'PPSS >= 8']
plt.xticks([0.5, 1, 1.5, 2], ticks, fontsize=16)
#plt.text(0.5, 1.05, 'N='+ str(len(df0)))
plt.text(0.5, 1.05, 'N='+ str(len(df1)))
plt.text(1, 1.05, 'N='+ str(len(df2)))
plt.text(1.5, 1.05, 'N='+ str(len(df3)))
plt.text(2, 1.05, 'N='+ str(len(df4)))

plt.savefig('/home/stijn/figures/figures_03_12_19_svg/peptidome_jaccard_spearman.svg')
plt.show()

'''
Figure 4B 'Individuals presenting similar microbial human gut peptides on their 
HLA have a more similar microbiome'
'''
list_df_E = [df1, df2, df3, df4]
Spearman_P = []
for i in range(len(list_df_E)):
    P_values_Spearman = []
    for j in range(len(list_df_E)):
        P_values_Spearman.append(stats.ks_2samp(list_df_E[i]['SpearmanR'],\
                                       list_df_E[j]['SpearmanR'], alternative='greater')[1])
    Spearman_P.append(P_values_Spearman)

fig, ax = plt.subplots()
im, cbar = heatmap(np.array(Spearman_P), ticks, ticks, ax=ax,
                   cmap="Reds", cbarlabel="P-value")
texts = annotate_heatmap(im, valfmt="{x:.6f}")

fig.tight_layout()
plt.title('Spearman two sample Kolmogorov-Smirnov\n item x-axis < item y-axis')
plt.savefig('/home/stijn/figures/figures_03_12_19_svg/heatmapKS_c1_jaccard_spearmanR.svg')
plt.show()

'''
Figure 4C 'Individuals presenting similar microbial human gut peptides on their 
HLA have a more similar microbiome'
'''
#weighted UniFrac
plt.figure(figsize=(10,8), dpi=100)
#plt.boxplot(df0['weighted_Unifrac'], positions=np.array([0.5]))
plt.boxplot(df1['weighted_Unifrac'], positions=np.array([0.5]))
plt.boxplot(df2['weighted_Unifrac'], positions=np.array([1]))
plt.boxplot(df3['weighted_Unifrac'], positions=np.array([1.5]))
plt.boxplot(df4['weighted_Unifrac'], positions=np.array([2]))
plt.xlim(0,2.4)
plt.ylabel('weighted_Unifrac', fontsize=16)
#plt.violinplot(df0['weighted_Unifrac'], positions=np.array([0.5]), showmeans=True)
plt.violinplot(df1['weighted_Unifrac'], positions=np.array([0.5]), showmeans=True)
plt.violinplot(df2['weighted_Unifrac'], positions=np.array([1]), showmeans=True)
plt.violinplot(df3['weighted_Unifrac'], positions=np.array([1.5]), showmeans=True)
plt.violinplot(df4['weighted_Unifrac'], positions=np.array([2]), showmeans=True)
ticks = ['PPSS < 4', 'PPSS 4-6',\
         'PPSS 6-8','PPSS >= 8']
plt.xticks([0.5, 1, 1.5, 2], ticks, fontsize=16)
#plt.text(0.5, 9.5, 'N='+ str(len(df0)))
plt.text(0.5, 9.5, 'N='+ str(len(df1)))
plt.text(1, 9.5, 'N='+ str(len(df2)))
plt.text(1.5, 9.5, 'N='+ str(len(df3)))
plt.text(2, 9.5, 'N='+ str(len(df4)))

plt.savefig('/home/stijn/figures/figures_03_12_19_svg/peptidome_jaccard_weighted_UniFrac.svg')
plt.show()


'''
Figure 4D 'Individuals presenting similar microbial human gut peptides on their 
HLA have a more similar microbiome'
'''

list_df_F = [df1, df2, df3, df4]
Spearman_P = []
for i in range(len(list_df_F)):
    P_values_Spearman = []
    for j in range(len(list_df_F)):
        P_values_Spearman.append(stats.ks_2samp(list_df_F[i]['weighted_Unifrac'],\
                                       list_df_F[j]['weighted_Unifrac'], alternative='greater')[1])
    Spearman_P.append(P_values_Spearman)

fig, ax = plt.subplots()
im, cbar = heatmap(np.array(Spearman_P), ticks, ticks, ax=ax,
                   cmap="Reds", cbarlabel="P-value")
texts = annotate_heatmap(im, valfmt="{x:.6f}")

fig.tight_layout()
plt.title('weighted Unifrac two sample Kolmogorov-Smirnov\n item x-axis < item y-axis')
plt.savefig('/home/stijn/figures/figures_03_12_19_svg/heatmapKS_c1_jaccard_wUniFrac.svg')
plt.show()

'''
Supplementary Figure 1A Jaccard peptide presentation between alleles. HLA-A
'''

# Heatmap ppss scores supplementary figure 7
#HLA-A
jaccard_2_A.keys()
alleles_A = []
for k in jaccard_2_A.keys():
    for i in k:
        if i not in alleles_A:
            alleles_A.append(i)
M_A = []
for i, k1 in enumerate(alleles_A):
    M_A.append([])
    for k2 in alleles_A:
        if (k1,k2) in jaccard_2_A.keys():
            M_A[i].append(jaccard_2_A[(k1,k2)])
        else:
            M_A[i].append(jaccard_2_A[(k2,k1)])
M_A_array = np.array(M_A)

fig, ax = plt.subplots()
im = ax.imshow(M_A_array, cmap='YlGn')
# We want to show all ticks...
ax.set_xticks(np.arange(len(alleles_A)))
ax.set_yticks(np.arange(len(alleles_A)))
# ... and label them with the respective list entries
ax.set_xticklabels(alleles_A)
ax.set_yticklabels(alleles_A)

plt.setp(ax.get_xticklabels(), rotation=45, ha="right",
         rotation_mode="anchor")

# Loop over data dimensions and create text annotations.
for i in range(len(alleles_A)):
    for j in range(len(alleles_A)):
        text = ax.text(j, i, round(M_A_array[i, j],3),
                       ha="center", va="center", color="black")

ax.set_title("PPSS HLA-A alleles")
fig.tight_layout()
fig.set_size_inches(10, 10, forward=True)
plt.savefig('/home/stijn/figures/figures_03_12_19_svg/PPSS_heatmap_HLAA.svg')
plt.show()

plt.colorbar(im)
plt.savefig('/home/stijn/figures/figures_03_12_19_svg/PPSS_colorbar_HLAA.svg')
plt.show()


'''
Supplementary Figure 1B Jaccard peptide presentation between alleles. HLA-B
'''


# HLA-B
alleles_B = []
for k in jaccard_2_B.keys():
    for i in k:
        if i not in alleles_B:
            alleles_B.append(i)
M_B = []
for i, k1 in enumerate(alleles_B):
    M_B.append([])
    for k2 in alleles_B:
        if (k1,k2) in jaccard_2_B.keys():
            M_B[i].append(jaccard_2_B[(k1,k2)])
        else:
            M_B[i].append(jaccard_2_B[(k2,k1)])
M_B_array = np.array(M_B)

fig, ax = plt.subplots()
im = ax.imshow(M_B_array, cmap='YlGn')
# We want to show all ticks...
ax.set_xticks(np.arange(len(alleles_B)))
ax.set_yticks(np.arange(len(alleles_B)))
# ... and label them with the respective list entries
ax.set_xticklabels(alleles_B)
ax.set_yticklabels(alleles_B)

plt.setp(ax.get_xticklabels(), rotation=45, ha="right",
         rotation_mode="anchor")

# Loop over data dimensions and create text annotations.
for i in range(len(alleles_B)):
    for j in range(len(alleles_B)):
        text = ax.text(j, i, round(M_B_array[i, j],3),
                       ha="center", va="center", color="black")

ax.set_title("PPSS HLA-B alleles")
fig.tight_layout()
fig.set_size_inches(10, 10, forward=True)
plt.savefig('/home/stijn/figures/figures_03_12_19_svg/PPSS_heatmap_HLAB.svg')
plt.show()

plt.colorbar(im)
plt.savefig('/home/stijn/figures/figures_03_12_19_svg/PPSS_colorbar_HLAB.svg')
plt.show()


'''
Supplementary Figure 1C Jaccard peptide presentation between alleles. HLA-C
'''


# HLA-C
alleles_C = []
for k in jaccard_2_C.keys():
    for i in k:
        if i not in alleles_C:
            alleles_C.append(i)
M_C = []
for i, k1 in enumerate(alleles_C):
    M_C.append([])
    for k2 in alleles_C:
        if (k1,k2) in jaccard_2_C.keys():
            M_C[i].append(jaccard_2_C[(k1,k2)])
        else:
            M_C[i].append(jaccard_2_C[(k2,k1)])
M_C_array = np.array(M_C)

fig, ax = plt.subplots()
im = ax.imshow(M_C_array, cmap='YlGn')
# We want to show all ticks...
ax.set_xticks(np.arange(len(alleles_C)))
ax.set_yticks(np.arange(len(alleles_C)))
# ... and label them with the respective list entries
ax.set_xticklabels(alleles_C)
ax.set_yticklabels(alleles_C)

plt.setp(ax.get_xticklabels(), rotation=45, ha="right",
         rotation_mode="anchor")

# Loop over data dimensions and create text annotations.
for i in range(len(alleles_C)):
    for j in range(len(alleles_C)):
        text = ax.text(j, i, round(M_C_array[i, j],3),
                       ha="center", va="center", color="black")

ax.set_title("PPSS HLA-C alleles")
fig.tight_layout()
fig.set_size_inches(10, 10, forward=True)
plt.savefig('/home/stijn/figures/figures_03_12_19_svg/PPSS_heatmap_HLAC.svg')
plt.show()

plt.colorbar(im)
plt.savefig('/home/stijn/figures/figures_03_12_19_svg/PPSS_colorbar_HLAC.svg')
plt.show()

'''
Supplementary Figure 2 Histogram of PPSS pairs of individuals.
'''
# Histogram PPSS
bina = []
c = 0
for i in range(0,49):
    bina.append(c)
    c = c + 0.25
n, bins, patches = plt.hist(data_jaccard['jaccard'], bina, density=False, alpha=0.75)
#plt.hist(data_jaccard['jaccard'], bins=bina, density=False, alpha=0.75)
plt.xlabel('PPSS', fontsize=16)
plt.ylabel('Counts', fontsize=16)


for i in range(0,16):
    print(i)
    patches[i].set_facecolor('dodgerblue')
for i in range(16,24):    
    patches[i].set_facecolor('darkorange')
for i in range(24, 32):
    patches[i].set_facecolor('g')
for i in range(32, len(patches)):
    patches[i].set_facecolor('red')
    
plt.grid(True)
plt.xlim(0, 12)
plt.ylim(0,650)
plt.savefig('/home/stijn/figures/figures_03_12_19_svg/PPSS_histogram.svg')
plt.show()



'''
Individual statistics in the main text

'''
#Correlation between richness and 16S libary size for Method statement
list_16S = []
list_richness = []
for k, v in bact_only_dd_selected_genus.items():
    list_16S.append((sum(bact_only_dd_selected_genus[k].values())))
    list_richness.append(alpha.margalef(list(bact_only_dd_selected_genus[k].values())))

stats.pearsonr(list_16S, list_richness)


#test if more reads map to HLA class I than class II.
list_HLA_reads_c1 = []
for k,v in HLA_profile_c1.items():
    for k2, v2 in v.items():
        for i in v2:
            if len(i) == 2:
                list_HLA_reads_c1.append(i[1])

list_HLA_reads_c2 = []
for k,v in HLA_profile_c2.items():
    for k2, v2 in v.items():
        for i in v2:
            if len(i) == 2:
                list_HLA_reads_c2.append(i[1])
stats.ks_2samp(list_HLA_reads_c1, list_HLA_reads_c2, alternative = 'less')