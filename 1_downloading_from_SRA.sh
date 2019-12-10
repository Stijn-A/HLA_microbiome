#!/bin/bash

# Script for downloading datasets from the Sequencing Read Archive (SRA)
# There are multiple ways of downloading; using prefetch, wget or esearch.
# I used prefetch for downloading most of the files and for some wget.


# Using prefetch (the tool provided by the SRA manual)
# prefetch downloads files at standard location ~/ncbi/public/sra/ because this
# is a backed up location on our systems and the files are large I directly move
# them to a non-backed up folder

while read p;
  do
    echo $p
    prefetch.2.9.2 $p
    mv -v ~/ncbi/public/sra/$p.sra ~/stijn2/1_SRA/location_SRA_files
  done < /home/stijn/Documents/3_lists/List_of_SRA_accession_numbers.txt

# Using wget
# use this from the location where you want your SRA files to be.

while read p;
  do
    echo $p
    wget ftp://ftp.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByRun/sra/${p:0:3}/${p:0:6}/${p:0:10}/${p}.sra
    #sleep $((5 + RANDOM % 15))s
  done < /home/stijn/Documents/3_lists/List_of_SRA_accession_numbers.txt

# Using esearch
# I didnt really use this one

while read p;
  do
    echo $p
    esearch -db sra -query $p | efetch --format runinfo | cut -d ',' -f 1 | grep SRR | xargs fastq-dump --split-files --bzip2
    #sleep $((5 + RANDOM % 15))s
  done < /home/stijn/Documents/3_lists/List_of_SRA_accession_numbers.txt
