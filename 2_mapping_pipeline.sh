#!/bin/bash

# This script contains the pipeline for unpacking the SRA files till the mapping
# to the SILVA and HLA allele databases
# Overview:
# Unpacking SRA (fastq-dump) > Quality Control (fastqc) > Trimming of reads (AdapterRemoval)
# > Quality Control (fastqc) > mapping of reads to reference database (BWA MEM)

# This script is build to run on one of the TBB mutants and for my personal directories
# Quality Control outputs are combined in text files to have an overview.
# and all the directories are created for the different files


# Moving the SRA files to the mutant
mkdir /linuxhome/tmp/stijn/4_SRA/
while read p;
  do
  echo $p
  cp -v /home/stijn/stijn2/1_SRA/1000-3400SRA/$p /linuxhome/tmp/stijn/4_SRA/
done < ~/Documents/3_lists/NA_download_2348/SRA_part2

echo "Mapping SSU Silva & Human Ref. genomes pipeline"

#creating the directories
# PE is Paired End, and SE is Single End
mkdir /linuxhome/tmp/stijn/SRA_fastq/
mkdir /linuxhome/tmp/stijn/SRA_fastq/PE/
mkdir /linuxhome/tmp/stijn/SRA_fastq/SE/
mkdir /linuxhome/tmp/stijn/SRA_fastQC/
mkdir /linuxhome/tmp/stijn/SRA_fastQC/untrimmed/
mkdir /linuxhome/tmp/stijn/SRA_fastQC/QC_reports/
mkdir /linuxhome/tmp/stijn/listsSRA/

# make headers for Quality Control files
echo 'f_name' > /linuxhome/tmp/stijn/SRA_fastQC/QC_reports/Sample_name_before_trimming.txt
echo 'amount_reads_BT' > /linuxhome/tmp/stijn/SRA_fastQC/QC_reports/number_of_reads_before_trimming.txt
echo 'poor_Q_BT' > /linuxhome/tmp/stijn/SRA_fastQC/QC_reports/Sequences_flagged_as_poor_quality_before_trimming.txt
echo 'GC_content_BT' > /linuxhome/tmp/stijn/SRA_fastQC/QC_reports/GC_content_before_trimming.txt
echo 'seq_lenght_BT' > /linuxhome/tmp/stijn/SRA_fastQC/QC_reports/sequence_length_before_trimming.txt


# Unpack SRA files into fastq files and do quality control on the fastq files
# settings fastq-dump from https://edwards.sdsu.edu/research/fastq-dump/
cd /linuxhome/tmp/stijn/4_SRA/
for fname in *.sra
	do
		base=${fname%.sra}
		echo "start fastq-dump for ${base}"
		fastq-dump.2.9.2 --gzip --skip-technical -F --readids --read-filter pass --dumpbase --split-3 --clip \
		 /linuxhome/tmp/stijn/4_SRA/${base}.sra --outdir /linuxhome/tmp/stijn/SRA_fastq/ #get fastq from sra file
		if [ -f /linuxhome/tmp/stijn/SRA_fastq/${base}_pass_1.fastq.gz ]; \
		then
			mv /linuxhome/tmp/stijn/SRA_fastq/${base}_pass_1.fastq.gz /linuxhome/tmp/stijn/SRA_fastq/PE/${base}_1.fastq.gz
			mv /linuxhome/tmp/stijn/SRA_fastq/${base}_pass_2.fastq.gz /linuxhome/tmp/stijn/SRA_fastq/PE/${base}_2.fastq.gz
			rm -rf /linuxhome/tmp/stijn/SRA_fastq/SE/${base}_pass.fastq.gz
			echo "${base} is a PE!"
			fastqc -t 20 -f fastq --extract -o /linuxhome/tmp/stijn/SRA_fastQC/untrimmed/ \
			/linuxhome/tmp/stijn/SRA_fastq/PE/${base}_1.fastq.gz
	    		fastqc -t 20 -f fastq --extract -o /linuxhome/tmp/stijn/SRA_fastQC/untrimmed/ \
			/linuxhome/tmp/stijn/SRA_fastq/PE/${base}_2.fastq.gz
		# extract information from fastqc files
		#file name
		    	sed -n '4p' /linuxhome/tmp/stijn/SRA_fastQC/untrimmed/${base}_1_fastqc/fastqc_data.txt \
			| awk '{print$2}'|tee -a /linuxhome/tmp/stijn/SRA_fastQC/QC_reports/Sample_name_before_trimming.txt
		    	sed -n '4p' /linuxhome/tmp/stijn/SRA_fastQC/untrimmed/${base}_2_fastqc/fastqc_data.txt \
			| awk '{print$2}'|tee -a /linuxhome/tmp/stijn/SRA_fastQC/QC_reports/Sample_name_before_trimming.txt
		#Number of reads before trimming
		    	sed -n '7p' /linuxhome/tmp/stijn/SRA_fastQC/untrimmed/${base}_1_fastqc/fastqc_data.txt \
			| awk '{print$3}' |tee -a /linuxhome/tmp/stijn/SRA_fastQC/QC_reports/number_of_reads_before_trimming.txt
	    		sed -n '7p' /linuxhome/tmp/stijn/SRA_fastQC/untrimmed/${base}_2_fastqc/fastqc_data.txt \
			| awk '{print$3}' |tee -a /linuxhome/tmp/stijn/SRA_fastQC/QC_reports/number_of_reads_before_trimming.txt
		#Amount of sequences flagged as poor quality
	    		sed -n '8p' /linuxhome/tmp/stijn/SRA_fastQC/untrimmed/${base}_1_fastqc/fastqc_data.txt \
			| awk '{print$6}' |tee -a /linuxhome/tmp/stijn/SRA_fastQC/QC_reports/Sequences_flagged_as_poor_quality_before_trimming.txt
	    		sed -n '8p' /linuxhome/tmp/stijn/SRA_fastQC/untrimmed/${base}_2_fastqc/fastqc_data.txt \
			| awk '{print$6}' |tee -a /linuxhome/tmp/stijn/SRA_fastQC/QC_reports/Sequences_flagged_as_poor_quality_before_trimming.txt
		#QC content before trimming
		    	sed -n '10p' /linuxhome/tmp/stijn/SRA_fastQC/untrimmed/${base}_1_fastqc/fastqc_data.txt \
			| awk '{print$2}'|tee -a /linuxhome/tmp/stijn/SRA_fastQC/QC_reports/GC_content_before_trimming.txt
		    	sed -n '10p' /linuxhome/tmp/stijn/SRA_fastQC/untrimmed/${base}_2_fastqc/fastqc_data.txt \
			| awk '{print$2}'|tee -a /linuxhome/tmp/stijn/SRA_fastQC/QC_reports/GC_content_before_trimming.txt
		#Sequence length before trimming
		    	sed -n '9p' /linuxhome/tmp/stijn/SRA_fastQC/untrimmed/${base}_1_fastqc/fastqc_data.txt \
			| awk '{print$3}' |tee -a /linuxhome/tmp/stijn/SRA_fastQC/QC_reports/sequence_length_before_trimming.txt
	    		sed -n '9p' /linuxhome/tmp/stijn/SRA_fastQC/untrimmed/${base}_2_fastqc/fastqc_data.txt \
			| awk '{print$3}' |tee -a /linuxhome/tmp/stijn/SRA_fastQC/QC_reports/sequence_length_before_trimming.txt
		#remove not zipped fastqc report
		    	rm -rf /linuxhome/tmp/stijn/SRA_fastQC/untrimmed/${base}_1_fastqc/
	    		rm -rf /linuxhome/tmp/stijn/SRA_fastQC/untrimmed/${base}_2_fastqc/
		else
			mv /linuxhome/tmp/stijn/SRA_fastq/${base}* /linuxhome/tmp/stijn/SRA_fastq/SE/${base}.fastq.gz
			echo "${base} is a SE!"
			fastqc -t 24 -f fastq --extract -o /linuxhome/tmp/stijn/SRA_fastQC/untrimmed/ \
			/linuxhome/tmp/stijn/SRA_fastq/SE/${base}.fastq.gz

		#file name
		    	sed -n '4p' /linuxhome/tmp/stijn/SRA_fastQC/untrimmed/${base}_fastqc/fastqc_data.txt \
			| awk '{print$2}'|tee -a /linuxhome/tmp/stijn/SRA_fastQC/QC_reports/Sample_name_before_trimming.txt
		#Number of reads before trimming
		    	sed -n '7p' /linuxhome/tmp/stijn/SRA_fastQC/untrimmed/${base}_fastqc/fastqc_data.txt \
			| awk '{print$3}' |tee -a /linuxhome/tmp/stijn/SRA_fastQC/QC_reports/number_of_reads_before_trimming.txt
		#Amount of sequences flagged as poor quality
		    	sed -n '8p' /linuxhome/tmp/stijn/SRA_fastQC/untrimmed/${base}_fastqc/fastqc_data.txt \
			| awk '{print$6}' |tee -a /linuxhome/tmp/stijn/SRA_fastQC/QC_reports/Sequences_flagged_as_poor_quality_before_trimming.txt
		#QC content before trimming
		    	sed -n '10p' /linuxhome/tmp/stijn/SRA_fastQC/untrimmed/${base}_fastqc/fastqc_data.txt \
			| awk '{print$2}'|tee -a /linuxhome/tmp/stijn/SRA_fastQC/QC_reports/GC_content_before_trimming.txt
		#Sequence length before trimming
		    	sed -n '9p' /linuxhome/tmp/stijn/SRA_fastQC/untrimmed/${base}_fastqc/fastqc_data.txt \
			| awk '{print$3}' |tee -a /linuxhome/tmp/stijn/SRA_fastQC/QC_reports/sequence_length_before_trimming.txt
		#remove not zipped fastqc report
		    	rm -rf /linuxhome/tmp/stijn/SRA_fastQC/untrimmed/${base}_fastqc/
		fi
		rm -rf /linuxhome/tmp/stijn/4_SRA/${base}.sra
	done

#make directories for trimmed fastq files
mkdir /linuxhome/tmp/stijn/SRA_trimmed/
mkdir /linuxhome/tmp/stijn/SRA_trimmed/PE/
mkdir /linuxhome/tmp/stijn/SRA_trimmed/SE/
mkdir /linuxhome/tmp/stijn/SRA_fastQC/trimmed/
#make headers for QC_data after trimming
echo 'f_name' > /linuxhome/tmp/stijn/SRA_fastQC/QC_reports/Sample_name_after_trimming.txt
echo 'amount_reads_AT' > /linuxhome/tmp/stijn/SRA_fastQC/QC_reports/number_of_reads_after_trimming.txt
echo 'poor_Q_AT' > /linuxhome/tmp/stijn/SRA_fastQC/QC_reports/Sequences_flagged_as_poor_quality_after_trimming.txt
echo 'GC_content_AT' > /linuxhome/tmp/stijn/SRA_fastQC/QC_reports/GC_content_after_trimming.txt
echo 'seq_lenght_AT' > /linuxhome/tmp/stijn/SRA_fastQC/QC_reports/sequence_length_after_trimming.txt

# trimming of fastq PE files using AdapterRemoval and afterwards again do quality control
cd /linuxhome/tmp/stijn/SRA_fastq/PE/
for fname in *_1.fastq.gz
	do
    		base=${fname%_1*}
		echo "start adapterremoval PE ${base}"
		AdapterRemoval --threads 24 --file1 ${base}_1.fastq.gz --file2 ${base}_2.fastq.gz --basename ${base} --trimns \
		--trimqualities --minquality 25 --gzip       #standard original settings: --minquality 25
    		mv ${base}.*.truncated.gz /linuxhome/tmp/stijn/SRA_trimmed/PE/
		fastqc -t 20 -f fastq --extract -o /linuxhome/tmp/stijn/SRA_fastQC/trimmed/ \
		/linuxhome/tmp/stijn/SRA_trimmed/PE/${base}.pair1.truncated.gz
		fastqc -t 20 -f fastq --extract -o /linuxhome/tmp/stijn/SRA_fastQC/trimmed/ \
		/linuxhome/tmp/stijn/SRA_trimmed/PE/${base}.pair2.truncated.gz
	# extract information from QC files.
	#file name
		sed -n '4p' /linuxhome/tmp/stijn/SRA_fastQC/trimmed/${base}.pair1.truncated_fastqc/fastqc_data.txt \
		| awk '{print$2}'|tee -a /linuxhome/tmp/stijn/SRA_fastQC/QC_reports/Sample_name_after_trimming.txt
		sed -n '4p' /linuxhome/tmp/stijn/SRA_fastQC/trimmed/${base}.pair2.truncated_fastqc/fastqc_data.txt \
		| awk '{print$2}'|tee -a /linuxhome/tmp/stijn/SRA_fastQC/QC_reports/Sample_name_after_trimming.txt
	#Number of reads before trimming
		sed -n '7p' /linuxhome/tmp/stijn/SRA_fastQC/trimmed/${base}.pair1.truncated_fastqc/fastqc_data.txt \
		| awk '{print$3}' |tee -a /linuxhome/tmp/stijn/SRA_fastQC/QC_reports/number_of_reads_after_trimming.txt
		sed -n '7p' /linuxhome/tmp/stijn/SRA_fastQC/trimmed/${base}.pair2.truncated_fastqc/fastqc_data.txt \
		| awk '{print$3}' |tee -a /linuxhome/tmp/stijn/SRA_fastQC/QC_reports/number_of_reads_after_trimming.txt
	#Amount of sequences flagged as poor quality
		sed -n '8p' /linuxhome/tmp/stijn/SRA_fastQC/trimmed/${base}.pair1.truncated_fastqc/fastqc_data.txt \
		| awk '{print$6}' |tee -a /linuxhome/tmp/stijn/SRA_fastQC/QC_reports/Sequences_flagged_as_poor_quality_after_trimming.txt
		sed -n '8p' /linuxhome/tmp/stijn/SRA_fastQC/trimmed/${base}.pair2.truncated_fastqc/fastqc_data.txt \
		| awk '{print$6}' |tee -a /linuxhome/tmp/stijn/SRA_fastQC/QC_reports/Sequences_flagged_as_poor_quality_after_trimming.txt
	#QC content before trimming
		sed -n '10p' /linuxhome/tmp/stijn/SRA_fastQC/trimmed/${base}.pair1.truncated_fastqc/fastqc_data.txt \
		| awk '{print$2}'|tee -a /linuxhome/tmp/stijn/SRA_fastQC/QC_reports/GC_content_after_trimming.txt
		sed -n '10p' /linuxhome/tmp/stijn/SRA_fastQC/trimmed/${base}.pair2.truncated_fastqc/fastqc_data.txt \
		| awk '{print$2}'|tee -a /linuxhome/tmp/stijn/SRA_fastQC/QC_reports/GC_content_after_trimming.txt
	#Sequence length before trimming
		sed -n '9p' /linuxhome/tmp/stijn/SRA_fastQC/trimmed/${base}.pair1.truncated_fastqc/fastqc_data.txt \
		| awk '{print$3}' |tee -a /linuxhome/tmp/stijn/SRA_fastQC/QC_reports/sequence_length_after_trimming.txt
		sed -n '9p' /linuxhome/tmp/stijn/SRA_fastQC/trimmed/${base}.pair2.truncated_fastqc/fastqc_data.txt \
		| awk '{print$3}' |tee -a /linuxhome/tmp/stijn/SRA_fastQC/QC_reports/sequence_length_after_trimming.txt
	#remove not zipped fastqc report
		rm -rf /linuxhome/tmp/stijn/SRA_fastQC/trimmed/${base}_1_fastqc/
		rm -rf /linuxhome/tmp/stijn/SRA_fastQC/trimmed/${base}_2_fastqc/
	done
#for SE
cd /linuxhome/tmp/stijn/SRA_fastq/SE/
for fname in *.fastq.gz
	do
    		base=${fname%.fastq*}
		echo "start adapterremoval SE ${base}"
    		AdapterRemoval --threads 24 --file1 ${base}.fastq.gz  --basename ${base} --trimns --trimqualities --minquality 25 --gzip
    		mv ${base}.truncated.gz /linuxhome/tmp/stijn/SRA_trimmed/SE/
		mv ${base}.discarded.gz /linuxhome/tmp/stijn/SRA_trimmed/SE/
		mv ${base}.settings /linuxhome/tmp/stijn/SRA_trimmed/SE/
	#run fastqc on trimmed SE reads
		fastqc -t 20 -f fastq --extract -o /linuxhome/tmp/stijn/SRA_fastQC/trimmed/ /linuxhome/tmp/stijn/SRA_trimmed/SE/${base}.truncated.gz
	#file name
		sed -n '4p' /linuxhome/tmp/stijn/SRA_fastQC/trimmed/${base}.truncated_fastqc/fastqc_data.txt \
		| awk '{print$2}'|tee -a /linuxhome/tmp/stijn/SRA_fastQC/QC_reports/Sample_name_after_trimming.txt
	#Number of reads before trimming
		sed -n '7p' /linuxhome/tmp/stijn/SRA_fastQC/trimmed/${base}.truncated_fastqc/fastqc_data.txt \
		| awk '{print$3}' |tee -a /linuxhome/tmp/stijn/SRA_fastQC/QC_reports/number_of_reads_after_trimming.txt
	#Amount of sequences flagged as poor quality
		sed -n '8p' /linuxhome/tmp/stijn/SRA_fastQC/trimmed/${base}.truncated_fastqc/fastqc_data.txt \
		| awk '{print$6}' |tee -a /linuxhome/tmp/stijn/SRA_fastQC/QC_reports/Sequences_flagged_as_poor_quality_after_trimming.txt
	#QC content before trimming
		sed -n '10p' /linuxhome/tmp/stijn/SRA_fastQC/trimmed/${base}.truncated_fastqc/fastqc_data.txt \
		| awk '{print$2}'|tee -a /linuxhome/tmp/stijn/SRA_fastQC/QC_reports/GC_content_after_trimming.txt
	#Sequence length before trimming
		sed -n '9p' /linuxhome/tmp/stijn/SRA_fastQC/trimmed/${base}.truncated_fastqc/fastqc_data.txt \
		| awk '{print$3}' |tee -a /linuxhome/tmp/stijn/SRA_fastQC/QC_reports/sequence_length_after_trimming.txt
	#remove not zipped fastqc report
		rm -rf /linuxhome/tmp/stijn/SRA_fastQC/trimmed/${base}_fastqc/
	done
#adding QC lists

a="$(wc -l /linuxhome/tmp/stijn/SRA_fastQC/QC_reports/Sample_name_before_trimming.txt | sed 's/ .*//')"
b="$(wc -l /linuxhome/tmp/stijn/SRA_fastQC/QC_reports/number_of_reads_before_trimming.txt | sed 's/ .*//')"
c="$(wc -l /linuxhome/tmp/stijn/SRA_fastQC/QC_reports/Sequences_flagged_as_poor_quality_before_trimming.txt | sed 's/ .*//')"
d="$(wc -l /linuxhome/tmp/stijn/SRA_fastQC/QC_reports/GC_content_before_trimming.txt | sed 's/ .*//')"
e="$(wc -l /linuxhome/tmp/stijn/SRA_fastQC/QC_reports/sequence_length_before_trimming.txt | sed 's/ .*//')"
echo "${a}"
echo "${b}"
echo "${c}"
echo "${d}"
echo "${e}"

# check if text files are of equal length otherwise something went wrong.
# if length is the same combine text files to one overview.
if test "${a}" -eq "${b}" && test "${a}" -eq "${c}" && test "${a}" -eq "${d}" && test "${a}" -eq "${e}";
	then
		echo 'QC lists have the same length'
		paste /linuxhome/tmp/stijn/SRA_fastQC/QC_reports/Sample_name_before_trimming.txt \
		/linuxhome/tmp/stijn/SRA_fastQC/QC_reports/number_of_reads_before_trimming.txt \
		/linuxhome/tmp/stijn/SRA_fastQC/QC_reports/Sequences_flagged_as_poor_quality_before_trimming.txt \
		/linuxhome/tmp/stijn/SRA_fastQC/QC_reports/GC_content_before_trimming.txt \
		/linuxhome/tmp/stijn/SRA_fastQC/QC_reports/sequence_length_before_trimming.txt \
		> /linuxhome/tmp/stijn/SRA_fastQC/QC_reports/QC_dataSRA_before_trimming.txt
		sort /linuxhome/tmp/stijn/SRA_fastQC/QC_reports/QC_dataSRA_before_trimming.txt \
		> /linuxhome/tmp/stijn/SRA_fastQC/QC_reports/QC_dataSRA_before_trimming.sorted.txt
	else
		echo 'QC LISTS (before) NOT THE SAME LENGTH WARNING!!!!'
fi
f="$(wc -l /linuxhome/tmp/stijn/SRA_fastQC/QC_reports/Sample_name_after_trimming.txt | sed 's/ .*//')"
g="$(wc -l /linuxhome/tmp/stijn/SRA_fastQC/QC_reports/number_of_reads_after_trimming.txt | sed 's/ .*//')"
h="$(wc -l /linuxhome/tmp/stijn/SRA_fastQC/QC_reports/Sequences_flagged_as_poor_quality_after_trimming.txt | sed 's/ .*//')"
i="$(wc -l /linuxhome/tmp/stijn/SRA_fastQC/QC_reports/GC_content_after_trimming.txt | sed 's/ .*//')"
j="$(wc -l /linuxhome/tmp/stijn/SRA_fastQC/QC_reports/sequence_length_after_trimming.txt | sed 's/ .*//')"
echo "${f}"
echo "${g}"
echo "${h}"
echo "${i}"
echo "${j}"

if test "${f}" -eq "${g}" && test "${f}" -eq "${h}" && test "${f}" -eq "${i}" && test "${f}" -eq "${j}";
	then
		echo 'QC lists have the same length'
		paste /linuxhome/tmp/stijn/SRA_fastQC/QC_reports/Sample_name_after_trimming.txt \
		/linuxhome/tmp/stijn/SRA_fastQC/QC_reports/number_of_reads_after_trimming.txt \
		/linuxhome/tmp/stijn/SRA_fastQC/QC_reports/Sequences_flagged_as_poor_quality_after_trimming.txt \
		/linuxhome/tmp/stijn/SRA_fastQC/QC_reports/GC_content_after_trimming.txt \
		/linuxhome/tmp/stijn/SRA_fastQC/QC_reports/sequence_length_after_trimming.txt \
		> /linuxhome/tmp/stijn/SRA_fastQC/QC_reports/QC_dataSRA_after_trimming.txt
		sort /linuxhome/tmp/stijn/SRA_fastQC/QC_reports/QC_dataSRA_after_trimming.txt \
		> /linuxhome/tmp/stijn/SRA_fastQC/QC_reports/QC_dataSRA_after_trimming.sorted.txt
	else
		echo 'QC LISTS (after) NOT THE SAME LENGTH WARNING!!!!'
fi

k="$(wc -l /linuxhome/tmp/stijn/SRA_fastQC/QC_reports/QC_dataSRA_before_trimming.sorted.txt | sed 's/ .*//')"
l="$(wc -l /linuxhome/tmp/stijn/SRA_fastQC/QC_reports/QC_dataSRA_after_trimming.sorted.txt | sed 's/ .*//')"
echo "${k}"
echo "${l}"
if test "${k}" -eq "${l}";
	then
		echo 'QC lists sorted.before & sorted.after have the same length'
		paste /linuxhome/tmp/stijn/SRA_fastQC/QC_reports/QC_dataSRA_before_trimming.sorted.txt \
		/linuxhome/tmp/stijn/SRA_fastQC/QC_reports/QC_dataSRA_after_trimming.sorted.txt \
		> /linuxhome/tmp/stijn/listsSRA/QC_dataSRA_All.txt
	else
		echo 'QC LISTS sorted.before & sorted.after NOT THE SAME LENGTH WARNING!!!!'
fi

#removing collapsed lists
rm -rf /linuxhome/tmp/stijn/SRA_fastQC/QC_reports/Sample_name_before_trimming.txt
rm -rf /linuxhome/tmp/stijn/SRA_fastQC/QC_reports/Sample_name_after_trimming.txt
rm -rf /linuxhome/tmp/stijn/SRA_fastQC/QC_reports/number_of_reads_before_trimming.txt
rm -rf /linuxhome/tmp/stijn/SRA_fastQC/QC_reports/number_of_reads_after_trimming.txt
rm -rf /linuxhome/tmp/stijn/SRA_fastQC/QC_reports/Sequences_flagged_as_poor_quality_before_trimming.txt
rm -rf /linuxhome/tmp/stijn/SRA_fastQC/QC_reports/Sequences_flagged_as_poor_quality_after_trimming.txt
rm -rf /linuxhome/tmp/stijn/SRA_fastQC/QC_reports/GC_content_before_trimming.txt
rm -rf /linuxhome/tmp/stijn/SRA_fastQC/QC_reports/GC_content_after_trimming.txt
rm -rf /linuxhome/tmp/stijn/SRA_fastQC/QC_reports/sequence_length_before_trimming.txt
rm -rf /linuxhome/tmp/stijn/SRA_fastQC/QC_reports/sequence_length_after_trimming.txt

# Make multiQC, I only used this step in the beginning of the project, so now
# the tool needs updating. This multiQC step is not essential.
#cd /linuxhome/tmp/stijn/SRA_fastQC/untrimmed/
#multiqc /linuxhome/tmp/stijn/SRA_fastQC/untrimmed/.
#mv /linuxhome/tmp/stijn/SRA_fastQC/untrimmed/multiqc_report.html /linuxhome/tmp/stijn/listsSRA/multiqc_report_untrimmedSRA.html

#cd /linuxhome/tmp/stijn/SRA_fastQC/trimmed/
#multiqc /linuxhome/tmp/stijn/SRA_fastQC/trimmed/.
#mv /linuxhome/tmp/stijn/SRA_fastQC/trimmed/multiqc_report.html /linuxhome/tmp/stijn/listsSRA/multiqc_report_trimmedSRA.html

#make directories for mapping
mkdir /linuxhome/tmp/stijn/bam_SRA/

#move Reference genomes
cp -rf ~/stijn2/8_Ref_Genomes/SILVA_HLA/ /linuxhome/tmp/stijn/  # skip if already moved

mkdir /linuxhome/tmp/stijn/unsortedbam/
#Mapping PE reads to reference genomes.
# before you can map, the reference genome needs to be indexed this needs to be done only once.
cd /linuxhome/tmp/stijn/SRA_trimmed/PE/
for fname in *.pair1.truncated.gz
	do
		base=${fname%.pair1.truncated*}
		echo "start mapping PE ${base}"
		cd /linuxhome/tmp/stijn/bam_SRA/
		bwa mem -t 24 /linuxhome/tmp/stijn/SILVA_HLA/SILVA_HLA.fasta /linuxhome/tmp/stijn/SRA_trimmed/PE/${base}.pair1.truncated.gz \
		/linuxhome/tmp/stijn/SRA_trimmed/PE/${base}.pair2.truncated.gz | samtools view -b -S - > /linuxhome/tmp/stijn/unsortedbam/${base}.bam
		samtools sort /linuxhome/tmp/stijn/unsortedbam/${base}.bam /linuxhome/tmp/stijn/bam_SRA/${base}.srt && samtools index /linuxhome/tmp/stijn/bam_SRA/${base}.srt.bam
		rm -rf /linuxhome/tmp/stijn/unsortedbam/${base}.bam
	done

#Mapping SE reads
cd /linuxhome/tmp/stijn/SRA_trimmed/SE/
for fname in *.truncated.gz
	do
		base=${fname%.truncated*}
		echo "start mapping SE ${base}"
		cd /linuxhome/tmp/stijn/bam_SRA/
		bwa mem -t 24 /linuxhome/tmp/stijn/SILVA_HLA/SILVA_HLA.fasta /linuxhome/tmp/stijn/SRA_trimmed/SE/${base}.truncated.gz \
		| samtools view -b -S - > /linuxhome/tmp/stijn/unsortedbam/${base}.bam
		samtools sort /linuxhome/tmp/stijn/unsortedbam/${base}.bam /linuxhome/tmp/stijn/bam_SRA/${base}.srt && samtools index /linuxhome/tmp/stijn/bam_SRA/${base}.srt.bam
		rm -rf /linuxhome/tmp/stijn/unsortedbam/${base}.bam
	done

#deleting .fastq files
rm -rf /linuxhome/tmp/stijn/SRA_fastq/
# only remove these last two if you are 100% sure no mistakes were made.
rm -rf /linuxhome/tmp/stijn/SRA_trimmed/
#deleting fastQC files
rm -rf /linuxhome/tmp/stijn/SRA_fastQC
