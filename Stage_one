#!/bin/bash

#Count the number of sequences in DNA.fa 
grep -c "^>" DNA.fa

#download dataset
wget https://raw.githubusercontent.com/HackBio-Internship/wale-home-tasks/main/DNA.fa
cat DNA.fa

#Write a one-line command in Bash to get the total A, T, G & C counts for all the sequences in the file above
echo -n "TGGGTTGATTCCACACCCCCGCCCGGCACCCGCGTCCGCGCCGTGGCCATCTACAAGCAGTCACAGCACATGACGGAGGTTGTGAGGCGCTGCCCCCACCATGAGCGCTGCTCAGATAGCGAT" | sed 's/\(.\)/\1\n/g'| sort | uniq -c

#OR
grep -Eo 'A|T|G|C' DNA.fa | sort | uniq -c | awk '{print $2": "$1}'

#Set up a conda (anaconda, miniconda or mini forge) environment on your terminal.
wget https://repo.anaconda.com/miniconda/Miniconda3-py39_4.12.0-Linux-x86_64.sh
chmod +x Miniconda3-py39_4.12.0-Linux-x86_64.sh
./Miniconda3-py39_4.12.0-Linux-x86_64.sh
 source ~/.bash.rc

 #Downloads some (>2) sample datasets from here: https://github.com/josoga2/yt-dataset/tree/main/dataset
wget https://github.com/josoga2/yt-dataset/blob/main/dataset/raw_reads/ACBarrie_R1.fastq.gz?raw=true/ -O ACBarrie_R1.fastq.gz
wget https://github.com/josoga2/yt-dataset/blob/main/dataset/raw_reads/ACBarrie_R2.fastq.gz?raw=true/ -O ACBarrie_R2.fastq.gz
wget https://github.com/josoga2/yt-dataset/blob/main/dataset/raw_reads/Alsen_R1.fastq.gz?raw=true/ -O Alsen_R1.fasq.gz
wget https://github.com/josoga2/yt-dataset/blob/main/dataset/raw_reads/Alsen_R2.fastq.gz?raw=true/ -O Alsen_R2.fasq.gz
wget https://github.com/josoga2/yt-dataset/blob/main/dataset/raw_reads/Baxter_R1.fastq.gz?raw=true/ -O Baxter_R1.fastq.gz
wget https://github.com/josoga2/yt-dataset/blob/main/dataset/raw_reads/Baxter_R2.fastq.gz?raw=true/ -O Baxter_R2.fastq.gz
wget https://github.com/josoga2/yt-dataset/blob/main/dataset/raw_reads/Chara_R1.fastq.gz?raw=true/ -O Chara_R1.fastq.gz
wget https://github.com/josoga2/yt-dataset/blob/main/dataset/raw_reads/Chara_R2.fastq.gz?raw=true/ -O Chara_R2.fastq.gz
wget https://github.com/josoga2/yt-dataset/blob/main/dataset/raw_reads/Drysdale_R1.fastq.gz?raw=true/ -O Drysdale_R1.fastq.gz
wget https://github.com/josoga2/yt-dataset/blob/main/dataset/raw_reads/Drysdale_R2.fastq.gz?raw=true/ -O Drysdale_R2.fastq.gz

#Create a folder called output
mkdir output

#Implement the three software on the downloaded files (sample datasets) and send all outputs to the output folder.

#For fastp
mkdir raw_reads
mv ACBarrie_R1.fastq.gz ACBarrie_R2.fastq.gz Alsen_R1.fastq.gz Alsen_R2.fastq.gz Baxter_R1.fastq.gz Baxter_R2.fastq.gz Chara_R1.fastq.gz Chara_R2.fastq.gz Drysdale_R1.fastq.gz Drysdale_R2.fastq.gz /home/damadegbite/raw_reads
cd raw_reads
wget https://github.com/josoga2/yt-dataset/blob/main/dataset/raw_reads/trim.sh
sudo apt install fastp
bash trim.sh
nano trim.sh
mv qc_reads/ trimmed_reads
cd trimmed_reads/
mv ACBarrie_fastq.html Alsen_fastq.html Baxter_fastq.html Chara_fastq.html Drysdale_fastq.html /home/damadegbite/output

#For fastqc
sudo apt install fastqc
ls
fastqc Alsen_R1.fasq.gz -O output/
fastqc Alsen_R2.fasq.gz -O output/
fastqc Baxter_R1.fastq.gz -O output/
fastqc Baxter_R2.fastq.gz -O output/
fastqc Chara_R1.fastq.gz -O output/
fastqc Chara_R2.fastq.gz -O output/
fastqc Drysdale_R1.fastq.gz -O output/
fastqc Drysdale_R2.fastq.gz -O output/

#For multiqc
sudo apt install multiqc
multiqc output/
