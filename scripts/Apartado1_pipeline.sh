#!/bin/bash

# Requiered software
# mamba install -c bioconda fastqc

#read -p "Enter the path where sequences are hosted: " ruta
#echo $ruta

# Quality of fastq files - FastQC software
# mkdir -p Apartado1/output/fastqc_results
# 
# for seqs in $(find Apartado1/input/ -name *.fastq)
# do
#     fastqc -o output/fastqc_results/ "$seqs"
# done

#FASTQ-SCREEN

# Trimmomatic pre-processing
# mamba install -c bioconda trimmomatic

# trimmomatic PE -phred33 -trimlog Apartado1/output/trimmomatic_results/trimlogFile.txt -summary Apartado1/output/trimmomatic_results/summary.txt \ 
#   Apartado1/input/SRR479052.chr21_1.fastq Apartado1/input/SRR479052.chr21_2.fastq \ 
#   Apartado1/output/trimmomatic_results/SRR479052.chr21_1_paired.fastq Apartado1/output/trimmomatic_results/SRR479052.chr21_1_unpaired.fastq \ 
#   Apartado1/output/trimmomatic_results/SRR479052.chr21_2_paired.fastq Apartado1/output/trimmomatic_results/SRR479052.chr21_2_unpaired.fastq \ 
#   ILLUMINACLIP:Apartado1/trimmomatic-0.39/adapters/TruSeq2-PE.fa:2:30:10 \
#   LEADING:20 TRAILING:20 SLIDINGWINDOW:4:20 MINLEN:50

# The LEADING and TRAILING options specify a minimum quality score of 20 for the first and last bases of the reads, respectively.
# The SLIDINGWINDOW option specifies a sliding window approach to quality trimming, where the average quality score of a window of 4 bases is calculated 
# and reads are trimmed from the 3' end if the average quality score falls below 20.
# Finally, the MINLEN option specifies the minimum length that a read must be after all trimming and filtering steps, which in this example is 50 bases. 
# You may need to adjust this value based on your specific needs.

# Indexing and Align wiht HISAT2

# mkdir -p Apartado1/output/hisat2/index \
#          Apartado1/output/hisat2/results \
#          Apartado1/output/hisat2/log

## Index with HISAT2

#hisat2-build --seed 123 -p 2 \
#             Apartado1/input/Homo_sapiens.GRCh38.dna.chromosome.21.fa \
#             Apartado1/output/hisat2/index/Homo_sapiens.GRCh38.dna.chromosome.21 > Apartado1/output/hisat2/log/hisat2_index.log


# Align HISAT2 (On going)
 for sid in $(cat Apartado1/output/hisat2/sample_id.txt); do \
    fw_path=$(find Apartado1/input -name "*$sid*_1*"); \
    rv_path=$(find Apartado1/input -name "*$sid*_2*"); \
    #echo "$fw_path and $rv_path"
    hisat2 --new-summary --summary-file Apartado1/output/hisat2/results/$sid.hisat2.summary \
            --rna-strandness RF --seed 123 --phred33 -p 2 -k 1 \
            -x Apartado1/output/hisat2/index/Homo_sapiens.GRCh38.dna.chromosome.21 \
            -1 $fw_path -2 $rv_path \
            -S Apartado1/output/hisat2/results/$sid.sam
    echo "\n"
 done

# STAR aligner
# Kallisto Aligner