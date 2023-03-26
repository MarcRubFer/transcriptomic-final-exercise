#!/bin/bash


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
## The SLIDINGWINDOW option specifies a sliding window approach to quality trimming, where the average quality score of a window of 4 bases is calculated 
## and reads are trimmed from the 3' end if the average quality score falls below 20.
## Finally, the MINLEN option specifies the minimum length that a read must be after all trimming and filtering steps, which in this example is 50 bases. 
## You may need to adjust this value based on your specific needs.
#
## Alinamiento de archivos paired
#hisat2 -x genome_index -1 muestra1_1_paired.fastq -2 muestra_1_2_paired.fastq -S muestra1_paired.sam
#
## Conversión de SAM a BAM
#samtools view -bS muestra1_paired.sam > muestra1_paired.bam
#
## Ordenamiento y indexación del archivo BAM
#samtools sort muestra1_paired.bam -o muestra1_paired_sorted.bam
#samtools index muestra1_paired_sorted.bam
#
## Alinamiento de archivos unpaired
#hisat2 -x genome_index -U muestra1_1_unpaired.fastq -U muestra1_2_unpaired.fastq -S muestra1_unpaired.sam
#
## Conversión de SAM a BAM
#samtools view -bS muestra1_unpaired.sam > muestra1_unpaired.bam
#
## Ordenamiento y indexación del archivo BAM
#samtools sort muestra1_unpaired.bam -o muestra1_unpaired_sorted.bam
#samtools index muestra1_unpaired_sorted.bam
#
## Fusión de archivos BAM
#samtools merge muestra1.bam muestra1_paired_sorted.bam muestra1_unpaired_sorted.bam
#
## Ordenamiento y indexación del archivo BAM fusionado
#samtools sort muestra1.bam -o muestra1_sorted.bam
#samtools index muestra1_sorted.bam
#