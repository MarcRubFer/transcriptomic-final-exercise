#!/bin/bash

# Requiered software

# mamba install -c bioconda fastqc
# mamba install -c bioconda hisat2
# mamba install -c bioconda samtools
# mamba install -c bioconda htseq
# pip install multiqc
# -----
# New for this pipeline
# mamba install -c agbiome bbtools


mkdir -p Apartado1/output_BBDuk

# FastQC control quality

mkdir -p Apartado1/output_BBDuk/fastqc_results/non-trimmed

for seqs in $(find Apartado1/input/ -name "*.fastq")
do
    fastqc "$seqs" -o Apartado1/output_BBDuk/fastqc_results/non-trimmed/ 2>>Apartado1/output_BBDuk/fastqc_results/non-trimmed/log_non-trimmed.txt
done

# Trimming, Index, Align and Counts

mkdir -p Apartado1/output_BBDuk/hisat2/index \
         Apartado1/output_BBDuk/hisat2/results \
         Apartado1/output_BBDuk/hisat2/log \
         Apartado1/output_BBDuk/htseq/results \
         Apartado1/output_BBDuk/htseq/log \
         Apartado1/output_BBDuk/BBDuk_results \
         Apartado1/output_BBDuk/fastqc_results/trimmed

## Index with HISAT2

hisat2-build --seed 123 -p 2 \
             Apartado1/input/Homo_sapiens.GRCh38.dna.chromosome.21.fa \
             Apartado1/output_BBDuk/hisat2/index/Homo_sapiens.GRCh38.dna.chromosome.21 > Apartado1/output_BBDuk/hisat2/log/hisat2_index.log

## Make sample_id.txt
find Apartado1/input/ -name '*.fastq' | xargs basename -s .fastq | cut -d'.' -f1 | uniq > Apartado1/output_BBDuk/hisat2/sample_id.txt


for sid in $(cat Apartado1/output_BBDuk/hisat2/sample_id.txt); 
do
    # Stablish path for fw/rv sequences
    fw_path=$(find Apartado1/input -name "*$sid*_1*")
    rv_path=$(find Apartado1/input -name "*$sid*_2*")
    
    # BBDuk adapter trimming 
    bbduk.sh in1=$fw_path in2=$rv_path \
            out1=Apartado1/output_BBDuk/BBDuk_results/"$sid"_1.trimmed.fastq out2=Apartado1/output_BBDuk/BBDuk_results/"$sid"_2.trimmed.fastq \
            ref=Apartado1/input/BBDuk_adapters/adapters.fa \
            ktrim=r k=23 mink=11 hdist=1 tpe tbo \

    #FASTQC Post-trimmed
    for seqs in $(find Apartado1/output_BBDuk/BBDuk_results -name "*.fastq")
    do
        fastqc "$seqs" -o Apartado1/output_BBDuk/fastqc_results/trimmed/ 2>>Apartado1/output_BBDuk/fastqc_results/trimmed/log_trimmed.txt
    done
    
    # Stablish path for fw/rv sequences
    fw_path=$(find Apartado1/output_BBDuk/BBDuk_results -name "*$sid*_1*")
    rv_path=$(find Apartado1/output_BBDuk/BBDuk_results -name "*$sid*_2*")
    
    # Create directory for each sample and align
    mkdir -p Apartado1/output_BBDuk/hisat2/results/$sid
    
    # Align paired sequences
    hisat2 --new-summary --summary-file Apartado1/output_BBDuk/hisat2/results/$sid/$sid.hisat2.summary \
            --seed 123 --phred33 -p 2 -k 1 \
            -x Apartado1/output_BBDuk/hisat2/index/Homo_sapiens.GRCh38.dna.chromosome.21 \
            -1 $fw_path -2 $rv_path \
            -S Apartado1/output_BBDuk/hisat2/results/$sid/$sid.sam
    
    samtools view -bS Apartado1/output_BBDuk/hisat2/results/$sid/$sid.sam > Apartado1/output_BBDuk/hisat2/results/$sid/$sid.bam
    samtools sort Apartado1/output_BBDuk/hisat2/results/$sid/$sid.bam -o Apartado1/output_BBDuk/hisat2/results/$sid/$sid.sorted.bam
    samtools index Apartado1/output_BBDuk/hisat2/results/$sid/$sid.sorted.bam


    
    #SAMTOOLS alignment statistics
    
    samtools flagstat Apartado1/output_BBDuk/hisat2/results/$sid/$sid.sorted.bam >> Apartado1/output_BBDuk/hisat2/results/$sid/$sid.samtools.summary
    
    #HTseq count
    htseq-count \
        --format=bam \
        --mode=intersection-nonempty \
        --minaqual=10 \
        --type=exon \
        --idattr=gene_id \
        --additional-attr=gene_name \
        Apartado1/output_BBDuk/hisat2/results/$sid/$sid.sorted.bam \
        Apartado1/input/Homo_sapiens.GRCh38.109.chr21.gtf \
        > Apartado1/output_BBDuk/htseq/results/$sid.gene_counts.txt 2> Apartado1/output_BBDuk/htseq/log/htseq_report.log
        
 done

 # MultiQC report

 mkdir -p Apartado1/output_BBDuk/multiqc_report

 multiqc -f -o Apartado1/output_BBDuk/multiqc_report/ Apartado1/output_BBDuk
