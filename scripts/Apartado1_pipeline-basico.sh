#!/bin/bash

# Requiered software
# mamba install -c bioconda fastqc
# mamba install -c bioconda hisat2
# mamba install -c bioconda samtools
# mamba install -c bioconda htseq

# FastQC software - Quality of fastq files 

mkdir -p Apartado1/output/fastqc_results

for seqs in $(find Apartado1/input/ -name "*.fastq")
do
    fastqc "$seqs" -o Apartado1/output/fastqc_results/ 2>>Apartado1/output/fastqc_results/log.txt
done

# Indexing and Align wiht HISAT2

mkdir -p Apartado1/output/hisat2/index \
         Apartado1/output/hisat2/results \
         Apartado1/output/hisat2/log \
         Apartado1/output/htseq/results \
         Apartado1/output/htseq/log \
         

## Index with HISAT2

hisat2-build --seed 123 -p 2 \
             Apartado1/input/Homo_sapiens.GRCh38.dna.chromosome.21.fa \
             Apartado1/output/hisat2/index/Homo_sapiens.GRCh38.dna.chromosome.21 > Apartado1/output/hisat2/log/hisat2_index.log

#Creación de archivo sample_id.txt
find Apartado1/input/ -name '*.fastq' | xargs basename -s .fastq | cut -d'.' -f1 | uniq > Apartado1/output/hisat2/sample_id.txt

# Align HISAT2 (On going)
# rna-strandness default(unstranded) - desconozco la orientación del experimento
 for sid in $(cat Apartado1/output/hisat2/sample_id.txt); 
 do 
        fw_path=$(find Apartado1/input -name "*$sid*_1*")
        rv_path=$(find Apartado1/input -name "*$sid*_2*")
        mkdir -p Apartado1/output/hisat2/results/$sid
        hisat2 --new-summary --summary-file Apartado1/output/hisat2/results/$sid/$sid.hisat2.summary \
                --seed 123 --phred33 -p 2 -k 1 \
                -x Apartado1/output/hisat2/index/Homo_sapiens.GRCh38.dna.chromosome.21 \
                -1 $fw_path -2 $rv_path \
                -S Apartado1/output/hisat2/results/$sid/$sid.sam

        echo -e "\nSamtools processing: SAM -> BAM / Sorting / Indexing"
        #1. conversion .sam to .bam
        samtools view -bS Apartado1/output/hisat2/results/$sid/$sid.sam > Apartado1/output/hisat2/results/$sid/$sid.bam
        #2. sorting
        samtools sort Apartado1/output/hisat2/results/$sid/$sid.bam -o Apartado1/output/hisat2/results/$sid/$sid.sorted.bam
        #3. indexing
        samtools index Apartado1/output/hisat2/results/$sid/$sid.sorted.bam

        echo -e "\n"

        #HTseq count
        htseq-count \
            --format=bam \
            --mode=intersection-nonempty \
            --minaqual=10 \
            --type=exon \
            --idattr=gene_id \
            --additional-attr=gene_name \
            Apartado1/output/hisat2/results/$sid/$sid.sorted.bam \
            Apartado1/input/Homo_sapiens.GRCh38.109.chr21.gtf \
            > Apartado1/output/htseq/results/"$sid".gene_counts.txt 2> Apartado1/output/htseq/log/htseq_report.log
        
 done
