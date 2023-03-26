#!/bin/bash

# mamba install -c bioconda picard  # Java 16o+ required - in conflict with fastqc (Java 11)
# mamba install -c bioconda htseq

mkdir -p Apartado1/output/picard/results \
        Apartado1/output/picard/log

for sid in $(cat Apartado1/output/hisat2/sample_id.txt); 
 do 
    mkdir -p Apartado1/output/picard/results/$sid
    #Picard
    picard MarkDuplicates \
        INPUT=Apartado1/output/hisat2/results/$sid/$sid.sorted.bam \
        OUTPUT=Apartado1/output/picard/results/$sid/$sid.sorted.marked.bam \
        METRICS_FILE=Apartado1/output/picard/log/$sid.metrics.txt \
        REMOVE_DUPLICATES=true \
        ASSUME_SORTED=true \
        VALIDATION_STRINGENCY=SILENT 
        2> Apartado1/output/picard/log/$sid.log

    #HTseq for Picard-samples
    htseq-count \
        --format=bam \
        --mode=intersection-nonempty \
        --minaqual=10 \
        --type=exon \
        --idattr=gene_id \
        --additional-attr=gene_name \
        Apartado1/output/picard/results/$sid/$sid.sorted.marked.bam \
        Apartado1/input/Homo_sapiens.GRCh38.109.chr21.gtf \
        > Apartado1/output/htseq/results/"$sid".marked.gene_counts.txt 2> Apartado1/output/htseq/log/htseq_marked_report.log

 done
