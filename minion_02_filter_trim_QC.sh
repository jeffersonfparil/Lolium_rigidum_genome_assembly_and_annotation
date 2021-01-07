#!/bin/bash
#######################################################
### Filter and trimm-off adapters from MinION reads ###
#######################################################

### Input:
### (1) MinION reads in compressed fastq format (minion.fastq.gz)

### Outputs:
### (1) Filtered reads (minion-filtered.fastq.gz)
### (2) Filtered and trimmed reads (minion-filtered-trimmed.fastq.gz)

### Parameters:
DIR=/data/Lolium_rigidum_ASSEMBLY/assembly_annotation_pipeline_tests_20210104/FASTQ/MINION/

### Navigate to the working directory
cd $DIR

### Filter-out reads with average PHRED score < 10 with nanofilt (PHRED threshold based on a qualitative look at the FastQC output)
time \
     gunzip -c minion.fastq.gz | \
     NanoFilt -q 10 | \
     gzip > minion-filtered.fastq.gz

### Trim-off adapters with porechop
time \
porechop \
     --threads 32 \
     --input minion-filtered.fastq.gz \
     --output minion-filtered-trimmed.fastq.gz

### Count the total number of bases sequenced
zcat minion-filtered-trimmed.fastq.gz | paste - - - - | cut -f 2 | tr -d '\n' | wc -c > minion-filtered-trimmed.base.count
### As of 2021-01-07 we have ~5 billion bases (5,369,872,027) sequenced which means ~2.68X depth of coverage
