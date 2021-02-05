#!/bin/bash
###########################################################
### Hybrid assembly using Illumina and MinION via HASLR ###
###########################################################

### Inputs:
### (1) Error-corrected Illumina reads in compressed fastq format (BayesHammer_output/*cor.fastq.gz)
### (1) MinION reads in compressed fastq format (BayesHammer_output/*cor.fastq.gz)

### Outputs:
### (1) 
### (2) 
### (3) 
### (4) 

### Parameters:
INPUT_DIR_ILLUMINA=/data/Lolium_rigidum_ASSEMBLY/assembly_annotation_pipeline_tests_20210104/FASTQ/ILLUMINA/BayesHammer_output/
INPUT_DIR_MINION=/data/Lolium_rigidum_ASSEMBLY/assembly_annotation_pipeline_tests_20210104/FASTQ/MINION/
HASLR=/data/Lolium_rigidum_ASSEMBLY/assembly_annotation_pipeline_tests_20210104/haslr/bin/haslr.py
OUTPUT_DIR=/data/Lolium_rigidum_ASSEMBLY/assembly_annotation_pipeline_tests_20210104/ASSEMBLY/

### Intialise and spades output folder
cd ${OUTPUT_DIR}
mkdir Lori_hh/

### Assemble
time \
${HASLR} \
--threads 12 \
--out ${OUTPUT_DIR}/Lori_hh/ \
--genome 2g \
--long ${INPUT_DIR_MINION}/minion-filtered-trimmed.fastq.gz \
--type nanopore \
--short ${INPUT_DIR_ILLUMINA%BayesHammer_output*}/LOL-WGS-1.0_combined_R1.fastq.00.0_0.cor.fastq.gz ${INPUT_DIR_ILLUMINA%BayesHammer_output*}/LOL-WGS-1.0_combined_R2.fastq.00.0_0.cor.fastq.gz
# --short ${INPUT_DIR_ILLUMINA}/LOL-WGS-1.0_combined_R1.fastq.00.0_0.cor.fastq.gz ${INPUT_DIR_ILLUMINA}/LOL-WGS-1.0_combined_R2.fastq.00.0_0.cor.fastq.gz
