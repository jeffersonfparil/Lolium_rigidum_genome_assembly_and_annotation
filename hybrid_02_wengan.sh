#!/bin/bash
#############################################################
### Hybrid assembly using Illumina with MinION via WENGAN ###
#############################################################

### Inputs:
### (1) Error-corrected Illumina reads in compressed fastq format (BayesHammer_output/*cor.fastq.gz)
### (1) MinION reads in compressed fastq format (BayesHammer_output/*cor.fastq.gz)

### Outputs:
### (1) Lori_hw.SPolished.asm.wengan.fasta
### (2) 
### (3) 
### (4) 

### Parameters:
INPUT_DIR_ILLUMINA=/data/Lolium_rigidum_ASSEMBLY/assembly_annotation_pipeline_tests_20210104/FASTQ/ILLUMINA/BayesHammer_output
INPUT_DIR_MINION=/data/Lolium_rigidum_ASSEMBLY/assembly_annotation_pipeline_tests_20210104/FASTQ/MINION
WENGAN=/data/Lolium_rigidum_ASSEMBLY/assembly_annotation_pipeline_tests_20210104/wengan-v0.2-bin-Linux/wengan.pl
OUTPUT_DIR=/data/Lolium_rigidum_ASSEMBLY/assembly_annotation_pipeline_tests_20210104/ASSEMBLY

### Intialise and spades output folder
cd ${OUTPUT_DIR}/
mkdir Lori_hw/

### Assemble: test on ssh_hose 20210217
time \
${WENGAN} \
    -x ontraw \
    -a M \
    -s ${INPUT_DIR_ILLUMINA}/Lrigidum_illumina_150bp_R1.fastq.gz,${INPUT_DIR_ILLUMINA}/Lrigidum_illumina_150bp_R2.fastq.gz \
    -l ${INPUT_DIR_MINION}/minion-filtered-trimmed.fastq.gz \
    -t 30 \
    -g 2000 \
    -p ${OUTPUT_DIR}/Lori_hw/Lori_hw
