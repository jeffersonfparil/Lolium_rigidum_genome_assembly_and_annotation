#!/bin/bash
#############################################################
### Assembly using Illumina reads via gatb-minia-pipeline ###
#############################################################

### Inputs:
### (1) Error-corrected Illumina reads in compressed fastq format (BayesHammer_output/*cor.fastq.gz)

### Outputs:
### (1) assembly.fasta
### (2) 

### Parameters:
INPUT_DIR=/data/Lolium_rigidum_ASSEMBLY/assembly_annotation_pipeline_tests_20210104/FASTQ/ILLUMINA/BayesHammer_output/
MINIA=/data/Lolium_rigidum_ASSEMBLY/assembly_annotation_pipeline_tests_20210104/gatb-minia-pipeline/gatb
OUTPUT_DIR=/data/Lolium_rigidum_ASSEMBLY/assembly_annotation_pipeline_tests_20210104/ASSEMBLY/

### Intialise and spades output folder
cd ${OUTPUT_DIR}
mkdir Lori_im/

### Execute
time \
${MINIA} \
    -1 ${INPUT_DIR}/Lrigidum_illumina_150bp_R1.fastq.gz \
    -2 ${INPUT_DIR}/Lrigidum_illumina_150bp_R2.fastq.gz \
    --max-memory 40000 \
    --nb-cores 12 \
    -o ${OUTPUT_DIR}/Lori_im/Lori_im

