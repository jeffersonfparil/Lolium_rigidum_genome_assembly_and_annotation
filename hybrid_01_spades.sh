#!/bin/bash
############################################################
### Hybrid assembly using Illumina and MinION via SPAdes ###
############################################################

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
SPADES=/data/Lolium_rigidum_ASSEMBLY/assembly_annotation_pipeline_tests_20210104/SPAdes-3.14.1-Linux/bin/spades.py
OUTPUT_DIR=/data/Lolium_rigidum_ASSEMBLY/assembly_annotation_pipeline_tests_20210104/ASSEMBLY/

### Intialise and spades output folder
cd ${OUTPUT_DIR}
mkdir HYBRID_SPADES_Lori/

### Hybrid assembly per pair of fastq read files
for i in 0 1.0 1.1 2 3 4 5
do
# i=0
mkdir ${OUTPUT_DIR}/HYBRID_SPADES_Lori/${i}/
time \
$SPADES \
    --only-assembler \
    --isolate \
    --threads 32 \
    --memory 280 \
    --pe1-1 ${INPUT_DIR_ILLUMINA}/LOL-WGS-${i}_combined_R1.fastq.00.0_0.cor.fastq.gz \
    --pe1-2 ${INPUT_DIR_ILLUMINA}/LOL-WGS-${i}_combined_R2.fastq.00.0_0.cor.fastq.gz \
    --nanopore ${INPUT_DIR_MINION}/minion-filtered-trimmed.fastq.gz \
    -o ${OUTPUT_DIR}/HYBRID_SPADES_Lori/${i}/
done
