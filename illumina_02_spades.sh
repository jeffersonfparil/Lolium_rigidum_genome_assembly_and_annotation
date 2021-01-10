#!/bin/bash
################################################
### Assembly using Illumina reads via SPAdes ###
################################################

### Inputs:
### (1) Error-corrected Illumina reads in compressed fastq format (*cor.fastq.gz)

### Outputs:
### (1) scaffolds.fasta - contains resulting scaffolds (recommended for use as resulting sequences)
### (2) contigs.fasta - contains resulting contigs
### (3) assembly_graph.gfa - contains SPAdes assembly graph and scaffolds paths in GFA 1.0 format
### (4) assembly_graph.fastg - contains SPAdes assembly graph in FASTG format


### Parameters:
INPUT_DIR=/data/Lolium_rigidum_ASSEMBLY/assembly_annotation_pipeline_tests_20210104/FASTQ/ILLUMINA/BayesHammer_output/
SPADES=/data/Lolium_rigidum_ASSEMBLY/assembly_annotation_pipeline_tests_20210104/SPAdes-3.14.1-Linux/bin/spades.py
OUTPUT_DIR=/data/Lolium_rigidum_ASSEMBLY/assembly_annotation_pipeline_tests_20210104/ASSMEBLY/

### Intialise and spades output folder
mkdir ${OUTPUT_DIR}/SPADES_Lori_i1/

### Assemble
time \
$SPADES \
    --only-assembler \
    --careful \
    --threads 32 \
    --memory 280 \
    --pe1-1 ${INPUT_DIR}/LOL-WGS-0_combined_R1.fastq.00.0_0.cor.fastq.gz --pe1-2 ${INPUT_DIR}/LOL-WGS-0_combined_R2.fastq.00.0_0.cor.fastq.gz \
    --pe2-1 ${INPUT_DIR}/LOL-WGS-1.0_combined_R1.fastq.00.0_0.cor.fastq.gz --pe2-2 ${INPUT_DIR}/LOL-WGS-1.0_combined_R2.fastq.00.0_0.cor.fastq.gz \
    --pe2-1 ${INPUT_DIR}/LOL-WGS-1.1_combined_R1.fastq.00.0_0.cor.fastq.gz --pe2-2 ${INPUT_DIR}/LOL-WGS-1.1_combined_R2.fastq.00.0_0.cor.fastq.gz \
    --pe3-1 ${INPUT_DIR}/LOL-WGS-2_combined_R1.fastq.00.0_0.cor.fastq.gz --pe3-2 ${INPUT_DIR}/LOL-WGS-2_combined_R2.fastq.00.0_0.cor.fastq.gz \
    --pe4-1 ${INPUT_DIR}/LOL-WGS-3_combined_R1.fastq.00.0_0.cor.fastq.gz --pe4-2 ${INPUT_DIR}/LOL-WGS-3_combined_R2.fastq.00.0_0.cor.fastq.gz \
    --pe5-1 ${INPUT_DIR}/LOL-WGS-4_combined_R1.fastq.00.0_0.cor.fastq.gz --pe5-2 ${INPUT_DIR}/LOL-WGS-4_combined_R2.fastq.00.0_0.cor.fastq.gz \
    --pe6-1 ${INPUT_DIR}/LOL-WGS-5_combined_R1.fastq.00.0_0.cor.fastq.gz --pe6-2 ${INPUT_DIR}/LOL-WGS-5_combined_R2.fastq.00.0_0.cor.fastq.gz \
    -o ${OUTPUT_DIR}/SPADES_Lori_i1/
