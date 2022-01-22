#!/bin/bash
############################################################
### Hybrid assembly using MinION with Illumina via HASLR ###
############################################################

### Inputs:
### (1) Error-corrected Illumina reads in compressed fastq format (BayesHammer_output/*cor.fastq.gz)
### (1) MinION reads in compressed fastq format (BayesHammer_output/*cor.fastq.gz)

### Outputs:
### (1) asm.final.fa
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
--threads 32 \
--out ${OUTPUT_DIR}/Lori_hh/ \
--genome 2g \
--long ${INPUT_DIR_MINION}/minion-filtered-trimmed.fastq.gz \
--type nanopore \
--short ${INPUT_DIR_ILLUMINA}/Lrigidum_illumina_150bp_R1.fastq.gz ${INPUT_DIR_ILLUMINA}/Lrigidum_illumina_150bp_R2.fastq.gz

### Rename
cp Lori_hh/asm_contigs*/asm.final.fa Lori_hh/Lori_hh.fasta

### Statistics
./assembly_statistics.sh Lori_hh/Lori_hh.fasta
