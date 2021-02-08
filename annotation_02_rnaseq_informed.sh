#!/bin/bash
##############################################
### Trinity de novo transcriptome assembly ###
##############################################

### Inputs:
### (1) RNAseq in compressed fastq (*.fastq.gz)
### (2) 

### Outputs:
### (1) Annotations per genome assembly per gene list (${assembly}.${species_gene_list}.gff)
### (2) 
### (3) 

### Parameters:
INPUT_DIR=/data/Lolium_rigidum_ASSEMBLY/assembly_annotation_pipeline_tests_20210104/ASSEMBLY
OUTPUT_DIR=/data/Lolium_rigidum_ASSEMBLY/assembly_annotation_pipeline_tests_20210104/ANNOTATION
AUGUSTUS=/data/Lolium_rigidum_ASSEMBLY/assembly_annotation_pipeline_tests_20210104/Augustus/bin/augustus

### Navigate to working directory
mkdir ${OUTPUT_DIR}/TRANSCRIPTOME/
OUTPUT_DIR=${OUTPUT_DIR}/TRANSCRIPTOME/
cd $OUTPUT_DIR

### Clean-up
