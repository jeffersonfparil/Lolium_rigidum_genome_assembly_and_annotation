#!/bin/bash
#############################################################
### Assembly using Illumina reads via gatb-minia-pipeline ###
#############################################################

### Inputs:
### (1) Error-corrected Illumina reads in compressed fastq format (BayesHammer_output/*cor.fastq.gz)

### Outputs:
### (1) Lori_im_final.contigs.fa
### (2) 

### Parameters:
INPUT_DIR=/data/Lolium_rigidum_ASSEMBLY/assembly_annotation_pipeline_tests_20210104/FASTQ/ILLUMINA/BayesHammer_output/
MINIA=/data/Lolium_rigidum_ASSEMBLY/assembly_annotation_pipeline_tests_20210104/gatb-minia-pipeline/gatb
OUTPUT_DIR=/data/Lolium_rigidum_ASSEMBLY/assembly_annotation_pipeline_tests_20210104/ASSEMBLY/

### Intialise and spades output folder
cd ${OUTPUT_DIR}
mkdir Lori_im/

### Execute (Use 75% of the RAM, e.g. from 280Gb to 210Gb; and 80% of the cores, e.g. from 32 cores to 25 cores)
### (use `--continue-scaffolding` flag to continue interrupted scaffolding step)
### IMPORTANT NOTE: CHANGE: `pysam.sort(bwa_output + ".bam", output_path)` into 
###                         `pysam.sort("-o", output_path + ".bam", bwa_output + ".bam")`
time \
${MINIA} \
    -1 ${INPUT_DIR}/Lrigidum_illumina_150bp_R1.fastq.gz \
    -2 ${INPUT_DIR}/Lrigidum_illumina_150bp_R2.fastq.gz \
    --max-memory 210000 \
    --nb-cores 25 \
    -o ${OUTPUT_DIR}/Lori_im/Lori_im

### Assess assembly
./assembly_statistics.sh ${OUTPUT_DIR}/Lori_im/Lori_im_final.contigs.fa
