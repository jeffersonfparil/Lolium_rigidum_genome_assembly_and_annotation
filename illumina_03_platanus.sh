#!/bin/bash
########################################################
### Assembly using Illumina reads via Platanus-allee ###
########################################################

### Inputs:
### (1) Error-corrected Illumina reads in compressed fastq format (BayesHammer_output/*cor.fastq.gz)

### Outputs:
### (1) out_allPhaseBlock.fa
### (2) out_consensusScaffolds.fa
### (3) out_consensusScaffoldComponent.bed 


### Parameters:
INPUT_DIR=/data/Lolium_rigidum_ASSEMBLY/assembly_annotation_pipeline_tests_20210104/FASTQ/ILLUMINA/BayesHammer_output/
PLATANUS=/data/Lolium_rigidum_ASSEMBLY/assembly_annotation_pipeline_tests_20210104/Platanus_allee_v2.2.2_Linux_x86_64/platanus_allee
OUTPUT_DIR=/data/Lolium_rigidum_ASSEMBLY/assembly_annotation_pipeline_tests_20210104/ASSEMBLY/

### Intialise and spades output folder
cd ${OUTPUT_DIR}
mkdir PLATANUS_Lori_i2/

time \
${PLATANUS} \
    assemble \
    -f ${INPUT_DIR}/Lrigidum_illumina_150bp_R1.fastq.gz ${INPUT_DIR}/Lrigidum_illumina_150bp_R2.fastq.gz \
    -k 32 \
    -t 32 \
    -m 200 \
    -o ${OUTPUT_DIR}/PLATANUS_Lori_i2/Lori_i2 \
    2> ${OUTPUT_DIR}/PLATANUS_Lori_i2/assembly.log

time \
${PLATANUS} \
    phase \
    -c ${OUTPUT_DIR}/PLATANUS_Lori_i2/Lori_i2_contig.fa \
    -IP1 ${INPUT_DIR}/Lrigidum_illumina_150bp_R1.fastq.gz ${INPUT_DIR}/Lrigidum_illumina_150bp_R2.fastq.gz \
    -t 12 \
    -o Lori_i2 \
    2>phase.log
mv Lori_i2_*.* ${OUTPUT_DIR}/PLATANUS_Lori_i2/
mv  Lori_i2_intermediateResults/ Lori_i2_intermediateResults_phasing/
mv  Lori_i2_intermediateResults_phasing/ ${OUTPUT_DIR}/PLATANUS_Lori_i2/
mv phase.log ${OUTPUT_DIR}/PLATANUS_Lori_i2/

platanus_allee \
    consensus \
    -c ${OUTPUT_DIR}/PLATANUS_Lori_i2/Lori_i2_primaryBubble.fa \
    -IP1 ${INPUT_DIR}/Lrigidum_illumina_150bp_R1.fastq.gz ${INPUT_DIR}/Lrigidum_illumina_150bp_R2.fastq.gz \
    2>consensus.log
