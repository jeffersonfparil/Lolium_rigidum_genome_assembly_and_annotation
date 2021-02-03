#!/bin/bash
########################################################
### Assembly using Illumina reads via Platanus-allee ###
########################################################

### Inputs:
### (1) Error-corrected Illumina reads in compressed fastq format (BayesHammer_output/*cor.fastq.gz)

### Outputs:
### (1) assembly.log
### (2) Lori_ip_32merFrq.tsv
### (3) Lori_ip_contig.fa


### Parameters:
INPUT_DIR=/data/Lolium_rigidum_ASSEMBLY/assembly_annotation_pipeline_tests_20210104/FASTQ/ILLUMINA/BayesHammer_output/
PLATANUS=/data/Lolium_rigidum_ASSEMBLY/assembly_annotation_pipeline_tests_20210104/Platanus_allee_v2.2.2_Linux_x86_64/platanus_allee
OUTPUT_DIR=/data/Lolium_rigidum_ASSEMBLY/assembly_annotation_pipeline_tests_20210104/ASSEMBLY/

### Intialise and spades output folder
cd ${OUTPUT_DIR}
mkdir Lori_ip/

### testing concatenated all reads
time \
${PLATANUS} \
    assemble \
    -f ${INPUT_DIR}/Lrigidum_illumina_150bp_R1.fastq.gz ${INPUT_DIR}/Lrigidum_illumina_150bp_R2.fastq.gz \
    -t 12 \
    -m 45 \
    -o ${OUTPUT_DIR}/Lori_ip/Lori_ip \
    2> ${OUTPUT_DIR}/Lori_ip/assembly.log








### Iterate across libraries
for i in 0 1.0 1.1 2 3 4 5
do
# i=1.0
mkdir ${i}/
cd ${OUTPUT_DIR}/Lori_ip/${i}/
# ################
# ### assemble ###
# ################
time \
${PLATANUS} \
    assemble \
    -f ${INPUT_DIR}/LOL-WGS-${i}_combined_R1.fastq.00.0_0.cor.fastq.gz ${INPUT_DIR}/LOL-WGS-${i}_combined_R2.fastq.00.0_0.cor.fastq.gz \
    -t 32 \
    -m 200 \
    -o ${OUTPUT_DIR}/Lori_ip/${i}/Lori_ip \
    2> ${OUTPUT_DIR}/Lori_ip/${i}/assembly.log
done
