#!/bin/bash
#####################################################################
### Quality check and Bayesian error correction of Illumina reads ###
#####################################################################

### Parameters:
DIR=/data/Lolium_rigidum_ASSEMBLY/jeff_lolium/FASTQ/ILLUMINA
FASTQC=/data/Lolium_rigidum_ASSEMBLY/jeff_lolium/FastQC/fastqc
SPADES=/data/Lolium_rigidum_ASSEMBLY/jeff_lolium/SPAdes-3.14.1-Linux/bin/spades.py

### Navigate to working directory
cd $DIR


### Quality check
time parallel ${FASTQC} {} ::: $(ls *.fastq.gz)
mkdir QC/
mv *.html QC/
mv *.zip QC/

### Bayesian error correction
mkdir ${DIR}/BayesHammer_out_test/
time \
${SPADES} \
    --only-error-correction \
    -1 ${DIR}/LOL-WGS2-2_combined_R1.fastq.gz \
    -2 ${DIR}/LOL-WGS2-2_combined_R2.fastq.gz \
    -o ${DIR}/BayesHammer_out_test/
### This test data needs at least 146 GB of RAM! 2020-12-29

