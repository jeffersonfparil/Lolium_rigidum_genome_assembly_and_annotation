#!/bin/bash
#####################################################################
### Quality check and Bayesian error correction of Illumina reads ###
#####################################################################

### Inputs:
### (1) Illumina reads in compressed fastq format (*.fastq.gz)
### Notes: - all reads from 2020.11 were renamed from LOL-WGS2-{1..4}.fastq.gz into simply LOL-WGS-{1..4}.fastq.gz
###        - LOL-WGS-*.fastq.gz from 2020.9 sequencing was renamed to LOL-WGS-0-*.fastq.gz

### Outputs:
### (1) Quality check html and zip files (*_fastqc.html and *_fastqc.zip)
### (2) Error corrected reads in compressed fastq format (*-corrected.fastq.gz)

### Parameters:
DIR=/data/Lolium_rigidum_ASSEMBLY/assembly_annotation_pipeline_tests_20210104/FASTQ/ILLUMINA
FASTQC=/data/Lolium_rigidum_ASSEMBLY/assembly_annotation_pipeline_tests_20210104/FastQC/fastqc
SPADES=/data/Lolium_rigidum_ASSEMBLY/assembly_annotation_pipeline_tests_20210104/SPAdes-3.14.1-Linux/bin/spades.py

### Navigate to working directory
cd $DIR

### Quality check
time parallel ${FASTQC} {} ::: $(ls *.fastq.gz)

### Bayesian error correction (may need to link python3 via: sudo ln -s /usr/bin/python3 /usr/bin/python)
time \
for i in 0 1 2 3 4
do
# i=1
mkdir ${DIR}/BayesHammer_output_WGS-${i}/
time \
$SPADES \
    --only-error-correction \
    --threads 32 \
    --memory 280 \
    -1 ${DIR}/LOL-WGS-${i}_combined_R1.fastq.gz \
    -2 ${DIR}/LOL-WGS-${i}_combined_R2.fastq.gz \
    -o ${DIR}/BayesHammer_output_WGS-${i}/
done

### Clean-up
mkdir QC/
mv *.html QC/
mv *.zip QC/

mkdir BayesHammer_output/
