#!/bin/bash
#####################################################################
### Quality check and Bayesian error correction of Illumina reads ###
#####################################################################

### Inputs:
### (1) Illumina reads in compressed fastq format (*.fastq.gz)
### Notes: - all reads from 2020.11 were renamed from LOL-WGS2-{2..5}.fastq.gz into simply LOL-WGS-{2..5}.fastq.gz
###        - LOL-WGS-*.fastq.gz from 2020.9 sequencing was split into 2 and renamed to LOL-WGS-{0..1}-*.fastq.gz
###        - LOL-WGS-1* was further dived into 2 because 280GB of RAM is insufficient, i.e. LOL-WGS-1.{0..1}*

### Outputs:
### (1) Quality check html and zip files (*_fastqc.html and *_fastqc.zip)
### (2) Error corrected reads in compressed fastq format (*cor.fastq.gz)

### Parameters:
DIR=/data/Lolium_rigidum_ASSEMBLY/assembly_annotation_pipeline_tests_20210104/FASTQ/ILLUMINA
FASTQC=/data/Lolium_rigidum_ASSEMBLY/assembly_annotation_pipeline_tests_20210104/FastQC/fastqc
SPADES=/data/Lolium_rigidum_ASSEMBLY/assembly_annotation_pipeline_tests_20210104/SPAdes-3.14.1-Linux/bin/spades.py

### Navigate to working directory
cd $DIR

### Quality check
time parallel ${FASTQC} {} ::: $(ls *.fastq.gz)
mkdir QC/
mv *.html QC/
mv *.zip QC/

### Bayesian error correction (may need to link python3 via: sudo ln -s /usr/bin/python3 /usr/bin/python)
time \
for i in 0 1.0 1.1 2 3 4 5
do
# i=1
mkdir BayesHammer_output_WGS-${i}/
time \
$SPADES \
    --only-error-correction \
    --threads 32 \
    --memory 280 \
    -1 LOL-WGS-${i}_combined_R1.fastq.gz \
    -2 LOL-WGS-${i}_combined_R2.fastq.gz \
    -o BayesHammer_output_WGS-${i}/
done

### Clean-up
mkdir BayesHammer_output/
mv BayesHammer_output_WGS-*/corrected/*.fastq.*.fastq.gz BayesHammer_output/ ### exclude the unpaired fastq.gz file from each read-pair

###@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
###@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
###@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
### Repeat LOL-WGS-1* but divide into 2 first so that we don;t run into "<jemalloc>: Error in malloc(): out of memory"
# time parallel gzip -d {} ::: $(ls *.fastq.gz)

# wc -l *.fastq

# head -n 220242640 LOL-WGS-1_combined_R1.fastq > 1.0-R1.fastq
# head -n 220242640 LOL-WGS-1_combined_R2.fastq > 1.0-R2.fastq

# tail -n 220242640 LOL-WGS-1_combined_R1.fastq > 1.1-R1.fastq
# tail -n 220242640 LOL-WGS-1_combined_R2.fastq > 1.1-R2.fastq

# gzip 1.0-R1.fastq
# gzip 1.0-R2.fastq
# gzip 1.1-R1.fastq
# gzip 1.1-R2.fastq

# mv 1.0-R1.fastq.gz LOL-WGS-1.0_combined_R1.fastq.gz
# mv 1.0-R2.fastq.gz LOL-WGS-1.0_combined_R2.fastq.gz
# mv 1.1-R1.fastq.gz LOL-WGS-1.1_combined_R1.fastq.gz
# mv 1.1-R2.fastq.gz LOL-WGS-1.1_combined_R2.fastq.gz
time \
for i in 1.0 1.1
do
# i=1.0
mkdir BayesHammer_output_WGS-${i}/
time \
$SPADES \
    --only-error-correction \
    --threads 32 \
    --memory 280 \
    -1 LOL-WGS-${i}_combined_R1.fastq.gz \
    -2 LOL-WGS-${i}_combined_R2.fastq.gz \
    -o BayesHammer_output_WGS-${i}/
done
