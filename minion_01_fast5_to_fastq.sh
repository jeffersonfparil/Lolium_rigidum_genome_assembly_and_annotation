#!/bin/bash
################################################################
### Basecalling MinION reads into fastq sequence information ###
################################################################

### Inputs:
### (1) MinION reads in fast5 format (*.fast5)
### Note: For subsequent MinION runs, modify the INPUT directory to point to the newly generated fast5 files

### Outputs:
### (1) MinION reads in compressed fastq format (minion.fastq.gz)
### (2) Quality check html and zip files (minion_fastqc.html and minion_fastqc.zip)

### Parameters:
FASTQC=/data/Lolium_rigidum_ASSEMBLY/assembly_annotation_pipeline_tests_20210104/FastQC/fastqc
INPUT=/data/Lolium_rigidum_ASSEMBLY/assembly_annotation_pipeline_tests_20210104/FAST5/
OUTPUT=/data/Lolium_rigidum_ASSEMBLY/assembly_annotation_pipeline_tests_20210104/FASTQ/MINION/

### Navigate to the output directory
cd $OUTPUT

### Basecalling with guppy ###~5 days 17 hours:::32cores:::~280Gb:::20210104
time \
guppy_basecaller \
    --input_path ${INPUT} \
    --save_path ${OUTPUT} \
    --flowcell FLO-MIN106 \
    --kit SQK-LSK109 \
    --cpu_threads_per_caller 32

### Concatenate all fastq into a single file and compress
cat ${OUTPUT}/*.fastq > ${OUTPUT}/minion.fastq

### Count the total number of bases sequenced
cat ${OUTPUT}/minion.fastq | paste - - - - | cut -f 2 | tr -d '\n' | wc -c > minion.base.count

### Compress
gzip ${OUTPUT}/minion.fastq

### Quality check
time ${FASTQC} -t 12 ${OUTPUT}/minion.fastq.gz

### Clean-up
mkdir guppy_output/
mv *.fastq guppy_output/
mv *.log guppy_output/
mv sequencing_* guppy_output/
