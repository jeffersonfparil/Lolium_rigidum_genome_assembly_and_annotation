#!/bin/bash
################################################################
### Basecalling minION reads into fastq sequence information ###
################################################################

### Parameters:
GUPPY=/data/Lolium_rigidum_ASSEMBLY/jeff_lolium/GUPPY/ont-guppy-cpu/bin/guppy_basecaller
FASTQC=/data/Lolium_rigidum_ASSEMBLY/jeff_lolium/FastQC/fastqc
INPUT=/data/Lolium_rigidum_ASSEMBLY/jeff_lolium/FAST5/
OUTPUT=/data/Lolium_rigidum_ASSEMBLY/jeff_lolium/FASTQ/MINION/

### Basecalling with guppy
time \
$GUPPY \
    --input_path ${INPUT} \
    --save_path ${OUTPUT} \
    --flowcell FLO-MIN106 \
    --kit SQK-LSK109 \
    --cpu_threads_per_caller 32

### Quality check
time \
FASTQC \
    ${OUTPUT}/*.fastq
