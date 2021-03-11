#!/bin/bash
####################################################
### De novo assembly of MinION reads with wtdbg2 ###
####################################################

### Input:
### (1) Filtered and trimmed MinION reads in compressed fastq format (minion-filtered-trimmed.fastq.gz)

### Outputs:
### (1) Assembly using wtdbg2 assembler in fasta format (Lori_mw.raw.fa)
### (2) Scaffold lengths in comma-separated format (Lori_m1-scaffold_stats.csv)
### (3) Assembly statistics in scalabel vector graphic format (Lori_m1-draft_genome_distribution_N50_L50.svg)

### Parameters:
INPUT_FASTQ_GZ=/data/Lolium_rigidum_ASSEMBLY/assembly_annotation_pipeline_tests_20210104/FASTQ/MINION/minion-filtered-trimmed.fastq.gz
OUTPUT_DIR=/data/Lolium_rigidum_ASSEMBLY/assembly_annotation_pipeline_tests_20210104/ASSEMBLY
WTDBG2=/data/Lolium_rigidum_ASSEMBLY/assembly_annotation_pipeline_tests_20210104/wtdbg2/wtdbg2
WTPOA_CNS=/data/Lolium_rigidum_ASSEMBLY/assembly_annotation_pipeline_tests_20210104/wtdbg2/wtpoa-cns
MINIMAP=/data/Lolium_rigidum_ASSEMBLY/assembly_annotation_pipeline_tests_20210104/minimap2/minimap2
PILON=/data/Lolium_rigidum_ASSEMBLY/assembly_annotation_pipeline_tests_20210104/pilon-1.24.jar

### Prepare output directory
cd ${OUTPUT_DIR}/
mkdir Lori_mw/

### De novo assembly with wtdbg2 (assemble ont read >=1G or preset3) ###~40minutes:::12cores:::47Gb:::20210104
time \
${WTDBG2} \
    -i ${INPUT_FASTQ_GZ} \
    -t 12 \
    -fo ${OUTPUT_DIR}/Lori_mw/Lori_mw

### Derive the consensus assembly in fasta format ###~3.5minutes:::12cores:::47Gb:::20210104
time \
${WTPOA_CNS} -t 32 \
             -i ${OUTPUT_DIR}/Lori_mw/Lori_mw.ctg.lay.gz \
             -fo ${OUTPUT_DIR}/Lori_mw/Lori_mw.raw.fa

### Improve the assembly with Pilon
### (1) map reads to the assembly
time \
${MINIMAP} -ax map-ont \
            ${OUTPUT_DIR}/Lori_mw/Lori_mw.raw.fa \
            ${INPUT_FASTQ_GZ} \
            -t 32 \
            > ${OUTPUT_DIR}/Lori_mw/Lori_mw.alignments.sam
### (2) filter, sort, compress, and index the alignments with samtools
time \
samtools view -q 20 -b ${OUTPUT_DIR}/Lori_mw/Lori_mw.alignments.sam | \
    samtools sort > ${OUTPUT_DIR}/Lori_mw/Lori_mw.alignments.bam
time \
samtools index ${OUTPUT_DIR}/Lori_mw/Lori_mw.alignments.bam
### (3) pilon error-correction
time \
java -Xmx100G -jar ${PILON} \
    --genome ${OUTPUT_DIR}/Lori_mw/Lori_mw.raw.fa \
    --unpaired ${OUTPUT_DIR}/Lori_mw/Lori_mw.alignments.bam \
    --outdir ${OUTPUT_DIR}/Lori_mw/ \
    --output Lori_mw.pilon \
    --tracks \
    --threads 32

