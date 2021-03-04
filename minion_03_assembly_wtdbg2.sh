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

### Assembly statistics
### FIX ME WITH SRC PATH??!?!?!? 
### ${path_to_src}/assembly_stats.sh ${OUTPUT_BASENAME}.raw.fa

### Miscellaneous
### polish consensus, not necessary if you want to polish the assemblies using other tools
# minimap2 -t16 \
#          -ax map-pb \
#          -r2k lolium5.raw.fa FASTQ/lolium5-porechoped-nanofilted.fastq | \
#          samtools sort -@4 > lolium5.bam
# samtools view -F0x900 lolium5.bam | \
#          wtdbg2/wtpoa-cns -t 16 \
#                           -d lolium5.raw.fa \
#                           -i - \
#                           -fo lolium5.cns.fa
# # Addtional polishment using short reads
# bwa index lolium5.cns.fa
# bwa mem -t 16 \
#         lolium5.cns.fa \
#         sr.1.fa \
#         sr.2.fa | \
#         samtools sort -O SAM | \
#         wtdbg2/wtpoa-cns -t 16 \
#                          -x sam-sr \
#                          -d lolium5.cns.fa \
#                          -i - \
#                          -fo lolium5.srp.fa