#!/bin/bash

### Assembly error correction with pilon or quiver/arrow

ASSEMBLY_FASTA=$1

### INPUT:
### (1) Genome assembly in fasta format

### OUTPUTS:
### (1) 

### Pilon
## install
DIR=/data/Lolium_rigidum_ASSEMBLY/assembly_annotation_pipeline_tests_20210104/
cd $DIR
# (1) install pilon
wget https://github.com/broadinstitute/pilon/releases/download/v1.24/pilon-1.24.jar
# (2) install minimap2
git clone https://github.com/lh3/minimap2
cd minimap2 && make
cd -
# (3) install bwa
sudo apt install bwa


### Genome assembly
ASSEMBLY=/data/Lolium_rigidum_ASSEMBLY/assembly_annotation_pipeline_tests_20210104/ASSEMBLY/Lori_mw/Lori_mw.raw.fa
### Navigate to working directory
ASSEMBLY_DIR=$(dirname $ASSEMBLY)
cd ${ASSEMBLY_DIR}
### Setup reference genome assembly name
ASSEMBLY_NAME=$(echo ${ASSEMBLY_DIR} | rev | cut -d'/' -f1 | rev)
### Index the assmebly to be used as reference for read alignment
bwa index -p ${ASSEMBLY_NAME} -a bwtsw ${ASSEMBLY}

## Align reads
# (1) minion reads
MINION_READS=/data/Lolium_rigidum_ASSEMBLY/assembly_annotation_pipeline_tests_20210104/FASTQ/MINION/minion-filtered-trimmed.fastq.gz
MINIMAP=/data/Lolium_rigidum_ASSEMBLY/assembly_annotation_pipeline_tests_20210104/minimap2/minimap2
time \
${MINIMAP} -ax map-ont \
    ${ASSEMBLY} \
    ${MINION_READS} \
    > ${ASSEMBLY_NAME}-MINION_READS.sam
time \
samtools view -q ${MAPQ} -b ${ASSEMBLY_NAME}-MINION_READS.sam | \
    samtools sort > ${ASSEMBLY_NAME}-MINION_READS.bam

# (2) illumina reads
time \
bwa mem ${FNAME_REF} ${FNAME_READ1} ${FNAME_READ2} | \
    samtools view -q ${MAPQ} -b | \
    samtools sort > ${FNAME_OUT}.bam


### Pilon error correction
PILON=/data/Lolium_rigidum_ASSEMBLY/assembly_annotation_pipeline_tests_20210104/pilon-1.24.jar
time \
java -Xmx16G -jar ${PILON} \
    --genome ${ASSEMBLY} \
    --frags ${ASSEMBLY_NAME}-MINION_READS.bam




### Quiver/Arrow
git clone https://github.com/PacificBiosciences/GenomicConsensus.git
