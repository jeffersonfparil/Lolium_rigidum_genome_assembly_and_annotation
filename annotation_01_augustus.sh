#!/bin/bash
#######################################
### Genome annotation with Augustus ###
#######################################

### Inputs:
### (1) Genome assemblies in fasta format (Lori_i1, Lori_i2, ... *.fasta or *.fa)
### (2) RNAseq in compressed fastq (*.fastq.gz)

### Outputs:
### (1) Annotations per genome assembly per gene list (${assembly}.${species_gene_list}.gff)
### (2) 
### (3) 

### Parameters:
INPUT_DIR=/data/Lolium_rigidum_ASSEMBLY/assembly_annotation_pipeline_tests_20210104/ASSEMBLY
OUTPUT_DIR=/data/Lolium_rigidum_ASSEMBLY/assembly_annotation_pipeline_tests_20210104/ANNOTATION
AUGUSTUS=/data/Lolium_rigidum_ASSEMBLY/assembly_annotation_pipeline_tests_20210104/Augustus/bin/augustus

### Navigate to working directory
cd $OUTPUT_DIR

### List of genome assemblies
ASSEMBLIES=$(ls ${INPUT_DIR}/ | grep "Lori_")

### List of gene lists we will be using to find homologs in the assemblies
GENE_LISTS=$(echo "rice maize arabidopsis rna")

### Annotation with rice, maize, and arabidopsis genes
### (1) split the genome assemblies by scaffold
echo 'from Bio import SeqIO
import pandas as pd
import sys
import os
fname_assembly = sys.argv[1]
# fname_assembly = "/data/Lolium_rigidum_ASSEMBLY/assembly_annotation_pipeline_tests_20210104/ASSEMBLY/Lori_m1/Lori_m1.raw.fa"
id_assembly = os.path.basename(fname_assembly).split(".")[0]
with open(fname_assembly, "rU") as genome:
    for scaffold in SeqIO.parse(genome, "fasta"):
        seq_string = SeqIO.FastaIO.as_fasta_2line(scaffold)
        f = open(id_assembly + "." + scaffold.id + ".fa", "w")
        f.write(seq_string)
        f.close()
' > split_assembly_by_scaffold.py
time \
parallel python3 split_assembly_by_scaffold.py {} ::: $(ls ${INPUT_DIR}/Lori_*/Lori_*.fa*)
### (2) run Augustus in parallel per scaffold per assembly per species genes
echo '#!/bin/bash
AUGUSTUS=$1
FASTA=$2
SPECIES=$3
# AUGUSTUS=/data/Lolium_rigidum_ASSEMBLY/assembly_annotation_pipeline_tests_20210104/Augustus/bin/augustus
# FASTA=ctg1.fa
# SPECIES=rice #SPECIES=maize #SPECIES=arabidopsis
${AUGUSTUS} \
    --species=${SPECIES} \
    --genemodel=partial \
    ${FASTA} \
    > ${FASTA}.${SPECIES}.gff
' > augustus_parallel.sh
chmod +x augustus_parallel.sh
time \
parallel ./augustus_parallel.sh ${AUGUSTUS} {1} {2} \
    ::: $(ls Lori_*.ctg*.fa) \
    ::: $(echo ${GENE_LISTS} | cut -d' ' -f1-3)

### Annotation with RNAseq data
# Read this:
# http://bioinf.uni-greifswald.de/augustus/binaries/readme.rnaseq.html


### Merge across scaffolds per genome assembly per gene list
echo '#!/bin/bash
assembly=$1
species_gene_list=$2
# assembly=Lori_m1
# species_gene_list=rice
f1=$(ls ${assembly}.*.${species_gene_list}.gff | head -n1)
line_number=$(echo $(grep -n "# ----- prediction" ${f1} | cut -d: -f1) - 1 | bc)
head -${line_number} ${f1} > ${assembly}.${species_gene_list}.gff
for f in $(ls ${assembly}.*.${species_gene_list}.gff)
do
    line_number=$(grep -n "# ----- prediction" ${f} | cut -d: -f1)
    tail -n+${line_number} ${f} >> ${assembly}.${species_gene_list}.gff
done
' > merge_gff_parallel.sh
chmod +x merge_gff_parallel.sh
time \
parallel ./merge_gff_parallel.sh {1} {2} \
    ::: ${ASSEMBLIES} \
    ::: ${GENE_LISTS}


### Clean-up
rm split_assembly_by_scaffold.py
rm augustus_parallel.sh
rm merge_gff_parallel.sh
rm Lori_*.*.fa ### remove fasta per scaffold
rm Lori_*.*.*.gff ### remove annotations per scaffold