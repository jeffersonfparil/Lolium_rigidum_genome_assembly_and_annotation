#!/bin/bash
####################################################
### De novo assembly of MinION reads with wtdbg2 ###
####################################################

### Input:
### (1) Filtered and trimmed MinION reads in compressed fastq format (minion-filtered-trimmed.fastq.gz)

### Outputs:
### (1) Assembly using wtdbg2 assembler in fasta format (Lori_m1.raw.fa)
### (2) Scaffold lengths in comma-separated format (Lori_m1-scaffold_stats.csv)
### (3) Assembly statistics in scalabel vector graphic format (Lori_m1-draft_genome_distribution_N50_L50.svg)

### Parameters:
INPUT_FASTQ_GZ=/data/Lolium_rigidum_ASSEMBLY/assembly_annotation_pipeline_tests_20210104/FASTQ/MINION/minion-filtered-trimmed.fastq.gz
OUTPUT_BASENAME=/data/Lolium_rigidum_ASSEMBLY/assembly_annotation_pipeline_tests_20210104/ASSEMBLY/Lori_m1
WTDBG2=/data/Lolium_rigidum_ASSEMBLY/assembly_annotation_pipeline_tests_20210104/wtdbg2/wtdbg2

### De novo assembly with wtdbg2 (assemble ont read >=1G or preset3) ###~40minutes:::12cores:::47Gb:::20210104
time \
$WTDBG2 \
    -p 19 \
    -AS 2 \
    -s 0.05 \
    -L 5000 \
    -i ${INPUT_FASTQ_GZ} \
    -t 32 \
    -fo ${OUTPUT_BASENAME}

### Derive the consensus assembly in fasta format ###~3.5minutes:::12cores:::47Gb:::20210104
time \
wtdbg2/wtpoa-cns -t 32 \
                 -i ${OUTPUT_BASENAME}.ctg.lay.gz \
                 -fo ${OUTPUT_BASENAME}.raw.fa

### Assembly statistics
### (1) count the total number of bases ###56,430,521:::20210112
grep "len" ${OUTPUT_BASENAME}.raw.fa | cut -d'=' -f2 | paste -s -d'+' | bc > ${OUTPUT_BASENAME}.base.count
### (2) count the length of each scaffold (generic script in python; a bash command like the one above will work faster but the use of the 'len' grep search key is not generic)
echo 'from Bio import SeqIO
import pandas as pd
import sys
fname_assembly = sys.argv[1]
# fname_assembly = "/data/Lolium_rigidum_ASSEMBLY/assembly_annotation_pipeline_tests_20210104/ASSEMBLY/Lori_m1.raw.fa"
id_len = []
with open(fname_assembly, "rU") as genome:
    for scaffold in SeqIO.parse(genome, "fasta"):
        id_len.append([str(scaffold.id), len(scaffold)])
out = pd.DataFrame(id_len)
out.columns = ["SCAFFOLD", "LENGTH"]
out = out.sort_values(["LENGTH"], ascending=[0])
out.to_csv(fname_assembly.split(".")[0] + "-scaffold_stats.csv", header=True, sep=",")
' > extract_scaffold_stats.py
python3 extract_scaffold_stats.py ${OUTPUT_BASENAME}.raw.fa
### (3) calculate N50, L50 and plot
echo '# generate histogram of scaffold sizes of the Lolium perenne reference genome
library(RColorBrewer)
library(stringr)
args = commandArgs(trailingOnly=TRUE)
fname_scaffold_stats = args[1]
# fname_scaffold_stats = "/data/Lolium_rigidum_ASSEMBLY/assembly_annotation_pipeline_tests_20210104/ASSEMBLY/Lori_m1-scaffold_stats.csv"
dat = read.csv(fname_scaffold_stats)
str(dat)
set1_colours = brewer.pal(9, "Set1")
pastel_colours = brewer.pal(9, "Pastel1")
min_size = min(dat$LENGTH)
max_size = max(dat$LENGTH)
mean_size = mean(dat$LENGTH)
genome_size = sum(dat$LENGTH)
dat = dat[order(dat$LENGTH, decreasing=FALSE),]
dat$CUMM_LENGTH = cumsum(dat$LENGTH)
dat$CUMM_COVER = dat$CUMM_LENGTH*100 / genome_size
svg(str_replace(fname_scaffold_stats, "-scaffold_stats.csv", "-draft_genome_distribution_N50_L50.svg"), width=15, height=8)
par(mfrow=c(1,2), mar=c(7,7,2,2))
hist(dat$LENGTH, col=pastel_colours[1], border=set1_colours[1], nclass=20, xlab="Scaffold Size (bp)", main="", las=2)
legend("right", legend=c(paste0("Total number of scaffolds = ", nrow(dat)),
						 paste0("Minimum scaffold length = ", round(min_size), " bp"),
						 paste0("Maximum scaffold length = ", round(max_size/1000), " kb"),
						 paste0("Average scaffold length = ", round(mean_size/1000), " kb")))
plot(dat$LENGTH, dat$CUMM_COVER, xlab="Scaffold Size (bp)", ylab="Assembly Covered (%)", xlim=c(0, max_size), ylim=c(0, 100), type="l", lty=1, lwd=3, col="Gray1")
N50 = min(dat$LENGTH[dat$CUMM_COVER>=50]) #contigs with size >= to this N50 value - cover 50% of the assembly
L50 = min( c(nrow(subset(dat, CUMM_COVER<=50)), nrow(subset(dat, CUMM_COVER>=50))) ) #the minimum number of contigs that cover 50% of the assembly
x1 = N50
y1 = 50
segments(x0=x1, x1=x1, y0=0, y1=y1, lty=2, lwd=2, col=set1_colours[1])
segments(x0=0, x1=x1, y0=50, y1=50, lty=2, lwd=2, col=set1_colours[1])
grid()
legend("right", legend=c(paste0("Assembly size = ", round(genome_size/1000000000, 2), " Gb"),
						 paste0("N50 = ", round(N50/1000), " kb"),
						 paste0("L50 = ", round(L50), " scaffolds")
						 ))
dev.off()
'  > extract_scaffold_stats.r
Rscript extract_scaffold_stats.r ${OUTPUT_BASENAME}-scaffold_stats.csv

### Clean-up
rm extract_scaffold_stats.*
mkdir Lori_m1/
mv Lori_m1.* Lori_m1/
mv Lori_m1-* Lori_m1/

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