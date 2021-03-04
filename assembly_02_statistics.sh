#!/bin/bash

### Assembly statistics
ASSEMBLY_FASTA=$1

### INPUT:
### (1) Genome assembly in fasta format

### OUTPUTS:
### (1) ${ASSEMBLY_FASTA}.base.count - total number of bases and excluding Ns
### (2) ${ASSEMBLY_FASTA}.scaffold_stats.csv - lengths per scaffold
### (3) ${ASSEMBLY_FASTA}.draft_genome_distribution_N50_L50.svg - svg image of scaffold distribution and and assembly statistics

### (1) count the total number of bases (line1: total; line2: excluding Ns)
grep -v ">" ${ASSEMBLY_FASTA} | sed ':a;N;$!ba;s/\n//g' | wc -c > ${ASSEMBLY_FASTA}.base.count
# cat ${ASSEMBLY_FASTA} | paste - - | cut -f2 | sed ':a;N;$!ba;s/\n//g' | wc -c > ${ASSEMBLY_FASTA}.base.count ### similarly with paste|cut
grep -v ">" ${ASSEMBLY_FASTA} | grep -v "N" | wc -c >> ${ASSEMBLY_FASTA}.base.count

### (2) count the length of each scaffold (generic script in python; a bash command like the one above will work faster but the use of the 'len' grep search key is not generic)
echo 'from Bio import SeqIO
import pandas as pd
import sys
fname_assembly = sys.argv[1]
# fname_assembly = "/data/Lolium_rigidum_ASSEMBLY/assembly_annotation_pipeline_tests_20210104/ASSEMBLY/Lori_m1.raw.fa"
id_len = []
with open(fname_assembly, "rU") as genome:
    for scaffold in SeqIO.parse(genome, "fasta"):
        id_len.append([str(scaffold.id), len(scaffold), len(str(scaffold.seq).replace("N", ""))])
out = pd.DataFrame(id_len)
out.columns = ["SCAFFOLD", "LENGTH", "LENGTH_NO_Ns"]
out = out.sort_values(["LENGTH"], ascending=[0])
out.to_csv(fname_assembly + ".scaffold_stats.csv", header=True, sep=",")
' > extract_scaffold_stats.py
python3 extract_scaffold_stats.py ${ASSEMBLY_FASTA} && \
	echo "Scaffold lengths: ${ASSEMBLY_FASTA}.scaffold_stats.csv" || \
	echo "Error. Install Bio and/or pandas via pip3 install BioPython and pip3 install pandas."

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
genome_size_no_Ns = sum(dat$LENGTH_NO_Ns)
dat = dat[order(dat$LENGTH, decreasing=FALSE),]
dat$CUMM_LENGTH = cumsum(dat$LENGTH)
dat$CUMM_COVER = dat$CUMM_LENGTH*100 / genome_size
svg(str_replace(fname_scaffold_stats, ".scaffold_stats.csv", ".draft_genome_distribution_N50_L50.svg"), width=15, height=8)
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
legend("right", legend=c(paste0("Assembly size = ", round(genome_size/1e9, 2), " Gb"),
						 paste0("Assembly size  (Ns excluded) = ", round(genome_size_no_Ns/1e9, 2), " Gb (", round((genome_size-genome_size_no_Ns)*100/genome_size, 4), "% Ns)"),
						 paste0("N50 = ", round(N50/1000), " kb"),
						 paste0("L50 = ", round(L50), " scaffolds")
						 ))
dev.off()
'  > extract_scaffold_stats.r
Rscript extract_scaffold_stats.r ${ASSEMBLY_FASTA}.scaffold_stats.csv && \
	echo "Assembly statistics: ${ASSEMBLY_FASTA%.*}.draft_genome_distribution_N50_L50.svg" || \
	echo "Error. Install RColorBrewer and/or stringr via > install.packages('RColorBrewer') and > install.packages('stringr')."

### Clean up
rm extract_scaffold_stats.*
