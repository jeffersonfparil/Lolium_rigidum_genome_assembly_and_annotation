### Count the number of each feature
grep -v "^#" Lolium_rigidum.gff | cut -f 3 | sort | uniq -c > Lolium_rigidum_annotation_counts.txt
### Extract the coordinates of genes
awk '{if ($3=="gene") print}' Lolium_rigidum.gff | cut -f4-5 > Lolium_rigidum_annotation_gene_coordinates.txt
### Find the mean length of gene models
echo 'dat = read.table("Lolium_rigidum_annotation_gene_coordinates.txt", header=FALSE)
MEAN_GENE_MODEL_LENGTH = mean(apply(dat, MARGIN=1, FUN=function(x){abs(diff(x))}))
print(MEAN_GENE_MODEL_LENGTH)
' > find_mean_gene_model_length.R
Rscript find_mean_gene_model_length.R
### Plot genome assembl diagram
time \
julia genome_statistics.jl


### K-mer analyis to estimate genome size
sudo apt install -y jellyfish
gunzip LOL-WGS_combined_R1.fastq.gz
gunzip LOL-WGS_combined_R2.fastq.gz
time \
for k in 17 19 21 25
do
jellyfish count \
    -m ${k} \
    -s 5G \
    -C \
    LOL-WGS_combined_R1.fastq \
    LOL-WGS_combined_R2.fastq \
    -t 32 \
    -o Lolium_rigidum-${k}mer.jf.out.tmp
jellyfish histo \
    Lolium_rigidum-${k}mer.jf.out.tmp \
    -o Lolium_rigidum-${k}mer.jf.out
rm Lolium_rigidum-${k}mer.jf.out.tmp
done

echo '
args = commandArgs(trailingOnly=TRUE)
# args = c("Lolium_rigidum-17mer.jf.out")
f = args[1]
dat = read.delim(f, header=FALSE, sep=" ")
dat = dat[5:(nrow(dat)-1), ] ### remove first 4 and last 1 outliers
idx = which(dat$V2==max(dat$V2))
kmer = dat$V1[idx]
estimated_size = (sum(dat$V1 * dat$V2) / kmer) / 1e9
write.table(estimated_size, file=gsub(".out", ".GBsize", f), col.names=FALSE, row.names=FALSE, quote=FALSE)
' > estimate_genome_size_in_GB.R

for f in $(ls *.mer.jf.out)
do
Rscript estimate_genome_size_in_GB.R ${f}
done
