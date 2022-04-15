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