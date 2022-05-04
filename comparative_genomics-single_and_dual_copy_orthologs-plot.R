### Plot figure 1

args = commandArgs(trailingOnly=TRUE)
# args = c("ORTHOGROUPS_SINGLE_GENE.NT.timetree.nex", "CONTRACTION_EXPANSION.txt", "orthogroups_summarised_gene_counts.csv", "orthogroups_gene_counts_families_go.out", ".4DTv")
fname_tree = args[1]
fname_conex = args[2]
fname_gene_groups = args[3]
fname_gene_counts = args[4]
extension_name_4DTv = args[5]

library(ape)
library(gplots)

par(mfrow=c(2,2))

### Tree: dendrogram
tree = read.nexus(fname_tree)
tree = ladderize(tree, right=FALSE)
par(mar=c(5,2,5,1))
plt = plot.phylo(tree, cex=1.2)
x_axis = round(seq(0, max(tree$edge.length), by=20))
axis(side=1, line=1.5, at=max(x_axis)-x_axis, lab=x_axis)
mtext(text="Million years ago", side=1, line=4.5, at=median(x_axis))

### Expansion / Contraction: middle area text
adj_frac = 0.04
conex = read.table(fname_conex, header=TRUE)
conex = conex[order(conex$Species), ]
tips.order = unlist(lapply(conex$Species, FUN=function(x) {
                            which(x == c("Arabidopsis_thaliana",
                                         "Oryza_sativa",
                                         "Sorghum_bicolor",
                                         "Zea_mays",
                                         "Secale_cereale",
                                         "Lolium_perenne",
                                         "Lolium_rigidum"))}))
conex$order = tips.order
conex = conex[order(conex$order), ]
conex$Expansion = formatC(conex$Expansion, format="d", big.mark=",")
conex$Contraction = formatC(conex$Contraction, format="d", big.mark=",")
conex_lab = paste0(conex$Expansion, " : ", conex$Contraction)
text(x=(plt$x.lim[2]-(adj_frac*plt$x.lim[2])), y=seq(plt$y.lim[1], plt$y.lim[2]), adj=0.5, lab=conex_lab, cex=1.2)
mtext(side=3, line=1, at=(plt$x.lim[2]-(adj_frac*plt$x.lim[2])), adj=0.5, text="Expansion : Contraction")

### Gene classifications: bar plot
gene_groups = read.csv(fname_gene_groups)
m = ncol(gene_groups)
gene_groups = gene_groups[order(gene_groups$Species), ]
gene_groups$order = tips.order
gene_groups = gene_groups[order(gene_groups$order), ]
X = t(as.matrix(gene_groups[, 3:m]))
rownames(X) = gsub("_", " ", colnames(gene_groups)[3:m])
colnames(X) = gene_groups$Species
colors = c("#a6cee3", "#1f78b4", "#b2df8a", "#33a02c")
par(mar=c(3.5,1,3.5,15))
barplot(X, col=colors, bord=NA, horiz=TRUE, yaxt="n", xaxt="n", xlim=c(0, signif(max(gene_groups$Total),0)),
        legend.text=TRUE, args.legend=list(x="bottomright", inset=c(-0.15, +0.05), bty="n", cex=1.2))
x_axis = seq(0, signif(max(gene_groups$Total),0), length=5)
axis(side=1, at=x_axis, lab=formatC(x_axis, format="d", big.mark=","))
mtext(text="Gene counts", side=1, line=3, at=median(x_axis))

### Venn diagram of shared gene families
gene_counts = read.delim(fname_gene_counts, header=TRUE)
X = gene_counts[, 1:(ncol(gene_counts)-4)]
X$Orthogroup = as.numeric(gsub("OG", "", X$Orthogroup))+1
X[,2:ncol(X)] = X[,2:ncol(X)] > 0
colnames(X) = gsub("_", " ", colnames(X))
par(mar=c(1,5,3,5))
venn(X[,c(3,4,5,6,8)]) ### picking only 5 species (maximum number of sets to draw a Venn diagram so far)

### Distribution of 4DTv (fraction of transverions among 4-fold degenerate codons - correlated with time from whole genome duplication using dual-copy paralogs and single-copy orthologs)
par(mar=c(5,5,5,5))
alpha = 0.05
FDTv_files = list.files(path=".", pattern=extension_name_4DTv)
for (f in FDTv_files){
    # f = FDTv_files[3]
    df = read.delim(f, header=FALSE, na.string="NA")
    x = df$V4
    hist(x)
}
