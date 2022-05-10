### Plot figure 2

### EXECUTION
# time \
# Rscript comparative_genomics-single_and_dual_copy_orthologs-plot.R \
#     /home/jeffersonfparil/Downloads/data/Lolium_genome_comparative_genomics/ORTHOGROUPS_SINGLE_GENE.NT.timetree.nex \
#     /home/jeffersonfparil/Downloads/data/Lolium_genome_comparative_genomics/CONTRACTION_EXPANSION.txt \
#     /home/jeffersonfparil/Downloads/data/Lolium_genome_comparative_genomics/orthogroups_summarised_gene_counts.csv \
#     /home/jeffersonfparil/Downloads/data/Lolium_genome_comparative_genomics/orthogroups_gene_counts_families_go.out \
#     /home/jeffersonfparil/Downloads/data/Lolium_genome_comparative_genomics \
#     .4DTv \
#     /home/jeffersonfparil/Downloads/data/Lolium_genome_comparative_genomics/ORTHOGROUPS_SINGLE_GENE.NT.4DTv \
#     test.svg

args = commandArgs(trailingOnly=TRUE)
# args = c("ORTHOGROUPS_SINGLE_GENE.NT.timetree.nex", "CONTRACTION_EXPANSION.txt", "orthogroups_summarised_gene_counts.csv", "orthogroups_gene_counts_families_go.out", "/home/jeffersonfparil/Downloads/data/Lolium_genome_comparative_genomics", ".4DTv", "ORTHOGROUPS_SINGLE_GENE.NT.4DTv", "test.svg")
fname_tree = args[1]
fname_conex = args[2]
fname_gene_groups = args[3]
fname_gene_counts = args[4]
dir_name_4DTv = args[5]
extension_name_4DTv = args[6]
fname_4DTv_singlecopy = args[7]
fname_svg_output = args[8]

library(ape)
library(gplots)

svg(fname_svg_output, width=12, height=9)

# par(mfrow=c(2,2))
layout(matrix(c(rep(1,times=6), rep(2,times=2), rep(3,times=6),
                rep(4,times=7), rep(5,times=7)), byrow=TRUE, nrow=2))

### Tree: dendrogram
tree = read.nexus(fname_tree)
tree = ladderize(tree, right=FALSE)
par(mar=c(5,2,5,0))
plt = plot.phylo(tree, cex=1.2)
x_axis = round(seq(0, max(tree$edge.length), by=20))
axis(side=1, line=1.5, at=max(x_axis)-x_axis, lab=x_axis)
mtext(text="Million years ago", side=1, line=4.5, at=median(x_axis))

### Expansion / Contraction: middle area text
par(mar=c(5,0,5,0))
plot(x=plt$x.lim, y=plt$y.lim, type="n", bty="n", xaxt="n", yaxt="n", xlab="", ylab="")
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
# text(x=(plt$x.lim[2]-(adj_frac*plt$x.lim[2])), y=seq(plt$y.lim[1], plt$y.lim[2]), adj=0.5, lab=conex_lab, cex=1.2)
text(x=median(plt$x.lim), y=seq(plt$y.lim[1], plt$y.lim[2]), adj=0.5, lab=conex_lab, cex=1.2)
mtext(side=3, line=1, at=median(plt$x.lim), adj=0.5, text="Expansion : Contraction")

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
par(mar=c(3.5, 1.0, 3.5, 8.0))
barplot(X, col=colors, bord=NA, horiz=TRUE, yaxt="n", xaxt="n", xlim=c(0, signif(max(gene_groups$Total),0)),
        legend.text=TRUE, args.legend=list(x="bottomright", inset=c(-0.25, +0.05), cex=1.2))
x_axis = seq(0, signif(max(gene_groups$Total),0), length=5)
axis(side=1, at=x_axis, lab=formatC(x_axis, format="d", big.mark=","))
mtext(text="Gene counts", side=1, line=3, at=median(x_axis))

### Venn diagram of shared gene families
gene_counts = read.delim(fname_gene_counts, header=TRUE)
X = gene_counts[, 1:(ncol(gene_counts)-4)]
X$Orthogroup = as.numeric(gsub("OG", "", X$Orthogroup))+1
X[,2:ncol(X)] = X[,2:ncol(X)] > 0
colnames(X) = gsub("_", "\n", colnames(X))
par(mar=c(2, 3, 3, 2))
species_to_include= c("Lolium\nrigidum", "Lolium\nperenne", "Oryza\nsativa", "Zea\nmays", "Secale\ncereale")
idx = colnames(X) %in% species_to_include
venn(X[, idx]) ### picking only 5 species (maximum number of sets to draw a Venn diagram so far)

### Distribution of 4DTv (fraction of transverions among 4-fold degenerate codons - correlated with time from whole genome duplication using dual-copy paralogs and single-copy orthologs)

###@@@ Extract within genome 4DTv (pairwise paralogs)
FDTv_files = file.path(dir_name_4DTv, list.files(path=dir_name_4DTv, pattern=extension_name_4DTv))
FDTv_files = FDTv_files[grep(fname_4DTv_singlecopy, FDTv_files, invert=TRUE)]
id = c()
x = c()
y = c()
for (f in FDTv_files){
    # f = FDTv_files[1]
    print(f)
    df = read.delim(f, header=FALSE, na.string="NA")
    df[is.na(df)] = 0
    d = density(df$V4)
    id = c(id, gsub("_", " ", rep(gsub(".4DTv", "", basename(f)), length(d$x))))
    x = c(x, d$x)
    # y = c(y, (d$y - min(d$y))/diff(range(d$y)))
    y = c(y, d$y)
}
###@@@ Extract across genomes 4DTv (pairwise single-copy orthologs) and append to x and y
dat = read.delim(fname_4DTv_singlecopy, header=TRUE)
dat$SPECIES_1 = as.factor(dat$SPECIES_1)
dat$SPECIES_2 = as.factor(dat$SPECIES_2)
for (species1 in levels(dat$SPECIES_1)){
    # species1 = levels(dat$SPECIES_1)[1]
    for (species2 in levels(dat$SPECIES_2)){
        # species2 = levels(dat$SPECIES_2)[2]
        if (species1 == species2){
            next
        }
        df = droplevels(dat[(dat$SPECIES_1==species1) & (dat$SPECIES_2==species2), ])
        df[is.na(df)] = 0
        d = density(df$X4DTv)
        id = c(id, gsub("_", " ", rep(paste0(species1, " X ", species2), length(d$x))))
        x = c(x, d$x)
        # y = c(y, (d$y - min(d$y))/diff(range(d$y)))
        y = c(y, d$y)
    }
}
df = data.frame(id=as.factor(id), x=x, y=y)

species_list = c("Lolium rigidum", "Lolium perenne", "Oryza sativa",
                 "Lolium rigidum X Lolium perenne",
                 "Lolium rigidum X Oryza sativa",
                 "Lolium rigidum X Secale cereale",
                 "Lolium rigidum X Zea mays")
df = droplevels(df[df$id %in% species_list, ])

par(mar=c(5, 5, 5, 2))
n = nlevels(df$id)
colours = c("#66c2a5", "#fc8d62", "#8da0cb", "#e78ac3", "#a6d854", "#ffd92f", "#e5c494", "#b3b3b3")
plot(0, 0, xlim=range(df$x), ylim=range(df$y), xlab="4DTv", ylab="Density", type="n")
for (i in 1:n){
    # i = 1
    id = levels(df$id)[i]
    subdf = df[df$id==id, ]
    lines(x=subdf$x, y=subdf$y, col=colours[i], lwd=2)
}
grid()
legend("topright", legend=levels(df$id), col=colours, lwd=2)


dev.off()