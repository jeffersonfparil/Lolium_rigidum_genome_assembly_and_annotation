args = commandsArgs(trailingOnly=TRUE)
# args = c("SOD-OG0016003.kaks.tmp", "0.001")
fname_input = args[1]
alpha = as.numeric(args[2])

dat = read.delim(fname_input, header=TRUE)
X = matrix(unlist(strsplit(as.character(dat$Sequence), "-")), ncol=12, byrow=TRUE)
Y = matrix(unlist(strsplit(X[,11], "\\(")), ncol=2, byrow=TRUE)
df = data.frame(
    Gene=X[,2],
    Orthogroup=X[,3],
    Species1=X[,1],
    Species_alignment1=X[,4],
    Species2=X[,8],
    Species_alignment2=Y[,1],
    Position_ini=as.numeric(Y[,2]),
    Position_fin=as.numeric(gsub("\\)", "", X[,12])),
    KaKs=dat$Ka.Ks,
    p=dat$P.Value.Fisher.
)
df$KaKs[is.na(df$KaKs)] = 0.0

alignments = unique(df$Species_alignment2)
n = ceiling(sqrt(length(alignments)))
m = ceiling(length(alignments)/n)

svg(paste0(fname_input, ".svg"), height=5*n, width=8*m)
par(mfrow=c(n,m))
for (aln in alignments) {
    # aln = unique(df$Species_alignment2)[1]
    subdf = df[df$Species_alignment2 == aln, ]
    plot(subdf$Position_ini, subdf$KaKs, type="l",
        main=paste0(subdf$Species1[1], "::", subdf$Species_alignment1[1], " vs ", subdf$Species2[1], "::", subdf$Species_alignment2[1]),
        sub=paste0(subdf$Gene[1], " (", subdf$Orthogroup[1], ")"),
        xlab="Position (bp)", ylab="Ka/Ks")
    grid()
    idx = c(1:nrow(subdf))[!is.na(subdf$KaKs) & (subdf$p <= alpha) & (subdf$KaKs > 1)]
    for (j in idx){
        text(x=subdf$Position_ini[j], y=subdf$KaKs[j], lab="*", col="red", cex=2)
    }
}
dev.off()