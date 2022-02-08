# Across multiple species, cluster gene families, align, estimate divergence times, and identify expanded and contracted gene families

**NOTE:** We need to have the gene annotation on-hand to perform these analyses.

## Cluster gene families with OrthoMCL
```{sh}
wget https://orthomcl.org/common/downloads/software/v2.0/orthomclSoftware-v2.0.9.tar.gz
tar -xvzf orthomclSoftware-v2.0.9.tar.gz
cd orthomclSoftware-v2.0.9/
```

## Align gene families across species with MAFFT
```{sh}
sudo apt install -y mafft
```

## Estimate divergence between species times using MCMCTREE and TimeTree.org fossil record estimates


## Identify expanded and contracted gene families for each species


**NOTE:** Naive cluster analysis using blast hits to build a tree with mafft

DIR=/data-weedomics-3
BLASTOUT=${DIR}/BLASTOUT-Overwatch-target_DXS_UniProt_Mesangiospermae.txt

```{R}
filename_blastout = "/data-weedomics-3/BLASTOUT-Overwatch-target_DXS_UniProt_Mesangiospermae.txt"
filename_blastout = "/data-weedomics-3/BLASTOUT-Glyphosate-target_EPSPS_UniProt_Mesangiospermae.txt"
dat = read.table(filename_blastout, header=FALSE)
colnames(dat) = c("qseqid", "staxids", "sstart", "send", "pident", "evalue", "qcovhsp", "bitscore", "stitle")

### Filter by query coverage and similarity
dat[(dat$qcovhsp >= 90) & (dat$pident>= 50), ]

```

```{julia}
filename_blastout = "/data-weedomics-3/BLASTOUT-Overwatch-target_DXS_UniProt_Mesangiospermae.txt"
filename_blastout = "/data-weedomics-3/BLASTOUT-Glyphosate-target_EPSPS_UniProt_Mesangiospermae.txt"
using CSV
using DataFrames
FILE = open(filename_blastout, "r")
dat = CSV.read(FILE, DataFrames.DataFrame, delim='\t', header=["qseqid", "staxids", "sstart", "send", "pident", "evalue", "qcovhsp", "bitscore", "stitle"])
close(FILE)

dat[(dat.qcovhsp .>= 90) .& (dat.pident .>= 50), :]

```


