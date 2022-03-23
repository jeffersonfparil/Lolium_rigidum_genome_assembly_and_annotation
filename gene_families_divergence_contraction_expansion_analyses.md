# Across multiple species, cluster gene families, align, estimate divergence times, and identify expanded and contracted gene families

## Set working directory
```{sh}
DIR=/data-weedomics-3
```

## Download genomes
Please refer to *Reference_genomes.md*

## Install GeMoMa (Gene Model Mapper)
```{sh}
wget http://www.jstacs.de/download.php?which=GeMoMa
mv 'download.php?which=GeMoMa' GeMoMa.zip
unzip GeMoMa.zip
java -jar GeMoMa-1.8.jar CLI -h
PATH=${PATH}:${DIR} ### for GeMoMa.jar
```

## Install mmseq for GeMoMa
```{sh}
wget https://github.com/soedinglab/MMseqs2/releases/download/13-45111/mmseqs-linux-avx2.tar.gz
tar -xvzf mmseqs-linux-avx2.tar.gz
mmseqs/bin/mmseqs -h
PATH=${PATH}:${DIR}/mmseqs/bin
```

## Download RefSeq version of the *Arabidopsis thaliana* and *Oryza sativa* reference genomes and gene annotations

We're using the Refseq data (i.e. non-GeneBank, GCF instead of GCA prefix) because we want to use the gff annotations and we don't want the hassle of converting the genebank gbff to gff

```{sh}
### Arabidopsis thaliana (TAIR10)
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/735/GCF_000001735.4_TAIR10.1/GCF_000001735.4_TAIR10.1_genomic.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/735/GCF_000001735.4_TAIR10.1/GCF_000001735.4_TAIR10.1_genomic.gff.gz
mv GCF_000001735.4_TAIR10.1_genomic.fna.gz Arabidopsis_thaliana.fasta.gz
mv GCF_000001735.4_TAIR10.1_genomic.gff.gz Arabidopsis_thaliana.gff.gz
### Oryza sativa (IRGSP-1.0)
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/001/433/935/GCF_001433935.1_IRGSP-1.0/GCF_001433935.1_IRGSP-1.0_genomic.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/001/433/935/GCF_001433935.1_IRGSP-1.0/GCF_001433935.1_IRGSP-1.0_genomic.gff.gz
mv GCF_001433935.1_IRGSP-1.0_genomic.fna.gz Oryza_sativa.fasta.gz
mv GCF_001433935.1_IRGSP-1.0_genomic.gff.gz Oryza_sativa.gff.gz
```

Identify genes across the genomes we want compare with *Gene Model Mapper*. For more information visit: [http://www.jstacs.de/index.php/GeMoMa#In_a_nutshell](http://www.jstacs.de/index.php/GeMoMa#In_a_nutshell)

Using the *Arabidopsis thaliana* gene annotations:

```{sh}
time \
for REF in Lolium_rigidum Lolium_perenne Oryza_sativa Zea_mays Secale_cereale Marchantia_polymorpha
do
echo ${REF}
java -jar -Xmx30G GeMoMa-1.8.jar CLI \
    GeMoMaPipeline \
    threads=15 \
    GeMoMa.Score=ReAlign \
    AnnotationFinalizer.r=NO \
    p=true pc=true pgr=true \
    o=true \
    t=${DIR}/${REF}/${REF}.fasta \
    i=Arabidopsis_thaliana \
    a=${DIR}/Arabidopsis_thaliana.gff.gz \
    g=${DIR}/Arabidopsis_thaliana.fasta.gz \
    outdir=${DIR}/GeMoMa_output_Arabidopsis_thaliana/${REF}
done

### WITH RNAseq data
# time \
# for REF in 
# java -jar -Xmx200G GeMoMa-1.8.jar CLI \
#     GeMoMaPipeline \
#     threads=31 \
#     GeMoMa.Score=ReAlign \
#     AnnotationFinalizer.r=NO \
#     p=true pc=true pgr=true \
#     o=true \
#     t=${GENOME} \
#     i=Arabidopsis_thaliana \
#     r=MAPPED \
#     ERE.m=${RNASEQ_BAM} \
#     a=${DIR}/Arabidopsis_thaliana.gff.gz \
#     g=${DIR}/Arabidopsis_thaliana.fasta.gz \
#     outdir=${DIR}/GEMOMA_ARABIDOPSIS_THALIANA_OUTPUT
```

Using the *Oryza sativa* gene annotations:

```{sh}
time \
for REF in Lolium_rigidum Lolium_perenne Arabidopsis_thaliana Zea_mays Secale_cereale Marchantia_polymorpha
do
echo ${REF}
java -jar -Xmx30G GeMoMa-1.8.jar CLI \
    GeMoMaPipeline \
    threads=15 \
    GeMoMa.Score=ReAlign \
    AnnotationFinalizer.r=NO \
    p=true pc=true pgr=true \
    o=true \
    t=${DIR}/${REF}/${REF}.fasta \
    i=Oryza_sativa \
    a=${DIR}/Oryza_sativa.gff.gz \
    g=${DIR}/Oryza_sativa.fasta.gz \
    outdir=${DIR}/GeMoMa_output_Oryza_sativa/${REF}
done
```

## Cluster gene families with OrthoMCL or Panther HMM
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

time \
parallel \
makeblastdb -in {}/{}.fasta \
            -title {} \
            -dbtype nucl ::: \
Lolium_rigidum \
Lolium_perenne \
Secale_cereale \
Zea_mays \
Oryza_sativa \
Arabidopsis_thaliana \
Marchantia_polymorpha

```{sh}
echo '#!/bin/bash
QUERY=$1
SPECIES=$2
echo "${SPECIES}----${QUERY}"
temp_name_1=$(basename $QUERY)
temp_name_2=${temp_name_1%.fasta*}
tblastn -db ${SPECIES}/${SPECIES}.fasta \
    -query ${QUERY} \
    -outfmt "6 qseqid staxids sstart send pident evalue qcovhsp bitscore stitle" \
    -out ${SPECIES}/BLASTOUT-${SPECIES}-${temp_name_2}.txt
' > tblastn_for_parallel_execution.sh
chmod +x tblastn_for_parallel_execution.sh

time \
parallel ./tblastn_for_parallel_execution.sh \
    {} \
    {} ::: $(ls ${DIR}/Lolium_rigidum_genome_assembly_and_annotation/misc/TSR_NTSR_etc_protein_sequences/*.fasta) \
       ::: Lolium_rigidum \
           Lolium_perenne \
           Secale_cereale \
           Zea_mays \
           Oryza_sativa \
           Arabidopsis_thaliana \
           Marchantia_polymorpha
```

```{julia}
using CSV
using DataFrames
using ProgressMeter

### Load the blast output
cd("/data-weedomics-3/GENE_FAMILIES")
vec_fnames_blastout = readdir()[match.(Regex("txt\$"), readdir()) .!= nothing]
@time for filename_blastout in vec_fnames_blastout
    # filename_blastout = "BLASTOUT-Glyphosate-target_EPSPS_UniProt_Mesangiospermae.txt"
    FILE = open(filename_blastout, "r")
    dat = CSV.read(FILE, DataFrames.DataFrame, delim='\t', header=["qseqid", "staxids", "sstart", "send", "pident", "evalue", "qcovhsp", "bitscore", "stitle"])
    close(FILE)

    ### Filter and sort
    subdat = dat[(dat.qcovhsp .>= 95) .& (dat.pident .>= 50), :]
    sort!(subdat, rev=false, [:sstart, :send])

    ### Find overlaps
    vec_overlaps = []
    @showprogress for i in 1:size(subdat, 1)
        # i=1
        x = subdat[i,:]
        for j in (i+1):size(subdat, 1)
            # j=2
            y = subdat[j,:]
            if (x.stitle == y.stitle)
                if (((x.send<=y.send) & (x.sstart>=y.sstart)) | 
                    ((x.send>=y.send) & (x.sstart>=y.sstart)) | 
                    ((x.send<=y.send) & (x.sstart<=y.sstart)) |
                    ((y.send<=x.send) & (y.sstart>=x.sstart)) | 
                    ((y.send>=x.send) & (y.sstart>=x.sstart)) | 
                    ((y.send<=x.send) & (y.sstart<=x.sstart)))
                    append!(vec_overlaps, j)
                end
            end
        end
    end
    unique!(vec_overlaps)

    ### Remove overlaps
    idx = map(x -> sum(x ∈ vec_overlaps)==0, collect(1:size(subdat,1)))
    subdat = subdat[idx, :]

    ### Save unique non-overlapping hits as csv
    filename_output = string(join(split(filename_blastout, '.')[1:(end-1)], '.'), "-UNIQUE_HITS.csv")
    OUT = open(filename_output, "w")
    CSV.write(OUT, subdat)
    close(OUT)
end

cd("/data-weedomics-3/GENE_FAMILIES")

using CSV, DataFrames
include("/home/jeffersonfparil/Documents/Lolium_rigidum_genome_assembly_and_annotation/genome_statistics.jl")

vec_fnames_blastout_uniques = readdir()[match.(Regex("-UNIQUE_HITS.csv\$"), readdir()) .!= nothing]
for f in vec_fnames_blastout_uniques
    # f = vec_fnames_blastout_uniques[1]
    gene_name = split(split(f, '-')[3], '_')[2]
    FILE = open(f, "r")
    dat = CSV.read(FILE, DataFrames.DataFrame, header=true)
    close(FILE)
end




```


