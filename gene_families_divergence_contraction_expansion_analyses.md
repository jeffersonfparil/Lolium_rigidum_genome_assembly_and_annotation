# Across multiple species, cluster gene families, align, estimate divergence times, and identify expanded and contracted gene families

## Set working directory
```{sh}
# DIR=/data-weedomics-3
DIR=/data/Lolium_rigidum_ASSEMBLY/COMPARATIVE_GENOMICS
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

## Install HMMER
```{sh}
sudo apt install hmmer
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

## Download PatherHMM library including 15,619 protein family HMMs and their GO terms
```{sh}
wget http://data.pantherdb.org/ftp/panther_library/current_release/PANTHER17.0_hmmscoring.tgz
tar -xvzf PANTHER17.0_hmmscoring.tgz
mv target/ PatherHMM_17.0/
cd PatherHMM_17.0/
wget http://data.pantherdb.org/ftp/hmm_classifications/current_release/PANTHER17.0_HMM_classifications
### Isolate family codes and names, i.e. exclude subfamily info
grep -v ':SF' PANTHER17.0_HMM_classifications > Panther17.0_HMM_familyIDs.txt
cd -
```

## Gene model mapping with GeMoMa

Identify genes across the genomes we want compare with *Gene Model Mapper*. For more information visit: [http://www.jstacs.de/index.php/GeMoMa#In_a_nutshell](http://www.jstacs.de/index.php/GeMoMa#In_a_nutshell)

Using the *Arabidopsis thaliana* gene annotations:

```{sh}
time \
for REF in Lolium_rigidum Lolium_perenne Arabidopsis_thaliana Oryza_sativa Zea_mays Secale_cereale Marchantia_polymorpha
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
    outdir=${DIR}/${REF}/GeMoMa_output_Arabidopsis_thaliana
done
```

Using the *Oryza sativa* gene annotations:

```{sh}
time \
for REF in Lolium_rigidum Lolium_perenne Arabidopsis_thaliana Oryza_sativa Zea_mays Secale_cereale Marchantia_polymorpha
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
    outdir=${DIR}/${REF}/GeMoMa_output_Oryza_sativa
done
```

## Classify GeMoMa-predicted proteins by gene families using PatherHMM database and HMMER
```{sh}
cd $DIR
# Define the location of the 15,619 protein family HMMs
DIR_PANTHER=${DIR}/PatherHMM_17.0/famlib/rel/PANTHER17.0_altVersion/hmmscoring/PANTHER17.0/books
GOT_PATHER=${DIR}/PatherHMM_17.0/PANTHER17.0_HMM_classifications
# Prepare parallelisable HMMER search script
echo '#!/bin/bash
PROTFA=$1
DIR_PANTHER=$2
d=$3
HMMODL=${DIR_PANTHER}/${d}/hmmer.hmm
OUTEMP=$(dirname ${PROTFA})/hhmer_gene_family_hits-${d}.tmp
hmmsearch -E 0.0001 --tblout ${OUTEMP} ${HMMODL} ${PROTFA}
' > hmmsearch_for_parallel_execution.sh
chmod +x hmmsearch_for_parallel_execution.sh
# Iteratively, for each GeMoMa-predicted protein sequences run hmmsearch in paralel for each PatherHMM protein family
for PROTFA in $(ls */GeMoMa_output_*/predicted_proteins.fasta)
do
    # PROTFA=Lolium_rigidum/GeMoMa_output_Arabidopsis_thaliana/predicted_proteins.fasta
    time \
    parallel ./hmmsearch_for_parallel_execution.sh \
        ${PROTFA} \
        ${DIR_PANTHER} \
        {} ::: $(ls $DIR_PANTHER)
    ### Concatenate hmmsearch output for each protein family into a single output file
    OUTPUT=$(dirname ${PROTFA})/hhmer_gene_family_hits.txt
    touch $OUTPUT
    for f in $(ls $(dirname ${PROTFA})/hhmer_gene_family_hits-*)
    do
        sed "/^#/d" ${f} | awk '{print $1,$3,$5}' >> ${OUTPUT}
    done
    sort ${OUTPUT} > ${OUTPUT}.tmp
    mv ${OUTPUT}.tmp ${OUTPUT}
done
```
## Generate a new annotation file with PantherHMM-derived gene families
```{julia}
using CSV, DataFrames, ProgressMeter
fname_input = ARGS[1]
fname_codes = ARGS[2]
fname_annot = ARGS[3]
# fname_input = "/data/Lolium_rigidum_ASSEMBLY/COMPARATIVE_GENOMICS/Lolium_rigidum/GeMoMa_output_Arabidopsis_thaliana/hhmer_gene_family_hits.txt"; fname_codes = "/data/Lolium_rigidum_ASSEMBLY/COMPARATIVE_GENOMICS/PatherHMM_17.0/Panther17.0_HMM_familyIDs.txt"; fname_annot = "/data/Lolium_rigidum_ASSEMBLY/COMPARATIVE_GENOMICS/Lolium_rigidum/GeMoMa_output_Arabidopsis_thaliana/final_annotation.gff"
# fname_input = "hhmer_gene_family_hits.txt"; fname_codes = "Panther17.0_HMM_familyIDs.txt"; fname_annot = "final_annotation.gff"
dat = CSV.read(open(fname_input), DataFrames.DataFrame, header=false)
### Find most likely gene family per predicted gene, i.e. remove duplicate family hits
gene_names = []
panther_codes = []
@showprogress for gene in levels(dat.Column1)
    # gene = levels(dat.Column1)[2]
    subset = dat[dat.Column1 .== gene, :]
    sort!(subset, :Column3, rev=false)
    push!(gene_names, subset[1,1])
    push!(panther_codes, split(subset[1,2], '.')[1])
end
### Identify PatherHMM gene families
gene_families = []
cod = CSV.read(open(fname_codes), DataFrames.DataFrame, header=false)
@showprogress for code in panther_codes
    # code = panther_codes[1]
    push!(gene_families, cod.Column2[code .== cod.Column1][1])
end
### Create a new annotation file with PatherHMM codes and gene family names
fname_output = "FINAL_ANNOTATION_PATHERHMM_GENE_FAMILIES.gff"
if dirname(fname_input) != ""
    fname_output = string(dirname(fname_input), "/", fname_output)
end
gff_input = open(fname_annot, "r")
gff_output = open(fname_output, "a")
seekend(gff_input); end_position = position(gff_input)
seekstart(gff_input)
p = Progress(end_position, 1)
while !eof(gff_input)
    line = readline(gff_input)
    if line[1] == '#'
        write(gff_output, string(line, '\n'))
        continue
    end
    vec_split = split(line, '\t')
    gene = replace(split(vec_split[end], ';')[1], "ID="=>"")
    idx = gene .== gene_names
    if sum(idx) == 1
        output_line = vec_split[1:(end-1)]
        push!(output_line, string("ID=", gene, ";",
                                        panther_codes[idx][1], ";", 
                                        gene_families[idx][1]))
        write(gff_output, string(join(output_line, '\t'), '\n'))
    end
    update!(p, position(gff_input))
end
close(gff_input)
close(gff_output)
## Misc
# gene_families[match.(Regex("ENOLPY"), gene_families) .!= nothing] ### EPSPS
# gene_families[match.(Regex("PSB"), gene_families) .!= nothing] ### psbA
# gene_families[match.(Regex("ACETO"), gene_families) .!= nothing] ### ALS
# gene_families[match.(Regex("ACETYL"), gene_families) .!= nothing] ### ACCase
# gene_families[match.(Regex("CYT"), gene_families) .!= nothing] ### NTSR cytochorme detox gene families
# gene_families[match.(Regex("UBIQUITIN"), gene_families) .!= nothing] ### use in assessing divergence between species?
```

Save the above Julia script as `generate_final_gff_annotation_file.jl` and generate final annotation files:
```{sh}
cd $DIR
PANTHER_CODES=${DIR}/PatherHMM_17.0/Panther17.0_HMM_familyIDs.txt
time \
parallel --link \
julia generate_final_gff_annotation_file.jl {1} ${PANTHER_CODES} {2} ::: \
    $(ls */GeMoMa_output_*/hhmer_gene_family_hits.txt) ::: \
    $(ls */GeMoMa_output_*/final_annotation.gff)
```


<!-- 
## Cluster gene families with OrthoMCL or Panther HMM
```{sh}
wget https://orthomcl.org/common/downloads/software/v2.0/orthomclSoftware-v2.0.9.tar.gz
tar -xvzf orthomclSoftware-v2.0.9.tar.gz
cd orthomclSoftware-v2.0.9/
``` -->

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
    gene_names = split(split(f, '-')[3], '_')[2]
    FILE = open(f, "r")
    dat = CSV.read(FILE, DataFrames.DataFrame, header=true)
    close(FILE)
end




```


