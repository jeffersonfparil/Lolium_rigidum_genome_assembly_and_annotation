# Comparative genomics

## Set working directories and the main orthogroup output filename (these variables and reiterated in their respective sections where they were used)
```sh}
DIR=/data/Lolium_rigidum_ASSEMBLY/COMPARATIVE_GENOMICS
DIR_ORTHOGROUPS=${DIR}/ORTHOGROUPS/OrthoFinder/Results_*
DIR_PANTHER=${DIR}/PantherHMM_17.0/famlib/rel/PANTHER17.0_altVersion/hmmscoring/PANTHER17.0/books
GOT_PATHER=${DIR}/PantherHMM_17.0/PANTHER17.0_HMM_classifications
MERGED_ORTHOGROUPS=${DIR}/ORTHOGROUPS/orthogroups.faa
ORTHOUT=${DIR}/ORTHOGROUPS/orthogroups_gene_counts_families_go.out
cd $DIR
```

## Download genomes, genome annotations, and predicted CDS sequences
```{sh}
### Arabidopsis thaliana
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/735/GCF_000001735.4_TAIR10.1/GCF_000001735.4_TAIR10.1_genomic.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/735/GCF_000001735.4_TAIR10.1/GCF_000001735.4_TAIR10.1_genomic.gff.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/735/GCF_000001735.4_TAIR10.1/GCF_000001735.4_TAIR10.1_cds_from_genomic.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/735/GCF_000001735.4_TAIR10.1/GCF_000001735.4_TAIR10.1_protein.faa.gz
gunzip -c GCF_000001735.4_TAIR10.1_genomic.fna.gz > Arabidopsis_thaliana.fasta
gunzip -c GCF_000001735.4_TAIR10.1_genomic.gff.gz > Arabidopsis_thaliana.gff
gunzip -c GCF_000001735.4_TAIR10.1_cds_from_genomic.fna.gz > Arabidopsis_thaliana.cds
gunzip -c GCF_000001735.4_TAIR10.1_protein.faa.gz > Arabidopsis_thaliana.faa
### Oryza sativa
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/001/433/935/GCF_001433935.1_IRGSP-1.0/GCF_001433935.1_IRGSP-1.0_genomic.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/001/433/935/GCF_001433935.1_IRGSP-1.0/GCF_001433935.1_IRGSP-1.0_genomic.gff.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/001/433/935/GCF_001433935.1_IRGSP-1.0/GCF_001433935.1_IRGSP-1.0_cds_from_genomic.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/001/433/935/GCF_001433935.1_IRGSP-1.0/GCF_001433935.1_IRGSP-1.0_protein.faa.gz
gunzip -c GCF_001433935.1_IRGSP-1.0_genomic.fna.gz > Oryza_sativa.fasta
gunzip -c GCF_001433935.1_IRGSP-1.0_genomic.gff.gz > Oryza_sativa.gff
gunzip -c GCF_001433935.1_IRGSP-1.0_cds_from_genomic.fna.gz > Oryza_sativa.cds
gunzip -c GCF_001433935.1_IRGSP-1.0_protein.faa.gz > Oryza_sativa.faa
### Zea mays
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/902/167/145/GCF_902167145.1_Zm-B73-REFERENCE-NAM-5.0/GCF_902167145.1_Zm-B73-REFERENCE-NAM-5.0_genomic.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/902/167/145/GCF_902167145.1_Zm-B73-REFERENCE-NAM-5.0/GCF_902167145.1_Zm-B73-REFERENCE-NAM-5.0_genomic.gff.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/902/167/145/GCF_902167145.1_Zm-B73-REFERENCE-NAM-5.0/GCF_902167145.1_Zm-B73-REFERENCE-NAM-5.0_cds_from_genomic.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/902/167/145/GCF_902167145.1_Zm-B73-REFERENCE-NAM-5.0/GCF_902167145.1_Zm-B73-REFERENCE-NAM-5.0_protein.faa.gz
gunzip -c GCF_902167145.1_Zm-B73-REFERENCE-NAM-5.0_genomic.fna.gz > Zea_mays.fasta
gunzip -c GCF_902167145.1_Zm-B73-REFERENCE-NAM-5.0_genomic.gff.gz > Zea_mays.gff
gunzip -c GCF_902167145.1_Zm-B73-REFERENCE-NAM-5.0_cds_from_genomic.fna.gz > Zea_mays.cds
gunzip -c GCF_902167145.1_Zm-B73-REFERENCE-NAM-5.0_protein.faa.gz > Zea_mays.faa
### Sorghum bicolor
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/003/195/GCF_000003195.3_Sorghum_bicolor_NCBIv3/GCF_000003195.3_Sorghum_bicolor_NCBIv3_genomic.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/003/195/GCF_000003195.3_Sorghum_bicolor_NCBIv3/GCF_000003195.3_Sorghum_bicolor_NCBIv3_genomic.gff.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/003/195/GCF_000003195.3_Sorghum_bicolor_NCBIv3/GCF_000003195.3_Sorghum_bicolor_NCBIv3_cds_from_genomic.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/003/195/GCF_000003195.3_Sorghum_bicolor_NCBIv3/GCF_000003195.3_Sorghum_bicolor_NCBIv3_protein.faa.gz
gunzip -c GCF_000003195.3_Sorghum_bicolor_NCBIv3_genomic.fna.gz > Sorghum_bicolor.fasta
gunzip -c GCF_000003195.3_Sorghum_bicolor_NCBIv3_genomic.gff.gz > Sorghum_bicolor.gff
gunzip -c GCF_000003195.3_Sorghum_bicolor_NCBIv3_cds_from_genomic.fna.gz > Sorghum_bicolor.cds
gunzip -c GCF_000003195.3_Sorghum_bicolor_NCBIv3_protein.faa.gz > Sorghum_bicolor.faa
### Lolium rigidum
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/022/539/505/GCF_022539505.1_APGP_CSIRO_Lrig_0.1/GCF_022539505.1_APGP_CSIRO_Lrig_0.1_genomic.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/022/539/505/GCF_022539505.1_APGP_CSIRO_Lrig_0.1/GCF_022539505.1_APGP_CSIRO_Lrig_0.1_genomic.gff.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/022/539/505/GCF_022539505.1_APGP_CSIRO_Lrig_0.1/GCF_022539505.1_APGP_CSIRO_Lrig_0.1_cds_from_genomic.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/022/539/505/GCF_022539505.1_APGP_CSIRO_Lrig_0.1/GCF_022539505.1_APGP_CSIRO_Lrig_0.1_protein.faa.gz
gunzip -c GCF_022539505.1_APGP_CSIRO_Lrig_0.1_genomic.fna.gz > Lolium_rigidum.fasta
gunzip -c GCF_022539505.1_APGP_CSIRO_Lrig_0.1_genomic.gff.gz > Lolium_rigidum.gff
gunzip -c GCF_022539505.1_APGP_CSIRO_Lrig_0.1_cds_from_genomic.fna.gz > Lolium_rigidum.cds
gunzip -c GCF_022539505.1_APGP_CSIRO_Lrig_0.1_protein.faa.gz > Lolium_rigidum.faa
### Lolium perenne
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/019/359/855/GCA_019359855.1_MPB_Lper_Kyuss_1697/GCA_019359855.1_MPB_Lper_Kyuss_1697_genomic.fna.gz
gunzip -c GCA_019359855.1_MPB_Lper_Kyuss_1697_genomic.fna.gz > Lolium_perenne.fasta
### Secale cereale
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/016/097/815/GCA_016097815.1_HAU_Weining_v1.0/GCA_016097815.1_HAU_Weining_v1.0_genomic.fna.gz
gunzip -c GCA_016097815.1_HAU_Weining_v1.0_genomic.fna.gz > Secale_cereale.fasta
### 
```

**Note:** The first column of the annotation files may not correspond to the chromosome IDs.

## Install GeMoMa (Gene Model Mapper)
```{sh}
wget http://www.jstacs.de/download.php?which=GeMoMa
mv 'download.php?which=GeMoMa' GeMoMa.zip
unzip GeMoMa.zip
java -jar GeMoMa-1.8.jar CLI -h
PATH=${PATH}:${DIR} ### for GeMoMa.jar
```

## Install mmseq for GeMoMa-based gene annotation mapping for unannotated *Lolium perenne* and *Secale cereale* genomes
```{sh}
wget https://github.com/soedinglab/MMseqs2/releases/download/13-45111/mmseqs-linux-avx2.tar.gz
tar -xvzf mmseqs-linux-avx2.tar.gz
cd mmseqs/bin/
./mmseqs -h
PATH=${PATH}:${DIR}/mmseqs/bin
cd -
```

## Install OrthoFinder for classifying genes into orthologs, and paralogs, as well as to build a tree for the analysis of gene family evolution
```{sh}
wget https://github.com/davidemms/OrthoFinder/releases/download/2.5.4/OrthoFinder.tar.gz
tar -xvzf OrthoFinder.tar.gz
cd OrthoFinder/
PATH=${PATH}:$(pwd)
cd -
```

## Install HMMER for mapping CDS to PantherHMM gene family models
```{sh}
sudo apt install hmmer
```

## Install CAFE to analyse gene family evolution
```{sh}
wget https://github.com/hahnlab/CAFE5/releases/download/v5.0/CAFE5-5.0.0.tar.gz
tar -xvzf CAFE5-5.0.0.tar.gz
cd CAFE5/
./configure
make
bin/cafe5 -h
PATH=${PATH}:$(pwd)/bin
cd -
```

## Install MACSE: Multiple Alignment of Coding SEquences Accounting for Frameshifts and Stop Codons
```{sh}
wget https://bioweb.supagro.inra.fr/macse/releases/macse_v2.06.jar
java -Xmx250G -jar macse_v2.06.jar -h
```

## Install KaKs_Calculator2.0 to assess signatures of selection
```{sh}
wget https://github.com/kullrich/kakscalculator2/archive/refs/tags/v2.0.1.tar.gz
tar -xvzf v2.0.1.tar.gz
cd kakscalculator2-2.0.1/src
make
PATH=${PATH}:$(pwd)
cd -
```

## Download PantherHMM library including 15,619 protein family HMMs and their GO terms
```{sh}
wget http://data.pantherdb.org/ftp/panther_library/current_release/PANTHER17.0_hmmscoring.tgz
tar -xvzf PANTHER17.0_hmmscoring.tgz
mv target/ PantherHMM_17.0/
cd PantherHMM_17.0/
wget http://data.pantherdb.org/ftp/hmm_classifications/current_release/PANTHER17.0_HMM_classifications
### Isolate family codes and names, i.e. exclude subfamily info
grep -v ':SF' PANTHER17.0_HMM_classifications > Panther17.0_HMM_familyIDs.txt
cd -
```

## Use GeMoMa to map our *Lolium rigidum* annotations into *Loliumm perenne* and *Secale cereale* genomes to extract CDS
```{sh}
time \
for REF in Lolium_perenne Secale_cereale
do
    echo ${REF}
    java -jar -Xmx280G GeMoMa-1.8.jar CLI \
        GeMoMaPipeline \
        threads=31 \
        GeMoMa.Score=ReAlign \
        AnnotationFinalizer.r=NO \
        p=true pc=true pgr=true \
        o=true \
        t=${DIR}/${REF}.fasta \
        i=Lolium_rigidum \
        a=${DIR}/Lolium_rigidum.gff \
        g=${DIR}/Lolium_rigidum.fasta \
        outdir=${DIR}/${REF}
    cp ${REF}/final_annotation.gff ${REF}.gff
    cp ${REF}/predicted_cds.fasta ${REF}.cds
    cp ${REF}/predicted_proteins.fasta ${REF}.faa
    rm -R ${REF}
done
```

## Use OrthoFinder to find orthologs and paralogs
```{sh}
mkdir ORTHOGROUPS/
cp *.faa ORTHOGROUPS/
for f in ORTHOGROUPS/*.faa
do
    # f=$(ls ORTHOGROUPS/*.faa | head -n1)
    fname=$(basename $f)
    species=${fname%.faa*}
    sed -i "s/^>/>$species|/g" $f
done
time \
orthofinder \
    -f ORTHOGROUPS/ \
    -t 32

DIR_ORTHOGROUPS=${DIR}/ORTHOGROUPS/OrthoFinder/Results_*/
```

## Install IQ-TREE for building trees with fossil root dates
```{sh}
sudo apt install libeigen3-dev
wget https://github.com/Cibiv/IQ-TREE/releases/download/v2.0.7/iqtree-2.0.7-Linux.tar.gz
tar -xvzf iqtree-2.0.7-Linux.tar.gz
PATH=${PATH}:$(pwd)/iqtree-2.0.7-Linux/bin
```

## Assign orthogroups into gene families
```{sh}
### Define the location of the 15,619 protein family HMMs
DIR_PANTHER=${DIR}/PantherHMM_17.0/famlib/rel/PANTHER17.0_altVersion/hmmscoring/PANTHER17.0/books
GOT_PATHER=${DIR}/PantherHMM_17.0/PANTHER17.0_HMM_classifications

### Iteratively, for each genome's predicted protein sequences run hmmsearch in paralel for each PantherHMM protein family

### Add Orthogroup name to each protein sequence name and merge so that we can be more efficient with hmmsearch
echo '#!/bin/bash
f=$1
ORTHOGROUP=$(basename ${f} | sed s/.fa$//g)
sed "s/^>/>${ORTHOGROUP}:/g" ${f} > ${ORTHOGROUP}.tmp
' > rename_sequences_with_orthogroup_ID.sh
chmod +x rename_sequences_with_orthogroup_ID.sh

### Split the large list of orthogroup filenames (too long for bash commands)
find ${DIR_ORTHOGROUPS}/Orthogroup_Sequences/ -name '*.fa' > \
    orthogroup_filenames.tmp
split -l 10000 orthogroup_filenames.tmp orthogroup_filenames-SPLIT-

### For each chunk of orthogroup filenames append the orthogroup names into the protein sequence names in parallel
time \
for F in $(ls orthogroup_filenames-SPLIT-*)
do
    echo $F
    parallel \
        ./rename_sequences_with_orthogroup_ID.sh \
        {} \
        ::: $(cat $F)
done

### Merge orthogroup sequences into a single fasta file
MERGED_ORTHOGROUPS=${DIR}/ORTHOGROUPS/orthogroups.faa
touch $MERGED_ORTHOGROUPS
time \
for f in $(ls | grep "^OG" | grep "tmp$")
do
    cat $f >> ${MERGED_ORTHOGROUPS}
    rm $f
done

### Clean-up
rm orthogroup_filenames*

### Prepare parallelisable HMMER search script
echo '#!/bin/bash
PROTFA=$1
DIR_PANTHER=$2
d=$3
HMMODL=${DIR_PANTHER}/${d}/hmmer.hmm
OUTEMP=${PROTFA}-hhmer_gene_family_hits-${d}.tmp
hmmsearch -E 0.0001 --tblout ${OUTEMP} ${HMMODL} ${PROTFA}
sed "/^#/d" ${OUTEMP} | awk @{print $1,$3,$5}@ > ${OUTEMP}.tmp
if [ $(cat ${OUTEMP}.tmp | wc -l) -eq 0 ]
then
    rm ${OUTEMP} ${OUTEMP}.tmp
else
    mv ${OUTEMP}.tmp ${OUTEMP}
fi
' | sed "s/@/'/g" > hmmsearch_for_parallel_execution.sh
chmod +x hmmsearch_for_parallel_execution.sh

### Find PantherHMM protein families for each orthogroup
time \
parallel \
./hmmsearch_for_parallel_execution.sh \
    ${MERGED_ORTHOGROUPS} \
    ${DIR_PANTHER} \
    {} \
    ::: $(ls $DIR_PANTHER)

### Concatenate hmmsearch output for each protein family into a single output file
PANTHER_ORTHOGROUPS=$(echo ${MERGED_ORTHOGROUPS} | sed s/.faa$//g).pthr
cat ${MERGED_ORTHOGROUPS}-hhmer_gene_family_hits-* > ${PANTHER_ORTHOGROUPS}
rm ${MERGED_ORTHOGROUPS}-hhmer_gene_family_hits-*

### Find the best fitting gene family to each unique sequence per orthogroup.
### This means that each orthogroup can have multiple gene families.
### Next, add family name and GO terms to each gene family.

grep "^>" ${MERGED_ORTHOGROUPS} | cut -d':' -f1 | sed 's/>//g' | sort | uniq > all_orthogroups.tmp

echo '
using CSV, DataFrames, ProgressMeter

fname_orthogroup_family_hits =  ARGS[1]
fname_family_GO =               ARGS[2]
fname_paralogs =             ARGS[3]
fname_orthogroup_gene_counts =  ARGS[4]
fname_unassigned_genes =        ARGS[5]
fname_output =                  ARGS[6]

# Load all orthogroup IDs
file = open(fname_paralogs, "r")
all_orthogroups = readlines(file)
close(file)
rm(fname_paralogs) # clean-up

# Load orthogroup hits
file = open(fname_orthogroup_family_hits, "r")
seekend(file); n = position(file); seekstart(file)
orthogroup = []
gene_name = []
family_ID = []
evalue = []
pb = Progress(n)
while !eof(file)
    line = split(readline(file), " "[1])
    seqName = split(line[1], ":"[1])
    speciesAndGeneName = split(seqName[2], "|"[1])
    geneName = join(split(speciesAndGeneName[2], "-"[1])[2:end], "-"[1])
    push!(orthogroup, seqName[1])
    push!(gene_name, geneName)
    push!(family_ID, replace(line[2], ".orig.30.pir"=>""))
    push!(evalue, parse(Float64, line[3]))
    update!(pb, position(file))
end
close(file)

# Load PantherHMM family descriptions
file = open(fname_family_GO, "r")
seekend(file); n = position(file); seekstart(file)
PTHR_family_ID = []
PTHR_family_name = []
PTHR_GO_term = []
pb = Progress(n)
while !eof(file)
    line = split(readline(file), "\t"[1])
    push!(PTHR_family_ID, line[1])
    push!(PTHR_family_name, line[2])
    push!(PTHR_GO_term, line[3])
    update!(pb, position(file))
end
close(file)

# Load gene counts orthogroup per species
df_counts = CSV.read(open(fname_orthogroup_gene_counts), DataFrames.DataFrame)

# Load list of unassigned genes, i.e. orthogroups with a single gene specific to each species
df_unassigned = CSV.read(open(fname_unassigned_genes), DataFrames.DataFrame)

# For each gene, set the family ID as the one with the lowest E-value
idx = sortperm(gene_name)
orthogroup = orthogroup[idx]
gene_name = gene_name[idx]
family_ID = family_ID[idx]
evalue = evalue[idx]

all_gene = unique(gene_name)
all_family_ID = []
all_orthogroup = []
i = 1
@showprogress for g in all_gene
    # g = all_gene[1]
    t = g == gene_name[i]
    while t == false
        i += 1
        t = g == gene_name[i]
    end
    f = []
    o = []
    e = []
    while t
        push!(f, family_ID[i])
        push!(o, orthogroup[i])
        push!(e, evalue[i])
        i += 1
        t = try
                g == gene_name[i]
            catch
                false
            end
    end
    push!(all_family_ID, f[e .== minimum(e)][1])
    push!(all_orthogroup, o[e .== minimum(e)][1])
end

# Identify the PantherHMM gene family names and GO terms
all_family = []
all_GO = []
@showprogress for ID in all_family_ID
    # ID = all_family_ID[1]
    idx = ID .== PTHR_family_ID
    if sum(idx) == 1
        push!(all_family, PTHR_family_name[idx][1])
        push!(all_GO, PTHR_GO_term[idx][1])
    else
        # for unclassified orthogroups
        push!(all_family, "UNKNOWN")
        push!(all_GO, "")
    end
end

# Summarise gene families per orthogroup
idx = sortperm(all_orthogroup)
all_orthogroup = all_orthogroup[idx]
all_family_ID = all_family_ID[idx]
all_family = all_family[idx]
all_GO = all_GO[idx]
final_orthogroup = unique(all_orthogroup)
sort!(final_orthogroup)
final_family_ID = []
final_family = []
final_GO = []
i = 1
@showprogress for o in final_orthogroup
    t = o == all_orthogroup[i]
    while t == false
        i += 1
        t = o == all_orthogroup[i]
    end
    fid = []
    fam = []
    fgo = []
    while t
        push!(fid, all_family_ID[i])
        push!(fam, all_family[i])
        push!(fgo, all_GO[i])
        i += 1
        t = try
                o == all_othogroup[i]
            catch
                false
            end
    end
    push!(final_family_ID, join(unique(fid), ";"[1]))
    push!(final_family, join(unique(fam), ";"[1]))
    push!(final_GO, join(unique(fgo), ";"[1]))
end

df_ID = DataFrames.DataFrame(Orthogroup=final_orthogroup,
                          Family_ID=final_family_ID,
                          Family=final_family,
                          GO=final_GO)

# Generate gene counts per species for the set of unassigned genes
df_append_unassigned = Int.(.!ismissing.(df_unassigned[:, 2:end]))
df_append_unassigned.Total = repeat([1], nrow(df_unassigned))
df_append_unassigned.Orthogroup = df_unassigned.Orthogroup

# Append unassigned orthogroups into the orthogroup gene counts dataframe
df_counts = vcat(df_counts, df_append_unassigned)

# Merge the orthogroup ID and orthogroup gene counts and save into a file
df = outerjoin(df_counts, df_ID, on=:Orthogroup)
CSV.write(open(fname_output, "w"), df, delim="\t"[1])
' > "orthogroup_classification_gene_family_GO_terms.jl"

time \
julia orthogroup_classification_gene_family_GO_terms.jl \
        ORTHOGROUPS/orthogroups.pthr \
        PantherHMM_17.0/Panther17.0_HMM_familyIDs.txt \
        all_orthogroups.tmp \
        ${DIR_ORTHOGROUPS}/Orthogroups/Orthogroups.GeneCount.tsv \
        ${DIR_ORTHOGROUPS}/Orthogroups/Orthogroups_UnassignedGenes.tsv \
        ORTHOGROUPS/orthogroups_gene_counts_families_go.out
```

## Preliminary assessment of the distribution of the genes, orthogroups and gene family classifications.
```{sh}
echo '
fname = ARGS[1]
fname_out = ARGS[2]
using CSV, DataFrames, ProgressMeter
file = open(fname, "r")
df = CSV.read(file, DataFrames.DataFrame, header=true)
close(file)
counts = df[:, 2:(end-4)]
species_names = names(counts)

# Unassigned gene
function count_unassigned_genes(counts, species_names)
    println("Count unassigned genes")
    idx = df.Total .== 1
    unassigned = []
    for species in species_names
        # species = species_names[1]
        push!(unassigned, sum(counts[idx, species_names .== species][:,1]))
    end
    return(unassigned)
end
@time unassigned = count_unassigned_genes(counts, species_names)

# Unique paralogs
function count_unique_paralogs(df, counts, species_names)
    println("Count genes belonging to unique paralogs for each species")
    unique_paralogs = []
    for species in species_names
        # species = species_names[1]
        idx = species_names .== species
        push!(unique_paralogs, sum((counts[:, idx] .== df.Total)[:,1]))
    end
    return(unique_paralogs)
end
@time unique_paralogs = count_unique_paralogs(df, counts, species_names)

# Single-copy gene orthologs
function count_single_copy_gene_orthologs(count, species_names)
    println("Count single-copy gene orthologs")
    idx = repeat([true], nrow(counts))
    for j in 1:ncol(counts)
        idx = idx .& (counts[:, j] .== 1)
    end
    single_copy_orthologs = repeat([sum(idx)], length(species_names))
    return(single_copy_orthologs)
end
@time single_copy_orthologs = count_single_copy_gene_orthologs(count, species_names)

# Multiple-copy orthologs
function count_multicopy_orthologs(counts, species_names)
    println("Count multi-copy orthologs")
    total_per_species = []
    for species in species_names
        push!(total_per_species, sum(counts[:, species.==species_names][:,1]))
    end
    multiple_orthologs = total_per_species - (unassigned + unique_paralogs + single_copy_orthologs)
    return(total_per_species, multiple_orthologs)
end
@time total_per_species, multiple_orthologs = count_multicopy_orthologs(counts, species_names)

# Merge gene count classifications per species and save
out = DataFrames.DataFrame(Species=species_names,
                           Total=total_per_species,
                           Multiple_Orthologs=multiple_orthologs,
                           Single_Copy_Orthologs=single_copy_orthologs,
                           Unique_Paralogs=unique_paralogs,
                           Unassigned_genes=unassigned
                          )
file = open(fname_out, "w")
CSV.write(file, out)
close(file)
' > count_genes_per_ortholog_paralog_classes.jl

time \
julia count_genes_per_ortholog_paralog_classes.jl \
        ORTHOGROUPS/orthogroups_gene_counts_families_go.out \
        ORTHOGROUPS/orthogroups_summarised_gene_counts.csv
```

## Locate paralogs in the Lolium rigidum genome for the Circos-like figure
```{sh}
ORT=${DIR_ORTHOGROUPS}/Orthogroups/Orthogroups.tsv
GFF=${DIR}/Lolium_rigidum.gff

echo '
using ProgressMeter
fname_paralogs = ARGS[1]
fname_annotations = ARGS[2]
fname_output = ARGS[3]
# fname_paralogs = "ORTHOGROUPS/OrthoFinder/Results_Apr05/Orthogroups/Orthogroups.tsv"
# fname_annotations = "Lolium_rigidum.gff"
# fname_output = "Lolium_rigidum.plg"

# Load gene names and start position into temproray vectors
function load_gene_names_and_coordinates(fname_annotations)
    temp_geneIDs = []
    temp_genes = []
    temp_chromosomes = []
    temp_positions = []
    file = open(fname_annotations, "r")
    seekend(file); n=position(file); seekstart(file)
    pb = Progress(n)
    while !eof(file)
        line = split(readline(file), "\t"[1])
        update!(pb, position(file))
        if line[1][1] .!= "#"[1]
            if line[3] == "CDS"

                if match(Regex("Name="), line[end]) != nothing
                    desc = split(line[end], ";"[1])
                    id = split(desc[match.(Regex("Dbxref=GeneID:"), desc) .!= nothing][1], ","[1])[1]
                    geneID = replace(id, "Dbxref=GeneID:" => "")
                    name = split(desc[match.(Regex("Name="), desc) .!= nothing][1], ","[1])[1]
                    gene = replace(name, "Name=" => "")
                    push!(temp_geneIDs, geneID)
                    push!(temp_genes, gene)
                    push!(temp_chromosomes, line[1])
                    push!(temp_positions, parse(Int, line[4]))
                end
            else
                continue
            end
        else
            continue
        end
    end
    close(file)
    @time idx = sortperm(temp_genes)
    temp_geneIDs = temp_geneIDs[idx]
    temp_genes = temp_genes[idx]
    temp_chromosomes = temp_chromosomes[idx]
    temp_positions = temp_positions[idx]

    # Use only one CDS/gene starting position
    genes = []
    chromosomes = []
    positions = []
    n = length(temp_genes)
    i = 1
    pb = Progress(n)
    while i < n
        id = temp_geneIDs[i]
        g = temp_genes[i]
        c = temp_chromosomes[i]
        p = temp_positions[i]
        push!(genes, g)
        push!(chromosomes, c)
        push!(positions, p)
        while (i < n) & ( (g == temp_genes[i]) | (id == temp_geneIDs[i]) )
            i += 1
        end
        update!(pb, i)
    end
    if genes[end] != temp_genes[n]
        push!(genes, temp_genes[n])
        push!(chromosomes, temp_chromosomes[n])
        push!(positions, temp_positions[n])
    end

    return(genes, chromosomes, positions)
end
genes, chromosomes, positions = load_gene_names_and_coordinates(fname_annotations)

# Add orthogroup labels, i.e. paralog classification
function add_paralogs(fname_paralogs, genes)
    paralogs = repeat(["None"], length(genes))
    file = open(fname_paralogs, "r")
    seekend(file); n=position(file); seekstart(file)
    pb = Progress(n)
    header = split(readline(file), "\t"[1])
    idx = header .== "Lolium_rigidum"
    while !eof(file)
        line = split(readline(file), "\t"[1])
        orthogroup = line[1]
        gene_names = replace.(split(line[idx][1], ", "), "Lolium_rigidum|" => "")
        for g in gene_names
            # g = gene_names[1]
            paralogs[genes .== g] .= orthogroup
        end
        update!(pb, position(file))
    end
    close(file)
    return(paralogs)
end
paralogs = add_paralogs(fname_paralogs, genes)

# Remove remove unclassified genes
idx = paralogs .!= "None"
paralogs = paralogs[idx]
genes = genes[idx]
chromosomes = chromosomes[idx]
positions = positions[idx]

# Save the list of gene names, chromosome, position, and ortholog info into a file
file = open(fname_output, "a")
for i in 1:length(genes)
    line = string(join([genes[i], chromosomes[i], positions[i], paralogs[i]], "\t"), "\n")
    write(file, line)
end
close(file)

# Save only the top 5 paralogs with the most genes
vec_p = unique(paralogs)
vec_p_counts = []
@showprogress for p in vec_p
    push!(vec_p_counts, sum(p .== paralogs))
end
idx = vec_p_counts .>= 10
vec_p = vec_p[idx]
file = open(replace(fname_output, ".plg"=>"-for_plotting.plg"), "a")
for i in 1:length(paralogs)
    if sum(paralogs[i] .== vec_p) > 0
        line = string(join([genes[i], chromosomes[i], positions[i], paralogs[i]], "\t"), "\n")
        write(file, line)
    end
end
close(file)
' > locate_paralogs.jl

time \
julia locate_paralogs.jl \
    ${ORT} \
    ${GFF} \
    ${GFF%.gff*}.plg

```

## What is the rate of gene family expansion and contraction in each species using CAFE?
```{sh}
ORTHOUT=${DIR}/ORTHOGROUPS/orthogroups_gene_counts_families_go.out
rev ${ORTHOUT} | cut -f5- | rev > col2_to_coln.tmp
awk -F'\t' '{print $(NF-1)}' ${ORTHOUT} > col1.tmp
paste -d'\t' col1.tmp col2_to_coln.tmp > counts.tmp
TREE=${DIR_ORTHOGROUPS}/Species_Tree/SpeciesTree_rooted.txt

# Run using the "Base" model where a single lambda (lambda = P(gene gain or gene loss per unit time)) is estimated.
time \
cafe5 \
    --infile counts.tmp \
    --tree ${TREE} \
    --cores 15 \
    --pvalue 0.01 \
    --output_prefix CAFE_Base_results

# Extract gene families significantly differentially mapped
echo $'#nexus\nbegin trees;' > Significant_trees.tre
grep "*" CAFE_Base_results/Base_asr.tre >> Significant_trees.tre
echo "end;" >> Significant_trees.tre

# Extract the number of expanded and contracted orthogroups (gene families)
cat CAFE_Base_results/Base_clade_results.txt

# Re-run with lambda_i ~ Gamma(alpha), for each of the i_th gene family category
time \
cafe5 \
    --infile counts.tmp \
    --tree ${TREE} \
    --n_gamma_cats 100 \
    --cores 15 \
    --pvalue 0.01 \
    --output_prefix CAFE_Gamma100_results

time \
cafe5 \
    --infile counts.tmp \
    --tree ${TREE} \
    --n_gamma_cats 10 \
    --cores 15 \
    --pvalue 0.01 \
    --output_prefix CAFE_Gamma10_results

time \
cafe5 \
    --infile counts.tmp \
    --tree ${TREE} \
    --n_gamma_cats 1000 \
    --cores 15 \
    --pvalue 0.01 \
    --output_prefix CAFE_Gamma1000_results

### Output
echo -e "Species\tExpansion\tContraction" > CONTRACTION_EXPANSION.txt
grep -v "^#" ${DIR}/CAFE_Gamma100_results/Gamma_clade_results.txt | \
    grep -v "^<" | \
    sed 's/<..>//g' | \
    sed 's/<.>//g' >> CONTRACTION_EXPANSION.txt
```

## Build tree using single-gene orthogroups

1. Identify single-copy orthogroups and their respective gene names across all 7 species:
```{sh}
awk '($2 == 1) && ($3 == 1) && ($4 == 1) && ($5 == 1) && ($6 == 1) && ($7 == 1) && ($8 == 1)'  $ORTHOUT | cut -f1 > single_gene_list.grep
grep -f single_gene_list.grep ${DIR_ORTHOGROUPS}/Orthogroups/Orthogroups.tsv > single_gene_list.geneNames
```

2. Extract the CDS of these genes (Outputs: ${ORTHONAME}.fasta [includes sequences from each species]):
```{sh}
echo '#!/bin/bash
i=$1
line=$(head -n${i} single_gene_list.geneNames | tail -n1)
ORTHONAME=$(echo $line | cut -d" " -f1)
for name in $(echo $line | cut -d" " -f2-)
do
    SPECIES=$(echo $name | cut -d"|" -f1)
    GENE_NAME=$(echo $name | cut -d"|" -f2)
    julia extract_sequence_using_name_query.jl \
        ${SPECIES}.cds \
        ${GENE_NAME} \
        ${ORTHONAME}-${SPECIES}.fasta \
        ${SPECIES} \
        false
done
if [ $(ls ${ORTHONAME}-*.fasta | wc -l) -eq 7 ]
then
    cat ${ORTHONAME}-*.fasta > ${ORTHONAME}.fasta
fi
rm ${ORTHONAME}-*.fasta
' > parallel_extract_single_gene_orthogroups.sh
chmod +x parallel_extract_single_gene_orthogroups.sh
time \
parallel \
./parallel_extract_single_gene_orthogroups.sh {} \
::: $(seq 1 $(cat single_gene_list.geneNames | wc -l))
```

3. Align CDS. NOTES: Set the final stop codons as "---", and internal stop codons as "NNN" so that PAML programs won't ask you to press enter to continue; also set frameshifts from "!" into "-". Also sort the alignment sequences, because MACSE jumbles them u for some reason with no option to return input order (Outputs: ${ORTHOLOG}.NT.cds [nucleotide alignments] and ${ORTHOLOG}.AA.prot [amino acid alignments])
```{sh}
echo '#!/bin/bash
f=$1
ORTHOLOG=${f%.fasta*}
# Align the CDS across species
java -Xmx8G \
    -jar macse_v2.06.jar \
    -prog alignSequences \
    -seq ${f} \
    -out_NT ${ORTHOLOG}.aligned.unsorted.cds.tmp \
    -out_AA ${ORTHOLOG}.aligned.unsorted.prot.tmp
# Convert stop codons and frameshifts as "---" for compatibility with downstream tools
java -Xmx8G \
    -jar macse_v2.06.jar \
    -prog exportAlignment \
    -align ${ORTHOLOG}.aligned.unsorted.cds.tmp \
    -codonForFinalStop --- \
    -codonForInternalStop NNN \
    -codonForExternalFS --- \
    -codonForInternalFS --- \
    -out_NT ${ORTHOLOG}.NT.cds \
    -out_AA ${ORTHOLOG}.AA.prot
# Clean-up
rm ${ORTHOLOG}*.tmp
' > parallel_align_cds.sh
chmod +x parallel_align_cds.sh
time \
parallel \
./parallel_align_cds.sh {} \
::: $(ls OG*.fasta)
```

4. Build the tree (Outputs: ORTHOGROUPS_SINGLE_GENE.[NT AA].timetree.nex)
```{sh}
time \
for TYPE in NT.cds AA.prot
do
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
### Extract sequences per species (Outputs: ${ORTHONAME}-${SPECIES}.fasta)
parallel \
julia extract_sequence_using_name_query.jl \
    {1}.${TYPE} \
    {2} \
    {1}-{2}.fasta \
    {1}-{2} \
    false \
::: $(ls *.${TYPE} | sed "s/.$TYPE//g") \
::: $(grep "^>" $(ls *.${TYPE} | head -n1) | sed 's/^>//g')

### Concatenate alignments per species (Outputs: ${SPECIES}.aln)
for SPECIES in $(grep "^>" $(ls *.${TYPE} | head -n1) | sed 's/^>//g')
do
    echo $SPECIES
    # Concatenate sequences
    echo ">${SPECIES}" > ${SPECIES}.aln.tmp
    cat *-${SPECIES}.fasta | sed '/^>/d' | sed -z 's/\n//g' >> ${SPECIES}.aln.tmp
    echo "" >> ${SPECIES}.aln.tmp
    julia reformat_fasta_sequence.jl ${SPECIES}.aln.tmp 50 ${SPECIES}-${TYPE%.*}.aln
    rm ${SPECIES}.aln.tmp
done

# Extract sequence lengths to build the sequence partitioning nexus file (Output: alignment_parition.nex)
SPECIES=$(grep "^>" $(ls *.${TYPE} | head -n1) | sed 's/^>//g' | head -n1)
echo '#nexus
begin sets;' > alignment_parition.nex
N0=0
for f in $(ls *-${SPECIES}.fasta)
do
    # f=$(ls *-${SPECIES}.fasta | head -n1)
    NAME=$(head -n1 $f | sed 's/>//g' | cut -d'-' -f1)
    N1=$(cat $f | sed '/^>/d' | sed -z 's/\n//g' | wc -c)
    START=$(echo "$N0 + 1" | bc)
    END=$(echo "$N0 + $N1" | bc)
    echo "charset $NAME = $START-$END;" >> alignment_parition.nex
    N0=$END
done
echo 'end;' >> alignment_parition.nex

### Concatenate species alignments (Output: ORTHOGROUPS_SINGLE_GENE.${TYPE%.*}.aln)
cat *-${TYPE%.*}.aln > ORTHOGROUPS_SINGLE_GENE.${TYPE%.*}.aln.tmp
mv ORTHOGROUPS_SINGLE_GENE.${TYPE%.*}.aln.tmp ORTHOGROUPS_SINGLE_GENE.${TYPE%.*}.aln
rm *-${TYPE%.*}.aln

# ### Convert codon and amino acid sequences from fasta into phylip format (Output: ${ORTHOLOG}.${TYPE%.*}.phylip)
# julia fasta_to_phylip.jl ORTHOGROUPS_SINGLE_GENE.${TYPE%.*}.aln

### Lookup divergence times between species (Output: dates.txt)
echo 'Arabidopsis_thaliana,Oryza_sativa     -160.00
Oryza_sativa,Zea_mays                        -50.00
Sorghum_bicolor,Zea_mays                     -12.19
Secale_cereale,Lolium_rigidum                -24.10
Lolium_perenne,Lolium_rigidum                 -1.65' > dates.txt

### Build tree
BOOTSTRAP_REPS=1000
THREADS=20
TIP_DATE=0
time \
iqtree2 \
    -s ORTHOGROUPS_SINGLE_GENE.${TYPE%.*}.aln \
    -p alignment_parition.nex \
    -B ${BOOTSTRAP_REPS} \
    -T ${THREADS} \
    --date dates.txt \
    --date-tip ${TIP_DATE} \
    --prefix ORTHOGROUPS_SINGLE_GENE.${TYPE%.*} \
    --redo
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
done

### Clean-up
rm OG*.fasta
rm OG*.NT.cds
rm OG*.AA.prot
for SPECIES in $(grep "^>" $(ls *.${TYPE} | head -n1) | sed 's/^>//g')
do
    rm *-${SPECIES}.fasta
done
```

## Estimate Ka/Ks (dN/dS) using pairs of paralogs per species (using maximum likelihood via "MS" or model selection method)
```{sh}
### Prepare parallelisable script
echo '#!/bin/bash
SPECIES=$1
COLUMN=$2
ORTHOGROUPS=$3
group=$4
# SPECIES=Lolium_rigidum
# COLUMN=4
# ORTHOGROUPS=${DIR_ORTHOGROUPS}/Orthogroups/Orthogroups.tsv
# group=$(cat orthogroups.tmp | head -n10 | tail -n1)
grep "^$group" $ORTHOGROUPS | \
    cut -f${COLUMN} | \
    sed "/^$/d" | \
    sed "s/$SPECIES|//g" > ${group}-gene_names.tmp
# set the first gene as the focal gene
focal_gene=$(cut -d"," -f1 ${group}-gene_names.tmp)
julia extract_sequence_using_name_query.jl \
    ${SPECIES}.cds \
    ${focal_gene} \
    ${group}.cds.tmp
# extract the rest of the sequences and concatenate into a single file together with the focal gene sequence
for gene in $(sed -z "s/, /\n/g" ${group}-gene_names.tmp | tail -n+2)
do
    # gene=$(sed -z "s/, /\n/g" ${group}-gene_names.tmp | tail -n+2 | head -n1)
    julia extract_sequence_using_name_query.jl \
        ${SPECIES}.cds \
        ${gene} \
        ${group}-Gi.cds.tmp
    cat ${group}-Gi.cds.tmp >> ${group}.cds.tmp
    rm ${group}-Gi.cds.tmp
done
# align
java -Xmx8G \
    -jar macse_v2.06.jar \
    -prog alignSequences \
    -seq ${group}.cds.tmp \
    -out_NT ${group}.NT.aln.tmp \
    -out_AA ${group}.AA.aln.tmp
rm ${group}.AA.aln.tmp
# pairwise Ka/Ks (dN/dS) estimation
touch ${group}.aln.tmp
g0=$(head -n2 ${group}.NT.aln.tmp | tail -n1)
for line in $(seq 4 2 $(cat ${group}.NT.aln.tmp | wc -l))
do
    echo ${group}-$(echo "($line / 2) - 1" | bc) >> ${group}.aln.tmp
    echo ${g0} >> ${group}.aln.tmp
    head -n${line} ${group}.NT.aln.tmp | tail -n1 >> ${group}.aln.tmp
    echo "" >> ${group}.aln.tmp
done
rm ${group}.NT.aln.tmp
KaKs_Calculator \
    -m MS \
    -i ${group}.aln.tmp \
    -o ${group}.kaks.out
# clean-up
rm ${group}*tmp
' > pairwise_within_species_KaKs_estimation.sh
chmod +x pairwise_within_species_KaKs_estimation.sh

### Run across gene families in parallel per species
time \
for SPECIES in $(head -n1 $ORTHOUT | cut -f2- | rev | cut -f5- | rev)
do
    # SPECIES=Lolium_rigidum
    COLUMN=$(head -n1 $ORTHOUT | sed -z "s/\t/\n/g" | grep -n $SPECIES - | cut -d":" -f1)
    awk -v col=${COLUMN} '($col > 9) && ($col <= 10)' $ORTHOUT | cut -f1 > orthogroups.tmp
    parallel ./pairwise_within_species_KaKs_estimation.sh \
        ${SPECIES} \
        ${COLUMN} \
        ${DIR_ORTHOGROUPS}/Orthogroups/Orthogroups.tsv \
        {} ::: $(cat orthogroups.tmp) ### Outputs: ${group}.kaks.out
    head -n1 $(ls *.kaks.out | head -n1) > ${SPECIES}.kaks
    for f in $(ls *.kaks.out)
    do
        ncol=$(head -n1 $f | awk '{print NF}')
        nrow=$(cat $f | wc -l)
        if [ $ncol -eq 22 ] && [ $nrow -gt 1 ]
        then
            tail -n+2 $f >> ${SPECIES}.kaks
        fi
    done
    rm *.kaks.out orthogroups.tmp
done



```


## Figure 2 plotting

```{R}
args = commandArgs(trailingOnly=TRUE)
# args = c("ORTHOGROUPS_SINGLE_GENE.NT.timetree.nex", "CONTRACTION_EXPANSION.txt", "ORTHOGROUPS/orthogroups_summarised_gene_counts.csv", "ORTHOGROUPS/orthogroups_gene_counts_families_go.out")
# args = c("ORTHOGROUPS_SINGLE_GENE.AA.timetree.nex", "CONTRACTION_EXPANSION.txt", "ORTHOGROUPS/orthogroups_summarised_gene_counts.csv", "ORTHOGROUPS/orthogroups_gene_counts_families_go.out")
fname_tree = args[1]
fname_conex = args[2]
fname_gene_groups = args[3]
fname_gene_counts = args[4]

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
conex$order = c(1, 6, 7, 2, 5, 3, 4)
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
gene_groups$order = c(1, 6, 7, 2, 5, 3, 4)
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


```


## Enrichment of stress-related genes: Do we have more ortholog members for herbicide and stress-related genes in Lolium rigidum compared with the other species?

Identify TSR and NTSR genes..
```{sh}

```


## dN/dS assessment: For the sress-related genes which are not more enriched, are there signs of selection?

## Phylogentic tree of stress-related genes: How did these stress-related gene which are under selection came about? 
