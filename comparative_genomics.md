# Comparative genomics

## Set working directories and the main orthogroup output filename (these variables and reiterated in their respective sections where they were used)
```sh}
DIR=/data/Lolium_rigidum_ASSEMBLY/COMPARATIVE_GENOMICS
DIR_ORTHOGROUPS=${DIR}/ORTHOGROUPS/OrthoFinder/Results_*
DIR_ORTHOGROUP_SEQS=${DIR}/ORTHOGROUPS/OrthoFinder/Results_*/Orthogroup_Sequences
DIR_PANTHER=${DIR}/PantherHMM_17.0/famlib/rel/PANTHER17.0_altVersion/hmmscoring/PANTHER17.0/books
GOT_PATHER=${DIR}/PantherHMM_17.0/PANTHER17.0_HMM_classifications
MERGED_ORTHOGROUPS=${DIR}/ORTHOGROUPS/orthogroups.faa
ORTHOUT=${DIR}/ORTHOGROUPS/orthogroups_gene_counts_families_go.out
TREE=${DIR_ORTHOGROUPS}/Species_Tree/SpeciesTree_rooted.txt
DIR_GENES=/data/Lolium_rigidum_ASSEMBLY/COMPARATIVE_GENOMICS/TSR_NTSR_GENES
cd $DIR
```

## Download genomes, genome annotations, and predicted CDS sequences
**Note:** The first column of the annotation files may not correspond to the chromosome IDs.
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
### Lolium perenne2
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/001/735/685/GCA_001735685.1_ASM173568v1/GCA_001735685.1_ASM173568v1_genomic.fna.gz
gunzip -c GCA_001735685.1_ASM173568v1_genomic.fna.gz > Lolium_perenne2.fasta
### Secale cereale
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/016/097/815/GCA_016097815.1_HAU_Weining_v1.0/GCA_016097815.1_HAU_Weining_v1.0_genomic.fna.gz
gunzip -c GCA_016097815.1_HAU_Weining_v1.0_genomic.fna.gz > Secale_cereale.fasta
```

## Install GeMoMa (Gene Model Mapper)
```{sh}
wget http://www.jstacs.de/download.php?which=GeMoMa
mv 'download.php?which=GeMoMa' GeMoMa.zip
unzip GeMoMa.zip
java -jar GeMoMa-1.8.jar CLI -h
PATH=${PATH}:${DIR} ### for GeMoMa.jar
```

## Install Java for GeMoMa (see below)
```{sh}
sudo apt install -y default-jre
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
./orthofinder -h
PATH=${PATH}:$(pwd)
cd -
```

## Install HMMER for mapping CDS to PantherHMM gene family models
```{sh}
sudo apt install hmmer
```

## Install CAFE5 to analyse gene family evolution
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
java -Xmx250G -jar macse_v2.06.jar -help
```

## Install IQ-TREE for building trees with fossil root dates
```{sh}
sudo apt install libeigen3-dev
wget https://github.com/Cibiv/IQ-TREE/releases/download/v2.0.7/iqtree-2.0.7-Linux.tar.gz
tar -xvzf iqtree-2.0.7-Linux.tar.gz
cd iqtree-2.0.7-Linux/bin
./iqtree2 -h
PATH=${PATH}:$(pwd)
cd -
```

## Install PAML (Phylogenetic Analysis by Maximum Likelihood) which includes MCMCTree for Bayesian phylogenetic analysis
```{sh}
wget http://abacus.gene.ucl.ac.uk/software/paml4.9j.tgz
tar -xvzf paml4.9j.tgz
cd paml4.9j/
PATH=${PATH}/bin
PATH=${PATH}/src
cd -
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

## Install R::ape
1. Install BLAS and LAPACK libraries for R::ape package
```{sh}
sudo apt install -y libblas-dev liblapack-dev
```

2. Install R::ape library
```{R}
install.packages("ape")
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
for REF in Lolium_perenne Lolium_perenne2 Secale_cereale
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
    cp ${DIR}/${REF}/final_annotation.gff ${REF}.gff
    cp ${DIR}/${REF}/predicted_cds.fasta ${REF}.cds
    cp ${DIR}/${REF}/predicted_proteins.fasta ${REF}.faa
    rm -R ${DIR}/${REF}
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
    -t 31

DIR_ORTHOGROUPS=${DIR}/ORTHOGROUPS/OrthoFinder/Results_*/
```

## Assign orthogroups into gene families
1. Define the location of the 15,619 protein family HMMs
```{sh}
DIR_PANTHER=${DIR}/PantherHMM_17.0/famlib/rel/PANTHER17.0_altVersion/hmmscoring/PANTHER17.0/books
GOT_PATHER=${DIR}/PantherHMM_17.0/PANTHER17.0_HMM_classifications
```

2. Iteratively, for each genome's predicted protein sequences run hmmsearch in paralel for each PantherHMM protein family
```{sh}
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
```

3. Find PantherHMM protein families for each orthogroup
```{sh}
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
```

4. Find the best fitting gene family to each unique sequence per orthogroup. This means that each orthogroup can have multiple gene families. Next, add family name and GO terms to each gene family.
```{sh}
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

## Infer gene family expansion and contraction in each species using CAFE5 (at alpha=1%)
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

### Output using n_gamma_cats=1,000
echo -e "Species\tExpansion\tContraction" > CONTRACTION_EXPANSION.txt
grep -v "^#" ${DIR}/CAFE_Gamma100_results/Gamma_clade_results.txt | \
    grep -v "^<" | \
    sed 's/<..>//g' | \
    sed 's/<.>//g' >> CONTRACTION_EXPANSION.txt
```

## Build tree using single-gene orthogroups
1. Identify single-copy orthogroups and their respective gene names across all 7 species:
```{sh}
awk '($2 == 1) && ($3 == 1) && ($4 == 1) && ($5 == 1) && ($6 == 1) && ($7 == 1) && ($8 == 1)' $ORTHOUT | cut -f1 > single_gene_list.grep
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

3. Align CDS (Outputs: ${ORTHOLOG}.NT.cds [nucleotide alignments] and ${ORTHOLOG}.AA.prot [amino acid alignments])
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
rm ${ORTHOLOG}*.tmp ${ORTHOLOG}.AA.prot # we are not using the amino acid sequences
' > parallel_align_cds.sh
chmod +x parallel_align_cds.sh
time \
parallel \
./parallel_align_cds.sh {} \
::: $(ls OG*.fasta)
```

4. Build the tree (Outputs: ORTHOGROUPS_SINGLE_GENE.NT.treefile [for PAML::codeml], ORTHOGROUPS_SINGLE_GENE.NT.timetree.nex [for plotting]; ORTHOGROUPS_SINGLE_GENE.NT.aln [for prepping PAML::codeml input]; and alignment_parition.NT.nex [for prepping PAML::codeml input])
```{sh}
TYPE=NT.cds
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
    for f in $(ls *-${SPECIES}.fasta)
    do
        sed '/^>/d' $f | sed -z 's/\n//g' >> ${SPECIES}.aln.tmp
    done
    echo "" >> ${SPECIES}.aln.tmp
    julia reformat_fasta_sequence.jl ${SPECIES}.aln.tmp 50 ${SPECIES}-${TYPE%.*}.aln
    rm ${SPECIES}.aln.tmp
done

# Extract sequence lengths to build the sequence partitioning nexus file (Output: alignment_parition.${TYPE%.*}.nex)
SPECIES=$(grep "^>" $(ls *.${TYPE} | head -n1) | sed 's/^>//g' | head -n1)
echo '#nexus
begin sets;' > alignment_parition.${TYPE%.*}.nex
N0=0
for f in $(ls *-${SPECIES}.fasta)
do
    # f=$(ls *-${SPECIES}.fasta | head -n1)
    NAME=$(head -n1 $f | sed 's/>//g' | cut -d'-' -f1)
    N1=$(cat $f | sed '/^>/d' | sed -z 's/\n//g' | wc -c)
    START=$(echo "$N0 + 1" | bc)
    END=$(echo "$N0 + $N1" | bc)
    echo "charset $NAME = $START-$END;" >> alignment_parition.${TYPE%.*}.nex
    N0=$END
done
echo 'end;' >> alignment_parition.${TYPE%.*}.nex

### Concatenate species alignments (Output: ORTHOGROUPS_SINGLE_GENE.${TYPE%.*}.aln)
cat *-${TYPE%.*}.aln > ORTHOGROUPS_SINGLE_GENE.${TYPE%.*}.aln.tmp
mv ORTHOGROUPS_SINGLE_GENE.${TYPE%.*}.aln.tmp ORTHOGROUPS_SINGLE_GENE.${TYPE%.*}.aln
rm *-${TYPE%.*}.aln

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
    -p alignment_parition.${TYPE%.*}.nex \
    -B ${BOOTSTRAP_REPS} \
    -T ${THREADS} \
    --date dates.txt \
    --date-tip ${TIP_DATE} \
    --prefix ORTHOGROUPS_SINGLE_GENE.${TYPE%.*} \
    --redo
```

5. **Additional**: Compute pairwise 4DTv
```{sh}
time \
parallel \
julia calculate_4DTv.jl {1} {1}.4DTv.tmp \
    ::: $(ls *.NT.cds)

```

6. Clean-up
```{sh}
rm OG*.fasta
rm OG*.NT.cds
rm single_gene_list.*
rm dates.txt
```

## Assess whole genome duplication (WGD) events using the distribution of four-fold degenerate sites (4DTv) across dual-copy paralogs within genomes and across sing-copy gene orthologs between pairs of species
1. Prepare script to extract CDS, and align in parallel
```{sh}
echo '#!/bin/bash
j=$1
# j=1017
line=$(head -n${j} dual_gene_list.geneNames | tail -n1)
ORTHONAME=$(echo $line | cut -d" " -f1)
# Extract CDS of these dual-copy genes
for name in $(echo $line | cut -d" " -f2- | sed -z "s/, /\n/g")
do
    # name=$(echo $line | cut -d" " -f2- | sed -z "s/, /\n/g" | head -n2 | tail -n1)
    SPECIES=$(echo $name | cut -d"|" -f1)
    GENE_NAME=$(echo $name | cut -d"|" -f2 | sed -z "s/\r//g")
    julia extract_sequence_using_name_query.jl \
        ${SPECIES}.cds \
        ${GENE_NAME} \
        ${ORTHONAME}-${GENE_NAME}.fasta \
        ${GENE_NAME} \
        false
done
# Concatenate the two sequences
if [ $(ls ${ORTHONAME}-*.fasta | wc -l) -eq 2 ]
then
    cat ${ORTHONAME}-*.fasta > ${ORTHONAME}.fasta
fi
# Align the CDS across species
java -Xmx8G \
    -jar macse_v2.06.jar \
    -prog alignSequences \
    -seq ${ORTHONAME}.fasta \
    -out_NT ${ORTHONAME}.aligned.unsorted.cds.tmp \
    -out_AA ${ORTHONAME}.aligned.unsorted.prot.tmp
# Convert stop codons and frameshifts as "---" for compatibility with downstream tools
java -Xmx8G \
    -jar macse_v2.06.jar \
    -prog exportAlignment \
    -align ${ORTHONAME}.aligned.unsorted.cds.tmp \
    -codonForFinalStop --- \
    -codonForInternalStop NNN \
    -codonForExternalFS --- \
    -codonForInternalFS --- \
    -out_NT ${ORTHONAME}.NT.cds \
    -out_AA ${ORTHONAME}.AA.prot
# Calculate 4DTv (i.e. the ratio of the number of 4-fold degenerate codons with trnasversion and the total number of 4-fold degenerate codons)
julia calculate_4DTv.jl ${ORTHONAME}.NT.cds ${ORTHONAME}.4DTv.tmp
# Calculate divergence time with KaKs_Calculator
grep "^>" ${ORTHONAME}.NT.cds | \
    sed "s/^>//g" | \
    sed -z "s/\n/:/g" | \
    sed -z "s/:$/\n/g" > ${ORTHONAME}.axt.tmp
grep -v "^>" ${ORTHONAME}.NT.cds >> ${ORTHONAME}.axt.tmp
KaKs_Calculator -i ${ORTHONAME}.axt.tmp \
                -m MA \
                -o ${ORTHONAME}.kaks.tmp
divergence_time=$(tail -n1 ${ORTHONAME}.kaks.tmp | cut -f16)
sed -i -z "s/\n/\t${divergence_time}\n/g" ${ORTHONAME}.4DTv.tmp
# Clean-up
rm ${ORTHONAME}*.fasta
rm ${ORTHONAME}.aligned.unsorted*.tmp ${ORTHONAME}.NT.cds ${ORTHONAME}.AA.prot
rm ${ORTHONAME}.axt.tmp ${ORTHONAME}.kaks.tmp
' > parallel_extract_dual_gene_orthogroups.sh
chmod +x parallel_extract_dual_gene_orthogroups.sh
```

2. Identify dual-copy paralogs per species, align, and estimate 4DTv
```{sh}
head -n1 ${ORTHOUT} | rev | cut -f5- | rev | cut -f2- | sed -z "s/\t/\n/g" > species_names.tmp
time \
for i in $(seq 1 $(cat species_names.tmp | wc -l))
do
    # i=1
    ### Extract species name
    SPECIES=$(head -n${i} species_names.tmp | tail -n1)
    idx=$(echo $i + 1 | bc)
    ### Exract names of orthogroups with 2 copies in the current species
    awk -v col="$idx" '($col == 2)' $ORTHOUT | cut -f1 > dual_gene_list.grep
    ### Extract names of the genes of these dual-copy orthogroups
    grep -f dual_gene_list.grep ${DIR_ORTHOGROUPS}/Orthogroups/Orthogroups.tsv | cut -f1,${idx} > dual_gene_list.geneNames
    ### Extract CDS, and align in parallel
    parallel \
    ./parallel_extract_dual_gene_orthogroups.sh {} \
    ::: $(seq 1 $(cat dual_gene_list.geneNames | wc -l))
    ### Concatenate 4DTv estimates
    cat *.4DTv.tmp > ${SPECIES}.4DTv
    ### Clean-up
    rm *.4DTv.tmp
    rm dual_gene_list.grep
    rm dual_gene_list.geneNames
done
rm species_names.tmp
```

## Identify herbicide TSR and NTSR genes
1. Download protein sequences of genes from UniProt (https://www.uniprot.org) (Outputs: ${GENE}.faa)
```{sh}
###############################
### SET WORKING DIRECTORIES ###
###############################
DIR_ORTHOGROUP_SEQS=${DIR}/ORTHOGROUPS/OrthoFinder/Results_*/Orthogroup_Sequences
DIR_GENES=/data/Lolium_rigidum_ASSEMBLY/COMPARATIVE_GENOMICS/TSR_NTSR_GENES
cd $DIR_GENES

#################
### CLETHODIM ###
#################
### TARGET: Acetyle Co-A carboxylase
echo 'acetyl-coenzyme a carboxylase AND organism:"Arabidopsis thaliana (Mouse-ear cress) [3702]" '
wget https://www.uniprot.org/uniprot/Q38970.fasta
wget https://www.uniprot.org/uniprot/Q9LD43.fasta
echo 'acetyl-coenzyme a carboxylase AND organism:"Oryza sativa (Rice) [4530]"'
wget https://www.uniprot.org/uniprot/P0C2Y2.fasta
echo 'acetyl-coenzyme a carboxylase AND organism:"Zea mays (Maize) [4577]"'
wget https://www.uniprot.org/uniprot/Q41743.fasta
wget https://www.uniprot.org/uniprot/A0A1D6J3Q3.fasta
echo 'acetyl-coenzyme a carboxylase AND organism:"Lolium rigidum (Annual ryegrass) [89674]"'
wget https://www.uniprot.org/uniprot/Q3HWZ8.fasta
wget https://www.uniprot.org/uniprot/A5JJU5.fasta
wget https://www.uniprot.org/uniprot/Q7XHK2.fasta
echo 'acetyl-coenzyme a carboxylase AND organism:"Lolium multiflorum (Italian ryegrass) (Lolium perenne subsp. multiflorum) [4521]"'
wget https://www.uniprot.org/uniprot/A0A4D6X5K8.fasta
wget https://www.uniprot.org/uniprot/Q8VWG1.fasta
cat *.fasta > ACCase.faa
rm *.fasta

##################
### GLYPHOSATE ###
##################
### TARGET: 5-enolpyruvylshikimate-3-phosphate synthase
echo 'enolpyruvylshikimate phosphate synthase AND organism:"Arabidopsis thaliana (Mouse-ear cress) [3702]"'
wget https://www.uniprot.org/uniprot/P05466.fasta
wget https://www.uniprot.org/uniprot/P57720.fasta
echo 'enolpyruvylshikimate phosphate synthase AND organism:"Oryza sativa (Rice) [4530]"'
echo "NONE"
echo 'enolpyruvylshikimate phosphate synthase AND organism:"Zea mays (Maize) [4577]"'
echo "NONE"
echo 'enolpyruvylshikimate phosphate synthase AND organism:"Lolium rigidum (Annual ryegrass) [89674]"'
wget https://www.uniprot.org/uniprot/F4YBG5.fasta
wget https://www.uniprot.org/uniprot/F4YBG4.fasta
wget https://www.uniprot.org/uniprot/B1P936.fasta
wget https://www.uniprot.org/uniprot/B1P932.fasta
wget https://www.uniprot.org/uniprot/Q2PT59.fasta
wget https://www.uniprot.org/uniprot/Q2PT49.fasta
echo 'enolpyruvylshikimate phosphate synthase AND organism:"Lolium multiflorum (Italian ryegrass) (Lolium perenne subsp. multiflorum) [4521]"'
wget https://www.uniprot.org/uniprot/A0A2D2AP30.fasta
wget https://www.uniprot.org/uniprot/A0A2D2AP54.fasta
cat *.fasta > EPSPS.faa
rm *.fasta

#################
### INTERCEPT ### Imazamox, and
################# Imazapyr
### TARGET: Acetolactate synthase a.k.a. acetohydroxy acid synthase
echo 'acetolactate synthase AND organism:"Arabidopsis thaliana (Mouse-ear cress) [3702]"'
wget https://www.uniprot.org/uniprot/A0A178VL64.fasta
wget https://www.uniprot.org/uniprot/P17597.fasta
echo 'acetolactate synthase AND organism:"Oryza sativa (Rice) [4530]"'
wget https://www.uniprot.org/uniprot/A0A5J6D4R6.fasta
wget https://www.uniprot.org/uniprot/Q01LD9.fasta
echo 'acetolactate synthase AND organism:"Zea mays (Maize) [4577]"'
wget https://www.uniprot.org/uniprot/Q41769.fasta
wget https://www.uniprot.org/uniprot/K7TWQ8.fasta
echo 'acetolactate synthase AND organism:"Lolium rigidum (Annual ryegrass) [89674]"'
wget https://www.uniprot.org/uniprot/A0A5B9T5W5.fasta
wget https://www.uniprot.org/uniprot/A7XBQ0.fasta
echo 'acetolactate synthase AND organism:"Lolium multiflorum (Italian ryegrass) (Lolium perenne subsp. multiflorum) [4521]"'
wget https://www.uniprot.org/uniprot/Q9FUD0.fasta
wget https://www.uniprot.org/uniprot/A0A2R4NC54.fasta
cat *.fasta > ALS.faa
rm *.fasta

###############
### LUXIMAX ###
############### Cinmethylin
### TARGET: acyl-acp thioesterase
echo 'acyl-acp thioesterase AND organism:"Arabidopsis thaliana (Mouse-ear cress) [3702]"'
wget https://www.uniprot.org/uniprot/Q42561.fasta
wget https://www.uniprot.org/uniprot/Q9SJE2.fasta
wget https://www.uniprot.org/uniprot/Q9SV64.fasta
wget https://www.uniprot.org/uniprot/Q8W583.fasta
wget https://www.uniprot.org/uniprot/F4HX80.fasta
echo 'acyl-acp thioesterase AND organism:"Oryza sativa (Rice) [4530]"'
echo "NONE!"
echo 'acyl-acp thioesterase AND organism:"Zea mays (Maize) [4577]"'
wget https://www.uniprot.org/uniprot/A0A3L6FUI9.fasta
wget https://www.uniprot.org/uniprot/A0A077D597.fasta
wget https://www.uniprot.org/uniprot/A0A3L6FVN3.fasta
wget https://www.uniprot.org/uniprot/A0A3L6E7J9.fasta
wget https://www.uniprot.org/uniprot/K7V747.fasta
echo 'acyl-acp thioesterase AND organism:"Lolium rigidum (Annual ryegrass) [89674]"'
echo "NONE!"
echo 'acyl-acp thioesterase AND organism:"Lolium multiflorum (Italian ryegrass) (Lolium perenne subsp. multiflorum) [4521]"'
echo "NONE!"
cat *.fasta > AACPT.faa
rm *.fasta

#################
### OVERWATCH ###
################# Bixlozone
### TARGET: deoxyxylulose 5-phosphate synthase (DXS)
echo 'deoxyxylulose 5-phosphate synthase AND organism:"Arabidopsis thaliana (Mouse-ear cress) [3702]"'
wget https://www.uniprot.org/uniprot/Q38854.fasta
wget https://www.uniprot.org/uniprot/Q8LAD0.fasta
wget https://www.uniprot.org/uniprot/Q9XFS9.fasta
wget https://www.uniprot.org/uniprot/Q0WUB4.fasta
echo 'deoxyxylulose 5-phosphate synthase AND organism:"Oryza sativa (Rice) [4530]"'
echo "NONE!"
echo 'deoxyxylulose 5-phosphate synthase AND organism:"Zea mays (Maize) [4577]"'
echo "NONE!"
echo 'deoxyxylulose 5-phosphate synthase AND organism:"Lolium rigidum (Annual ryegrass) [89674]"'
echo "NONE!"
echo 'deoxyxylulose 5-phosphate synthase AND organism:"Lolium multiflorum (Italian ryegrass) (Lolium perenne subsp. multiflorum) [4521]"'
echo "NONE!"
cat *.fasta > DXS.faa
rm *.fasta

#################
### BOXERGOLD ### S-metolachlor, and
################# Prosulfocarb
##############
### SAKURA ###
############## Pyroxasulfone
#################
### TRIALLATE ###
#################
### TARGET: fatty acid elongase (FAE)
echo 'fatty acid elongase AND organism:"Arabidopsis thaliana (Mouse-ear cress) [3702]"'
wget https://www.uniprot.org/uniprot/Q38860.fasta
wget https://www.uniprot.org/uniprot/V9Z5W3.fasta
wget https://www.uniprot.org/uniprot/Q5XEP9.fasta
wget https://www.uniprot.org/uniprot/Q9MAM3.fasta
wget https://www.uniprot.org/uniprot/Q570B4.fasta
echo 'fatty acid elongase AND organism:"Oryza sativa (Rice) [4530]"'
echo "NONE!"
echo 'fatty acid elongase AND organism:"Zea mays (Maize) [4577]"'
wget https://www.uniprot.org/uniprot/A0A1D6MV48.fasta
wget https://www.uniprot.org/uniprot/A0A1D6L0Y4.fasta
wget https://www.uniprot.org/uniprot/A0A3L6FEW8.fasta
wget https://www.uniprot.org/uniprot/Q6A4M2.fasta
wget https://www.uniprot.org/uniprot/A0A1D6P2P8.fasta
echo 'fatty acid elongase AND organism:"Lolium rigidum (Annual ryegrass) [89674]"'
echo "NONE!"
echo 'fatty acid elongase AND organism:"Lolium multiflorum (Italian ryegrass) (Lolium perenne subsp. multiflorum) [4521]"'
echo "NONE!"
cat *.fasta > FAE.faa
rm *.fasta

################
### ATRAZINE ###
################
### TARGET: Photosystem II protein D1
echo 'Photosystem II protein D1 AND organism:"Arabidopsis thaliana (Mouse-ear cress) [3702]"'
wget https://www.uniprot.org/uniprot/P83755.fasta
wget https://www.uniprot.org/uniprot/A0A1B1W4S7.fasta
echo 'Photosystem II protein D1 AND organism:"Oryza sativa (Rice) [4530]"'
wget https://www.uniprot.org/uniprot/P0C432.fasta
wget https://www.uniprot.org/uniprot/A0A0K0LK42.fasta
echo 'Photosystem II protein D1 AND organism:"Zea mays (Maize) [4577]"'
wget https://www.uniprot.org/uniprot/P48183.fasta
wget https://www.uniprot.org/uniprot/A0A3L6DET6.fasta
echo 'Photosystem II protein D1 AND organism:"Lolium rigidum (Annual ryegrass) [89674]"'
echo "NONE!"
echo 'Photosystem II protein D1 AND organism:"Lolium multiflorum (Italian ryegrass) (Lolium perenne subsp. multiflorum) [4521]"'
wget https://www.uniprot.org/uniprot/B8Q684.fasta
cat *.fasta > psbA.faa
rm *.fasta

################
### PARAQUAT ###
################
############################
### DETOXIFICATION GENES ###
############################

##################################
### 1.) superoxide dismutase (SOD)
echo 'superoxide dismutase AND organism:"Arabidopsis thaliana (Mouse-ear cress) [3702]"'
wget https://www.uniprot.org/uniprot/P24704.fasta
wget https://www.uniprot.org/uniprot/O78310.fasta
echo 'superoxide dismutase AND organism:"Oryza sativa (Rice) [4530]"'
wget https://www.uniprot.org/uniprot/Q43803.fasta
wget https://www.uniprot.org/uniprot/Q01JW6.fasta
echo 'superoxide dismutase AND organism:"Zea mays (Maize) [4577]"'
wget https://www.uniprot.org/uniprot/P09233.fasta
wget https://www.uniprot.org/uniprot/P23346.fasta
echo 'superoxide dismutase AND organism:"Lolium rigidum (Annual ryegrass) [89674]"'
echo "NONE!"
echo 'superoxide dismutase AND organism:"Lolium multiflorum (Italian ryegrass) (Lolium perenne subsp. multiflorum) [4521]"'
echo "NONE!"
cat *.fasta > SOD.faa
rm *.fasta

##################################
### 2.) ascorbate peroxidase (APX)
echo 'ascorbate peroxidase AND organism:"Arabidopsis thaliana (Mouse-ear cress) [3702]"'
wget https://www.uniprot.org/uniprot/Q42593.fasta
wget https://www.uniprot.org/uniprot/Q42564.fasta
echo 'ascorbate peroxidase AND organism:"Oryza sativa (Rice) [4530]"'
wget https://www.uniprot.org/uniprot/Q9ST81.fasta
wget https://www.uniprot.org/uniprot/Q25AM6.fasta
echo 'ascorbate peroxidase AND organism:"Zea mays (Maize) [4577]"'
wget https://www.uniprot.org/uniprot/B4G031.fasta
wget https://www.uniprot.org/uniprot/C0HEE6.fasta
echo 'ascorbate peroxidase AND organism:"Lolium rigidum (Annual ryegrass) [89674]"'
echo "NONE!"
echo 'ascorbate peroxidase AND organism:"Lolium multiflorum (Italian ryegrass) (Lolium perenne subsp. multiflorum) [4521]"'
echo "NONE!"
cat *.fasta > APX.faa
rm *.fasta

#######################################
### 3.) glutathione S-transferase (GST)
echo 'glutathione S-transferase AND organism:"Arabidopsis thaliana (Mouse-ear cress) [3702]"'
wget https://www.uniprot.org/uniprot/Q96266.fasta
wget https://www.uniprot.org/uniprot/Q9FRL8.fasta
echo 'glutathione S-transferase AND organism:"Oryza sativa (Rice) [4530]"'
echo "NONE!"
echo 'glutathione S-transferase AND organism:"Zea mays (Maize) [4577]"'
wget https://www.uniprot.org/uniprot/P12653.fasta
wget https://www.uniprot.org/uniprot/P46420.fasta
wget https://www.uniprot.org/uniprot/P04907.fasta
wget https://www.uniprot.org/uniprot/Q9ZP62.fasta
wget https://www.uniprot.org/uniprot/C0PG07.fasta
wget https://www.uniprot.org/uniprot/Q9FQC6.fasta
wget https://www.uniprot.org/uniprot/A0A3L6E9E2.fasta
wget https://www.uniprot.org/uniprot/A0A1D6LVU5.fasta
echo 'glutathione S-transferase AND organism:"Lolium rigidum (Annual ryegrass) [89674]"'
echo "NONE!"
echo 'glutathione S-transferase AND organism:"Lolium multiflorum (Italian ryegrass) (Lolium perenne subsp. multiflorum) [4521]"'
echo "NONE!"
cat *.fasta > GST.faa
rm *.fasta

#############################################
### 4.) monodehydroascorbate reductase (MDAR)
echo 'monodehydroascorbate reductase AND organism:"Arabidopsis thaliana (Mouse-ear cress) [3702]"'
wget https://www.uniprot.org/uniprot/P92947.fasta
wget https://www.uniprot.org/uniprot/Q9LK94.fasta
wget https://www.uniprot.org/uniprot/Q9LFA3.fasta
wget https://www.uniprot.org/uniprot/F4J849.fasta
echo 'monodehydroascorbate reductase AND organism:"Oryza sativa (Rice) [4530]"'
echo "NONE!"
echo 'monodehydroascorbate reductase AND organism:"Zea mays (Maize) [4577]"'
wget https://www.uniprot.org/uniprot/C4J4E4.fasta
wget https://www.uniprot.org/uniprot/B8A028.fasta
wget https://www.uniprot.org/uniprot/A0A3L6G1K4.fasta
wget https://www.uniprot.org/uniprot/B4FQK0.fasta
wget https://www.uniprot.org/uniprot/A0A3L6ESY8.fasta
wget https://www.uniprot.org/uniprot/C0P4M0.fasta
echo 'monodehydroascorbate reductase AND organism:"Lolium rigidum (Annual ryegrass) [89674]"'
echo "NONE!"
echo 'monodehydroascorbate reductase AND organism:"Lolium multiflorum (Italian ryegrass) (Lolium perenne subsp. multiflorum) [4521]"'
echo "NONE!"
cat *.fasta > MDAR.faa
rm *.fasta

####################################
### 5.) glutathione peroxidase (GPX)
echo 'glutathione peroxidase AND organism:"Arabidopsis thaliana (Mouse-ear cress) [3702]"'
wget https://www.uniprot.org/uniprot/O22850.fasta
wget https://www.uniprot.org/uniprot/Q9LYB4.fasta
echo 'glutathione peroxidase AND organism:"Oryza sativa (Rice) [4530]"'
wget https://www.uniprot.org/uniprot/Q259Q9.fasta
wget https://www.uniprot.org/uniprot/Q9FEV2.fasta
wget https://www.uniprot.org/uniprot/Q8L8G3.fasta
echo 'glutathione peroxidase AND organism:"Zea mays (Maize) [4577]"'
wget https://www.uniprot.org/uniprot/Q6JAH6.fasta
wget https://www.uniprot.org/uniprot/B4FRF0.fasta
wget https://www.uniprot.org/uniprot/B6T5N2.fasta
wget https://www.uniprot.org/uniprot/A0A1D6K2F9.fasta
wget https://www.uniprot.org/uniprot/A0A3L6DIL3.fasta
echo 'glutathione peroxidase AND organism:"Lolium rigidum (Annual ryegrass) [89674]"'
echo "NONE!"
echo 'glutathione peroxidase AND organism:"Lolium multiflorum (Italian ryegrass) (Lolium perenne subsp. multiflorum) [4521]"'
echo "NONE!"
cat *.fasta > GPX.faa
rm *.fasta

################################
### 7.) cytochrome P450 (CYP450)
echo 'cytochrome P450 AND organism:"Arabidopsis thaliana (Mouse-ear cress) [3702]"'
wget https://www.uniprot.org/uniprot/O65782.fasta
wget https://www.uniprot.org/uniprot/Q9SMP5.fasta
echo 'cytochrome P450 AND organism:"Oryza sativa (Rice) [4530]"'
wget https://www.uniprot.org/uniprot/Q712G4.fasta
wget https://www.uniprot.org/uniprot/Q8S3F7.fasta
echo 'cytochrome P450 AND organism:"Zea mays (Maize) [4577]"'
wget https://www.uniprot.org/uniprot/Q43246.fasta
wget https://www.uniprot.org/uniprot/A0A1D6IZ40.fasta
echo 'cytochrome P450 AND organism:"Lolium rigidum (Annual ryegrass) [89674]"'
wget https://www.uniprot.org/uniprot/A5LFV4.fasta
wget https://www.uniprot.org/uniprot/Q6F4F2.fasta
wget https://www.uniprot.org/uniprot/Q6F4F4.fasta
wget https://www.uniprot.org/uniprot/A0A858NLU3.fasta
echo 'cytochrome P450 AND organism:"Lolium multiflorum (Italian ryegrass) (Lolium perenne subsp. multiflorum) [4521]"'
echo "NONE!"
cat *.fasta > CYP450.faa
rm *.fasta

##############################################
### 8.) ATP-binding cassette transporter (ABC)
echo 'ABC transporter AND organism:"Arabidopsis thaliana (Mouse-ear cress) [3702]"'
wget https://www.uniprot.org/uniprot/Q9LVM1.fasta
wget https://www.uniprot.org/uniprot/Q94FB9.fasta
echo 'ABC transporter AND organism:"Oryza sativa (Rice) [4530]"'
wget https://www.uniprot.org/uniprot/A0A650FN06.fasta
wget https://www.uniprot.org/uniprot/Q94HL6.fasta
wget https://www.uniprot.org/uniprot/Q01HZ4.fasta
wget https://www.uniprot.org/uniprot/Q25AJ5.fasta
echo 'ABC transporter AND organism:"Zea mays (Maize) [4577]"'
wget https://www.uniprot.org/uniprot/A7KVC2.fasta
wget https://www.uniprot.org/uniprot/C0P683.fasta
wget https://www.uniprot.org/uniprot/A0A317YHC9.fasta
wget https://www.uniprot.org/uniprot/A0A1D6E8S2.fasta
echo 'ABC transporter AND organism:"Lolium rigidum (Annual ryegrass) [89674]"'
echo "NONE!"
echo 'ABC transporter AND organism:"Lolium multiflorum (Italian ryegrass) (Lolium perenne subsp. multiflorum) [4521]"'
echo "NONE!"
cat *.fasta > ABC.faa
rm *.fasta
```

2. Generate BLAST database for each orthogroup (Outputs: ${ORTHOGROUP}.*)
```{sh}
echo '#!/bin/bash
f=$1
# f=$(find $DIR_ORTHOGROUP_SEQS -name "*fa" | head -n1)
makeblastdb \
    -in $f \
    -dbtype prot
' > prepare_blastdb_per_orthogroup.sh
chmod +x prepare_blastdb_per_orthogroup.sh

find $DIR_ORTHOGROUP_SEQS -name "*.fa" > orthogroup_faa_list.tmp
split --additional-suffix=".tmp" --lines 10000 orthogroup_faa_list.tmp

time \
for f in $(ls x*.tmp)
do
parallel ./prepare_blastdb_per_orthogroup.sh {} ::: $(cat $f)
done
```

3. Blastp (Outputs: ${GENE}-${ORTHOGROUP}.blastout)
```{sh}
echo '#!/bin/bash
GENE=$1
f=$2
# GENE=EPSPS
# f=$(find $DIR_ORTHOGROUP_SEQS -name "*.fa" | head -n10 | tail -n1)
ortname=$(basename $f)
blastp \
    -db $f \
    -query ${GENE}.faa \
    -out ${GENE}-${ortname%.fa*}.blastout \
    -evalue 1e-10 \
    -max_hsps 1 \
    -outfmt "6 qseqid sseqid slen qlen length evalue bitscore qcovs qcovhsp"
if [ $(cat ${GENE}-${ortname%.fa*}.blastout | wc -l) -eq 0 ]
then
rm ${GENE}-${ortname%.fa*}.blastout
fi
' > blastp_and_remove_no_hits.sh
chmod +x blastp_and_remove_no_hits.sh

time \
for f in $(ls x*.tmp)
do
parallel ./blastp_and_remove_no_hits.sh \
    {1} \
    {2} \
    ::: $(ls *.faa | sed 's/.faa$//g') \
    ::: $(cat $f)
done

rm *.tmp
```

## Enrichment of stress-related genes: Do we have more ortholog members for herbicide and stress-related genes in Lolium rigidum compared with the other species?
1. Extract orthologs per gene (Outputs: ${GENE}.ortho)
```{sh}
echo 'args = commandArgs(trailingOnly=TRUE)
# args = c("EPSPS")
gene = args[1]
fnames = list.files()
fnames = fnames[grep(paste0(gene, "-"), fnames)]
fnames = fnames[grep("blastout$", fnames)]
orthoout = c()
for (f in fnames){
    # f = fnames[1]
    orthogroup = unlist(strsplit(unlist(strsplit(f, "-"))[2], "[.]"))[1]
    df = read.delim(f, header=FALSE)
    if (mean(df[,ncol(df)]) >= 50){
        orthoout = c(orthoout, orthogroup)
    }
}
out = data.frame(gene=rep(gene, times=length(orthoout)), orthogroup=orthoout)
write.table(out, file=paste0(gene, ".ortho"), col.names=FALSE, row.names=FALSE, quote=FALSE, sep="\t")
' > extract_orthogroups_per_gene.R

time \
parallel Rscript extract_orthogroups_per_gene.R \
    {} ::: $(ls *.blastout | cut -d"-" -f1 | sort | uniq)

mkdir BLASTOUT/
mv *.blastout BLASTOUT/
```

2. Infer gene family expansion and contraction (Outpus: ${GENE}.conex)
```{sh}
### Input files
ORTHOUT=${DIR}/ORTHOGROUPS/orthogroups_gene_counts_families_go.out
TREE=${DIR_ORTHOGROUPS}/Species_Tree/SpeciesTree_rooted.txt

### Extract TSR and NTSR gene families (Output: ${GENE}.orthocounts)
rev ${ORTHOUT} | cut -f5- | rev > coln.tmp
awk -F'\t' '{print $(NF-1)}' ${ORTHOUT} > col1.tmp
paste -d'\t' col1.tmp coln.tmp > counts.tmp

### Infer expansion and contraction of gene families
time \
for GENE in $(ls *.ortho | sed 's/.ortho//g')
do
    # GENE=$(ls *.ortho | sed 's/.ortho//g' | head -n1)
    ### Extract gene family counts
    cat ${GENE}.ortho | cut -f2 > ${GENE}.ortho.tmp
    head -n1 counts.tmp > ${GENE}.orthocounts
    grep -f ${GENE}.ortho.tmp counts.tmp >> ${GENE}.orthocounts
    ### Run with lambda_i ~ Gamma(alpha), for each of the i_th gene family category (Output: ${GENE}_CAFE_Gamma100_results/)
    cafe5 \
        --infile ${GENE}.orthocounts \
        --tree ${TREE} \
        --n_gamma_cats 100 \
        --cores 31 \
        --pvalue 0.01 \
        --output_prefix ${GENE}_CAFE_Gamma100_results
    ### Output(Output: ${GENE}.conex)
    echo -e "Species\tExpansion\tContraction" > ${GENE}.conex
    grep -v "^#" ${GENE}_CAFE_Gamma100_results/Gamma_clade_results.txt | \
        grep -v "^<" | \
        sed 's/<..>//g' | \
        sed 's/<.>//g' >> ${GENE}.conex
done

### Clean-up
rm *.tmp
```

## dN/dS assessment: For the sress-related genes which are not more enriched (well we're including everything not just the ones that were not enriched just to be thorough), are there signs of selection?
Note: we use "gene" to refer to TSR and NTSR genes, and alignments even genes again for the genes within orthogroups within TSR/NTSR genes per species. Apologies for any misunderstandings.

1. Extract CDS per gene (i.e. all orthologs and paralogs within blast-hit orthologs) (Outputs: ${species}-${gene}-${ortho}.cds)
```{sh}
### Extract species names and number of species
head -n1 ${DIR_ORTHOGROUPS}/Orthogroups/Orthogroups.tsv | cut -f2- > \
    species_names.tmp
NSPECIES=$(awk -F"\t" "{print NF}" species_names.tmp)

### Extract gene names
echo '#!/bin/bash
f=$1
DIR_ORTHOGROUPS=$2
NSPECIES=$3
# f=$(ls *.ortho | head -n13 | tail -n1)
gene=$(echo ${f%.ortho*})
for ortho in $(cut -f2 $f)
do
    # ortho=$(cut -f2 $f | head -n1 | tail -n1)
    grep "$ortho" ${DIR_ORTHOGROUPS}/Orthogroups/Orthogroups.tsv > \
        ${gene}-${ortho}-list_gene_names.tmp
    for i in $(seq 1 $NSPECIES)
    do
        # i=1
        species=$(cut -f$i species_names.tmp | sed -z "s/\r//g" | sed -z "s/\n//g")
        col=$(echo $i + 1 | bc)
        cut -f${col} ${gene}-${ortho}-list_gene_names.tmp | \
            sed -z "s/, /\n/g" | \
            sed "s/$species|//g" > \
            ${species}-${gene}-${ortho}-list_gene_names.tmp
        if [ $(cat ${species}-${gene}-${ortho}-list_gene_names.tmp | wc -l) -lt 2 ]
        then
            rm ${species}-${gene}-${ortho}-list_gene_names.tmp
        fi
    done
    rm  ${gene}-${ortho}-list_gene_names.tmp
done
' > extract_gene_names.sh
chmod +x extract_gene_names.sh
time \
parallel ./extract_gene_names.sh \
    {} \
    ${DIR_ORTHOGROUPS} \
    ${NSPECIES} \
    ::: $(ls *.ortho)

### Extract gene sequences
echo '#!/bin/bash
f=$1
DIR=$2
species=$(echo $f | cut -d"-" -f1)
gene=$(echo $f | cut -d"-" -f2)
ortho=$(echo $f | cut -d"-" -f3)
for query in $(cat $f)
do
    julia ${DIR}/extract_sequence_using_name_query.jl \
                    ${DIR}/${species}.cds \
                    ${query} \
                    ${species}-${gene}-${ortho}-${query}.cds.tmp \
                    ${species}-${gene}-${ortho}-${query} \
                    false
done
cat ${species}-${gene}-${ortho}-*.cds.tmp > ${species}-${gene}-${ortho}.cds
rm ${species}-${gene}-${ortho}-*.cds.tmp
rm $f
' > extract_sequences_in_parallel.sh
chmod +x extract_sequences_in_parallel.sh
time \
parallel ./extract_sequences_in_parallel.sh \
    {} \
    ${DIR} \
    ::: $(ls *-list_gene_names.tmp)
```

2. Merge per orthogropup prior to alignment (Outputs: ${gene}-${ortho}.cds)
```{sh}
for ortho in $(ls *.cds | cut -d"-" -f3 | cut -d"." -f1 | sort | uniq)
do
    # ortho=$(ls *.cds | cut -d"-" -f3 | cut -d"." -f1 | sort | uniq | head -n13 | tail -n1)
    fname_output=$(ls *-${ortho}.cds | cut -d"-" -f2 | sort | uniq | sed -z "s/\n/-/g")${ortho}.cds.tmp
    cat *-${ortho}.cds > ${fname_output}
    rm *-${ortho}.cds
    mv ${fname_output} ${fname_output%.tmp*}
done
```

3. Align CDS per orthogroup per gene (Outputs: ${gene}-${ortho}.aln)
```{sh}
echo '#!/bin/bash
f=$1
ext=$2
DIR=$3
# f=Zea_mays-SOD-OG0026358.cds; ext=cds
java -Xmx8G \
    -jar ${DIR}/macse_v2.06.jar \
    -prog alignSequences \
    -seq ${f} \
    -out_NT ${f%.${ext}*}.aligned.unsorted.cds.tmp \
    -out_AA ${f%.${ext}*}.aligned.unsorted.prot.tmp
# Convert stop codons and frameshifts as "---" for compatibility with downstream tools
java -Xmx8G \
    -jar ${DIR}/macse_v2.06.jar \
    -prog exportAlignment \
    -align ${f%.${ext}*}.aligned.unsorted.cds.tmp \
    -codonForFinalStop --- \
    -codonForInternalStop NNN \
    -codonForExternalFS --- \
    -codonForInternalFS --- \
    -out_NT ${f%.${ext}*}.aln \
    -out_AA ${f%.${ext}*}.AA.prot
' > align_in_parallel.sh
chmod +x align_in_parallel.sh
time \
parallel ./align_in_parallel.sh {} cds ${DIR} \
    ::: $(ls *.cds)

rm *.tmp *.AA.prot
```

4. Remove alignments (also the corresponding cds) without Lolium rigidum genes (Outputs: ${gene}-${ortho}.aln)
```{sh}
for f in $(ls *.aln)
do
    # f=$(ls *.aln | head -n1 | tail -n1)
    x=$(grep "^>Lolium_rigidum" $f | wc -l)
    if [ $x -eq 0 ]
    then
        rm $f
        rm ${f%.aln*}.cds
    fi
done
```

5. Create pairwise cds alignments using the first Lolium rigidum alignment as the focal alignment per orthogroup per gene (Outputs: ${gene}-${ortho}.aln.pw)
```{sh}
echo '#!/bin/bash
f=$1
# f=$(ls *.aln | head -n1 | tail -n1)
focal_aln=$(grep -A1 "^>Lolium_rigidum" $f | head -n1 | sed "s/^>//g")
grep -A1 "^>Lolium_rigidum" $f | head -n2 | tail -n1 > ${f}.FOCAL_SEQ.tmp
touch ${f}.pw.tmp
### temp file without the focal alignment
sed -e "/$focal_aln/,+1d" $f > ${f}.tmp
for line in $(seq 4 2 $(cat ${f}.tmp | wc -l))
do
    # line=4
    curr_aln_name=$(head -n$(echo $line -1 | bc) ${f}.tmp | tail -n1 | sed "s/>//g")
    name=${focal_aln}--:--${curr_aln_name}
    echo $name >> ${f}.pw.tmp                           ### alignment pair name
    cat ${f}.FOCAL_SEQ.tmp >> ${f}.pw.tmp               ### focal alignment
    head -n${line} ${f}.tmp | tail -n1 >> ${f}.pw.tmp   ### current alignment
    echo "" >> ${f}.pw.tmp
done
### Clean-up
mv ${f}.pw.tmp ${f}.pw
rm ${f}.FOCAL_SEQ.tmp ${f}.tmp
### Remove single alignments
if [ $(cat ${f}.pw | wc -l) -eq 0 ]
then
    rm ${f}.pw
fi
' > prepare_pairwise_alignments_with_Lolium_rigidum_as_focus_in_parallel.sh
chmod +x prepare_pairwise_alignments_with_Lolium_rigidum_as_focus_in_parallel.sh
time \
parallel ./prepare_pairwise_alignments_with_Lolium_rigidum_as_focus_in_parallel.sh \
    {} ::: $(ls *.aln)
```

6. KaKs_calculator2 for pairwise orthogroup gene comparisons
```{sh}
time \
parallel \
KaKs_Calculator \
    -m MS \
    -i {1} \
    -o {1}.kaks.tmp \
    ::: $(ls *.aln.pw)
```

7. Find kaks file with significant (p<= 0.001) Ka/Ks > 1.0
**NOTE**: If we find high dN/dS between a pair of sequences in Lolium rigidum while the focal alignment if not that different from those of other species,then we have to change the focal alignment to that gene,and re-run KaKs_calculator!
```{sh}
echo 'args = commandArgs(trailingOnly=TRUE)
# args = c("GPX-OG0000728.aln.pw.kaks.tmp")
f = args[1]
p = 0.001
dNdS = 1.00
dat = read.delim(f, header=TRUE)
idx = (dat$P.Value.Fisher. < p) & (dat$Ka.Ks > dNdS)
if (sum(idx, na.rm=TRUE) > 0){
    print(f)
}
' > find_signs_ofsignificant_selection.R
for f in $(ls *.aln.pw.kaks.tmp)
do
    Rscript find_signs_ofsignificant_selection.R $f
done
```

## @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ ##
##               MISCELLANEOUS             ##
## @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ ##

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





### TESTING KaKs_Calucator2.0 with sliding windows
```{sh}
echo '#!/bin/bash
f=$1
julia split_alignment_pairs.jl \
    ${f} \
    60 \
    30 \
    ${f}.windows.tmp

KaKs_Calculator \
    -m MS \
    -i ${f}.windows.tmp \
    -o ${f%.aln.pw*}.kaks.tmp

Rscript plot_KaKs_across_windows.R \
    ${f%.aln.pw*}.kaks.tmp \
    0.001
' > KaKs_per_window_and_plot_in_parallel.sh
chmod +x KaKs_per_window_and_plot_in_parallel.sh
time \
parallel ./KaKs_per_window_and_plot_in_parallel.sh \
    {} ::: $(ls *.aln.pw)


```


### OR OR OR SIMPLY USE PAML::codeml on these small datasets, i.e. per TSR/NTSR gene per orthogroup
```{sh}
echo '
**************
*** INPUTS ***
**************
      seqfile = test.aln   * sequence data file name
     treefile = ../ORTHOGROUPS_SINGLE_GENE.NT.treefile * tree structure file name
**************
*** OUPUTS ***
**************
      outfile = test.codeml   * main result file
        noisy = 3                                   * 0,1,2,3: how much rubbish on the screen
      verbose = 1                                   * 1: detailed output, 0: concise output
      runmode = 0                                   * 0: user tree; 1: semi-automatic; 2: automatic; 3: StepwiseAddition; (4,5):PerturbationNNI
********************************
*** PROPERTIES OF THE INPUTS ***
********************************
        ndata = 1                                   * number of datasets
      seqtype = 1                                   * 1:codons; 2:AAs; 3:codons-->AAs 
    CodonFreq = 2                                   * 0:1/61 each, 1:F1X4, 2:F3X4, 3:codon table, 4:F1x4MG, 5:F3x4MG, 6:FMutSel0, 7:FMutSel
    aaDist = 0                                      * 0:equal, +:geometric; -:linear, 1-6:G1974,Miyata,c,p,v,a
********************************
*** CODON SUBSTITUTION MODEL ***
********************************
        model = 0                                   * 0: JC69, 1: K80 (free-ratios model to detect), 2: F81, 3: F84, 4: HKY85, 5: T92, 6: TN93, 7: GTR (REV), 8: UNREST (also see: https://github.com/ddarriba/modeltest/wiki/Models-of-Evolution)
      NSsites = 0 1 2 7 8                           * 0: M0 (one ratio), 1: M1a (neutral), 2: M2a (selection), ...
        icode = 0                                   * 0:universal code, 1:mammalian mt, ...
        Mgene = 0                                   * only for combined sequence data files, i.e. with option G in the sequence file: 0:rates, 1:separate; 2:diff pi, 3:diff k&w, 4:all diff; set as 0 if G option was not used
********************************************
*** TRANSITION / TRANSVERSION RATE RATIO ***
********************************************
    fix_kappa = 0                                   * 0: estimate kappa, 1: fix kappa, 2: kappa for branches
        kappa = 2                                   * initial or fixed kappa
*********************************************************
*** dN/dS: NONSYNONYNOUS / SYNONYNOUS VARIATION RATIO ***
*********************************************************
    fix_omega = 0                                   * 0: estimate omega, 1: fix omega
        omega = 0.5                                 * initial or fixed omega
***************************************************************************************
*** GAMMA DISTRIBUTION SHAPE PARAMETER FOR VARIABLE EVOLUTIONARY RATES ACROSS SITES ***
***************************************************************************************
    fix_alpha = 1                                   * 0: estimate alpha; 1: fix alpha
        alpha = 0.0                                   * initial or fixed alpha or is equal 0:infinity (constant rate)
       Malpha = 0                                   * 0: one alpha, 1: different alphas for genes
        ncatG = 3                                  * number of categories in the dG, AdG, or nparK models of rates
**********************
*** CLOCK SETTINGS ***
**********************
        clock = 0                                   * 0:no clock, 1:global clock; 2:local clock; 3:CombinedAnalysis
        getSE = 1
 RateAncestor = 0
*********************
*** MISCELLANEOUS ***
*********************
   Small_Diff = .1e-6
    cleandata = 1
       method = 0
  fix_blength = 0                                  * 0: ignore, -1: random, 1: initial, 2: fixed
' > test.ctl
time codeml test.ctl



```

## Evolutionary tree of stress-related genes: How did these stress-related gene which are under selection came about? 
