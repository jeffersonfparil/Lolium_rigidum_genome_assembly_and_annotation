# Comparative genomics

## Set working directory
```sh}
DIR=/data/Lolium_rigidum_ASSEMBLY/COMPARATIVE_GENOMICS
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

## Install HMMER for mapping CDS to PatherHMM gene family models
```{sh}
sudo apt install hmmer
```

## Install MACSE: Multiple Alignment of Coding SEquences Accounting for Frameshifts and Stop Codons
```{sh}
wget https://bioweb.supagro.inra.fr/macse/releases/macse_v2.06.jar
java -Xmx250G -jar macse_v2.06.jar -h
```

## Install Gblocks to remove poorly aligned sequences
```{sh}
wget http://molevol.cmima.csic.es/castresana/Gblocks/Gblocks_Linux64_0.91b.tar.Z
tar -xvzf Gblocks_Linux64_0.91b.tar.Z
cd Gblocks_0.91b/
PATH=${PATH}:$(pwd)
cd -
```

## Install RaxML-ng for building trees
```{sh}
wget https://github.com/amkozlov/raxml-ng/releases/download/1.1.0/raxml-ng_v1.1.0_linux_x86_64.zip
unzip raxml-ng_v1.1.0_linux_x86_64.zip -d raxml-ng_v1.1.0
cd raxml-ng_v1.1.0/
PATH=${PATH}:$(pwd)
cd -
```

## Install Clann for merging trees
```{sh}
wget https://github.com/ChrisCreevey/clann/archive/refs/tags/v4.2.4.tar.gz
tar -xvzf v4.2.4.tar.gz
rm v4.2.4.tar.gz
cd clann-4.2.4/
./configure
make
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

## Classify predicted proteins by gene families using PatherHMM database and HMMER
```{sh}
### Define the location of the 15,619 protein family HMMs
DIR_PANTHER=${DIR}/PatherHMM_17.0/famlib/rel/PANTHER17.0_altVersion/hmmscoring/PANTHER17.0/books
GOT_PATHER=${DIR}/PatherHMM_17.0/PANTHER17.0_HMM_classifications

### Prepare parallelisable HMMER search script
echo '#!/bin/bash
PROTFA=$1
DIR_PANTHER=$2
d=$3
HMMODL=${DIR_PANTHER}/${d}/hmmer.hmm
OUTEMP=${PROTFA}-hhmer_gene_family_hits-${d}.tmp
hmmsearch -E 0.0001 --tblout ${OUTEMP} ${HMMODL} ${PROTFA}
' > hmmsearch_for_parallel_execution.sh
chmod +x hmmsearch_for_parallel_execution.sh

### Iteratively, for each genome's predicted protein sequences run hmmsearch in paralel for each PatherHMM protein family
for PROTFA in $(ls *.faa)
do
    # PROTFA=Lolium_rigidum.faa
    echo ${PROTFA}
    time \
    parallel \
    ./hmmsearch_for_parallel_execution.sh \
        ${PROTFA} \
        ${DIR_PANTHER} \
        {} ::: $(ls $DIR_PANTHER)
    ### Concatenate hmmsearch output for each protein family into a single output file
    CONCAT=${PROTFA}-hhmer_gene_family_hits.txt
    touch $CONCAT
    for f in $(ls ${PROTFA}-hhmer_gene_family_hits-*)
    do
        sed "/^#/d" ${f} | awk '{print $1,$3,$5}' >> ${CONCAT}
        rm $f
    done
done

### Find the best fitting protein family (lowest E-value) for each predicted protein
echo 'using CSV, DataFrames, ProgressMeter
fname_input = ARGS[1]
fname_codes = ARGS[2]
fname_output = string(split(basename(fname_input), "."[1])[1], ".pthr")
if dirname(fname_input) != ""
    fname_output = string(dirname(fname_input), "/", fname_output)
end
# fname_input = "Arabidopsis_thaliana.faa-hhmer_gene_family_hits.txt"; fname_codes = "PatherHMM_17.0/Panther17.0_HMM_familyIDs.txt"
dat = CSV.read(open(fname_input, "r"), DataFrames.DataFrame, header=false)
### Find most likely gene family per predicted gene, i.e. first row has the best E-value
gene_names = []
panther_codes = []
@showprogress for gene in levels(dat.Column1)
    # gene = levels(dat.Column1)[2]
    subset = dat[dat.Column1 .== gene, :]
    sort!(subset, :Column3, rev=false)
    push!(gene_names, subset[1,1])
    push!(panther_codes, split(subset[1,2], "."[1])[1])
end
### Identify PatherHMM gene families
gene_families = []
cod = CSV.read(open(fname_codes), DataFrames.DataFrame, header=false)
@showprogress for code in panther_codes
    # code = panther_codes[1]
    push!(gene_families, cod.Column2[code .== cod.Column1][1])
end
### Save gene names and their corresponding PantherHMM protein family classifications
out = DataFrames.DataFrame(GENE_NAME=gene_names, PANTHER_CODE=panther_codes, GENE_FAMILY=gene_families)
CSV.write(open(fname_output, "w"), out, delim="\t"[1], header=false)
' > extrac_PantherHMM_protein_families_per_predicted_protein.jl

### Parallel execution to find best fitting protein family
PANTHER_CODES=${DIR}/PatherHMM_17.0/Panther17.0_HMM_familyIDs.txt
time \
parallel \
julia extrac_PantherHMM_protein_families_per_predicted_protein.jl \
    {} \
    ${PANTHER_CODES} ::: $(ls *.faa-hhmer_gene_family_hits.txt)
```

## Estimate divergence between species times using MCMCTREE and TimeTree.org fossil record estimates

Prepare Rscript to find single-copy gene families:
```{find_single_copy_gene_families.R}
args = commandArgs(trailingOnly=TRUE)
dat = read.table(args[1], header=FALSE)
frq = as.data.frame(table(dat$V1))
write.table(frq$Var1[frq$Freq==1], file=args[1], row.names=FALSE, col.names=FALSE, quote=FALSE)
```

Prepare Rscript to find the common single-copy gene families across the genomes:
```{find_common_single_copy_gene_families_across_genomes.R}
args = commandArgs(trailingOnly=TRUE)
fnames_input = args[1:(length(args)-1)]
fname_output = args[length(args)]
for (f in fnames_input){
    dat = read.table(f, header=FALSE)
    if (exists("MERGED")==FALSE){
        MERGED = dat
    } else {
        MERGED = merge(MERGED, dat, by="V1")
    }
}
write.table(MERGED, file=fname_output, row.names=FALSE, col.names=FALSE, quote=FALSE)
```

Prepare a tree plotting Rscript:
```{draw_tree.R}
args = commandArgs(trailingOnly=TRUE)
myTree = ape::read.tree(args[1])
svg(args[2], width=20, height=5)
plot(myTree)
dev.off()
```

Create a simple sequence extractor in julia:
```{extract_sequence_using_name_query.jl}
using ProgressMeter
fasta_input = ARGS[1]
sequence_name_query = ARGS[2]
# fasta_input = "Secale_cereale.cds"
# sequence_name_query = "Lolium_rigidum_rna-XM_047205051.1_R0"
# fasta_output = ""
fasta_output = try 
                    ARGS[3]
                catch
                    ""
                end
new_sequence_name = try 
                    ARGS[4]
                catch
                    ""
                end
add_gene_coordinates = try 
                    parse(Bool, ARGS[5])
                catch
                    true
                end
if fasta_output == ""
    fasta_output = "temp.fasta"
end
file_input = open(fasta_input, "r")
seekend(file_input); n = position(file_input)
seekstart(file_input)
pb = Progress(n)
while !eof(file_input)
    line = readline(file_input)
    if line[1] == '>'
        while match(Regex(sequence_name_query), line) != nothing
            @show line
            file_output = open(fasta_output, "a")
            if (new_sequence_name != "")
                vec_line = split(line, " ")
                line = string(">", new_sequence_name)
            end
            write(file_output, string(line, '\n'))
            line = readline(file_input)
            while line[1] != '>'
                write(file_output, line)
                line = readline(file_input)
                update!(pb, position(file_input))
            end
            write(file_output, '\n')
            close(file_output)
        end
    end
    update!(pb, position(file_input))
end
close(file_input)
```

Prepare julia script to convert fasta into phylip format:
```{fasta_to_phylip.jl}
filename_input = ARGS[1]
filename_output = string(join(split(filename_input, '.')[1:(end-1)], '.'), ".phylip")
# count number of sequences
function COUNT_SEQUENCES(filename_input)
    FILE = open(filename_input, "r")
    count_sequences = 0
    while !eof(FILE)
        if readline(FILE)[1] == '>'
        count_sequences += 1
        end
    end
    close(FILE)
    return(count_sequences)
end
count_sequences = COUNT_SEQUENCES(filename_input)
# count sequence length
function SEQUENCE_LENGTH(filename_input)
    FILE = open(filename_input, "r")
    _ = readline(FILE)
    line = readline(FILE)
    sequence_length = length(line)
    while line[1] != '>'
        line = readline(FILE)
        sequence_length = sequence_length + length(line)
    end
    sequence_length = sequence_length - length(filename_input = ARGS[1]
filename_output = string(join(split(filename_input, '.')[1:(end-1)], '.'), ".phylip")
# count number of sequences
function COUNT_SEQUENCES(filename_input)
    FILE = open(filename_input, "r")
    count_sequences = 0
    while !eof(FILE)
        if readline(FILE)[1] == '>'
        count_sequences += 1
        end
    end
    close(FILE)
    return(count_sequences)
end
count_sequences = COUNT_SEQUENCES(filename_input)
# count sequence length
function SEQUENCE_LENGTH(filename_input)
    FILE = open(filename_input, "r")
    _ = readline(FILE)
    line = readline(FILE)
    sequence_length = length(line)
    while line[1] != '>'
        line = readline(FILE)
        sequence_length = sequence_length + length(line)
    end
    sequence_length = sequence_length - length(line)
    close(FILE)
    return(sequence_length)
end
sequence_length = SEQUENCE_LENGTH(filename_input)
# output phylip format
function CONVERT_FASTA_TO_PHYLIP(filename_input, filename_output, count_sequences, sequence_length)
    FILE_INPUT = open(filename_input, "r")
    FILE_OUTPUT = open(filename_output, "a")
    write(FILE_OUTPUT, string(count_sequences, " ", sequence_length, '\n'))
    while !eof(FILE_INPUT)
        line = string(readline(FILE_INPUT), '\n')
        if line[1] == '>'
            line = string(line[2:(end-1)], "  ")
        end
        write(FILE_OUTPUT, line)
    end
    close(FILE_INPUT)
    close(FILE_OUTPUT)
    return(0)
end
CONVERT_FASTA_TO_PHYLIP(filename_input, filename_output, count_sequences, sequence_length)

    write(FILE_OUTPUT, string(count_sequences, " ", sequence_length, '\n'))
    while !eof(FILE_INPUT)
        line = string(readline(FILE_INPUT), '\n')
        if line[1] == '>'
            line = string(line[2:(end-1)], "  ")
        end
        write(FILE_OUTPUT, line)
    end
    close(FILE_INPUT)
    close(FILE_OUTPUT)
    return(0)
end
CONVERT_FASTA_TO_PHYLIP(filename_input, filename_output, count_sequences, sequence_length)
```

Identify single-copy gene families:
**Note:** The genes for each protein family can be different for each species. We are assuming that the gene in each of these single-copy protein families represent its corresponding family solely and completely.
```{sh}
for PTHR in $(ls *.pthr)
do
    cut -f2 ${PTHR} > ${PTHR%.pthr*}.tmp
    Rscript find_single_copy_gene_families.R \
        ${PTHR%.pthr*}.tmp
done
Rscript find_common_single_copy_gene_families_across_genomes.R \
    $(ls *.tmp) \
    SINGLE_COPY_GENE_FAMILIES.pthr
rm *.tmp
```

Extract gene names belonging to the single-copy protein families:
```{sh}
time \
for SPECIES in $(ls *.pthr | grep -v "SINGLE_COPY_GENE_FAMILIES" | sed 's/.pthr//g')
do
    for FAMILY in $(cat SINGLE_COPY_GENE_FAMILIES.pthr)
    do
        grep ${FAMILY} ${SPECIES}.pthr | cut -f1 >> ${SPECIES}-SINGLE_COPY_GENE_NAMES.tmp
    done
done
```

Extract the CDS sequences from the predicted protein sequences using the gene names:
```{sh}
### Create parallelisable bash script
echo '#!/bin/bash
CDS=$1
QUERY=$2
SPECIES=$(basename $CDS); SPECIES=${SPECIES%.cds*}
PANTHER_GENE_FAMILY=$(grep $QUERY ${SPECIES}.pthr | cut -f2)
julia \
extract_sequence_using_name_query.jl \
    ${CDS} \
    ${QUERY} \
    ${CDS}-${QUERY}-SINGLE_COPY_SEQUENCE.tmp \
    ${SPECIES}-${PANTHER_GENE_FAMILY} \
    false
' > extract_sequence_using_name_query_PARALLEL.sh
chmod +x extract_sequence_using_name_query_PARALLEL.sh

### Extract the CDS of the single-copy protein family genes
time \
for SPECIES in $(ls *.cds | grep -v "SINGLE_COPY_GENE_FAMILIES" | sed 's/.cds//g')
do
    parallel ./extract_sequence_using_name_query_PARALLEL.sh \
            ${SPECIES}.cds \
            {} ::: $(cat ${SPECIES}-SINGLE_COPY_GENE_NAMES.tmp)
done
cat *-SINGLE_COPY_SEQUENCE.tmp > SINGLE_COPY_GENE_FAMILIES.cds
rm *.tmp
```

Use `MACSE` to align the CDS into codons and protein sequences, build the trees for each ortholog with `RaxML-ng`, convert fasta alignemnts into phylip format for `PAML`, and merge all trees into a single multi-tree file `temp_ALL_TREES.trees`:
```{sh}
### Generate parallelisable MACSE alignement script
### NOTES: Set the final stop codons as "---", and internal stop codons as "NNN" so that PAML programs won't ask you to press enter to continue; also set frameshifts from "!" into "-"
###        Also sort the alignment sequences, because MACSE jumbles them u for some reason with no option to return input order
echo '#!/bin/bash
ORTHOLOG=$1
# Extract CDS of the gene
cat SINGLE_COPY_GENE_FAMILIES.cds | \
    grep -A1 ${ORTHOLOG} | \
    sed "/^--$/d" | \
    sed "s/-${ORTHOLOG}//g" > ${ORTHOLOG}.cds.tmp
# Align the CDS across species
java -Xmx8G \
     -jar macse_v2.06.jar \
     -prog alignSequences \
     -seq ${ORTHOLOG}.cds.tmp \
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
     -out_NT temp_${ORTHOLOG}.NT.cds \
     -out_AA temp_${ORTHOLOG}.AA.prot
# Clean-up
rm ${ORTHOLOG}.cds.tmp
rm ${ORTHOLOG}.aligned.unsorted.cds.tmp
rm ${ORTHOLOG}.aligned.unsorted.prot.tmp
# Build tree
raxml-ng --all \
         --msa temp_${ORTHOLOG}.NT.cds \
         --model GTR+G8+F \
         --outgroup Arabidopsis_thaliana \
         --prefix temp_${ORTHOLOG} \
         --threads 1 \
         --seed 42069 \
         --tree pars{25},rand{25}
# Convert codon and amino acid sequences from fasta into phylip format
julia fasta_to_phylip.jl \
      temp_${ORTHOLOG}.NT.cds ### output: temp_${ORTHOLOG}.NT.phylip
julia fasta_to_phylip.jl \
      temp_${ORTHOLOG}.AA.prot ### output: temp_${ORTHOLOG}.AA.phylip
# Sort the sequences alphabetically by species
mv temp_${ORTHOLOG}.NT.phylip temp_${ORTHOLOG}.NT.phylip.tmp
mv temp_${ORTHOLOG}.AA.phylip temp_${ORTHOLOG}.AA.phylip.tmp
sort temp_${ORTHOLOG}.NT.phylip.tmp > temp_${ORTHOLOG}.NT.phylip
sort temp_${ORTHOLOG}.AA.phylip.tmp > temp_${ORTHOLOG}.AA.phylip
' > macse_RaxML-ng_PARALLEL.sh
chmod +x macse_RaxML-ng_PARALLEL.sh

### Parallel execution
time \
parallel \
./macse_RaxML-ng_PARALLEL.sh {} ::: $(cat SINGLE_COPY_GENE_FAMILIES.pthr)

### Merge phylip alignements and tree files
cat temp_*.NT.phylip > SINGLE_COPY_GENE_FAMILIES.phylip
cat temp_*.bestTree  > SINGLE_COPY_GENE_FAMILIES.trees
rm temp_*
```

Add number of species and trees into `SINGLE_COPY_GENE_FAMILIES.trees`:
```{sh}
mv SINGLE_COPY_GENE_FAMILIES.trees SINGLE_COPY_GENE_FAMILIES.trees.bk
echo "7 $(cat SINGLE_COPY_GENE_FAMILIES.trees.bk | wc -l)" > SINGLE_COPY_GENE_FAMILIES.trees
cat SINGLE_COPY_GENE_FAMILIES.trees.bk >> SINGLE_COPY_GENE_FAMILIES.trees
rm SINGLE_COPY_GENE_FAMILIES.trees.bk
```

<!-- Using `timetree.org` tree:
```{sh}
echo "7 1
(Marchantia_polymorpha:532.29191200,
    (
        (
            (Oryza_sativa:49.60000000,
                (Secale_cereale:24.10373857,
                    (Lolium_rigidum:1.65000000, Lolium_perenne:1.65000000) #1 'A': 22.45373857 '@0.0165'
                ) #2 'B': 25.49626143 '@0.241'
            ) #3 'C': 0.00000000 '@0.5', Zea_mays:49.60000000
        ) #4 'D':110.87712725 '@0.5', Arabidopsis_thaliana:160.47712725
    ) #5 'E':371.81478475 '@1.6'
)'@5.32';
" > SINGLE_COPY_GENE_FAMILIES.tree
``` -->



Estimate the parammeters for the gamma distribution prior for the substitution rates in gene using the CODEML control file `FIND_GENE_SUBS_RATE.ctl`:
```{sh}
echo '
**************
*** INPUTS ***
**************
      seqfile = SINGLE_COPY_GENE_FAMILIES.phylip    * sequence data file name
     treefile = SINGLE_COPY_GENE_FAMILIES.trees     * tree structure file name
**************
*** OUPUTS ***
**************
      outfile = SINGLE_COPY_GENE_FAMILIES.codeml    * main result file
        noisy = 0                                   * 0,1,2,3: how much rubbish on the screen
      verbose = 0                                   * 1: detailed output, 0: concise output
      runmode = 0                                   * 0: user tree; 1: semi-automatic; 2: automatic; 3: StepwiseAddition; (4,5):PerturbationNNI
********************************
*** PROPERTIES OF THE INPUTS ***
********************************
        ndata = 222                                 * number of datasets
      seqtype = 1                                   * 1:codons; 2:AAs; 3:codons-->AAs 
    CodonFreq = 2                                   * 0:1/61 each, 1:F1X4, 2:F3X4, 3:codon table, 4:F1x4MG, 5:F3x4MG, 6:FMutSel0, 7:FMutSel
********************************
*** CODON SUBSTITUTION MODEL ***
********************************
        model = 0                                   * 0: JC69, 1: K80, 2: F81, 3: F84, 4: HKY85, 5: T92, 6: TN93, 7: GTR (REV), 8: UNREST (also see: https://github.com/ddarriba/modeltest/wiki/Models-of-Evolution)
      NSsites = 0                                   * 0: M0 (one ratio), 1: M1a (neutral), 2: M2a (selection), ...
        icode = 0                                   * 0:universal code, 1:mammalian mt, ...
        Mgene = 0                                   * only for combined sequence data files, i.e. with option G in the sequence file: 0:rates, 1:separate; 2:diff pi, 3:diff k&w, 4:all diff; set as 0 if G option was not used
********************************************
*** TRANSITION / TRANSVERSION RATE RATIO ***
********************************************
    fix_kappa = 0                                   * 0: estimate kappa, 1: fix kappa, 2: kappa for branches
        kappa = 1.6                                 * initial or fixed kappa
*********************************************************
*** dN/dS: NONSYNONYNOUS / SYNONYNOUS VARIATION RATIO ***
*********************************************************
    fix_omega = 0                                   * 0: estimate omega, 1: fix omega
        omega = .9                                  * initial or fixed omega
******************************************
*** GAMMA DISTRIBUTION SHAPE PARAMETER ***
******************************************
    fix_alpha = 0                                   * 0: estimate alpha; 1: fix alpha
        alpha = 1                                   * initial or fixed alpha or is equal 0:infinity (constant rate)
       Malpha = 0                                   * 0: one alpha, 1: different alphas for genes
        ncatG = 10                                  * number of categories in the dG, AdG, or nparK models of rates
**********************
*** CLOCK SETTINGS ***
**********************
        clock = 1                                   * 0:no clock, 1:global clock; 2:local clock; 3:CombinedAnalysis
        getSE = 0
 RateAncestor = 0
*********************
*** MISCELLANEOUS ***
*********************
   Small_Diff = .1e-6
    cleandata = 1
       method = 0
  fix_blength = 0                                  * 0: ignore, -1: random, 1: initial, 2: fixed
' > SINGLE_COPY_GENE_FAMILIES-CODEML.ctl

time codeml SINGLE_COPY_GENE_FAMILIES-CODEML.ctl


```



Prepare MCMCTree control or script file: `SINGLE_COPY_GENE_FAMILIES.ctl`:
```{MCMCTree-ctl}
          seed = 42069
       seqfile = SINGLE_COPY_GENE_FAMILIES.phylip
      treefile = SINGLE_COPY_GENE_FAMILIES.tree
      mcmcfile = SINGLE_COPY_GENE_FAMILIES-MCMCTREE.mcmc
       outfile = SINGLE_COPY_GENE_FAMILIES-MCMCTREE.out

         ndata = 222
       seqtype = 2    * 0: nucleotides; 1:codons; 2:AAs
       usedata = 2    * 0: no data; 1:seq like; 2:normal approximation; 3:out.BV (in.BV)
         clock = 1    * 1: global clock; 2: independent rates; 3: correlated rates
       RootAge = '<5.33'  * safe constraint on root age, used if no fossil for root.
    
       runmode = 2
         model = 0    * 0:JC69, 1:K80, 2:F81, 3:F84, 4:HKY85
         alpha = 0    * alpha for gamma rates at sites
         ncatG = 5    * No. categories in discrete gamma

     cleandata = 0    * remove sites with ambiguity data (1:yes, 0:no)?

       BDparas = 2 2 0.1  * birth, death, sampling
   kappa_gamma = 6 2      * gamma prior for kappa
   alpha_gamma = 1 1      * gamma prior for alpha

   rgene_gamma = 2 20 1   * gammaDir prior for rate for genes
  sigma2_gamma = 1 10 1   * gammaDir prior for sigma^2     (for clock=2 or 3)

      finetune = 1: .1 .1 .1 .1 .1 .1 * auto (0 or 1): times, musigma2, rates, mixing, paras, FossilErr

         print = 1   * 0: no mcmc sample; 1: everything except branch rates 2: everything
        burnin = 2000
      sampfreq = 10
       nsample = 20000
```

Run `MCMCTree`:
```{sh}
time \
mcmctree SINGLE_COPY_GENE_FAMILIES.ctl
```







## Herbicide resistance genes
echo -e "CYTOCHROME P450\tNTSR
GLUTATHIONE S-TRANSFERASE\tNTSR
ACETYL-COENZYME A CARBOXYLASE CARBOXYL TRANSFERASE\tACCase
ACETOLACTATE SYNTHASE\tAcetolactate synthase
AROM/DEHYDROQUINATE SYNTHASE\tEPSP
FATTY ACID SYNTHASE SUBUNIT BETA\tFatty acid biosynthesis
" > TSR_NTSR_PATHERHMM.list
for GENE in $(cut -f1 TSR_NTSR_PATHERHMM.list)
do
    echo $GENE
    for SPECIES in $(ls *.pthr | grep -v "SINGLE_COPY" | sed 's/.pthr//g')
    do
        COUNT=$(grep '$GENE' ${SPECIES}.pthr | wc -l)
        echo ${SPECIES}-${GENE}-${COUNT}
    done
done
