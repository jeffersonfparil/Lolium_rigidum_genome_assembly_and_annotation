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

## Install RaxML-ng
```{sh}
wget https://github.com/amkozlov/raxml-ng/releases/download/1.1.0/raxml-ng_v1.1.0_linux_x86_64.zip
unzip raxml-ng_v1.1.0_linux_x86_64.zip -d raxml-ng_v1.1.0
cd raxml-ng_v1.1.0/
PATH=${PATH}:$(pwd)
cd -
```

## Install Clann for merging trees generated using each ortholog into a single supertree
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

## Install the ape R package:
```{R}
install.packages("ape")
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
        rm $f
    done
    sort ${OUTPUT} > ${OUTPUT}.tmp
    mv ${OUTPUT}.tmp ${OUTPUT}
done
```

## Generate a new annotation file with PantherHMM-derived gene families

Save the following Julia script as `generate_final_gff_annotation_file.jl`:
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

Generate the final annotation files (`FINAL_ANNOTATION_PATHERHMM_GENE_FAMILIES.gff`) in parallel:
```{sh}
cd $DIR
PANTHER_CODES=${DIR}/PatherHMM_17.0/Panther17.0_HMM_familyIDs.txt
time \
parallel --link \
julia generate_final_gff_annotation_file.jl {1} ${PANTHER_CODES} {2} ::: \
    $(ls */GeMoMa_output_*/hhmer_gene_family_hits.txt) ::: \
    $(ls */GeMoMa_output_*/final_annotation.gff)
```

## Estimate divergence between species times using MCMCTREE and TimeTree.org fossil record estimates

Prepare Rscript to find single-copy gene families `find_single_copy_gene_families.R`:
```{R}
args = commandArgs(trailingOnly=TRUE)
dat = read.table(args[1], header=FALSE)
frq = as.data.frame(table(dat$V1))
frq = frq[order(frq$Var1), ]
write.table(frq$Var1[frq$Freq==1], file=args[1], row.names=FALSE, col.names=FALSE, quote=FALSE)
```

Prepare Rscript to find the common single-copy gene families across the genomes `find_common_single_copy_gene_families_across_genomes.R`:
```{R}
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

Using `R::ape` prepare the tree plotting Rscript `draw_tree.R`:
```{R}
args = commandArgs(trailingOnly=TRUE)
myTree = ape::read.tree(args[1])
svg(args[2], width=20, height=5)
plot(myTree)
dev.off()
```

Create a fast and simple sequence extractor in julia: `extract_sequence_using_name_query.jl`:
```{julia}
using ProgressMeter
fasta_input = ARGS[1]
sequence_name_query = ARGS[2]
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
# fasta_input = "/data-weedomics-3/Arabidopsis_thaliana/GeMoMa_output_Arabidopsis_thaliana/predicted_proteins.fasta"
# sequence_name_query = "Arabidopsis_thaliana_gene-AT1G01010"
# fasta_input = "/data-weedomics-3/Oryza_sativa/GeMoMa_output_Oryza_sativa/predicted_proteins.fasta"
# sequence_name_query = "Oryza_sativa_rna-XM_015766610.2_R0"
# fasta_output = ""
if fasta_output == ""
    fasta_output = "/data-weedomics-3/temp.fasta"
end
file_input = open(fasta_input, "r")
seekend(file_input); n = position(file_input)
seekstart(file_input)
pb = Progress(n)
while !eof(file_input)
    line = readline(file_input)
    if line[1] == '>'
        while match(Regex(sequence_name_query), line) != nothing
            file_output = open(fasta_output, "a")
            if (new_sequence_name != "")
                vec_line = split(line, " ")
                chr = vec_line[match.(Regex("chr"), vec_line) .!= nothing][1]
                interval = vec_line[match.(Regex("interval"), vec_line) .!= nothing][1]
                if add_gene_coordinates
                    line = string(">", new_sequence_name, "-", chr, "-", interval)
                else
                    line = string(">", new_sequence_name)
                end
            end
            write(file_output, string(line, '\n'))
            line = readline(file_input)
            while line[1] != '>'
                write(file_output, string(line, '\n'))
                line = readline(file_input)
                update!(pb, position(file_input))
            end
            close(file_output)
        end
    end
    update!(pb, position(file_input))
end
close(file_input)
```

Prepare julia script to convert fasta into phylip format `fasta_to_phylip.jl`:
```{julia}
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
```

Identify single-copy gene families:
```{sh}
cd $DIR
REF="Arabidopsis_thaliana"
# REF="Oryza_sativa" ### Which reference genome to use?
# wc -l */GeMoMa_output_*/FINAL*.gff ### Arabidopsis thaliana generates consistently less mapped annotations than rice.
for GENOME in Lolium_rigidum Lolium_perenne Arabidopsis_thaliana Oryza_sativa Zea_mays Secale_cereale Marchantia_polymorpha
do
    sed '/^#/d' ${GENOME}/GeMoMa_output_${REF}/FINAL_ANNOTATION_PATHERHMM_GENE_FAMILIES.gff | \
        cut -f9 | \
        cut -d';' -f2 > temp_gene_families-${REF}-${GENOME}.txt
    Rscript find_single_copy_gene_families.R temp_gene_families-${REF}-${GENOME}.txt
done
Rscript find_common_single_copy_gene_families_across_genomes.R $(ls temp_*.txt) SINGLE_COPY_GENE_FAMILIES-${REF}.txt
rm temp_*.txt
```

Extract gene names belonging to the single-copy gene families from the annotation files:
```{sh}
for GENOME in Lolium_rigidum Lolium_perenne Arabidopsis_thaliana Oryza_sativa Zea_mays Secale_cereale Marchantia_polymorpha
do
    for gene in $(cat SINGLE_COPY_GENE_FAMILIES-${REF}.txt)
    do
        grep ${gene} ${GENOME}/GeMoMa_output_${REF}/FINAL_ANNOTATION_PATHERHMM_GENE_FAMILIES.gff | \
                sed '/^#/d' | \
                cut -f9 | \
                cut -d';' -f1 | \
                sed 's/ID=//g' >> temp_gene_names-${GENOME}-${REF}.txt
    done
done
```

Extract the CDS sequences from the predicted protein sequences using the gene names:
```{sh}
# Create parallelisable bash script
echo '#!/bin/bash
DIR=$1
GENOME=$2
QUERY=$3
PANTHER_GENE_FAMILY=$4
REF=$5
julia \
extract_sequence_using_name_query.jl \
    ${DIR}/${GENOME}/GeMoMa_output_${REF}/predicted_cds.fasta \
    ${QUERY} \
    $(echo temp_${GENOME}-GENE_${QUERY} | sed -z "s/ /_/g").fasta \
    ${GENOME}-${PANTHER_GENE_FAMILY} \
    false
' > extract_sequence_using_name_query_PARALLEL.sh
chmod +x extract_sequence_using_name_query_PARALLEL.sh

### Extract sequences for each species in parallel
for GENOME in Lolium_rigidum Lolium_perenne Arabidopsis_thaliana Oryza_sativa Zea_mays Secale_cereale Marchantia_polymorpha
do
    parallel --link \
    ./extract_sequence_using_name_query_PARALLEL.sh \
        $DIR \
        ${GENOME} \
        {1} \
        {2} \
        ${REF} \
        ::: $(cat temp_gene_names-${GENOME}-${REF}.txt) \
        ::: $(cat SINGLE_COPY_GENE_FAMILIES-${REF}.txt)
    cat temp_${GENOME}-GENE_*.fasta > ${GENOME}-SINGLE_COPY_GENE_FAMILY_Athaliana_genes.fasta
    rm temp_${GENOME}-GENE_*.fasta
done
rm temp_gene_names-*-${REF}.txt
```

Use `MACSE` to align the CDS into codons and protein sequences, build the trees for each ortholog with `RaxML-ng`, convert fasta alignemnts into phylip format for `PAML`, and merge all trees into a single multi-tree file `temp_ALL_TREES.trees`:
```{sh}
### Generate parallelisable MACSE alignement script
### NOTES: Set the final stop codons as "---", and internal stop codons as "NNN" so that PAML programs won't ask you to press enter to continue; also set frameshifts from "!" into "-"
###        Also sort the alignment sequences, because MACSE jumbles them u for some reason with no option to return input order
echo '#!/bin/bash
ORTHOLOG=$1
cat *-SINGLE_COPY_GENE_FAMILY_Athaliana_genes.fasta | \
    grep -A1 ${ORTHOLOG} | \
    sed "/^--$/d" | \
    sed "s/-${ORTHOLOG}//g" > temp_${ORTHOLOG}.fasta
java -Xmx8G \
     -jar macse_v2.06.jar \
     -prog alignSequences \
     -seq temp_${ORTHOLOG}.fasta \
     -out_NT temp_${ORTHOLOG}-ALIGNED-NT.tmp \
     -out_AA temp_${ORTHOLOG}-ALIGNED-AA.tmp
java -Xmx8G \
     -jar macse_v2.06.jar \
     -prog exportAlignment \
     -align temp_${ORTHOLOG}-ALIGNED-NT.tmp \
     -codonForFinalStop --- \
     -codonForInternalStop NNN \
     -codonForExternalFS --- \
     -codonForInternalFS --- \
     -out_NT temp_${ORTHOLOG}-ALIGNED-NT.fasta \
     -out_AA temp_${ORTHOLOG}-ALIGNED-AA.fasta
raxml-ng --all \
         --msa temp_${ORTHOLOG}-ALIGNED-NT.fasta \
         --model GTR+G8+F \
         --outgroup Marchantia_polymorpha \
         --prefix temp_${ORTHOLOG} \
         --threads 1 \
         --seed 42069 \
         --tree pars{25},rand{25}
julia fasta_to_phylip.jl \
      temp_${ORTHOLOG}-ALIGNED-NT.fasta ### outputs: temp_${ORTHOLOG}-ALIGNED-NT.phylip
julia fasta_to_phylip.jl \
      temp_${ORTHOLOG}-ALIGNED-AA.fasta ### outputs: temp_${ORTHOLOG}-ALIGNED-AA.phylip
mv temp_${ORTHOLOG}-ALIGNED-NT.phylip temp_${ORTHOLOG}-ALIGNED-NT.phylip.tmp
mv temp_${ORTHOLOG}-ALIGNED-AA.phylip temp_${ORTHOLOG}-ALIGNED-AA.phylip.tmp
sort temp_${ORTHOLOG}-ALIGNED-NT.phylip.tmp > temp_${ORTHOLOG}-ALIGNED-NT.phylip
sort temp_${ORTHOLOG}-ALIGNED-AA.phylip.tmp > temp_${ORTHOLOG}-ALIGNED-AA.phylip
' > macse_RaxML-ng_PARALLEL.sh
chmod +x macse_RaxML-ng_PARALLEL.sh

### Parallel execution
time \
parallel \
./macse_RaxML-ng_PARALLEL.sh {} ::: $(cat SINGLE_COPY_GENE_FAMILIES-${REF}.txt)

### Merge phylip alignements and tree files
cat temp_*-ALIGNED-NT.phylip > SINGLE_COPY_GENE_FAMILIES.phylip
cat temp_*.bestTree       > SINGLE_COPY_GENE_FAMILIES.trees
rm temp_*.*
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




## Align TSR and NTSR gene families across species with MAFFT


Extract protein sequences using `extract_sequence_using_name_query.jl`:
```{sh}
# Create parallelisable bash script
echo '#!/bin/bash
DIR=$1
GENOME=$2
QUERY=$3
GENE_FAMILY_REGEX=$4
REF=$5
julia \
extract_sequence_using_name_query.jl \
    ${DIR}/${GENOME}/GeMoMa_output_${REF}/predicted_proteins.fasta \
    ${QUERY} \
    $(echo temp_${GENOME}-${GENE_FAMILY_REGEX}-${QUERY} | sed -z "s/ /_/g").fasta \
    $(echo ${GENOME}-${GENE_FAMILY_REGEX}-${QUERY} | sed -z "s/ /_/g")
' > extract_sequence_using_name_query_PARALLEL.sh
chmod +x extract_sequence_using_name_query_PARALLEL.sh

# Execute in parallel and merge
REF="Arabidopsis_thaliana"
# REF="Oryza_sativa"
GENE_FAMILY_REGEX="ACETOLACTATE SYNTHASE"
GENE_FAMILY=$(echo ${GENE_FAMILY_REGEX} | sed -z 's/ /_/g')
OUTPUT=${GENE_FAMILY}-GeMoMa_${REF}_genes-PatherHMM-predicted.fasta
time \
parallel \
./extract_sequence_using_name_query_PARALLEL.sh \
    ${DIR} \
    {1} \
    {2} \
    ${GENE_FAMILY} \
    ${REF} \
    ::: Lolium_rigidum Lolium_perenne Arabidopsis_thaliana Oryza_sativa Zea_mays Secale_cereale Marchantia_polymorpha \
    ::: $(grep "$GENE_FAMILY_REGEX" */GeMoMa_output_${REF}/FINAL*.gff  | cut -f9 | cut -d';' -f1 | sed 's/ID=//g' | sort | uniq)

cat temp_*.fasta > ${OUTPUT}
rm temp_*.fasta
```

Align sequences with `MAFFT`, build the phylogenetic tree with `FastTreeDbl`, and plot the tree:
```{sh}
mafft --maxiterate 1000 --localpair ${OUTPUT} > ${OUTPUT%.fasta*}-ALIGNED.fasta
FastTreeDbl ${OUTPUT%.fasta*}-ALIGNED.fasta > ${OUTPUT%.fasta*}-ALIGNED.tree
Rscript draw_tree.R \
    ${OUTPUT%.fasta*}-ALIGNED.tree \
    ${OUTPUT%.fasta*}-TREE.svg
```

## Identify expanded and contracted gene families for each species


