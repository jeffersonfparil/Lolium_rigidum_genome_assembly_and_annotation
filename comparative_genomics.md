# Comparative genomics

## Set working directories, executables, and output files
```sh}
DIR=/data/Lolium_rigidum_ASSEMBLY/COMPARATIVE_GENOMICS
PATH=${PATH}:${DIR}/OrthoFinder
PATH=${PATH}:${DIR}/CAFE5/bin
MACSE=${DIR}/MACSE/macse_v2.06.jar
PATH=${PATH}:${DIR}/iqtree-2.0.7-Linux/bin
PATH=${PATH}:${DIR}/paml4.9j/bin
PATH=${PATH}:${DIR}/paml4.9j/src
PATH=${PATH}:${DIR}/kakscalculator2-2.0.1/src

DIR_ORTHOGROUPS=${DIR}/ORTHOGROUPS/OrthoFinder/Results_May11
DIR_ORTHOGROUP_SEQS=${DIR}/ORTHOGROUPS/OrthoFinder/Results_*/Orthogroup_Sequences
DIR_PANTHER=${DIR}/PantherHMM_17.0/famlib/rel/PANTHER17.0_altVersion/hmmscoring/PANTHER17.0/books
GOT_PATHER=${DIR}/PantherHMM_17.0/PANTHER17.0_HMM_classifications
MERGED_ORTHOGROUPS=${DIR}/ORTHOGROUPS/orthogroups.faa
ORTHOUT=${DIR}/ORTHOGROUPS/orthogroups_gene_counts_families_go.out
TREE=${DIR_ORTHOGROUPS}/Species_Tree/SpeciesTree_rooted.txt
DIR_ORTHOGROUP_SEQS=${DIR}/ORTHOGROUPS/OrthoFinder/Results_May11/Orthogroup_Sequences
DIR_GENES=${DIR}/TSR_NTSR_GENES
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
# wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/019/359/855/GCA_019359855.1_MPB_Lper_Kyuss_1697/GCA_019359855.1_MPB_Lper_Kyuss_1697_genomic.fna.gz
# gunzip -c GCA_019359855.1_MPB_Lper_Kyuss_1697_genomic.fna.gz > Lolium_perenne.fasta
wget https://de.cyverse.org/anon-files//iplant/home/shared/commons_repo/curated/Copetti_Kyuss_assembly_annotation_March_2021/Kyuss_1697_assembly.fa
wget https://de.cyverse.org/anon-files//iplant/home/shared/commons_repo/curated/Copetti_Kyuss_assembly_annotation_March_2021/Kyuss_1697_KYUS.gff
wget https://de.cyverse.org/anon-files//iplant/home/shared/commons_repo/curated/Copetti_Kyuss_assembly_annotation_March_2021/Kyuss_1697_KYUS_CDS.fa
wget https://de.cyverse.org/anon-files//iplant/home/shared/commons_repo/curated/Copetti_Kyuss_assembly_annotation_March_2021/Kyuss_1697_KYUS_proteins.fa
mv Kyuss_1697_assembly.fa Lolium_perenne.fasta
mv Kyuss_1697_KYUS.gff Lolium_perenne.gff
mv Kyuss_1697_KYUS_CDS.fa Lolium_perenne.cds
mv Kyuss_1697_KYUS_proteins.fa Lolium_perenne.faa
### Secale cereale
wget https://doi.ipk-gatersleben.de/DOI/b9a5ca69-6263-4e6e-a400-f669bee5e92c/3d38d912-c8d7-4587-be08-8b22b7e9f3b2/1/DOWNLOAD
mv DOWNLOAD Secale_cereale.fasta
wget https://doi.ipk-gatersleben.de/DOI/8afb3971-b5e1-4748-8f0e-1b929ba73248/43ce7b0a-f4c6-4996-82bf-9ba90b9ec6b0/1/DOWNLOAD
mv DOWNLOAD Secale_cereale.gff
wget https://doi.ipk-gatersleben.de/DOI/8afb3971-b5e1-4748-8f0e-1b929ba73248/ecfaef55-5b60-4b1f-9c6f-88b3f0f602ad/1/DOWNLOAD
mv DOWNLOAD Secale_cereale.cds
wget https://doi.ipk-gatersleben.de/DOI/8afb3971-b5e1-4748-8f0e-1b929ba73248/01369868-8f23-4a21-834f-113b1a9d922d/1/DOWNLOAD
mv DOWNLOAD Secale_cereale.faa
### Clean-up
rm GCF*.gz GCA*.gz
```

## Install R and Julia
```{sh}
sudo apt install -y r-base julia
```

## Install OrthoFinder for classifying genes into orthologs, and paralogs, as well as to build a tree for the analysis of gene family evolution
```{sh}
wget https://github.com/davidemms/OrthoFinder/releases/download/2.5.4/OrthoFinder.tar.gz
tar -xvzf OrthoFinder.tar.gz
cd OrthoFinder/
./orthofinder -h
PATH=${PATH}:${DIR}/OrthoFinder
cd -
rm OrthoFinder.tar.gz
```

## Install HMMER for mapping CDS to PantherHMM gene family models
```{sh}
sudo apt install -y hmmer
```

## Install CAFE5 to analyse gene family evolution
```{sh}
wget https://github.com/hahnlab/CAFE5/releases/download/v5.0/CAFE5-5.0.0.tar.gz
tar -xvzf CAFE5-5.0.0.tar.gz
cd CAFE5/
./configure
make
PATH=${PATH}:${DIR}/CAFE5/bin
cafe5 -h
cd -
rm CAFE5-5.0.0.tar.gz
```

## Install MACSE: Multiple Alignment of Coding SEquences Accounting for Frameshifts and Stop Codons
```{sh}
mkdir MACSE/
cd MACSE/
wget https://bioweb.supagro.inra.fr/macse/releases/macse_v2.06.jar
MACSE=${DIR}/MACSE/macse_v2.06.jar
java -Xmx250G -jar ${MACSE} -help
cd -
```

## Install IQ-TREE for building trees with fossil root dates
```{sh}
sudo apt install libeigen3-dev
wget https://github.com/Cibiv/IQ-TREE/releases/download/v2.0.7/iqtree-2.0.7-Linux.tar.gz
tar -xvzf iqtree-2.0.7-Linux.tar.gz
PATH=${PATH}:${DIR}/iqtree-2.0.7-Linux/bin
iqtree2 -h
rm iqtree-2.0.7-Linux.tar.gz
```

## Install PAML (Phylogenetic Analysis by Maximum Likelihood) which includes MCMCTree for Bayesian phylogenetic analysis
```{sh}
wget http://abacus.gene.ucl.ac.uk/software/paml4.9j.tgz
tar -xvzf paml4.9j.tgz
PATH=${PATH}:${DIR}/paml4.9j/bin
PATH=${PATH}:${DIR}/paml4.9j/src
rm paml4.9j.tgz
```

## Install KaKs_Calculator2.0 to assess signatures of selection
```{sh}
wget https://github.com/kullrich/kakscalculator2/archive/refs/tags/v2.0.1.tar.gz
tar -xvzf v2.0.1.tar.gz
cd kakscalculator2-2.0.1/src
make
PATH=${PATH}:${DIR}/kakscalculator2-2.0.1/src
KaKs_Calculator -h
cd -
rm v2.0.1.tar.gz
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

## Install Julia packages: DataFrames, CSV, and ProgressMeter
```{julia}
using Pkg
Pkg.add(["DataFrames", "CSV", "ProgressMeter"])
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
rm PANTHER17.0_hmmscoring.tgz
cd -
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

DIR_ORTHOGROUPS=${DIR}/ORTHOGROUPS/OrthoFinder/Results_May11/
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

time \
cafe5 \
    --infile counts.tmp \
    --tree ${TREE} \
    --n_gamma_cats 1000 \
    --cores 32 \
    --pvalue 0.01 \
    --output_prefix CAFE_Gamma1000_results

### Output using n_gamma_cats=1,000
echo -e "Species\tExpansion\tContraction" > CONTRACTION_EXPANSION.txt
grep -v "^#" ${DIR}/CAFE_Gamma1000_results/Gamma_clade_results.txt | \
    grep -v "^<" | \
    sed 's/<..>//g' | \
    sed 's/<.>//g' >> CONTRACTION_EXPANSION.txt
```

## GO term enrichment analysis of contracted and expanded gene families
```{sh}
ORTHOUT=${DIR}/ORTHOGROUPS/orthogroups_gene_counts_families_go.out
n=$(head -n1 ${DIR}/CAFE_Gamma1000_results/Gamma_change.tab | sed -z "s/\t/\n/g" | grep -n "Lolium_rigidum" | cut -d":" -f1)
cut -f1,${n} ${DIR}/CAFE_Gamma1000_results/Gamma_change.tab | grep -v "+0" | grep "+" | cut -f1 > expanded_orthogroups_for_grep.tmp
cut -f1,${n} ${DIR}/CAFE_Gamma1000_results/Gamma_change.tab | grep -v "+" | cut -f1 > contracted_orthogroups_for_grep.tmp

grep -wf expanded_orthogroups_for_grep.tmp $ORTHOUT | cut -f$(head -n1 $ORTHOUT | awk '{printf NF-2}') | grep -v "^$" > expanded_orthogroups.pthr.tmp
grep -wf contracted_orthogroups_for_grep.tmp $ORTHOUT | cut -f$(head -n1 $ORTHOUT | awk '{printf NF-2}') | grep -v "^$" > contracted_orthogroups.pthr.tmp
wget http://data.pantherdb.org/ftp/sequence_classifications/current_release/PANTHER_Sequence_Classification_files/PTHR17.0_rice
grep -wf expanded_orthogroups.pthr.tmp PTHR17.0_rice | cut -f3 > expanded_orthogroups.forgo
grep -wf contracted_orthogroups.pthr.tmp PTHR17.0_rice | cut -f3 > contracted_orthogroups.forgo
### Then paste with Arabidopsis thalian list into http://geneontology.org/
### Set to "Biological process"
### Set the reference gene list to "Oryza sativa"
### Set the correction to "Bonferron correction" and Relaunch analysis
### Sort by decreasing "Fold Enrichment"
### Export as table
```

### expanded_orthogroups.goea
```{sh}
echo '
Analysis Type:	PANTHER Overrepresentation Test (Released 20220202)
Annotation Version and Release Date:	GO Ontology database DOI:  10.5281/zenodo.6399963 Released 2022-03-22
Analyzed List:	upload_1 (Oryza sativa)
Reference List:	Oryza sativa (all genes in database)
Test Type:	FISHER
Correction:	BONFERRONI
Bonferroni count:	2108
GO biological process complete	Oryza sativa - REFLIST (43658)	upload_1 (16219)	upload_1 (expected)	upload_1 (over/under)	upload_1 (fold Enrichment)	upload_1 (P-value)
xylan biosynthetic process (GO:0045492)	39	39	14.49	+	2.69	3.70E-02
export across plasma membrane (GO:0140115)	62	61	23.03	+	2.65	2.52E-04
xyloglucan metabolic process (GO:0010411)	53	51	19.69	+	2.59	4.54E-03
hydrogen peroxide catabolic process (GO:0042744)	144	137	53.50	+	2.56	2.84E-11
hydrogen peroxide metabolic process (GO:0042743)	144	137	53.50	+	2.56	2.84E-11
xenobiotic transport (GO:0042908)	50	47	18.58	+	2.53	1.93E-02
reactive oxygen species metabolic process (GO:0072593)	159	149	59.07	+	2.52	4.18E-12
lignin metabolic process (GO:0009808)	55	51	20.43	+	2.50	9.99E-03
hormone metabolic process (GO:0042445)	93	85	34.55	+	2.46	9.98E-06
oligopeptide transport (GO:0006857)	51	46	18.95	+	2.43	3.83E-02
peptide transport (GO:0015833)	51	46	18.95	+	2.43	3.83E-02
cellular response to chemical stress (GO:0062197)	58	52	21.55	+	2.41	1.54E-02
auxin-activated signaling pathway (GO:0009734)	124	111	46.07	+	2.41	9.37E-08
microtubule-based movement (GO:0007018)	56	50	20.80	+	2.40	2.94E-02
hemicellulose metabolic process (GO:0010410)	112	100	41.61	+	2.40	9.96E-07
plant-type secondary cell wall biogenesis (GO:0009834)	56	50	20.80	+	2.40	2.94E-02
cellular response to ethylene stimulus (GO:0071369)	64	57	23.78	+	2.40	5.83E-03
phosphorelay signal transduction system (GO:0000160)	117	104	43.47	+	2.39	5.18E-07
ethylene-activated signaling pathway (GO:0009873)	63	56	23.40	+	2.39	7.91E-03
endocytosis (GO:0006897)	71	63	26.38	+	2.39	1.69E-03
cytokinin-activated signaling pathway (GO:0009736)	62	55	23.03	+	2.39	1.08E-02
protein targeting to membrane (GO:0006612)	77	68	28.61	+	2.38	9.38E-04
cellular component macromolecule biosynthetic process (GO:0070589)	77	68	28.61	+	2.38	9.38E-04
cell wall macromolecule biosynthetic process (GO:0044038)	77	68	28.61	+	2.38	9.38E-04
cellular response to auxin stimulus (GO:0071365)	126	111	46.81	+	2.37	1.88E-07
amide transport (GO:0042886)	66	58	24.52	+	2.37	7.58E-03
cell wall polysaccharide biosynthetic process (GO:0070592)	74	65	27.49	+	2.36	1.55E-03
hormone biosynthetic process (GO:0042446)	64	56	23.78	+	2.36	9.42E-03
cellular response to cytokinin stimulus (GO:0071368)	63	55	23.40	+	2.35	1.26E-02
response to water (GO:0009415)	70	61	26.01	+	2.35	5.32E-03
cell wall biogenesis (GO:0042546)	168	146	62.41	+	2.34	5.22E-10
response to water deprivation (GO:0009414)	68	59	25.26	+	2.34	6.48E-03
hormone-mediated signaling pathway (GO:0009755)	410	355	152.32	+	2.33	4.66E-27
regulation of hormone levels (GO:0010817)	134	116	49.78	+	2.33	1.44E-07
pectin metabolic process (GO:0045488)	89	77	33.06	+	2.33	2.15E-04
galacturonan metabolic process (GO:0010393)	89	77	33.06	+	2.33	2.15E-04
cell wall organization (GO:0071555)	314	271	116.65	+	2.32	3.89E-20
cell wall macromolecule metabolic process (GO:0044036)	146	126	54.24	+	2.32	3.19E-08
pectin catabolic process (GO:0045490)	65	56	24.15	+	2.32	1.59E-02
cellular response to hormone stimulus (GO:0032870)	418	360	155.29	+	2.32	3.43E-27
potassium ion transmembrane transport (GO:0071805)	64	55	23.78	+	2.31	2.18E-02
external encapsulating structure organization (GO:0045229)	340	291	126.31	+	2.30	2.37E-21
response to acid chemical (GO:0001101)	76	65	28.23	+	2.30	3.07E-03
cellular response to endogenous stimulus (GO:0071495)	425	363	157.89	+	2.30	5.23E-27
steroid metabolic process (GO:0008202)	102	87	37.89	+	2.30	5.26E-05
carbohydrate transmembrane transport (GO:0034219)	88	75	32.69	+	2.29	6.68E-04
response to cytokinin (GO:0009735)	67	57	24.89	+	2.29	1.39E-02
cellulose metabolic process (GO:0030243)	73	62	27.12	+	2.29	7.66E-03
response to oomycetes (GO:0002239)	79	67	29.35	+	2.28	2.85E-03
defense response to oomycetes (GO:0002229)	79	67	29.35	+	2.28	2.85E-03
detoxification (GO:0098754)	276	234	102.53	+	2.28	1.89E-16
peptidyl-serine phosphorylation (GO:0018105)	84	71	31.21	+	2.28	1.45E-03
response to ethylene (GO:0009723)	84	71	31.21	+	2.28	1.45E-03
steroid biosynthetic process (GO:0006694)	71	60	26.38	+	2.27	9.45E-03
transmembrane receptor protein serine/threonine kinase signaling pathway (GO:0007178)	77	65	28.61	+	2.27	5.32E-03
enzyme-linked receptor protein signaling pathway (GO:0007167)	77	65	28.61	+	2.27	5.32E-03
cell wall organization or biogenesis (GO:0071554)	415	350	154.17	+	2.27	2.52E-25
regulation of cyclin-dependent protein serine/threonine kinase activity (GO:0000079)	70	59	26.01	+	2.27	1.26E-02
protein autophosphorylation (GO:0046777)	70	59	26.01	+	2.27	1.26E-02
xylan metabolic process (GO:0045491)	63	53	23.40	+	2.26	4.52E-02
plant-type cell wall organization or biogenesis (GO:0071669)	156	131	57.95	+	2.26	4.51E-08
regulation of protein serine/threonine kinase activity (GO:0071900)	81	68	30.09	+	2.26	2.49E-03
cellular oxidant detoxification (GO:0098869)	218	183	80.99	+	2.26	3.67E-12
sterol metabolic process (GO:0016125)	87	73	32.32	+	2.26	1.36E-03
response to toxic substance (GO:0009636)	295	247	109.59	+	2.25	6.08E-17
cell wall polysaccharide metabolic process (GO:0010383)	122	102	45.32	+	2.25	8.12E-06
peptidyl-serine modification (GO:0018209)	85	71	31.58	+	2.25	1.70E-03
SCF-dependent proteasomal ubiquitin-dependent protein catabolic process (GO:0031146)	91	76	33.81	+	2.25	9.43E-04
cellular response to toxic substance (GO:0097237)	234	195	86.93	+	2.24	8.33E-13
cellular detoxification (GO:1990748)	234	195	86.93	+	2.24	8.33E-13
regulation of cyclin-dependent protein kinase activity (GO:1904029)	72	60	26.75	+	2.24	1.57E-02
potassium ion transport (GO:0006813)	70	58	26.01	+	2.23	1.99E-02
cellular response to organic substance (GO:0071310)	483	399	179.44	+	2.22	6.81E-28
polysaccharide biosynthetic process (GO:0000271)	177	146	65.76	+	2.22	7.88E-09
cellular response to chemical stimulus (GO:0070887)	758	624	281.60	+	2.22	1.32E-44
carbohydrate transport (GO:0008643)	113	93	41.98	+	2.22	5.51E-05
positive regulation of mRNA metabolic process (GO:1903313)	67	55	24.89	+	2.21	4.71E-02
positive regulation of cellular catabolic process (GO:0031331)	111	91	41.24	+	2.21	9.63E-05
cellular response to abscisic acid stimulus (GO:0071215)	83	68	30.83	+	2.21	4.67E-03
cellular response to alcohol (GO:0097306)	83	68	30.83	+	2.21	4.67E-03
secondary metabolic process (GO:0019748)	176	144	65.38	+	2.20	1.57E-08
plant-type cell wall biogenesis (GO:0009832)	108	88	40.12	+	2.19	2.32E-04
protein dephosphorylation (GO:0006470)	187	152	69.47	+	2.19	4.72E-09
positive regulation of catabolic process (GO:0009896)	128	104	47.55	+	2.19	1.57E-05
response to cold (GO:0009409)	90	73	33.44	+	2.18	2.86E-03
abscisic acid-activated signaling pathway (GO:0009738)	79	64	29.35	+	2.18	1.49E-02
cellular polysaccharide biosynthetic process (GO:0033692)	158	128	58.70	+	2.18	2.91E-07
phyllome development (GO:0048827)	100	81	37.15	+	2.18	7.46E-04
response to hormone (GO:0009725)	613	495	227.73	+	2.17	2.02E-33
polysaccharide metabolic process (GO:0005976)	395	318	146.74	+	2.17	1.74E-20
response to endogenous stimulus (GO:0009719)	620	498	230.33	+	2.16	3.97E-33
cellular response to lipid (GO:0071396)	161	129	59.81	+	2.16	4.00E-07
cellular response to oxygen-containing compound (GO:1901701)	206	165	76.53	+	2.16	1.39E-09
response to auxin (GO:0009733)	215	172	79.87	+	2.15	4.89E-10
export from cell (GO:0140352)	154	123	57.21	+	2.15	1.35E-06
response to lipid (GO:0033993)	250	199	92.88	+	2.14	9.90E-12
protein deubiquitination (GO:0016579)	97	77	36.04	+	2.14	2.66E-03
polysaccharide catabolic process (GO:0000272)	188	149	69.84	+	2.13	3.61E-08
peptidyl-proline modification (GO:0018208)	80	63	29.72	+	2.12	2.58E-02
beta-glucan metabolic process (GO:0051273)	84	66	31.21	+	2.11	1.73E-02
phenylpropanoid metabolic process (GO:0009698)	112	88	41.61	+	2.11	5.42E-04
cellular polysaccharide metabolic process (GO:0044264)	259	203	96.22	+	2.11	1.61E-11
response to abscisic acid (GO:0009737)	151	118	56.10	+	2.10	7.00E-06
regulation of protein kinase activity (GO:0045859)	95	74	35.29	+	2.10	7.00E-03
response to oxidative stress (GO:0006979)	275	213	102.16	+	2.08	7.40E-12
intracellular signal transduction (GO:0035556)	352	272	130.77	+	2.08	1.27E-15
regulation of protein phosphorylation (GO:0001932)	101	78	37.52	+	2.08	5.76E-03
glucan metabolic process (GO:0044042)	202	156	75.04	+	2.08	4.32E-08
protein localization to membrane (GO:0072657)	144	111	53.50	+	2.07	3.51E-05
cellular glucan metabolic process (GO:0006073)	196	151	72.81	+	2.07	1.11E-07
establishment of protein localization to membrane (GO:0090150)	117	90	43.47	+	2.07	8.39E-04
response to chemical (GO:0042221)	1189	914	441.71	+	2.07	2.32E-57
metal ion transport (GO:0030001)	224	172	83.22	+	2.07	5.27E-09
positive regulation of transcription by RNA polymerase II (GO:0045944)	155	119	57.58	+	2.07	1.45E-05
response to organic substance (GO:0010033)	757	580	281.23	+	2.06	5.29E-35
peptidyl-tyrosine modification (GO:0018212)	94	72	34.92	+	2.06	1.37E-02
peptidyl-tyrosine phosphorylation (GO:0018108)	94	72	34.92	+	2.06	1.37E-02
cellular carbohydrate biosynthetic process (GO:0034637)	205	157	76.16	+	2.06	5.66E-08
secondary metabolite biosynthetic process (GO:0044550)	97	74	36.04	+	2.05	1.22E-02
regulation of kinase activity (GO:0043549)	97	74	36.04	+	2.05	1.22E-02
regulation of protein modification process (GO:0031399)	168	128	62.41	+	2.05	4.99E-06
response to alcohol (GO:0097305)	155	118	57.58	+	2.05	2.11E-05
response to oxygen-containing compound (GO:1901700)	428	325	159.00	+	2.04	3.63E-18
signal transduction (GO:0007165)	1097	832	407.54	+	2.04	4.04E-50
regulation of cellular catabolic process (GO:0031329)	132	100	49.04	+	2.04	3.39E-04
regulation of phosphorylation (GO:0042325)	103	78	38.26	+	2.04	7.07E-03
inorganic anion transport (GO:0015698)	98	74	36.41	+	2.03	1.38E-02
chloroplast organization (GO:0009658)	133	100	49.41	+	2.02	3.83E-04
organelle localization (GO:0051640)	124	93	46.07	+	2.02	1.10E-03
plant organ development (GO:0099402)	255	191	94.73	+	2.02	1.97E-09
response to bacterium (GO:0009617)	159	119	59.07	+	2.01	4.27E-05
signaling (GO:0023052)	1112	832	413.11	+	2.01	1.83E-48
fatty acid biosynthetic process (GO:0006633)	127	95	47.18	+	2.01	9.65E-04
glucan biosynthetic process (GO:0009250)	115	86	42.72	+	2.01	4.44E-03
vesicle organization (GO:0016050)	99	74	36.78	+	2.01	2.10E-02
positive regulation of RNA metabolic process (GO:0051254)	386	288	143.40	+	2.01	3.73E-15
shoot system development (GO:0048367)	265	197	98.45	+	2.00	1.09E-09
positive regulation of protein metabolic process (GO:0051247)	132	98	49.04	+	2.00	1.00E-03
regulation of catabolic process (GO:0009894)	155	115	57.58	+	2.00	8.96E-05
positive regulation of cellular protein metabolic process (GO:0032270)	115	85	42.72	+	1.99	6.28E-03
regulation of transferase activity (GO:0051338)	111	82	41.24	+	1.99	9.30E-03
anatomical structure morphogenesis (GO:0009653)	196	144	72.81	+	1.98	2.86E-06
positive regulation of nucleobase-containing compound metabolic process (GO:0045935)	399	293	148.23	+	1.98	7.56E-15
carbohydrate catabolic process (GO:0016052)	293	215	108.85	+	1.98	2.91E-10
regulation of phosphate metabolic process (GO:0019220)	154	113	57.21	+	1.98	1.69E-04
regulation of phosphorus metabolic process (GO:0051174)	154	113	57.21	+	1.98	1.69E-04
root system development (GO:0022622)	120	88	44.58	+	1.97	4.66E-03
root development (GO:0048364)	120	88	44.58	+	1.97	4.66E-03
homeostatic process (GO:0042592)	303	222	112.56	+	1.97	1.11E-10
positive regulation of RNA biosynthetic process (GO:1902680)	317	232	117.77	+	1.97	2.97E-11
positive regulation of transcription, DNA-templated (GO:0045893)	317	232	117.77	+	1.97	2.97E-11
positive regulation of nucleic acid-templated transcription (GO:1903508)	317	232	117.77	+	1.97	2.97E-11
response to osmotic stress (GO:0006970)	115	84	42.72	+	1.97	9.03E-03
anion transport (GO:0006820)	204	149	75.79	+	1.97	1.88E-06
positive regulation of cellular biosynthetic process (GO:0031328)	348	254	129.28	+	1.96	2.42E-12
positive regulation of macromolecule biosynthetic process (GO:0010557)	340	248	126.31	+	1.96	5.27E-12
positive regulation of biosynthetic process (GO:0009891)	350	255	130.03	+	1.96	2.82E-12
positive regulation of nitrogen compound metabolic process (GO:0051173)	522	380	193.92	+	1.96	2.36E-19
glycosylation (GO:0070085)	208	151	77.27	+	1.95	1.83E-06
chemical homeostasis (GO:0048878)	270	196	100.31	+	1.95	6.34E-09
positive regulation of macromolecule metabolic process (GO:0010604)	528	383	196.15	+	1.95	2.60E-19
positive regulation of cellular metabolic process (GO:0031325)	515	373	191.32	+	1.95	1.09E-18
positive regulation of metabolic process (GO:0009893)	544	393	202.10	+	1.94	1.11E-19
microtubule-based process (GO:0007017)	187	135	69.47	+	1.94	2.30E-05
transmembrane transport (GO:0055085)	1365	984	507.10	+	1.94	7.93E-53
defense response to bacterium (GO:0042742)	136	98	50.52	+	1.94	2.76E-03
carbohydrate biosynthetic process (GO:0016051)	278	200	103.28	+	1.94	5.97E-09
response to inorganic substance (GO:0010035)	174	125	64.64	+	1.93	9.57E-05
regulation of biological quality (GO:0065008)	585	420	217.33	+	1.93	9.29E-21
cellular homeostasis (GO:0019725)	159	114	59.07	+	1.93	3.49E-04
protein modification by small protein removal (GO:0070646)	120	86	44.58	+	1.93	1.31E-02
response to temperature stimulus (GO:0009266)	194	139	72.07	+	1.93	1.94E-05
cellular carbohydrate metabolic process (GO:0044262)	427	305	158.63	+	1.92	2.94E-14
cell communication (GO:0007154)	1274	909	473.29	+	1.92	3.91E-47
response to abiotic stimulus (GO:0009628)	597	425	221.79	+	1.92	1.55E-20
fatty acid metabolic process (GO:0006631)	194	138	72.07	+	1.91	2.70E-05
carbohydrate metabolic process (GO:0005975)	1039	739	385.99	+	1.91	1.86E-37
negative regulation of RNA biosynthetic process (GO:1902679)	114	81	42.35	+	1.91	3.31E-02
negative regulation of nucleic acid-templated transcription (GO:1903507)	114	81	42.35	+	1.91	3.31E-02
localization within membrane (GO:0051668)	162	115	60.18	+	1.91	4.36E-04
cell division (GO:0051301)	196	139	72.81	+	1.91	3.13E-05
negative regulation of transcription, DNA-templated (GO:0045892)	113	80	41.98	+	1.91	3.18E-02
positive regulation of cellular process (GO:0048522)	578	408	214.73	+	1.90	4.37E-19
cellular response to stimulus (GO:0051716)	1912	1346	710.31	+	1.89	1.79E-69
negative regulation of catalytic activity (GO:0043086)	229	161	85.07	+	1.89	2.94E-06
negative regulation of molecular function (GO:0044092)	229	161	85.07	+	1.89	2.94E-06
regulation of RNA metabolic process (GO:0051252)	2135	1495	793.16	+	1.88	1.13E-76
regulation of RNA biosynthetic process (GO:2001141)	2028	1420	753.40	+	1.88	1.71E-72
regulation of transcription, DNA-templated (GO:0006355)	2028	1420	753.40	+	1.88	1.71E-72
regulation of nucleic acid-templated transcription (GO:1903506)	2028	1420	753.40	+	1.88	1.71E-72
dephosphorylation (GO:0016311)	374	261	138.94	+	1.88	4.10E-11
system development (GO:0048731)	492	343	182.78	+	1.88	3.37E-15
cellular chemical homeostasis (GO:0055082)	125	87	46.44	+	1.87	2.58E-02
inorganic cation transmembrane transport (GO:0098662)	299	208	111.08	+	1.87	2.45E-08
sulfur compound metabolic process (GO:0006790)	305	212	113.31	+	1.87	1.40E-08
protein folding (GO:0006457)	229	159	85.07	+	1.87	7.78E-06
ion homeostasis (GO:0050801)	193	134	71.70	+	1.87	1.30E-04
response to biotic stimulus (GO:0009607)	346	240	128.54	+	1.87	8.33E-10
organic substance transport (GO:0071702)	1246	864	462.89	+	1.87	2.82E-41
anion transmembrane transport (GO:0098656)	127	88	47.18	+	1.87	2.97E-02
lipid catabolic process (GO:0016042)	140	97	52.01	+	1.87	9.79E-03
response to external biotic stimulus (GO:0043207)	319	221	118.51	+	1.86	6.85E-09
response to other organism (GO:0051707)	319	221	118.51	+	1.86	6.85E-09
ion transmembrane transport (GO:0034220)	529	366	196.52	+	1.86	5.59E-16
regulation of nucleobase-containing compound metabolic process (GO:0019219)	2192	1515	814.33	+	1.86	4.36E-75
establishment of localization (GO:0051234)	2519	1741	935.81	+	1.86	1.76E-87
transport (GO:0006810)	2485	1717	923.18	+	1.86	4.59E-86
peptidyl-amino acid modification (GO:0018193)	469	323	174.23	+	1.85	1.29E-13
negative regulation of RNA metabolic process (GO:0051253)	122	84	45.32	+	1.85	4.19E-02
organic acid transport (GO:0015849)	138	95	51.27	+	1.85	1.22E-02
localization (GO:0051179)	2618	1802	972.59	+	1.85	8.59E-90
response to external stimulus (GO:0009605)	455	313	169.03	+	1.85	3.66E-13
organic hydroxy compound metabolic process (GO:1901615)	285	196	105.88	+	1.85	1.78E-07
regulation of macromolecule biosynthetic process (GO:0010556)	2237	1536	831.05	+	1.85	9.59E-75
inorganic ion transmembrane transport (GO:0098660)	347	238	128.91	+	1.85	2.29E-09
inorganic ion homeostasis (GO:0098771)	178	122	66.13	+	1.84	8.97E-04
plastid organization (GO:0009657)	165	113	61.30	+	1.84	1.98E-03
regulation of cellular biosynthetic process (GO:0031326)	2258	1546	838.85	+	1.84	9.31E-75
response to stimulus (GO:0050896)	3698	2528	1373.81	+	1.84	1.22E-128
regulation of biosynthetic process (GO:0009889)	2265	1548	841.45	+	1.84	1.95E-74
nitrogen compound transport (GO:0071705)	1005	686	373.36	+	1.84	1.25E-30
cell surface receptor signaling pathway (GO:0007166)	261	178	96.96	+	1.84	2.42E-06
regulation of nitrogen compound metabolic process (GO:0051171)	2644	1798	982.25	+	1.83	1.66E-86
positive regulation of biological process (GO:0048518)	687	467	255.22	+	1.83	7.63E-20
protein transport (GO:0015031)	689	468	255.96	+	1.83	8.64E-20
ion transport (GO:0006811)	631	428	234.42	+	1.83	6.52E-18
regulation of gene expression (GO:0010468)	2484	1683	922.81	+	1.82	1.33E-79
regulation of molecular function (GO:0065009)	527	357	195.78	+	1.82	1.56E-14
regulation of cellular process (GO:0050794)	3867	2619	1436.60	+	1.82	1.96E-130
regulation of primary metabolic process (GO:0080090)	2669	1807	991.54	+	1.82	7.19E-86
cation homeostasis (GO:0055080)	173	117	64.27	+	1.82	2.31E-03
regulation of cellular metabolic process (GO:0031323)	2690	1814	999.34	+	1.82	3.36E-85
biological regulation (GO:0065007)	4867	3282	1808.10	+	1.82	1.61E-167
establishment of protein localization (GO:0045184)	697	470	258.94	+	1.82	1.86E-19
flower development (GO:0009908)	141	95	52.38	+	1.81	2.70E-02
regulation of catalytic activity (GO:0050790)	518	349	192.44	+	1.81	5.98E-14
regulation of macromolecule metabolic process (GO:0060255)	2862	1928	1063.24	+	1.81	5.44E-91
reproductive shoot system development (GO:0090567)	147	99	54.61	+	1.81	1.55E-02
protein localization (GO:0008104)	767	516	284.94	+	1.81	1.94E-21
response to light stimulus (GO:0009416)	213	143	79.13	+	1.81	1.76E-04
regulation of biological process (GO:0050789)	4286	2869	1592.25	+	1.80	3.98E-140
vesicle-mediated transport (GO:0016192)	493	330	183.15	+	1.80	9.89E-13
regulation of proteolysis (GO:0030162)	139	93	51.64	+	1.80	4.48E-02
cation transmembrane transport (GO:0098655)	362	242	134.48	+	1.80	1.02E-08
regulation of metabolic process (GO:0019222)	2943	1967	1093.33	+	1.80	7.40E-91
multicellular organism development (GO:0007275)	651	435	241.85	+	1.80	2.09E-17
protein phosphorylation (GO:0006468)	1515	1011	562.82	+	1.80	1.37E-43
protein glycosylation (GO:0006486)	147	98	54.61	+	1.79	2.81E-02
macromolecule glycosylation (GO:0043413)	147	98	54.61	+	1.79	2.81E-02
catabolic process (GO:0009056)	1640	1093	609.26	+	1.79	2.73E-47
cellular protein localization (GO:0034613)	553	367	205.44	+	1.79	5.06E-14
cellular macromolecule localization (GO:0070727)	555	368	206.18	+	1.78	5.59E-14
glycoprotein biosynthetic process (GO:0009101)	148	98	54.98	+	1.78	3.00E-02
defense response to other organism (GO:0098542)	275	182	102.16	+	1.78	7.33E-06
cation transport (GO:0006812)	387	256	143.77	+	1.78	4.98E-09
glycoprotein metabolic process (GO:0009100)	159	105	59.07	+	1.78	1.64E-02
monocarboxylic acid metabolic process (GO:0032787)	391	258	145.26	+	1.78	4.62E-09
reproductive structure development (GO:0048608)	291	192	108.11	+	1.78	2.97E-06
regulation of cellular protein metabolic process (GO:0032268)	501	330	186.12	+	1.77	5.05E-12
regulation of protein metabolic process (GO:0051246)	527	347	195.78	+	1.77	8.07E-13
organic hydroxy compound biosynthetic process (GO:1901617)	155	102	57.58	+	1.77	2.42E-02
protein modification process (GO:0036211)	3176	2090	1179.89	+	1.77	8.72E-93
cellular protein modification process (GO:0006464)	3176	2090	1179.89	+	1.77	8.72E-93
reproductive system development (GO:0061458)	292	192	108.48	+	1.77	4.16E-06
regulation of hydrolase activity (GO:0051336)	146	96	54.24	+	1.77	4.87E-02
regulation of transcription by RNA polymerase II (GO:0006357)	516	339	191.69	+	1.77	2.57E-12
phosphorylation (GO:0016310)	1841	1209	683.93	+	1.77	1.92E-50
intracellular protein transport (GO:0006886)	477	313	177.21	+	1.77	3.60E-11
monosaccharide metabolic process (GO:0005996)	157	103	58.33	+	1.77	2.66E-02
monocarboxylic acid biosynthetic process (GO:0072330)	174	114	64.64	+	1.76	8.39E-03
response to radiation (GO:0009314)	220	144	81.73	+	1.76	4.74E-04
macromolecule modification (GO:0043412)	3495	2286	1298.40	+	1.76	1.14E-100
lipid metabolic process (GO:0006629)	961	628	357.01	+	1.76	2.74E-24
protein metabolic process (GO:0019538)	4589	2991	1704.82	+	1.75	1.15E-135
nucleobase-containing compound transport (GO:0015931)	152	99	56.47	+	1.75	4.97E-02
macromolecule localization (GO:0033036)	954	621	354.41	+	1.75	1.11E-23
protein targeting (GO:0006605)	226	147	83.96	+	1.75	4.81E-04
cellular protein metabolic process (GO:0044267)	4088	2657	1518.70	+	1.75	3.22E-117
cellular localization (GO:0051641)	819	532	304.26	+	1.75	7.85E-20
protein modification by small protein conjugation or removal (GO:0070647)	915	594	339.92	+	1.75	2.85E-22
anatomical structure development (GO:0048856)	768	497	285.31	+	1.74	4.54E-18
regulation of developmental process (GO:0050793)	223	144	82.84	+	1.74	9.74E-04
cellular catabolic process (GO:0044248)	1349	871	501.16	+	1.74	2.58E-33
developmental process (GO:0032502)	899	580	333.98	+	1.74	3.05E-21
RNA splicing, via transesterification reactions (GO:0000375)	217	140	80.62	+	1.74	1.30E-03
cellular macromolecule biosynthetic process (GO:0034645)	951	613	353.30	+	1.74	1.69E-22
biological_process (GO:0008150)	17397	11206	6463.01	+	1.73	0.00E00
organonitrogen compound metabolic process (GO:1901564)	5580	3592	2072.98	+	1.73	2.54E-162
proteasomal protein catabolic process (GO:0010498)	275	177	102.16	+	1.73	5.62E-05
organic substance catabolic process (GO:1901575)	1434	922	532.73	+	1.73	4.93E-35
proteasome-mediated ubiquitin-dependent protein catabolic process (GO:0043161)	246	158	91.39	+	1.73	3.11E-04
RNA splicing, via transesterification reactions with bulged adenosine as nucleophile (GO:0000377)	215	138	79.87	+	1.73	2.06E-03
macromolecule catabolic process (GO:0009057)	952	611	353.67	+	1.73	4.14E-22
post-embryonic development (GO:0009791)	381	244	141.54	+	1.72	1.58E-07
protein modification by small protein conjugation (GO:0032446)	794	508	294.97	+	1.72	9.05E-18
cellular process (GO:0009987)	12820	8196	4762.65	+	1.72	0.00E00
establishment of localization in cell (GO:0051649)	703	449	261.17	+	1.72	2.26E-15
proteolysis (GO:0006508)	1219	778	452.86	+	1.72	2.62E-28
phosphorus metabolic process (GO:0006793)	2618	1663	972.59	+	1.71	3.13E-64
methylation (GO:0032259)	334	212	124.08	+	1.71	4.59E-06
biological process involved in interspecies interaction between organisms (GO:0044419)	353	224	131.14	+	1.71	1.60E-06
developmental process involved in reproduction (GO:0003006)	320	203	118.88	+	1.71	1.25E-05
metabolic process (GO:0008152)	10316	6541	3832.41	+	1.71	0.00E00
cellular macromolecule metabolic process (GO:0044260)	5139	3257	1909.14	+	1.71	2.90E-137
phosphate-containing compound metabolic process (GO:0006796)	2581	1635	958.84	+	1.71	1.96E-62
protein ubiquitination (GO:0016567)	755	478	280.48	+	1.70	5.13E-16
negative regulation of macromolecule biosynthetic process (GO:0010558)	228	144	84.70	+	1.70	2.20E-03
defense response (GO:0006952)	1040	656	386.36	+	1.70	1.58E-22
regulation of cell cycle (GO:0051726)	211	133	78.39	+	1.70	6.92E-03
organic substance metabolic process (GO:0071704)	9447	5946	3509.57	+	1.69	3.63E-289
organic acid biosynthetic process (GO:0016053)	429	270	159.37	+	1.69	5.35E-08
response to stress (GO:0006950)	2142	1346	795.76	+	1.69	4.95E-49
RNA splicing (GO:0008380)	261	164	96.96	+	1.69	4.81E-04
lipid biosynthetic process (GO:0008610)	548	344	203.58	+	1.69	1.51E-10
intracellular transport (GO:0046907)	671	421	249.28	+	1.69	2.25E-13
primary metabolic process (GO:0044238)	8914	5588	3311.56	+	1.69	4.46E-262
peptide metabolic process (GO:0006518)	715	448	265.62	+	1.69	2.67E-14
mitotic cell cycle (GO:0000278)	222	139	82.47	+	1.69	5.00E-03
post-transcriptional regulation of gene expression (GO:0010608)	270	169	100.31	+	1.68	4.11E-04
macromolecule biosynthetic process (GO:0009059)	1321	826	490.75	+	1.68	3.92E-28
mitotic cell cycle process (GO:1903047)	184	115	68.36	+	1.68	4.22E-02
organic acid metabolic process (GO:0006082)	883	551	328.04	+	1.68	7.99E-18
macromolecule metabolic process (GO:0043170)	6945	4333	2580.08	+	1.68	7.09E-185
cellular metabolic process (GO:0044237)	8971	5597	3332.74	+	1.68	1.33E-258
carboxylic acid biosynthetic process (GO:0046394)	405	252	150.46	+	1.67	6.13E-07
biosynthetic process (GO:0009058)	2747	1705	1020.51	+	1.67	6.32E-61
mRNA splicing, via spliceosome (GO:0000398)	195	121	72.44	+	1.67	3.04E-02
carboxylic acid metabolic process (GO:0019752)	840	519	312.06	+	1.66	4.65E-16
negative regulation of cellular biosynthetic process (GO:0031327)	238	147	88.42	+	1.66	4.46E-03
oxoacid metabolic process (GO:0043436)	851	524	316.15	+	1.66	5.56E-16
aromatic compound catabolic process (GO:0019439)	242	149	89.90	+	1.66	4.08E-03
cytoskeleton organization (GO:0007010)	223	137	82.84	+	1.65	1.18E-02
cellular biosynthetic process (GO:0044249)	2505	1536	930.61	+	1.65	5.53E-52
organic substance biosynthetic process (GO:1901576)	2612	1601	970.36	+	1.65	2.13E-54
cellular amide metabolic process (GO:0043603)	826	506	306.86	+	1.65	4.67E-15
cellular macromolecule catabolic process (GO:0044265)	793	485	294.60	+	1.65	2.85E-14
organic cyclic compound catabolic process (GO:1901361)	247	151	91.76	+	1.65	4.96E-03
cellular lipid metabolic process (GO:0044255)	684	418	254.11	+	1.64	5.51E-12
small molecule biosynthetic process (GO:0044283)	578	353	214.73	+	1.64	9.01E-10
nitrogen compound metabolic process (GO:0006807)	7413	4525	2753.94	+	1.64	2.58E-181
negative regulation of biosynthetic process (GO:0009890)	241	147	89.53	+	1.64	8.52E-03
cellular component organization or biogenesis (GO:0071840)	2322	1411	862.63	+	1.64	9.91E-46
cellular component organization (GO:0016043)	1998	1211	742.26	+	1.63	3.09E-38
generation of precursor metabolites and energy (GO:0006091)	440	266	163.46	+	1.63	1.50E-06
protein catabolic process (GO:0030163)	622	376	231.07	+	1.63	4.89E-10
proteolysis involved in cellular protein catabolic process (GO:0051603)	612	369	227.36	+	1.62	9.31E-10
membrane organization (GO:0061024)	224	135	83.22	+	1.62	2.69E-02
small molecule metabolic process (GO:0044281)	1436	865	533.48	+	1.62	8.51E-26
cellular protein catabolic process (GO:0044257)	613	369	227.73	+	1.62	1.23E-09
amide biosynthetic process (GO:0043604)	662	397	245.93	+	1.61	2.26E-10
RNA modification (GO:0009451)	302	181	112.19	+	1.61	1.49E-03
organonitrogen compound catabolic process (GO:1901565)	786	471	292.00	+	1.61	1.25E-12
multicellular organismal process (GO:0032501)	830	496	308.35	+	1.61	3.11E-13
translation (GO:0006412)	584	348	216.96	+	1.60	1.38E-08
regulation of response to stimulus (GO:0048583)	306	182	113.68	+	1.60	1.72E-03
peptide biosynthetic process (GO:0043043)	596	353	221.41	+	1.59	1.65E-08
cellular developmental process (GO:0048869)	235	139	87.30	+	1.59	4.04E-02
carbohydrate derivative metabolic process (GO:1901135)	623	368	231.45	+	1.59	7.66E-09
organonitrogen compound biosynthetic process (GO:1901566)	1372	810	509.70	+	1.59	5.04E-22
organic cyclic compound biosynthetic process (GO:1901362)	868	510	322.46	+	1.58	9.03E-13
aromatic compound biosynthetic process (GO:0019438)	769	450	285.68	+	1.58	8.87E-11
nucleic acid-templated transcription (GO:0097659)	282	165	104.76	+	1.57	1.10E-02
modification-dependent macromolecule catabolic process (GO:0043632)	547	320	203.21	+	1.57	4.82E-07
modification-dependent protein catabolic process (GO:0019941)	537	313	199.50	+	1.57	9.04E-07
protein localization to organelle (GO:0033365)	292	170	108.48	+	1.57	1.19E-02
negative regulation of nitrogen compound metabolic process (GO:0051172)	368	214	136.71	+	1.57	6.20E-04
transcription, DNA-templated (GO:0006351)	276	160	102.53	+	1.56	2.32E-02
mRNA metabolic process (GO:0016071)	474	274	176.09	+	1.56	2.06E-05
negative regulation of macromolecule metabolic process (GO:0010605)	525	303	195.04	+	1.55	3.94E-06
cellular component biogenesis (GO:0044085)	1116	644	414.60	+	1.55	2.22E-15
cellular amino acid metabolic process (GO:0006520)	403	232	149.71	+	1.55	3.64E-04
RNA biosynthetic process (GO:0032774)	287	165	106.62	+	1.55	2.58E-02
negative regulation of biological process (GO:0048519)	707	406	262.65	+	1.55	1.07E-08
negative regulation of cellular process (GO:0048523)	507	291	188.35	+	1.54	1.14E-05
cellular nitrogen compound biosynthetic process (GO:0044271)	1298	742	482.21	+	1.54	2.60E-17
gene expression (GO:0010467)	1725	984	640.84	+	1.54	1.25E-23
cellular response to stress (GO:0033554)	639	364	237.39	+	1.53	3.22E-07
ubiquitin-dependent protein catabolic process (GO:0006511)	513	292	190.58	+	1.53	2.11E-05
ribonucleoprotein complex biogenesis (GO:0022613)	450	256	167.18	+	1.53	1.91E-04
negative regulation of cellular metabolic process (GO:0031324)	378	215	140.43	+	1.53	2.08E-03
negative regulation of metabolic process (GO:0009892)	541	306	200.98	+	1.52	1.54E-05
ribosome biogenesis (GO:0042254)	358	202	133.00	+	1.52	6.97E-03
carbohydrate derivative biosynthetic process (GO:1901137)	392	221	145.63	+	1.52	2.48E-03
RNA processing (GO:0006396)	805	452	299.06	+	1.51	7.60E-09
RNA metabolic process (GO:0016070)	1458	817	541.65	+	1.51	1.28E-17
organelle organization (GO:0006996)	1279	714	475.15	+	1.50	9.19E-15
mRNA processing (GO:0006397)	341	190	126.68	+	1.50	2.50E-02
heterocycle biosynthetic process (GO:0018130)	711	392	264.14	+	1.48	1.13E-06
nucleobase-containing compound biosynthetic process (GO:0034654)	541	298	200.98	+	1.48	1.60E-04
organic cyclic compound metabolic process (GO:1901360)	2840	1545	1055.06	+	1.46	1.85E-31
cellular aromatic compound metabolic process (GO:0006725)	2745	1489	1019.77	+	1.46	1.15E-29
ncRNA processing (GO:0034470)	428	232	159.00	+	1.46	1.13E-02
cell cycle (GO:0007049)	399	216	148.23	+	1.46	2.39E-02
cellular nitrogen compound metabolic process (GO:0034641)	3203	1733	1189.92	+	1.46	8.15E-35
ncRNA metabolic process (GO:0034660)	518	279	192.44	+	1.45	1.60E-03
organophosphate metabolic process (GO:0019637)	616	325	228.84	+	1.42	9.77E-04
reproduction (GO:0000003)	576	303	213.98	+	1.42	2.89E-03
heterocycle metabolic process (GO:0046483)	2613	1364	970.73	+	1.41	9.14E-22
nucleobase-containing compound metabolic process (GO:0006139)	2375	1224	882.32	+	1.39	9.29E-18
nucleic acid metabolic process (GO:0090304)	1995	1022	741.14	+	1.38	7.61E-14
reproductive process (GO:0022414)	556	284	206.55	+	1.37	3.45E-02
cellular component assembly (GO:0022607)	670	336	248.91	+	1.35	1.86E-02
Unclassified (UNCLASSIFIED)	26261	5013	9755.99	-	.51	0.00E00
' > expanded_orthogroups.goea
```

### contracted_orthogroups.goea
```{sh}
echo '
Analysis Type:	PANTHER Overrepresentation Test (Released 20220202)
Annotation Version and Release Date:	GO Ontology database DOI:  10.5281/zenodo.6399963 Released 2022-03-22
Analyzed List:	upload_1 (Oryza sativa)
Reference List:	Oryza sativa (all genes in database)
Test Type:	FISHER
Correction:	BONFERRONI
Bonferroni count:	2108
GO biological process complete	Oryza sativa - REFLIST (43658)	upload_1 (1638)	upload_1 (expected)	upload_1 (over/under)	upload_1 (fold Enrichment)	upload_1 (P-value)
cellular manganese ion homeostasis (GO:0030026)	9	9	.34	+	26.65	7.89E-06
sulfation (GO:0051923)	27	26	1.01	+	25.67	2.16E-20
negative regulation of cell population proliferation (GO:0008285)	19	18	.71	+	25.25	1.98E-13
cellulose catabolic process (GO:0030245)	27	25	1.01	+	24.68	2.98E-19
beta-glucan catabolic process (GO:0051275)	27	25	1.01	+	24.68	2.98E-19
manganese ion homeostasis (GO:0055071)	10	9	.38	+	23.99	1.45E-05
monosaccharide transmembrane transport (GO:0015749)	29	25	1.09	+	22.98	9.78E-19
glucosamine-containing compound catabolic process (GO:1901072)	19	16	.71	+	22.44	3.56E-11
chitin catabolic process (GO:0006032)	19	16	.71	+	22.44	3.56E-11
chitin metabolic process (GO:0006030)	19	16	.71	+	22.44	3.56E-11
aminoglycan catabolic process (GO:0006026)	19	16	.71	+	22.44	3.56E-11
amino sugar catabolic process (GO:0046348)	19	16	.71	+	22.44	3.56E-11
polyketide biosynthetic process (GO:0030639)	31	26	1.16	+	22.35	2.36E-19
polyketide metabolic process (GO:0030638)	31	26	1.16	+	22.35	2.36E-19
intracellular sequestering of iron ion (GO:0006880)	11	9	.41	+	21.81	2.55E-05
iron ion transmembrane transport (GO:0034755)	15	12	.56	+	21.32	1.06E-07
nucleoside monophosphate phosphorylation (GO:0046940)	16	12	.60	+	19.99	1.80E-07
sequestering of iron ion (GO:0097577)	12	9	.45	+	19.99	4.33E-05
sequestering of metal ion (GO:0051238)	12	9	.45	+	19.99	4.33E-05
mitotic cell cycle phase transition (GO:0044772)	59	41	2.21	+	18.52	2.55E-29
cell cycle phase transition (GO:0044770)	59	41	2.21	+	18.52	2.55E-29
glucosamine-containing compound metabolic process (GO:1901071)	24	16	.90	+	17.77	4.64E-10
plasmodesmata-mediated intercellular transport (GO:0010497)	11	7	.41	+	16.96	3.77E-03
intercellular transport (GO:0010496)	11	7	.41	+	16.96	3.77E-03
regulation of cell population proliferation (GO:0042127)	30	18	1.13	+	15.99	5.60E-11
regulation of cyclin-dependent protein serine/threonine kinase activity (GO:0000079)	70	41	2.63	+	15.61	3.79E-27
regulation of cyclin-dependent protein kinase activity (GO:1904029)	72	41	2.70	+	15.18	8.75E-27
flavonoid biosynthetic process (GO:0009813)	46	26	1.73	+	15.06	3.16E-16
flavonoid metabolic process (GO:0009812)	46	26	1.73	+	15.06	3.16E-16
aminoglycan metabolic process (GO:0006022)	31	16	1.16	+	13.76	8.73E-09
cytokinin-activated signaling pathway (GO:0009736)	62	32	2.33	+	13.76	1.72E-19
cellular polysaccharide catabolic process (GO:0044247)	49	25	1.84	+	13.60	1.02E-14
cellular response to cytokinin stimulus (GO:0071368)	63	32	2.36	+	13.54	2.51E-19
regulation of protein serine/threonine kinase activity (GO:0071900)	81	41	3.04	+	13.49	2.98E-25
carbohydrate transmembrane transport (GO:0034219)	88	43	3.30	+	13.02	4.31E-26
response to cytokinin (GO:0009735)	67	32	2.51	+	12.73	1.07E-18
magnesium ion transmembrane transport (GO:1903830)	21	10	.79	+	12.69	1.74E-04
regulation of root development (GO:2000280)	15	7	.56	+	12.44	1.78E-02
glucan catabolic process (GO:0009251)	55	25	2.06	+	12.12	8.54E-14
amino sugar metabolic process (GO:0006040)	36	16	1.35	+	11.85	5.08E-08
regulation of protein kinase activity (GO:0045859)	95	41	3.56	+	11.50	3.80E-23
regulation of kinase activity (GO:0043549)	97	41	3.64	+	11.27	7.19E-23
magnesium ion transport (GO:0015693)	24	10	.90	+	11.11	4.65E-04
regulation of protein phosphorylation (GO:0001932)	101	41	3.79	+	10.82	2.48E-22
iron ion transport (GO:0006826)	30	12	1.13	+	10.66	4.09E-05
cell wall macromolecule catabolic process (GO:0016998)	40	16	1.50	+	10.66	1.78E-07
regulation of phosphorylation (GO:0042325)	103	41	3.86	+	10.61	4.54E-22
cellular iron ion homeostasis (GO:0006879)	23	9	.86	+	10.43	2.88E-03
carbohydrate transport (GO:0008643)	113	43	4.24	+	10.14	1.32E-22
regulation of transferase activity (GO:0051338)	111	41	4.16	+	9.84	4.55E-21
cellular carbohydrate catabolic process (GO:0044275)	70	25	2.63	+	9.52	7.76E-12
beta-glucan metabolic process (GO:0051273)	84	30	3.15	+	9.52	1.40E-14
cellulose metabolic process (GO:0030243)	73	26	2.74	+	9.49	2.31E-12
RNA modification (GO:0009451)	302	100	11.33	+	8.83	8.88E-51
carbohydrate derivative catabolic process (GO:1901136)	50	16	1.88	+	8.53	2.60E-06
manganese ion transport (GO:0006828)	29	9	1.09	+	8.27	1.38E-02
manganese ion transmembrane transport (GO:0071421)	29	9	1.09	+	8.27	1.38E-02
maintenance of location in cell (GO:0051651)	30	9	1.13	+	8.00	1.73E-02
sterol biosynthetic process (GO:0016126)	41	12	1.54	+	7.80	6.83E-04
phosphorelay signal transduction system (GO:0000160)	117	32	4.39	+	7.29	6.99E-13
secondary metabolite biosynthetic process (GO:0044550)	97	26	3.64	+	7.14	6.07E-10
regulation of phosphate metabolic process (GO:0019220)	154	41	5.78	+	7.10	1.14E-16
regulation of phosphorus metabolic process (GO:0051174)	154	41	5.78	+	7.10	1.14E-16
regulation of protein modification process (GO:0031399)	168	41	6.30	+	6.50	1.66E-15
maturation of LSU-rRNA (GO:0000470)	47	11	1.76	+	6.24	1.34E-02
mitotic cell cycle process (GO:1903047)	184	41	6.90	+	5.94	2.68E-14
polysaccharide catabolic process (GO:0000272)	188	41	7.05	+	5.81	5.15E-14
cell division (GO:0051301)	196	42	7.35	+	5.71	3.65E-14
regulation of cell cycle (GO:0051726)	211	41	7.92	+	5.18	1.68E-12
mitotic cell cycle (GO:0000278)	222	43	8.33	+	5.16	3.70E-13
steroid metabolic process (GO:0008202)	102	17	3.83	+	4.44	2.99E-03
recognition of pollen (GO:0048544)	108	17	4.05	+	4.20	5.98E-03
cell recognition (GO:0008037)	108	17	4.05	+	4.20	5.98E-03
cellular glucan metabolic process (GO:0006073)	196	30	7.35	+	4.08	2.06E-06
pollen-pistil interaction (GO:0009875)	112	17	4.20	+	4.05	9.25E-03
glucan metabolic process (GO:0044042)	202	30	7.58	+	3.96	3.87E-06
secondary metabolic process (GO:0019748)	176	26	6.60	+	3.94	5.01E-05
carbohydrate catabolic process (GO:0016052)	293	41	10.99	+	3.73	2.39E-08
cell cycle process (GO:0022402)	325	43	12.19	+	3.53	3.89E-08
multi-multicellular organism process (GO:0044706)	139	18	5.22	+	3.45	3.51E-02
pollination (GO:0009856)	139	18	5.22	+	3.45	3.51E-02
defense response to other organism (GO:0098542)	275	35	10.32	+	3.39	7.06E-06
nucleic acid phosphodiester bond hydrolysis (GO:0090305)	385	47	14.44	+	3.25	5.61E-08
polysaccharide metabolic process (GO:0005976)	395	47	14.82	+	3.17	1.23E-07
cellular polysaccharide metabolic process (GO:0044264)	259	30	9.72	+	3.09	5.80E-04
cell cycle (GO:0007049)	399	46	14.97	+	3.07	5.02E-07
response to external biotic stimulus (GO:0043207)	319	35	11.97	+	2.92	2.05E-04
response to other organism (GO:0051707)	319	35	11.97	+	2.92	2.05E-04
carbohydrate metabolic process (GO:0005975)	1039	112	38.98	+	2.87	6.71E-18
cell wall organization or biogenesis (GO:0071554)	415	44	15.57	+	2.83	1.25E-05
sulfur compound metabolic process (GO:0006790)	305	32	11.44	+	2.80	1.74E-03
protein phosphorylation (GO:0006468)	1515	155	56.84	+	2.73	1.55E-23
response to biotic stimulus (GO:0009607)	346	35	12.98	+	2.70	1.19E-03
biological process involved in interspecies interaction between organisms (GO:0044419)	353	35	13.24	+	2.64	1.82E-03
hormone-mediated signaling pathway (GO:0009755)	410	39	15.38	+	2.54	1.28E-03
phosphorylation (GO:0016310)	1841	173	69.07	+	2.50	8.96E-23
cellular response to hormone stimulus (GO:0032870)	418	39	15.68	+	2.49	1.89E-03
cellular response to endogenous stimulus (GO:0071495)	425	39	15.95	+	2.45	3.97E-03
regulation of molecular function (GO:0065009)	527	48	19.77	+	2.43	2.27E-04
intracellular signal transduction (GO:0035556)	352	32	13.21	+	2.42	2.84E-02
regulation of catalytic activity (GO:0050790)	518	47	19.43	+	2.42	3.39E-04
nucleic acid metabolic process (GO:0090304)	1995	181	74.85	+	2.42	2.58E-22
RNA metabolic process (GO:0016070)	1458	128	54.70	+	2.34	7.82E-14
cellular carbohydrate metabolic process (GO:0044262)	427	37	16.02	+	2.31	1.75E-02
regulation of cellular protein metabolic process (GO:0032268)	501	43	18.80	+	2.29	4.72E-03
response to external stimulus (GO:0009605)	455	38	17.07	+	2.23	3.12E-02
nucleobase-containing compound metabolic process (GO:0006139)	2375	198	89.11	+	2.22	1.83E-20
cellular response to organic substance (GO:0071310)	483	40	18.12	+	2.21	3.12E-02
macromolecule modification (GO:0043412)	3495	286	131.13	+	2.18	2.67E-30
regulation of protein metabolic process (GO:0051246)	527	43	19.77	+	2.17	1.53E-02
transmembrane transport (GO:0055085)	1365	105	51.21	+	2.05	1.10E-07
heterocycle metabolic process (GO:0046483)	2613	198	98.04	+	2.02	4.48E-16
organic cyclic compound metabolic process (GO:1901360)	2840	214	106.55	+	2.01	2.30E-17
phosphate-containing compound metabolic process (GO:0006796)	2581	189	96.84	+	1.95	8.47E-14
cellular aromatic compound metabolic process (GO:0006725)	2745	200	102.99	+	1.94	1.24E-14
phosphorus metabolic process (GO:0006793)	2618	189	98.22	+	1.92	3.37E-13
cellular nitrogen compound metabolic process (GO:0034641)	3203	224	120.17	+	1.86	9.30E-15
metabolic process (GO:0008152)	10316	670	387.04	+	1.73	2.10E-48
primary metabolic process (GO:0044238)	8914	576	334.44	+	1.72	2.25E-38
macromolecule metabolic process (GO:0043170)	6945	445	260.57	+	1.71	2.28E-26
organic substance metabolic process (GO:0071704)	9447	605	354.44	+	1.71	5.61E-40
biological_process (GO:0008150)	17397	1042	652.72	+	1.60	2.87E-77
nitrogen compound metabolic process (GO:0006807)	7413	437	278.13	+	1.57	1.52E-18
cellular metabolic process (GO:0044237)	8971	527	336.58	+	1.57	9.07E-24
protein modification process (GO:0036211)	3176	186	119.16	+	1.56	1.47E-05
cellular protein modification process (GO:0006464)	3176	186	119.16	+	1.56	1.47E-05
cellular process (GO:0009987)	12820	749	480.99	+	1.56	2.67E-39
cellular macromolecule metabolic process (GO:0044260)	5139	295	192.81	+	1.53	1.25E-09
cellular protein metabolic process (GO:0044267)	4088	219	153.38	+	1.43	5.07E-04
Unclassified (UNCLASSIFIED)	26261	596	985.28	-	.60	0.00E00
organic acid metabolic process (GO:0006082)	883	11	33.13	-	.33	3.69E-02
carboxylic acid metabolic process (GO:0019752)	840	10	31.52	-	.32	3.49E-02
oxoacid metabolic process (GO:0043436)	851	10	31.93	-	.31	2.58E-02
establishment of protein localization (GO:0045184)	697	7	26.15	-	.27	4.99E-02
protein localization (GO:0008104)	767	7	28.78	-	.24	7.57E-03
intracellular transport (GO:0046907)	671	6	25.18	-	.24	2.23E-02
establishment of localization in cell (GO:0051649)	703	6	26.38	-	.23	1.22E-02
nitrogen compound transport (GO:0071705)	1005	8	37.71	-	.21	2.46E-05
proteolysis (GO:0006508)	1219	9	45.74	-	.20	2.10E-07
positive regulation of biological process (GO:0048518)	687	5	25.78	-	.19	3.45E-03
proteolysis involved in cellular protein catabolic process (GO:0051603)	612	4	22.96	-	.17	7.22E-03
cellular protein catabolic process (GO:0044257)	613	4	23.00	-	.17	7.26E-03
protein catabolic process (GO:0030163)	622	4	23.34	-	.17	5.10E-03
positive regulation of macromolecule metabolic process (GO:0010604)	528	3	19.81	-	.15	2.12E-02
positive regulation of metabolic process (GO:0009893)	544	3	20.41	-	.15	1.01E-02
positive regulation of cellular metabolic process (GO:0031325)	515	2	19.32	-	.10	4.67E-03
positive regulation of nitrogen compound metabolic process (GO:0051173)	522	2	19.58	-	.10	3.15E-03
positive regulation of cellular process (GO:0048522)	578	2	21.69	-	.09	5.16E-04
negative regulation of cellular metabolic process (GO:0031324)	378	1	14.18	-	.07	4.71E-02
positive regulation of RNA metabolic process (GO:0051254)	386	1	14.48	-	.07	3.23E-02
positive regulation of nucleobase-containing compound metabolic process (GO:0045935)	399	1	14.97	-	.07	2.23E-02
vesicle-mediated transport (GO:0016192)	493	1	18.50	-	.05	7.64E-04
' > contracted_orthogroups.goea
```

### Findings: Stress response and xenobiotic metabolism genes are significantly expanded and significantly enriched
Among the top 20 significantly enriched GO terms are:
- xylan biosynthetic process (most enriched) - integral to the integrity of cell walls and plays a role in defence against herbivory and mechanical stress
- hydrogen peroxide metabolic/catabolic process
- xenobiotic transport
- reactive oxygen species metabolic process
- cellular response to chemical stress

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
MACSE=$2
ORTHOLOG=${f%.fasta*}
# Align the CDS across species
java -Xmx8G \
    -jar ${MACSE} \
    -prog alignSequences \
    -seq ${f} \
    -out_NT ${ORTHOLOG}.aligned.unsorted.cds.tmp \
    -out_AA ${ORTHOLOG}.aligned.unsorted.prot.tmp
# Convert stop codons and frameshifts as "---" for compatibility with downstream tools
java -Xmx8G \
    -jar ${MACSE} \
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
./parallel_align_cds.sh {} ${MACSE} \
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

5. **Additional**: Compute pairwise 4DTv (Output: ORTHOGROUPS_SINGLE_GENE.NT.4DTv)
```{sh}
### Compute the transversion rate among 4-fold degenerate sites (Output: ${ORTHOLOG}.NT.cds.4DTv.tmp)
time \
parallel \
julia calculate_4DTv.jl {1} {1}.4DTv.tmp \
    ::: $(ls *.NT.cds)
### Concatenate pairwise 4DTv among species across single-copy orthogroups
echo -e "ORTHOGROUP\tSPECIES_1\tSPECIES_2\tn4D_sites\tnTv4D_sites\t4DTv" > ORTHOGROUPS_SINGLE_GENE.${TYPE%.*}.4DTv
for f in $(ls *.cds.4DTv.tmp)
do
    # f=$(ls *.cds.4DTv.tmp | head -n10 | tail -n1)
    n=$(cat $f | wc -l)
    printf "${f%.${TYPE}*}\n%.0s" $(seq 1 $n) > col1.tmp
    sed -z "s/ /\t/g" $f | sed -z "s/:/\t/g"  > col2n.tmp
    paste col1.tmp col2n.tmp >> ORTHOGROUPS_SINGLE_GENE.${TYPE%.*}.4DTv
done
```

6. Clean-up
```{sh}
rm OG*.fasta
rm OG*.NT.cds
rm single_gene_list.*
rm dates.txt
rm *.tmp
```

## Assess whole genome duplication (WGD) events 
We will use the distribution of four-fold degenerate sites (4DTv) across multi-copy paralogs within genomes and across sing-copy gene orthologs between pairs of species
1. Prepare script to extract CDS, align, calculate 4DTv (pairwise), and divergence time (pairwise) in parallel
```{sh}
echo '#!/bin/bash
### NOTE: The file: dual_gene_list.geneNames contains the gene names of one species. The gene names are in the second column (space-delimited), and each gene is comma-space-delimited
j=$1
MACSE=$2
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
# Concatenate the sequences
if [ $(ls ${ORTHONAME}-*.fasta | wc -l) -gt 1 ]
then
    cat ${ORTHONAME}-*.fasta > ${ORTHONAME}.fasta
fi
# Align the CDS across species
java -Xmx8G \
    -jar ${MACSE} \
    -prog alignSequences \
    -seq ${ORTHONAME}.fasta \
    -out_NT ${ORTHONAME}.aligned.unsorted.cds.tmp \
    -out_AA ${ORTHONAME}.aligned.unsorted.prot.tmp
# Convert stop codons and frameshifts as "---" for compatibility with downstream tools
java -Xmx8G \
    -jar ${MACSE} \
    -prog exportAlignment \
    -align ${ORTHONAME}.aligned.unsorted.cds.tmp \
    -codonForFinalStop --- \
    -codonForInternalStop NNN \
    -codonForExternalFS --- \
    -codonForInternalFS --- \
    -out_NT ${ORTHONAME}.NT.cds \
    -out_AA ${ORTHONAME}.AA.prot
# Calculate 4DTv (i.e. the ratio of the number of 4-fold degenerate codons with transversion and the total number of 4-fold degenerate codons)
julia calculate_4DTv.jl ${ORTHONAME}.NT.cds ${ORTHONAME}.4DTv.tmp
##### # Calculate divergence time with KaKs_Calculator
##### grep "^>" ${ORTHONAME}.NT.cds | \
#####     sed "s/^>//g" | \
#####     sed -z "s/\n/:/g" | \
#####     sed -z "s/:$/\n/g" > ${ORTHONAME}.axt.tmp
##### grep -v "^>" ${ORTHONAME}.NT.cds >> ${ORTHONAME}.axt.tmp
##### KaKs_Calculator -i ${ORTHONAME}.axt.tmp \
#####                 -m MA \
#####                 -o ${ORTHONAME}.kaks.tmp
##### divergence_time=$(tail -n1 ${ORTHONAME}.kaks.tmp | cut -f16)
##### sed -i -z "s/\n/\t${divergence_time}\n/g" ${ORTHONAME}.4DTv.tmp
# Clean-up
rm ${ORTHONAME}*.fasta
rm ${ORTHONAME}.aligned.unsorted*.tmp ${ORTHONAME}.NT.cds ${ORTHONAME}.AA.prot
##### rm ${ORTHONAME}.axt.tmp ${ORTHONAME}.kaks.tmp
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
    awk -v col="$idx" '($col >= 2) && ($col <= 5)' $ORTHOUT | cut -f1 > dual_gene_list.grep
    ### Extract names of the genes of these dual-copy orthogroups
    grep -f dual_gene_list.grep ${DIR_ORTHOGROUPS}/Orthogroups/Orthogroups.tsv | cut -f1,${idx} > dual_gene_list.geneNames
    ### Extract CDS, and align in parallel
    parallel \
    ./parallel_extract_dual_gene_orthogroups.sh {} ${MACSE} \
    ::: $(seq 1 $(cat dual_gene_list.geneNames | wc -l))
    ### Concatenate 4DTv estimates
    cat *.4DTv.tmp > ${SPECIES}.4DTv
    ### Clean-up
    rm *.4DTv.tmp
    rm dual_gene_list.grep
    rm dual_gene_list.geneNames
done
### Clean-up
rm species_names.tmp
rm *.tmp
```

## Identify herbicide TSR and NTSR genes
1. Download protein sequences of genes from UniProt (https://www.uniprot.org) (Outputs: ${GENE}.faa)
```{sh}
###############################
### SET WORKING DIRECTORIES ###
###############################
DIR_ORTHOGROUP_SEQS=${DIR}/ORTHOGROUPS/OrthoFinder/Results_*/Orthogroup_Sequences
DIR_GENES=${DIR}/TSR_NTSR_GENES
mkdir $DIR_GENES
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

## TESTING KaKs_Calucator2.0 with sliding windows
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

## OR OR OR SIMPLY USE PAML::codeml on these small datasets, i.e. per TSR/NTSR gene per orthogroup
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
