#!/bin/bash
#######################################
### Genome annotation with Augustus ###
#######################################

### Inputs:
### (1) Genome assemblies in fasta format (Lori_i1, Lori_i2, ... *.fasta or *.fa)

### Outputs:
### (1) Annotations per genome assembly per gene list (${assembly}.${species_gene_list}.gff)
### (2) 
### (3) 

### Parameters:
INPUT_DIR=/data/Lolium_rigidum_ASSEMBLY/assembly_annotation_pipeline_tests_20210104/ASSEMBLY
OUTPUT_DIR=/data/Lolium_rigidum_ASSEMBLY/assembly_annotation_pipeline_tests_20210104/ANNOTATION
AUGUSTUS=/data/Lolium_rigidum_ASSEMBLY/assembly_annotation_pipeline_tests_20210104/Augustus/bin/augustus

### Navigate to working directory
cd $OUTPUT_DIR

### List of genome assemblies
ASSEMBLIES=$(ls ${INPUT_DIR}/ | grep "Lori_" | sed 's/.fasta//g')

### List of gene lists we will be using to find homologs in the assemblies
GENE_LISTS=$(echo "rice maize arabidopsis")

### Annotation with rice, maize, and arabidopsis genes
### (1) split the genome assemblies by scaffold
echo 'from Bio import SeqIO
import pandas as pd
import sys
import os
fname_assembly = sys.argv[1]
# fname_assembly = "/data/Lolium_rigidum_ASSEMBLY/assembly_annotation_pipeline_tests_20210104/ASSEMBLY/Lori_m1/Lori_m1.raw.fa"
id_assembly = os.path.basename(fname_assembly).split(".")[0]
with open(fname_assembly, "rU") as genome:
    for scaffold in SeqIO.parse(genome, "fasta"):
        seq_string = SeqIO.FastaIO.as_fasta_2line(scaffold)
        f = open(id_assembly + "." + scaffold.id + ".fa", "w")
        f.write(seq_string)
        f.close()
' > split_assembly_by_scaffold.py
time \
parallel python3 split_assembly_by_scaffold.py {} ::: $(ls ${INPUT_DIR}/Lori_*.fasta)
### (2) run Augustus in parallel per scaffold per assembly per species genes
echo '#!/bin/bash
AUGUSTUS=$1
FASTA=$2
SPECIES=$3
# AUGUSTUS=/data/Lolium_rigidum_ASSEMBLY/assembly_annotation_pipeline_tests_20210104/Augustus/bin/augustus
# FASTA=ctg1.fa
# SPECIES=rice #SPECIES=maize #SPECIES=arabidopsis
${AUGUSTUS} \
    --species=${SPECIES} \
    --genemodel=partial \
    ${FASTA} \
    > ${FASTA}.${SPECIES}.gff
' > augustus_parallel.sh
chmod +x augustus_parallel.sh
time \
parallel ./augustus_parallel.sh ${AUGUSTUS} {1} {2} \
    ::: $(ls Lori_*.fa) \
    ::: $(echo ${GENE_LISTS} | cut -d' ' -f1-3)

### Merge across scaffolds per genome assembly per gene list
echo '#!/bin/bash
assembly=$1
species_gene_list=$2
# ### test
# assembly=Lori_hw
# species_gene_list=rice
f1=$(ls ${assembly}.*.${species_gene_list}.gff | head -n1)
line_number=$(echo $(grep -n "# ----- prediction" ${f1} | cut -d: -f1) - 1 | bc)
if [ ${line_number} -eq "-1" ]
then
    touch ${assembly}.${species_gene_list}.gff
else
    head -${line_number} ${f1} > ${assembly}.${species_gene_list}.gff
fi
for f in $(ls ${assembly}.*.${species_gene_list}.gff)
do
    line_number=$(grep -n "# ----- prediction" ${f} | cut -d: -f1)
    n_match=$(grep -n "# ----- prediction" ${f} | wc -l)
    # echo $line_number
    if [ ${n_match} -ne 0 ]
    then
        tail -n+${line_number} ${f} >> ${assembly}.${species_gene_list}.gff
    fi
done
' > merge_gff_parallel.sh
chmod +x merge_gff_parallel.sh
time \
parallel ./merge_gff_parallel.sh {1} {2} \
    ::: ${ASSEMBLIES} \
    ::: ${GENE_LISTS}

### Move per scaffold and per gene list sequence and annotations into a separate folder
mkdir INDIVIDUAL_AUGUSTUS_PER_SCAFFOLD_OUTPUT/
ls | grep "Lori_h" | grep "\.fa$" | xargs -I {} mv {} INDIVIDUAL_AUGUSTUS_PER_SCAFFOLD_OUTPUT/
ls | grep "Lori_h" | grep '.fa.' | grep "gff$" | xargs -I {} mv {} INDIVIDUAL_AUGUSTUS_PER_SCAFFOLD_OUTPUT/

### Clean-up
rm split_assembly_by_scaffold.py
rm augustus_parallel.sh
rm merge_gff_parallel.sh



### TESTING BRAKER2

### install dependencies
DIR=/data/Lolium_rigidum_ASSEMBLY/assembly_annotation_pipeline_tests_20210104
cd $DIR
### (1) GeneMark-EX
wget http://topaz.gatech.edu/GeneMark/tmp/GMtool_mNfY6/gmes_linux_64.tar.gz ### dowload GeneMar-EX
wget http://topaz.gatech.edu/GeneMark/tmp/GMtool_AWk_2/gm_key_64.gz ### Download key
tar -xvzf gmes_linux_64.tar.gz ### decompress the software
gunzip -d gm_key_64.gz; mv gm_key_64 .gm_key ### decompress and rename the key
cd gmes_linux_64/
sudo apt install -y cpanminus ### install perl module installer cpanm
sudo cpanm Hash::Merge ### install Hash::Merge perl module
sudo cpanm MCE::Mutex ### install MCE::Mutex perl module
sudo cpanm Math::Utils ### install MCE::Mutex perl module
./check_install.bash ### check installation of GeneMark-EX
echo "export GENEMARK_PATH=${DIR}/gmes_linux_64/" >> ~/.bashrc ### add to path
cd -

### (2) Augustus
git clone https://github.com/Gaius-Augustus/Augustus.git
# install required packages
sudo apt update
sudo apt install -y build-essential wget git autoconf
# install dependencies for AUGUSTUS comparative gene prediction mode (CGP)
sudo apt install -y libgsl-dev libboost-all-dev libsuitesparse-dev liblpsolve55-dev
sudo apt install -y libsqlite3-dev libmysql++-dev
# install dependencies for the optional support of gzip compressed input files
sudo apt install -y libboost-iostreams-dev zlib1g-dev
# install dependencies for bam2hints and filterBam 
sudo apt install -y libbamtools-dev
# install additional dependencies for bam2wig
sudo apt install -y samtools libhts-dev
# install additional dependencies for homGeneMapping and utrrnaseq
sudo apt install -y libboost-all-dev
# install additional dependencies for scripts
sudo apt install -y cdbfasta diamond-aligner libfile-which-perl libparallel-forkmanager-perl libyaml-perl libdbd-mysql-perl
sudo apt install -y --no-install-recommends python3-biopython
# # install HTSLib from source (even after install all of the above making Augustus still spit out error because HTSlib is not installed)
# wget https://github.com/samtools/htslib/releases/download/1.3.2/htslib-1.3.2.tar.bz2 -O htslib.tar.bz2
# tar -xjvf htslib.tar.bz2
# cd htslib-*/
# make
# sudo make install
# cd -
# make install Augustus but first comment-out bam2wig make which causes problems with HTSlib not being located when it is actually installed
cd Augustus/
### if you encounter any compilation errors at sudo make install then execute the commented-out script below:
# cp auxprogs/Makefile auxprogs/Makefile.bk
# sed -i 's/	cd bam2wig; make/	#cd bam2wig; make/g' auxprogs/Makefile
# make clean
make
sudo make install
# test
make unit_test
### add Augustus to path
echo "export AUGUSTUS_CONFIG_PATH=${DIR}/Augustus/config/" >> ~/.bashrc
echo "export AUGUSTUS_BIN_PATH=${DIR}/Augustus/bin/" >> ~/.bashrc
echo "export AUGUSTUS_SCRIPTS_PATH=${DIR}/Augustus/scripts/" >> ~/.bashrc
cd -

### (3) Bamtools
sudo apt install -y bamtools

### (4) ncbi-blast+
sudo apt install -y ncbi-blast+

### (5) ProtHint
git clone https://github.com/gatech-genemark/ProtHint.git
cd ProtHint/
### add ProtHint to path
echo "export PROTHINT_PATH=${DIR}/ProtHint/bin/"  >> ~/.bashrc
