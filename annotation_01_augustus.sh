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

# ### List of genome assemblies
# ASSEMBLIES=$(ls ${INPUT_DIR}/ | grep "Lori_" | sed 's/.fasta//g')

# ### List of gene lists we will be using to find homologs in the assemblies
# GENE_LISTS=$(echo "rice maize arabidopsis")

# ### Annotation with rice, maize, and arabidopsis genes
# ### (1) split the genome assemblies by scaffold
# echo 'from Bio import SeqIO
# import pandas as pd
# import sys
# import os
# fname_assembly = sys.argv[1]
# # fname_assembly = "/data/Lolium_rigidum_ASSEMBLY/assembly_annotation_pipeline_tests_20210104/ASSEMBLY/Lori_m1/Lori_m1.raw.fa"
# id_assembly = os.path.basename(fname_assembly).split(".")[0]
# with open(fname_assembly, "rU") as genome:
#     for scaffold in SeqIO.parse(genome, "fasta"):
#         seq_string = SeqIO.FastaIO.as_fasta_2line(scaffold)
#         f = open(id_assembly + "." + scaffold.id + ".fa", "w")
#         f.write(seq_string)
#         f.close()
# ' > split_assembly_by_scaffold.py
# time \
# parallel python3 split_assembly_by_scaffold.py {} ::: $(ls ${INPUT_DIR}/Lori_*.fasta)
# ### (2) run Augustus in parallel per scaffold per assembly per species genes
# echo '#!/bin/bash
# AUGUSTUS=$1
# FASTA=$2
# SPECIES=$3
# # AUGUSTUS=/data/Lolium_rigidum_ASSEMBLY/assembly_annotation_pipeline_tests_20210104/Augustus/bin/augustus
# # FASTA=ctg1.fa
# # SPECIES=rice #SPECIES=maize #SPECIES=arabidopsis
# ${AUGUSTUS} \
#     --species=${SPECIES} \
#     --genemodel=partial \
#     ${FASTA} \
#     > ${FASTA}.${SPECIES}.gff
# ' > augustus_parallel.sh
# chmod +x augustus_parallel.sh
# time \
# parallel ./augustus_parallel.sh ${AUGUSTUS} {1} {2} \
#     ::: $(ls Lori_*.fa) \
#     ::: $(echo ${GENE_LISTS} | cut -d' ' -f1-3)

# ### Merge across scaffolds per genome assembly per gene list
# echo '#!/bin/bash
# assembly=$1
# species_gene_list=$2
# # ### test
# # assembly=Lori_hw
# # species_gene_list=rice
# f1=$(ls ${assembly}.*.${species_gene_list}.gff | head -n1)
# line_number=$(echo $(grep -n "# ----- prediction" ${f1} | cut -d: -f1) - 1 | bc)
# if [ ${line_number} -eq "-1" ]
# then
#     touch ${assembly}.${species_gene_list}.gff
# else
#     head -${line_number} ${f1} > ${assembly}.${species_gene_list}.gff
# fi
# for f in $(ls ${assembly}.*.${species_gene_list}.gff)
# do
#     line_number=$(grep -n "# ----- prediction" ${f} | cut -d: -f1)
#     n_match=$(grep -n "# ----- prediction" ${f} | wc -l)
#     # echo $line_number
#     if [ ${n_match} -ne 0 ]
#     then
#         tail -n+${line_number} ${f} >> ${assembly}.${species_gene_list}.gff
#     fi
# done
# ' > merge_gff_parallel.sh
# chmod +x merge_gff_parallel.sh
# time \
# parallel ./merge_gff_parallel.sh {1} {2} \
#     ::: ${ASSEMBLIES} \
#     ::: ${GENE_LISTS}

# ### Move per scaffold and per gene list sequence and annotations into a separate folder
# mkdir INDIVIDUAL_AUGUSTUS_PER_SCAFFOLD_OUTPUT/
# ls | grep "Lori_h" | grep "\.fa$" | xargs -I {} mv {} INDIVIDUAL_AUGUSTUS_PER_SCAFFOLD_OUTPUT/
# ls | grep "Lori_h" | grep '.fa.' | grep "gff$" | xargs -I {} mv {} INDIVIDUAL_AUGUSTUS_PER_SCAFFOLD_OUTPUT/

# ### Clean-up
# rm split_assembly_by_scaffold.py
# rm augustus_parallel.sh
# rm merge_gff_parallel.sh


#######################
###                 ###
### TESTING BRAKER2 ###
###                 ###
#######################

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

### (5) Download/install ProtHint
sudo cpanm threads ### install threads perl module
sudo cpanm YAML ### install YAML perl module
sudo cpanm Thread::Queue ### install Thread::Queue perl module
git clone https://github.com/gatech-genemark/ProtHint.git
cd ProtHint/
### add ProtHint to path
echo "export PROTHINT_PATH=${DIR}/ProtHint/bin/"  >> ~/.bashrc
# ### TEST
# cd /data/Lolium_rigidum_ASSEMBLY/assembly_annotation_pipeline_tests_20210104/ASSEMBLY/
# ### download viridiplantae protein sequences
# wget https://v100.orthodb.org/download/odb10_plants_fasta.tar.gz
# ### ProtHint test
# PROTHINT=/data/Lolium_rigidum_ASSEMBLY/assembly_annotation_pipeline_tests_20210104/ProtHint/bin/prothint.py
# ${PROTHINT} Lori_hh.fasta odb10_plants_fasta.tar.gz
cd -

### (6) Download/install STAR executable
git clone https://github.com/alexdobin/STAR.git
STAR/bin/Linux_x86_64_static/STAR -h

### (7) Download/install GeneMark
wget http://topaz.gatech.edu/GeneMark/tmp/GMtool_mNfY6/gmes_linux_64.tar.gz ### dowload GeneMar-EX
tar -xvzf gmes_linux_64.tar.gz; rm gmes_linux_64.tar.gz ### decompress the software
cd gmes_linux_64/
wget http://topaz.gatech.edu/GeneMark/tmp/GMtool_AWk_2/gm_key_64.gz ### Download key
gunzip -d gm_key_64.gz; mv gm_key_64 .gm_key ### decompress and rename the key
sudo apt install -y cpanminus ### install perl module installer cpanm
sudo cpanm Hash::Merge ### install Hash::Merge perl module
sudo cpanm MCE::Mutex ### install MCE::Mutex perl module
sudo cpanm Math::Utils ### install MCE::Mutex perl module
./check_install.bash ### check installation of GeneMark-EX
cp -R * /data/Lolium_rigidum_ASSEMBLY/assembly_annotation_pipeline_tests_20210104/ProtHint/dependencies/GeneMarkES/
cp .gm_key /data/Lolium_rigidum_ASSEMBLY/assembly_annotation_pipeline_tests_20210104/ProtHint/dependencies/GeneMarkES/
cd -

### (8) Download/install BRAKER2
git clone https://github.com/Gaius-Augustus/BRAKER.git
sudo cpanm File::HomeDir
sudo cpanm Scalar::Util::Numeric


###@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
###@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
###@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

### step 0: navigate to working directory
cd /data/Lolium_rigidum_ASSEMBLY/assembly_annotation_pipeline_tests_20210104/ASSEMBLY

###step 1: Download OrthoDB protein sequencies and gene list (to identify the proteins) and concatenate the individual protein fastas into a single fasta file
# wget https://v100.orthodb.org/download/odb10_plants_fasta.tar.gz
# tar -xvzf odb10_plants_fasta.tar.gz
# cat plants/Rawdata/*.fs > odb10_plants.fasta
# rm -R plants/ odb10_plants_fasta.tar.gz
### or use all proteins
wget https://v101.orthodb.org/download/odb10v1_all_fasta.tab.gz
gunzip -d odb10v1_all_fasta.tab.gz
mv odb10v1_all_fasta.tab odb10v1_all.fasta

wget https://v101.orthodb.org/download/odb10v1_genes.tab.gz
gunzip -d odb10v1_genes.tab.gz


# ### step 2: ProtHint (generate gff annotations using Viridiplantae protein sequences from OrthoDB: )
# PROTHINT=/data/Lolium_rigidum_ASSEMBLY/assembly_annotation_pipeline_tests_20210104/ProtHint/bin/prothint.py
# time ${PROTHINT} Lori_hh.fasta odb10_plants.fasta
# ### outputs:
# ### (1) prothint.gff - all hints including introns, starts and stops
# ### (2) evidence.gff  high-confidence subset of prothint.gff suitable for GeneMark-EP
# ### (3) prothint-augustus.gff - BRAKER- and AUGUSTUS-compatible format

# ### step 3: GeneMark-EP+ (generate gtf annotations)
# GENEMARK_EPP=/data/Lolium_rigidum_ASSEMBLY/assembly_annotation_pipeline_tests_20210104/gmes_linux_64/gmes_petap.pl
# time \
# ${GENEMARK_EPP} \
#     --EP prothint.gff \
#     --evidence evidence.gff \
#     --seq Lori_hh.fasta \
#     --soft_mask 1000 \
#     --cores 12 \
#     --verbose
# ### output:
# ### (1) genemark.gtf
# sed 's/ from/_from/g' genemark.gtf | sed 's/ to/_to/g' > genemark_col1_fixed.gtf

### step 4: STAR RNAseq alignment (generate transcript alignment bam files)
STAR=/data/Lolium_rigidum_ASSEMBLY/assembly_annotation_pipeline_tests_20210104/STAR/bin/Linux_x86_64_static/STAR
### prepare reference genome indices
time \
${STAR} --runMode genomeGenerate \
        --genomeDir /data/Lolium_rigidum_ASSEMBLY/assembly_annotation_pipeline_tests_20210104/ASSEMBLY/ \
        --genomeFastaFiles /data/Lolium_rigidum_ASSEMBLY/assembly_annotation_pipeline_tests_20210104/ASSEMBLY/Lori_hh.fasta \
        --genomeSAindexNbases 13 \
        --runThreadN 12
### align
time \
${STAR} --genomeDir /data/Lolium_rigidum_ASSEMBLY/assembly_annotation_pipeline_tests_20210104/ASSEMBLY/ \
        --readFilesIn \
            /data/Lolium_rigidum_ASSEMBLY/assembly_annotation_pipeline_tests_20210104/FASTQ/ILLUMINA/RNAseq/INFLO-1_combined_R1.fastq.gz \
            /data/Lolium_rigidum_ASSEMBLY/assembly_annotation_pipeline_tests_20210104/FASTQ/ILLUMINA/RNAseq/INFLO-1_combined_R2.fastq.gz \
        --readFilesCommand zcat \
        --runThreadN 12 \
        --outFileNamePrefix /data/Lolium_rigidum_ASSEMBLY/assembly_annotation_pipeline_tests_20210104/ASSEMBLY/Lori_hh_RNAseq
### samtools sort and compression
MAPQ=20
time \
samtools view -q ${MAPQ} -b Lori_hh_RNAseqAligned.out.sam | samtools sort > Lori_hh_RNAseq.bam

### step 5: BRAKER pipeline D
BRAKER2=/data/Lolium_rigidum_ASSEMBLY/assembly_annotation_pipeline_tests_20210104/BRAKER/scripts/braker.pl
time \
${BRAKER2} --genome Lori_hh.fasta \
           --prot_seq odb10v1_all.fasta \
           --bam Lori_hh_RNAseq.bam \
           --etpmode \
           --softmasking \
           --cores 12


