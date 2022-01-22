#!/bin/bash

### Assembly error correction with pilon or quiver/arrow

ASSEMBLY_FASTA=$1

### INPUT:
### (1) Genome assembly in fasta format

### OUTPUTS:
### (1) 

### Pilon
## install
DIR=/data/Lolium_rigidum_ASSEMBLY/assembly_annotation_pipeline_tests_20210104/
cd $DIR
# (1) install dependencies
sudo apt-get update; sudo apt-get install build-essential software-properties-common -y;
sudo add-apt-repository ppa:ubuntu-toolchain-r/test -y; sudo apt update; 
sudo apt install gcc-snapshot -y; sudo apt update
sudo apt install gcc-8 g++-8 -y; 
sudo update-alternatives --install /usr/bin/gcc gcc /usr/bin/gcc-8 60 --slave /usr/bin/g++ g++ /usr/bin/g++-8
sudo apt install cmake
# (2) install HTSLIB
sudo -i
cd /usr/bin
wget https://github.com/samtools/htslib/releases/download/1.10/htslib-1.10.tar.bz2
tar -vxjf htslib-1.10.tar.bz2
cd htslib-1.10
make
logout
export PATH="$PATH:/usr/bin/htslib-1.10"
source ~/.profile
# (3) install kmc for suk
cd $DIR
git clone https://github.com/refresh-bio/KMC
cd KMC/
make
export PATH="$PATH:$(pwd)"
source ~/.profile
# (4) install hypo
cd $DIR
git clone --recursive https://github.com/kensung-lab/hypo hypo
cd hypo
chmod +x install_deps.sh
./install_deps.sh
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Release -Doptimise_for_native=ON ..
make -j 8
# (5) install minimap2
git clone https://github.com/lh3/minimap2
cd minimap2 && make
cd -
# (6) install bwa
sudo apt install bwa





### Genome assembly
# ASSEMBLY=/data/Lolium_rigidum_ASSEMBLY/assembly_annotation_pipeline_tests_20210104/ASSEMBLY/Lori_mw/Lori_mw.raw.fa
ASSEMBLY=/data/Lolium_rigidum_ASSEMBLY/assembly_annotation_pipeline_tests_20210104/ASSEMBLY/Lori_mw/Lori_mw.fasta
### Navigate to working directory
ASSEMBLY_DIR=$(dirname $ASSEMBLY)
cd ${ASSEMBLY_DIR}
### Setup reference genome assembly name
ASSEMBLY_NAME=$(echo ${ASSEMBLY_DIR} | rev | cut -d'/' -f1 | rev)
### Index the assmebly to be used as reference for read alignment
bwa index -p ${ASSEMBLY_NAME} -a bwtsw ${ASSEMBLY}

## Align reads
# (1) illumina reads

# (2) minion reads
MINION_READS=/data/Lolium_rigidum_ASSEMBLY/assembly_annotation_pipeline_tests_20210104/FASTQ/MINION/minion-filtered-trimmed.fastq.gz
MINIMAP=/data/Lolium_rigidum_ASSEMBLY/assembly_annotation_pipeline_tests_20210104/minimap2/minimap2
time \
${MINIMAP} -ax map-ont \
    ${ASSEMBLY} \
    ${MINION_READS} \
    > ${ASSEMBLY_NAME}-MINION_READS.sam
MAPQ=20
time \
samtools view -q ${MAPQ} -b ${ASSEMBLY_NAME}-MINION_READS.sam | \
    samtools sort > ${ASSEMBLY_NAME}-MINION_READS.bam
time \
samtools index ${ASSEMBLY_NAME}-MINION_READS.bam



### HYPO error correction
