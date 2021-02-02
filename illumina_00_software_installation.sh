#!/bin/bash
################################################
### Software installation for Illumina reads ###
################################################

### Parameters:
DIR=/data/Lolium_rigidum_ASSEMBLY/assembly_annotation_pipeline_tests_20210104/

### Navigate to working directory
cd $DIR


### Download FastQC
wget http://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.7.zip
unzip fastqc*
cd FastQC
chmod +x fastqc
cd -

### Download SPAdes (includes BayesHammer)
wget http://cab.spbu.ru/files/release3.14.1/SPAdes-3.14.1-Linux.tar.gz
tar -xzf SPAdes-3.14.1-Linux.tar.gz
cd SPAdes-3.14.1-Linux/bin/
python3 spades.py -h
cd -

### Install Platanus-allee
wget -O Platanus_allee_v2.2.2_Linux_x86_64.tgz http://platanus.bio.titech.ac.jp/?ddownload=431
tar -xzf Platanus_allee_v2.2.2_Linux_x86_64.tgz
cd Platanus_allee_v2.2.2_Linux_x86_64/
./platanus_allee
cd -

### Install Biopython, and pandas
pip3 install biopython
pip3 install pandas

### Install php and mummer for MACS meta-assembler
sudo apt install -y php
wget https://sourceforge.net/projects/mummer/files/mummer/3.23/MUMmer3.23.tar.gz
tar -xzf MUMmer3.23.tar.gz
cd MUMmer3.23/
make install
cd -
git clone https://github.com/bioinfomaticsCSU/MAC.git

### Clean-up
rm fastqc_v0.11.7.zip
rm SPAdes-3.14.1-Linux.tar.gz
rm Platanus_allee_v2.2.2_Linux_x86_64.tgz
rm MUMmer3.23.tar.gz