#!/bin/bash
################################################
### Software installation for Illumina reads ###
################################################

### Parameters:
DIR=/data/Lolium_rigidum_ASSEMBLY/assembly_annotation_pipeline_tests_20210104/

### Navigate to working directory
cd $DIR


### Download FastQC
wget http://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.9.zip
unzip fastqc*
cd FastQC
chmod +x fastqc
cd -

### Install Minia specifically gatb-minia-pipeline
###### install dependency, i.e. bwa
sudo apt install -y bwa
###### clone and make test
git clone https://github.com/GATB/gatb-minia-pipeline.git
cd gatb-minia-pipeline/
make test ### may get python2 errors just remove duplicate libraries and reinstall and/or install missing libraries with pip
# pip3 install scipy mathstats pysam

### Install BESST for scaffolding via conda
###### install Miniconda
# wget https://repo.anaconda.com/miniconda/Miniconda2-latest-Linux-x86_64.sh
# chmod +x Miniconda2-latest-Linux-x86_64.sh ### respond appropriately to the prompts
# rm Miniconda2-latest-Linux-x86_64.sh
###### activate conda environment and install BESST via pip
source activate 
sudo pip install BESST

### Clean-up
rm fastqc_v0.11.7.zip
