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

### minia specifically gatb-minia-pipeline
###### install dependency, i.e. bwa
sudo apt install -y bwa
###### clone and make test
git clone https://github.com/GATB/gatb-minia-pipeline.git
cd gatb-minia-pipeline/
make test ### may get python 2 errors just remove duplicate libraries and reinstall and/or install missing libraries with pip

### Clean-up
rm fastqc_v0.11.7.zip
