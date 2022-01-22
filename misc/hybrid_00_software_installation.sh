#!/bin/bash
###################################################
### Software installation for hybrid assemblies ###
###################################################

### Parameters:
DIR=/data/Lolium_rigidum_ASSEMBLY/assembly_annotation_pipeline_tests_20210104/

### Navigate to working directory
cd $DIR

### Install HASLR
git clone https://github.com/vpc-ccg/haslr.git
cd haslr
make
cd -

### Install wengan
wget https://github.com/adigenova/wengan/releases/download/v0.2/wengan-v0.2-bin-Linux.tar.gz
tar -xvzf wengan-v0.2-bin-Linux.tar.gz
rm wengan-v0.2-bin-Linux.tar.gz
cd wengan-v0.2-bin-Linux/
chmod +x wengan.pl
cd -
