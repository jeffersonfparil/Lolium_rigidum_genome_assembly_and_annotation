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
