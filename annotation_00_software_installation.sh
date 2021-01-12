#!/bin/bash
###################################################
### Software installation for genome annotation ###
###################################################

### Parameters:
DIR=/data/Lolium_rigidum_ASSEMBLY/assembly_annotation_pipeline_tests_20210104/

### Navigate to working directory
cd $DIR

### Install Augustus
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
cp auxprogs/Makefile auxprogs/Makefile.bk
sed -i 's/	cd bam2wig; make/	#cd bam2wig; make/g' auxprogs/Makefile
make ### execute 'make clean' first to restart, then 'make' again
sudo make install
# test
make unit_test

