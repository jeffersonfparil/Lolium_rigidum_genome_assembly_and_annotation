#!/bin/bash
##############################################
### Software installation for MinION reads ###
##############################################

### Parameters:
DIR=/data/Lolium_rigidum_ASSEMBLY/assembly_annotation_pipeline_tests_20210104/

### Navigate to working directory
cd $DIR

### Install guppy on ubuntu
sudo apt-get update
sudo apt-get install wget lsb-release
PLATFORM=$(lsb_release -cs)
wget -O- https://mirror.oxfordnanoportal.com/apt/ont-repo.pub | sudo apt-key add -
echo "deb http://mirror.oxfordnanoportal.com/apt ${PLATFORM}-stable non-free" | sudo tee /etc/apt/sources.list.d/nanoporetech.sources.list
sudo apt-get update
sudo apt-get install ont-guppy ### cpu+gpu version; or sudo apt-get install ont-guppy-cpu for the cpu-only version

### Install nanofilt using pip3 (assumes python3)
pip3 install nanofilt
pip3 install nanofilt --upgrade

### Install porechop
git clone https://github.com/rrwick/Porechop.git
cd Porechop
python3 setup.py install --user
porechop -h
cd -

### Install wtdbg2
git clone https://github.com/ruanjue/wtdbg2
cd wtdbg2 && make
#quick start with wtdbg2.pl
./wtdbg2.pl -h
cd -

### Clean-up
rm canu-1.9.Linux-amd64.tar.xz
