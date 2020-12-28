#!/bin/bash
################################################
### Software installation for Illumina reads ###
################################################

### Parameters:
DIR=/data/Lolium_rigidum_ASSEMBLY/jeff_lolium/

### Navigate to working directory
cd $DIR

### Download SPAdes (includes BayesHammer)
wget http://cab.spbu.ru/files/release3.14.1/SPAdes-3.14.1-Linux.tar.gz
tar -xzf SPAdes-3.14.1-Linux.tar.gz
cd SPAdes-3.14.1-Linux/bin/
python3 spades.py -h

### Install Platanus-allee
wget -O Platanus_allee_v2.2.2_Linux_x86_64.tgz http://platanus.bio.titech.ac.jp/?ddownload=431
tar -xzf Platanus_allee_v2.2.2_Linux_x86_64.tgz
cd Platanus_allee_v2.2.2_Linux_x86_64/
./platanus_allee
