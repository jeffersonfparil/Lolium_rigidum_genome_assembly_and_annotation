#!/bin/bash
### Install Ubuntu system pre-requisites
sudo apt-get install zlib1g zlib1g-dev libexpat1-dev libpng-dev libgd-dev build-essential git software-properties-common python make

### Download and install JBrowse
curl -L -O https://github.com/GMOD/jbrowse/releases/download/1.16.11-release/JBrowse-1.16.11.zip
unzip JBrowse-1.16.11.zip
sudo mv JBrowse-1.16.11 /var/www/html/jbrowse
cd /var/www/html
sudo chown `whoami` jbrowse
cd jbrowse
./setup.sh

### Set-up the genome and annotation track
cd /var/www/html/jbrowse/
mkdir data/
cd data/
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/022/539/505/GCF_022539505.1_APGP_CSIRO_Lrig_0.1/GCF_022539505.1_APGP_CSIRO_Lrig_0.1_genomic.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/022/539/505/GCF_022539505.1_APGP_CSIRO_Lrig_0.1/GCF_022539505.1_APGP_CSIRO_Lrig_0.1_genomic.gff.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/022/539/505/GCF_022539505.1_APGP_CSIRO_Lrig_0.1/GCF_022539505.1_APGP_CSIRO_Lrig_0.1_cds_from_genomic.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/022/539/505/GCF_022539505.1_APGP_CSIRO_Lrig_0.1/GCF_022539505.1_APGP_CSIRO_Lrig_0.1_protein.faa.gz
gunzip -c GCF_022539505.1_APGP_CSIRO_Lrig_0.1_genomic.fna.gz > Lolium_rigidum.fasta
gunzip -c GCF_022539505.1_APGP_CSIRO_Lrig_0.1_genomic.gff.gz > Lolium_rigidum.gff
gunzip -c GCF_022539505.1_APGP_CSIRO_Lrig_0.1_cds_from_genomic.fna.gz > Lolium_rigidum.cds
gunzip -c GCF_022539505.1_APGP_CSIRO_Lrig_0.1_protein.faa.gz > Lolium_rigidum.faa
rm GCF_022539505.1_APGP_CSIRO_Lrig_0.1_*
### NOTE: Symbolic links inside data/ also works if disk space is limited in / and heaps of space in some mounted drive ;-P
