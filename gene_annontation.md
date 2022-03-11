# Gaisu-Augustus/BRAKER pipeline D

![](misc/braker2_pipeline_D.png)

## Working directory
```{sh}
DIR=/data/Lolium_rigidum_ASSEMBLY/ANNOTATION
```

## Install dependencies

1. Perl modules
```{sh}
sudo cpan Hash::Merge MCE::Mutex Math::Utils Parallel:ForkManager ### GeneMark-EX dependencies
sudo cpan threads YAML Thread::Queue ### ProtHint dependencies
sudo cpan File::Spec::Functions List::Util Module::Load::Conditional POSIX File::HomeDir ### Braker dependencies
sudo cpan Scalar::Util::Numeric List::MoreUtils ### Braker dependencies
```

2. GeneMark-EX
Download **GeneMark-ES/ET/EP** manually from (http://exon.gatech.edu/GeneMark/license_download.cgi)[http://exon.gatech.edu/GeneMark/license_download.cgi]. Enter the credentials being required. You will need to download the software and its corresponding key.
```{sh}
tar -xvzf gmes_linux_64.tar.gz ### decompress the software
gunzip -d gm_key_64.gz; mv gm_key_64 gmes_linux_64/.gm_key ### decompress, rename, set as hidden, and move to the GeneMark-EX directory
cd gmes_linux_64/
./check_install.bash ### check installation of GeneMark-EX
cd -
```

3. Augustus
```{sh}
wget https://github.com/Gaius-Augustus/Augustus/releases/download/v3.4.0/augustus-3.4.0.tar.gz
tar -xvzf augustus-3.4.0.tar.gz
sudo apt install libboost-iostreams-dev zlib1g-dev libgsl-dev libboost-all-dev libsuitesparse-dev liblpsolve55-dev \
                 libsqlite3-dev libmysql++-dev \
                 libbamtools-dev libboost-all-dev libboost-all-dev \
                 libhts-dev
cd augustus-3.4.0/
make
bin/augustus --species=help
auxprogs/bam2hints/bam2hints -h
sudo make install
cd -    
```

4. Python 3
```{sh}
sudo apt install python3.8.10
```

5. Samtools and Bamtools
```{sh}
sudo apt install samtools bamtools
```

6. NCBI+
```{sh}
sudo apt install ncbi-blast+
```

7. ProtHint
```{sh}
wget https://github.com/gatech-genemark/ProtHint/releases/download/v2.6.0/ProtHint-2.6.0.tar.gz
tar -xvzf ProtHint-2.6.0.tar.gz
cd ProtHint-2.6.0
bin/prothint.py -h
cd -
```

8. Biopython
```{sh}
sudo pip3 install biopython
```

9. cdbfasta
```{sh}
sudo apt install cdbfasta
```

10. GenomeThreader
```{sh}
wget https://genomethreader.org/distributions/gth-1.7.3-Linux_x86_64-64bit.tar.gz
tar -xvzf gth-1.7.3-Linux_x86_64-64bit.tar.gz
```

11. Exonorate
```{sh}
sudo apt install exonerate
```

12. GUSHR
```{sh}
sudo apt install openjdk-8-jdk
git clone https://github.com/Gaius-Augustus/GUSHR.git
```

13. MakHub
```{sh}
wget https://github.com/Gaius-Augustus/MakeHub/archive/refs/tags/1.0.6.tar.gz
```

14. Install Star transcriptome aligner
```{sh}
sudo apt install -y rna-star
```

## Install BRAKER
```{sh}
wget https://github.com/Gaius-Augustus/BRAKER/archive/refs/tags/v2.1.6.tar.gz
tar -xvzf v2.1.6.tar.gz
```

## Download Viridiplantae protein database
```{sh}
wget https://v100.orthodb.org/download/odb10_plants_fasta.tar.gz
tar -xvzf odb10_plants_fasta.tar.gz
cat plants/Rawdata/* > plant_proteins.fasta
```

## Genome assembly and RNAseq data
```{sh}
GENOME=/data/Lolium_rigidum_ASSEMBLY/GENOME_ASSEMBLY/APGP_CSIRO_Lrig_flye-racon-polca-allhic-juicebox_v0.1n_clean1.fasta
DIR_RAW_RNASEQ=/data/Lolium_rigidum_ASSEMBLY/TRANSCRIPTOME_ASSEMBLY/raw_reads
cd $DIR_RAW_RNASEQ
cat *_R1.fastq > ${DIR}/RNAseq_R1.fastq
cat *_R2.fastq > ${DIR}/RNAseq_R2.fastq
cd -
RNASEQ_BAM=RNAseqAlignedSorted.bam
```

## Align transcriptome to the genome
```{sh}
### Prepare genome assembly for alignment
time \
STAR --runMode genomeGenerate \
    --genomeDir $(dirname ${GENOME}) \
    --genomeFastaFiles ${GENOME} \
    --genomeSAindexNbases 13 \
    --runThreadN 31

### Align
time \
STAR --genomeDir $(dirname ${GENOME}) \
    --readFilesIn \
    ${DIR}/RNAseq_R1.fastq \
    ${DIR}/RNAseq_R2.fastq \
    --runThreadN 31 \
    --outFileNamePrefix RNAseq

### Sort and compress
time \
samtools view -bS RNAseqAligned.out.sam | \
    samtools sort > ${RNASEQ_BAM}
```

## Run

Quote from BRAKER README.md:
"Even though BRAKER supports the combination of RNA-Seq and protein data within the BRAKER pipeline, we strongly recommend to run BRAKER twice (1x with RNA-Seq only, 1x with protein data only) and subsequently combine the results of both runs with TSEBRA, the Transcript Selector for BRAKER (https://github.com/Gaius-Augustus/TSEBRA). You find more information on TSEBRA at https://www.biorxiv.org/content/10.1101/2021.06.07.447316v1"

Input file paths:
```{sh}
GENOME=/data/Lolium_rigidum_ASSEMBLY/GENOME_ASSEMBLY/APGP_CSIRO_Lrig_flye-racon-polca-allhic-juicebox_v0.1n_clean1.fasta
RNASEQ_BAM=/data/Lolium_rigidum_ASSEMBLY/ANNOTATION/RNAseqAlignedSorted.bam
PROTEIN=/data/Lolium_rigidum_ASSEMBLY/ANNOTATION/plant_proteins.fasta
```

Set-up paths:
```{sh}
DIR=/data/Lolium_rigidum_ASSEMBLY/ANNOTATION
export GENEMARK_PATH=${DIR}/gmes_linux_64
export AUGUSTUS_CONFIG_PATH=${DIR}/augustus-3.4.0/config
export AUGUSTUS_BIN_PATH=${DIR}/augustus-3.4.0/bin
export AUGUSTUS_SCRIPTS_PATH=${DIR}/augustus-3.4.0/scripts
export PROTHINT_PATH=${DIR}/ProtHint-2.6.0/bin
PATH=${PATH}:${DIR}/gth-1.7.3-Linux_x86_64-64bit/bin
PATH=${PATH}:${DIR}/BRAKER-2.1.6/scripts
```

Step 1 of 3: BRAKER run using RNAseq data (1,261 minutes)
```{sh}
time \
braker.pl \
    --species=Lolium_rigidum_RNAseq \
    --genome=${GENOME} \
    --bam=${RNASEQ_BAM} \
    --cores 32
```

Check log
```{sh}
bat braker/errors/new_species.stderr
bat braker/braker.log
```

Step 2 of 3: BRAKER run using protein database information
```{sh}
time \
braker.pl \
    --species=Lolium_rigidum_Viridiprot \
    --genome=${GENOME} \
    --prot_seq=${PROTEIN} \
    --epmode \
    --cores 32
```
Step 3 of 3: TSEBRA
```{sh}
```