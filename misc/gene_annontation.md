# Gaia-Augustus/BRAKER pipeline D + GeMoMa

![](misc/braker2_pipeline_D.png)

## Working directory
```{sh}
DIR=/data/Lolium_rigidum_ASSEMBLY/ANNOTATION
cd $DIR
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

15. Java
```{sh}
sudo apt install default-jre
```

16. Bedtools
```{sh}
sudo apt install bedtools
```

## Install BRAKER
```{sh}
wget https://github.com/Gaius-Augustus/BRAKER/archive/refs/tags/v2.1.6.tar.gz
tar -xvzf v2.1.6.tar.gz
```

## Install TSEBRA
```{sh}
git clone https://github.com/Gaius-Augustus/TSEBRA
cd TSEBRA/bin
./tsebra.py -h
cd -
```


## Download Viridiplantae protein database
```{sh}
wget https://v100.orthodb.org/download/odb10_plants_fasta.tar.gz
tar -xvzf odb10_plants_fasta.tar.gz
cat plants/Rawdata/* > plant_proteins.fasta
```

## Define genome assembly and RNAseq data
```{sh}
GENOME=/data/Lolium_rigidum_ASSEMBLY/GENOME_ASSEMBLY/APGP_CSIRO_Lrig_flye-racon-polca-allhic-juicebox_v0.1n_clean1.fasta
DIR_RAW_RNASEQ=/data/Lolium_rigidum_ASSEMBLY/TRANSCRIPTOME_ASSEMBLY/rRNAdepleted_reads/
cd $DIR_RAW_RNASEQ
zcat *_paired_*.fq.1.gz > ${DIR}/RNAseq_R1.fastq.gz
zcat *_paired_*.fq.2.gz > ${DIR}/RNAseq_R2.fastq.gz
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
    ${DIR}/RNAseq_R1.fastq.gz \
    ${DIR}/RNAseq_R2.fastq.gz \
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

### Input file paths:
```{sh}
GENOME=/data/Lolium_rigidum_ASSEMBLY/GENOME_ASSEMBLY/APGP_CSIRO_Lrig_flye-racon-polca-allhic-juicebox_v0.1n_clean1.fasta
RNASEQ_BAM=/data/Lolium_rigidum_ASSEMBLY/ANNOTATION/RNAseqAlignedSorted.bam
PROTEIN=/data/Lolium_rigidum_ASSEMBLY/ANNOTATION/plant_proteins.fasta
```

### Set-up paths:
```{sh}
DIR=/data/Lolium_rigidum_ASSEMBLY/ANNOTATION
export GENEMARK_PATH=${DIR}/gmes_linux_64
export AUGUSTUS_CONFIG_PATH=${DIR}/augustus-3.4.0/config
export AUGUSTUS_BIN_PATH=${DIR}/augustus-3.4.0/bin
export AUGUSTUS_SCRIPTS_PATH=${DIR}/augustus-3.4.0/scripts
export PROTHINT_PATH=${DIR}/ProtHint-2.6.0/bin
PATH=${PATH}:${DIR}/gth-1.7.3-Linux_x86_64-64bit/bin
PATH=${PATH}:${DIR}/BRAKER-2.1.6/scripts
PATH=${PATH}:${GENEMARK_PATH}
PATH=${PATH}:${DIR}/TSEBRA/bin
```

### Step 1 of 4: BRAKER run using RNAseq data (1,261 minutes)
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

Rename braker output folder
```{sh}
mv braker braker_RNAseq
```

### Step 2 of 4: BRAKER run using protein database information (Can be run in parallel with Step 1)

But first, run ProtHint separately so we can finish the algorithm before the heat death of the universe using the option: `--maxProteinsPerSeed 5`:
```{sh}
time \
${PROTHINT_PATH}/prothint.py \
    ${GENOME} \
    ${PROTEIN} \
    --maxProteinsPerSeed 5 \
    --threads 30

mkdir PROTHINT_plant_proteins
mv prothint.gff PROTHINT_plant_proteins
mv top_chains.gff PROTHINT_plant_proteins
mv evidence.gff PROTHINT_plant_proteins
mv prothint_augustus.gff PROTHINT_plant_proteins
PROTHINT_AUGUSTUS_HINTS=${DIR}/PROTHINT_plant_proteins/prothint_augustus.gff
```

Now run braker with prothint-generated hints:
```{sh}
time \
braker.pl \
    --species=Lolium_rigidum_Viridiprot \
    --genome=${GENOME} \
    --hints=${PROTHINT_AUGUSTUS_HINTS} \
    --cores 32

bat braker/braker.log
```

Rename braker output folder
```{sh}
mv braker braker_proteins
```


### Step 3 of 4: TSEBRA
```{sh}
time \
tsebra.py \
    -g ${DIR}/braker_RNAseq/augustus.hints.gtf,${DIR}/braker_proteins/augustus.hints.gtf \
    -c ${DIR}/TSEBRA/config/default.cfg \
    -e ${DIR}/braker_RNAseq/hintsfile.gff,${DIR}/braker_proteins/hintsfile.gff \
    -o ${DIR}/BRAKER_OUTPUT.gtf
```

### Step 4 of 4: BLASTX predicted genes against the plant protein database
```{sh}
# Using the fixed fasta file, i.e. non-terminal lines have the same number of characters per line
GENOME_FIXED_CHARNUM_PER_LINE=/data/Lolium/Quantitative_Genetics/02_FASTQ/REFERENCE/Reference.fasta
# Filter and modify name field in ${DIR}/BRAKER_OUTPUT.gtf
grep -v "start_codon\|stop_codon" ${DIR}/BRAKER_OUTPUT.gtf > ${DIR}/BRAKER_OUTPUT-FILTERED.tmp
awk -F'\t' -v OFS='\t' '$3=$9' ${DIR}/BRAKER_OUTPUT-FILTERED.tmp > ${DIR}/BRAKER_OUTPUT-FILTERED.gtf

time \
bedtools getfasta \
    -fi ${GENOME_FIXED_CHARNUM_PER_LINE} \
    -bed ${DIR}/BRAKER_OUTPUT-FILTERED.gtf \
    -name+ \
    -fo ${DIR}/BRAKER_OUTPUT.fasta

sed -i 's/ "/-/g' ${DIR}/BRAKER_OUTPUT.fasta
sed -i 's/"; /-/g' ${DIR}/BRAKER_OUTPUT.fasta
sed -i 's/";::/::/g' ${DIR}/BRAKER_OUTPUT.fasta

time \
makeblastdb \
    -in ${PROTEIN} \
    -dbtype prot

time \
blastx -db ${PROTEIN} \
    -query ${DIR}/BRAKER_OUTPUT.fasta \
    -outfmt "6 qseqid staxids sstart send pident evalue qcovhsp bitscore stitle" \
    -max_target_seqs 5 \
    -num_threads 31 \
    -out ${DIR}/BRAKER_OUTPUT.blastout
```
