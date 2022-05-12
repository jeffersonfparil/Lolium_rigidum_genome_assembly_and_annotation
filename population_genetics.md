# Population genetics analysis of 60 Lolium rigidum population from SE Australi

## Set-up working directories
```{sh}
DIR=/data-weedomics-1/Lolium_rigidum_population_genetics_analysis
```

## Install npstat
```{sh}
sudo apt install samtools gsl-bin libgsl0-dev
git clone https://github.com/lucaferretti/npstat.git
cd npstat
make
PATH=${PATH}:$(pwd)
cd -
```

## Prepare the reference genome sequence and annotation
```{sh}
### Extract each of the 7 chromosomes
cd REFERENCE/
for i in $(seq 1 7)
do
    # i=1
    CHR="Chromosome${i}"
    julia extract_sequence_using_name_query.jl \
        Reference.fasta \
        ${CHR} \
        ${CHR}.fasta.tmp \
        ${CHR} \
        false
    julia reformat_fasta_sequence.jl \
        ${CHR}.fasta.tmp \
        50 \
        ${CHR}.fasta
    rm ${CHR}.fasta.tmp
done

### Extract genome annotation file, split by chromosome and rename to be same as the fasta's i.e. Chromosome${1..7}
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/022/539/505/GCF_022539505.1_APGP_CSIRO_Lrig_0.1/GCF_022539505.1_APGP_CSIRO_Lrig_0.1_genomic.gff.gz
gunzip -c GCF_022539505.1_APGP_CSIRO_Lrig_0.1_genomic.gff.gz > Reference.gff

grep "^NC" Reference.gff | cut -f1 | sort | uniq > chr_id.tmp ### assumes we get the 7 chromosomes hence 7 lines here!!!
for i in $(seq 1 7)
do
    # i=1
    id=$(head -n${i} chr_id.tmp | tail -n1)
    CHR=Chromosome${i}
    grep ${id} Reference.gff | sed "s/${id}/${CHR}/g" > ${CHR}.gff
done

### Clean-up
rm chr_id.tmp
cd -
```

## Run NPSTAT
```{sh}
echo '#!/bin/bash
i=$1
BAM=$2
# i=3
# BAM=BAM/ACC09.bam
CHR=Chromosome${i}
samtools index ${BAM}
samtools mpileup \
    -r ${CHR} \
    ${BAM} > \
    ${BAM%.bam*}-${CHR}.pileup
npstat \
    -n 42 \
    -l 1000000 \
    -outgroup REFERENCE/${CHR}.fasta \
    -annot REFERENCE/${CHR}.gff \
     ${BAM%.bam*}-${CHR}.pileup
' > parallel_pileup_npstat.sh
chmod +x parallel_pileup_npstat.sh
time \
parallel \
./parallel_pileup_npstat.sh \
    {1} {2} \
    ::: $(seq 1 7) \
    ::: $(find ${DIR}/BAM/ -name '*.bam')
```

## Population genetics analyses
```{R}

```
