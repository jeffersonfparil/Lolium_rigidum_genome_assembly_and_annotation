# Population genetics analysis of 60 Lolium rigidum population from SE Australia
- Here we derive site frequency spectrum (SFS) summary statistics, Watterson's theta, Tajima's D, Fst, etc...

## Set-up working directories
```{sh}
DIR=/data-weedomics-1/Lolium_rigidum_population_genetics_analysis
SRC=${DIR}/Lolium_rigidum_genome_assembly_and_annotation
mkdir ${DIR}/REFERENCE
mkdir ${DIR}/FASTQ
mkdir ${DIR}/BAM
mkdir ${DIR}/PILEUP
mkdir ${DIR}/NPSTAT
```

## Clone analysis repository
```{sh}
git clone https://github.com/jeffersonfparil/Lolium_rigidum_genome_assembly_and_annotation.git
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

## Install Popoolation2
```{sh}
wget https://sourceforge.net/projects/popoolation2/files/popoolation2_1201.zip
unzip popoolation2_1201.zip
rm popoolation2_1201.zip
```

## Prepare the reference genome sequence and annotation
```{sh}
cd ${DIR}/REFERENCE
### Download reference genome and annotation files
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/022/539/505/GCF_022539505.1_APGP_CSIRO_Lrig_0.1/GCF_022539505.1_APGP_CSIRO_Lrig_0.1_genomic.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/022/539/505/GCF_022539505.1_APGP_CSIRO_Lrig_0.1/GCF_022539505.1_APGP_CSIRO_Lrig_0.1_genomic.gff.gz
gunzip -c GCF_022539505.1_APGP_CSIRO_Lrig_0.1_genomic.fna.gz > Lolium_rigidum.fasta
gunzip -c GCF_022539505.1_APGP_CSIRO_Lrig_0.1_genomic.gff.gz > Lolium_rigidum.gff
rm GCF_022539505.1_APGP_CSIRO_Lrig_0.1_genomic.fna.gz
rm GCF_022539505.1_APGP_CSIRO_Lrig_0.1_genomic.gff.gz

### Extract each of the 7 chromosomes
echo '#!/bin/bash
i=$1
SRC=$2
CHR="chromosome ${i}"
NEW_SEQ_NAME=$(echo $CHR | sed "s/ /_/g")
FNAME_TMP=${NEW_SEQ_NAME}.fasta.tmp ### NO SPACES PLEASE!
FNAME=${FNAME_TMP%.tmp*}
julia ${SRC}/extract_sequence_using_name_query.jl \
    Lolium_rigidum.fasta \
    ${CHR} \
    ${FNAME_TMP} \
    ${NEW_SEQ_NAME} \
    false
julia ${SRC}/reformat_fasta_sequence.jl \
    ${FNAME_TMP} \
    50 \
    ${FNAME}
rm ${FNAME_TMP}
' > extract_chrom_and_rename.sh
chmod +x extract_chrom_and_rename.sh

time \
parallel ./extract_chrom_and_rename.sh \
    {} ${SRC} ::: $(seq 1 7)


### Extract genome annotation file, split by chromosome and rename to be same as the fasta's i.e. Chromosome${1..7}
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/022/539/505/GCF_022539505.1_APGP_CSIRO_Lrig_0.1/GCF_022539505.1_APGP_CSIRO_Lrig_0.1_genomic.gff.gz
gunzip -c GCF_022539505.1_APGP_CSIRO_Lrig_0.1_genomic.gff.gz > Lolium_rigidum.gff

grep "^NC" Lolium_rigidum.gff | cut -f1 | sort | uniq > chr_id.tmp ### assumes we get the 7 chromosomes hence 7 lines here!!!
for i in $(seq 1 7)
do
    # i=1
    id=$(head -n${i} chr_id.tmp | tail -n1)
    CHR=Chromosome${i}
    grep ${id} Lolium_rigidum.gff | sed "s/${id}/${CHR}/g" > ${CHR}.gff
done

### Clean-up
rm chr_id.tmp
cd ${DIR}
```





## Run NPSTAT to estimate Watterson's theta, Tajima's D, Fay and Wu's 
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

## Clean-up
```{sh}
mkdir ${DIR}/PILEUP
mkdir ${DIR}/NPSTAT
mv ${DIR}/BAM/*.pileup ${DIR}/PILEUP
mv ${DIR}/BAM/*.stats ${DIR}/NPSTAT
```

## Merge NPSTAT output
```{sh}
echo -e "Population\tChromosome" > col1_to_2.tmp
head -n1 $(ls ${DIR}/NPSTAT/*.stats | head -n1) > col3_to_n.tmp
paste col1_to_2.tmp col3_to_n.tmp > Lolium_rigidum.npstat
for f in $(ls ${DIR}/NPSTAT/*.stats)
do
    # f=$(ls BAM/*.stats | head -n13 | tail -n1)
    b=$(basename $f)
    pop=${b%-*}
    pop=$(echo $b | cut -d"-" -f1)
    chr=$(echo $b | cut -d"-" -f2 | cut -d"." -f1)
    printf "$pop\n%.s" $(seq 1 $(cat $f | wc -l)) > pop.tmp
    printf "$chr\n%.s" $(seq 1 $(cat $f | wc -l)) > chr.tmp
    paste pop.tmp chr.tmp $f > merged.tmp
    tail -n+2 merged.tmp >> Lolium_rigidum.npstat
done
rm *.tmp
```

## Run Popoolation2 to estimate
```{sh}

perl <popoolation2-path>/fst-sliding.pl --input p1_p2.sync --output p1_p2_w500.fst --min-count 6 --min-coverage 50 --max-coverage 200 --min-covered-fraction 1 --window-size 500 --step-size 500 --pool-size 500
```


## Population genetics analyses
```{R}
dat = read.delim("Lolium_rigidum.npstat", header=T)
vec_pops = unique(dat$Population)
vec_chrs = unique(dat$Chromosome)
vec_resp = colnames(dat)[5:ncol(dat)]

window_length = 1e+6

i = 34
j = 3

l = length(vec_resp)
n = ceiling(sqrt(l))
m = ceiling(l / n)
par(mfrow=c(n, m), mar=c(5, 5, 2, 1))
for (k in 1:l){
    # k = 9
    pop = dat$Population[i]
    chr = dat$Chromosome[j]
    res = vec_resp[k]

    idx = (dat$Population==pop) & (dat$Chromosome==chr)
    df = droplevels(dat[idx, ])

    x = df$window
    y = eval(parse(text=paste0("df$", res)))
    y = eval(parse(text=paste0("df$", res)))
    plot(x, y, type="l", las=2, xlab="", ylab="", xaxt="n", main=res)
    x_min = min(x)
    x_max = max(x)
    pos = round(seq(x_min, x_max, length=5))
    axis(side=1, at=pos, lab=pos)
    mtext(side=1, text=paste0(chr, ": window (Mb)"), cex=0.6, padj=3.5)
    grid()
}

```
