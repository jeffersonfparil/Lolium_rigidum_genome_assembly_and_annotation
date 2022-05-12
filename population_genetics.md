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

## Merge NPSTAT output
```{sh}
echo -e "Population\tChromosome" > col1_to_2.tmp
head -n1 $(ls BAM/*.stats | head -n1) > col3_to_n.tmp
paste col1_to_2.tmp col3_to_n.tmp > Lolium_rigidum.popgen
for f in $(ls BAM/*.stats)
do
    # f=$(ls BAM/*.stats | head -n13 | tail -n1)
    b=$(basename $f)
    pop=${b%-*}
    pop=$(echo $b | cut -d"-" -f1)
    chr=$(echo $b | cut -d"-" -f2 | cut -d"." -f1)
    printf "$pop\n%.s" $(seq 1 $(cat $f | wc -l)) > pop.tmp
    printf "$chr\n%.s" $(seq 1 $(cat $f | wc -l)) > chr.tmp
    paste pop.tmp chr.tmp $f > merged.tmp
    tail -n+2 merged.tmp >> Lolium_rigidum.popgen
done
rm *.tmp
```

## Population genetics analyses
```{R}
dat = read.delim("Lolium_rigidum.popgen", header=T)
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
