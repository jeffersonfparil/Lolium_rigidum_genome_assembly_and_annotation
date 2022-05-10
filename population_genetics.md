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

## Generate pileup files: 1 pileup : 1 chromosome and run npstat on each
```{sh}
echo '#!/bin/bash
i=$1
BAM=$2
# i=3
CHR=Chromosome${i}
samtools index ${BAM}
samtools mpileup \
    -r ${CHR} \
    ${BAM} > \
    ${BAM%.bam*}-${CHR}.pileup
npstat \
    -n 42 \
    -l 1000000 \
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