# Why are we mapping less Lolium rigidum Pool-seq reads into the new genome that the old Byrne's Lolium perenne genome?

## Install Mummer
```{sh}
wget https://sourceforge.net/projects/mummer/files/mummer/3.23/MUMmer3.23.tar.gz
tar -xvzf MUMmer3.23.tar.gz
cd MUMmer3.23/
sed -z -i 's/#ifdef SIXTYFOURBITS/#define SIXTYFOURBITS\n#ifdef SIXTYFOURBITS/g' src/kurtz/libbasedir/types.h
make

wget https://github.com/mummer4/mummer/releases/download/v4.0.0rc1/mummer-4.0.0rc1.tar.gz
tar -xvzf mummer-4.0.0rc1.tar.gz
cd mummer-4.0.0rc1/
./configure
make
```

## Run
```{sh}
# MUMMER=/data-weedomics-3/MUMmer3.23/mummer
# REF_LORI=/data-weedomics-3/APGP_CSIRO_Lrig_flye-racon-polca-allhic-juicebox_v0.1n_clean1.fasta
# REF_LOPE=/data-weedomics-3/Lolium_perenne/lope_V1.0.fasta
MUMMER=/data/Lolium/Quantitative_Genetics/mummer-4.0.0rc1/mummer
REF_LORI=/data/Lolium/Quantitative_Genetics/02_FASTQ/REFERENCE/Reference.fasta
REF_LOPE=/data/Lolium/Quantitative_Genetics/02_FASTQ/REFERENCE/BK_Lperenne/Reference.fasta
time \
${MUMMER} \
    -thread 20 \
    ${REF_LORI} ${REF_LOPE} > Lope_on_Lori.mums
```

<!-- ## Plot output
```{sh}
# MUMMERPLOT=/data-weedomics-3/MUMmer3.23/mummerplot
MUMMERPLOT=/data/Lolium/Quantitative_Genetics/mummer-4.0.0rc1/mummerplot
time \
${MUMMERPLOT} \
    --png \
    --prefix=Lope_on_Lori \
    Lope_on_Lori.mums

gnuplot Lope_on_Lori.gp
eog Lope_on_Lori.ps

```
 -->

## Run nucmer
```{sh}
NUCMER=/data/Lolium/Quantitative_Genetics/mummer-4.0.0rc1/nucmer
${NUCMER} \
    -p Lope_on_Lori \
    -t 20 \
    ${REF_LORI} ${REF_LOPE}

```