#!/bin/bash
### Filter-out reads with average PHRED score < 10 with nanofilt
### e.g.
cd /data/Lolium_rigidum_ASSEMBLY
echo '#!/bin/bash
     FNAME=$1
     NanoFilt -q 10 ${FNAME} > ${FNAME%.fastq*}-nanofilted.fastq
     ' > nanofilt_for_parallel.sh
chmod +x nanofilt_for_parallel.sh
time \
parallel ./nanofilt_for_parallel.sh {} ::: $(ls FAST5/lol_full_protocol/*.fastq)
cat FAST5/lol_full_protocol/*-nanofilted.fastq > FASTQ/lolium5.fastq ### About half of the reads were discarded! So maybe realx -q from 10 to 5?
rm nanofilt_for_parallel.sh

### Remove  adapters with porechop
### e.g.
cd /data/Lolium_rigidum_ASSEMBLY
time \
porechop --threads 32 --input FASTQ/lolium5.fastq --output FASTQ/lolium5-porechoped.fastq
