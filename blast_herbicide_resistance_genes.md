# BLAST herbicide resistance genes

## Set our working directory directory, reference genome assembly, and TSR & NTSR protein query sequences
```{sh}
DIR=/data-weedomics-3
REF=${DIR}/APGP_CSIRO_Lrig_flye-racon-polca-allhic-juicebox_v0.1n_clean1.fasta
DIR_GENES=${DIR}/Lolium_rigidum_genome_assembly_and_annotation/misc/TSR_NTSR_etc_protein_sequences
cd $DIR
```

## Download the protein query sequences from our github repo, and install BLAST from NCBI
```{sh}
git clone https://github.com/jeffersonfparil/Lolium_rigidum_genome_assembly_and_annotation.git
### Copy the Lolium rigidum genome assembly
sudo apt install ncbi-blast+
```

## Generate the BLAST database from our genome assembly
```{sh}
time \
makeblastdb -in ${REF} \
            -title "Lolium_rigidum" \
            -dbtype nucl
```

## BLAST using protein sequences as query and nucleotide sequence database (our genome assembly)
```{sh}
echo '#!/bin/bash
query=$1
REF=$2
DIR=$3
echo $query
temp_name_1=$(basename $query)
temp_name_2=${temp_name_1%.fasta*}
tblastn -db ${REF} \
    -query ${query} \
    -outfmt "6 qseqid staxids pident evalue qcovhsp bitscore stitle" \
    -out ${DIR}/BLASTOUT-${temp_name_2}.txt
' > tblastn_for_parallel_execution.sh
chmod +x tblastn_for_parallel_execution.sh
time parallel -j 14 ./tblastn_for_parallel_execution.sh {} ${REF} ${DIR} ::: $(find ${DIR_GENES}/*.fasta)
```

*Note:* Atrazine, clethodim, and paraquat will take a very very very long time to finsh > 1 week - I had to stop it manually.
