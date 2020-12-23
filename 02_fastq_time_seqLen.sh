#!/bin/bash
######################################################################
### counting the length of reads and noting the time of sequencing ###
######################################################################
### iterate across sequencing experiments
for PREFIX in lolium1 lolium2
do
    # PREFIX=lolium1
    ### setwd
    DIR=/mnt/Lolium_rigidum_ASSEMBLY/Im_bored_20200702/FASTQ/${PREFIX}
    ### initialize output file
    ### headerless: ${time}/t${seq_len}
    touch ${PREFIX}.length.time
    ### iterate across fastq files
    for f in $(find ${DIR}/*.fastq)
    do
        # f=/mnt/Lolium_rigidum_ASSEMBLY/Im_bored_20200702/FASTQ/lolium1/fastq_runid_d86a9d393d89555bee944f5bb8c98b4f44ebf601_0_0.fastq
        echo $f
        ### count the number of lines and loci in the fastq file (the sequence line is assumed to occupy a single line)
        nLines=$(cat $f | wc -l)
        nLoci=$(echo $nLines / 4 | bc)
        ### iterate across sequence reads
        for i in $(seq 1 $nLoci)
        do
            ### extract the line number for the meta data to extract the tiome of sequencing as well as the sequence line itself
            idx_meta=$(echo "( (${i}-1) * 4 ) + 1" | bc)
            idx_seq=$(echo "( (${i}-1) * 4 ) + 2" | bc)
            ### extract the time and sequence length
            metaData=$(sed "${idx_meta}q;d" $f)
            time=$(echo "$metaData" | cut -d" " -f6 | sed 's/start_time=//g')
            sequence=$(sed "${idx_seq}q;d" $f)
            seq_len=$(echo $sequence | wc -c)
            ### write-out time and sequence length into the output file
            echo -e "$time\t$seq_len" >> ${PREFIX}.length.time
        done # ${i}
    done # ${f}
done # ${PREFIX}

