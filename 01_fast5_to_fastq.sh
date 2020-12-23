#!/bin/bash
for i in 1 2
do
    /mnt/Lolium_rigidum_ASSEMBLY/GUPPY/ont-guppy-cpu/bin/guppy_basecaller \
            --input_path /mnt/Lolium_rigidum_ASSEMBLY/Im_bored_20200702/FAST5/lolium${i} \
            --save_path /mnt/Lolium_rigidum_ASSEMBLY/Im_bored_20200702/FASTQ/lolium${i} \
            --flowcell FLO-MIN106 \
            --kit SQK-LSK109 \
            --cpu_threads_per_caller 12
done

