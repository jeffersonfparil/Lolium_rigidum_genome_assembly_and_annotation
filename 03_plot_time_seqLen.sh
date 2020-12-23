#!/bin/bash
##########################################
### read length and a function of time ###
##########################################
time \
for f in $(ls *.length.time)
do
	echo $f
	Rscript 03_plot_time_seqLen.R $f
done
