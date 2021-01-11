#!/bin/bash
#######################################################
### Filter and trimm-off adapters from MinION reads ###
#######################################################

### Input:
### (1) MinION reads in compressed fastq format (minion.fastq.gz)

### Outputs:
### (1) Filtered reads (minion-filtered.fastq.gz)
### (2) Filtered and trimmed reads (minion-filtered-trimmed.fastq.gz)
### (3) Number of bases sequenced after filtering and trimming (minion-filtered-trimmed.base.count)
### (4) Histogram of read length distribution (minion-filtered-trimmed.readlen.hist.svg)

### Parameters:
DIR=/data/Lolium_rigidum_ASSEMBLY/assembly_annotation_pipeline_tests_20210104/FASTQ/MINION/

### Navigate to the working directory
cd $DIR

### Filter-out reads with average PHRED score < 10 with nanofilt (PHRED threshold based on a qualitative look at the FastQC output)
time \
     gunzip -c minion.fastq.gz | \
     NanoFilt -q 10 | \
     gzip > minion-filtered.fastq.gz

### Trim-off adapters with porechop
time \
porechop \
     --threads 32 \
     --input minion-filtered.fastq.gz \
     --output minion-filtered-trimmed.fastq.gz

### Count the total number of bases sequenced
zcat minion-filtered-trimmed.fastq.gz | paste - - - - | cut -f 2 | tr -d '\n' | wc -c > minion-filtered-trimmed.base.count
### As of 2021-01-07 we have ~5 billion bases (5,369,872,027) sequenced which means ~2.68X depth of coverage

### Plot histogram of read lengths
echo '
using FASTX
using CodecZlib
using Plots
vec_lengths = []
try
     for record in FASTX.FASTQ.Reader(CodecZlib.GzipDecompressorStream(open("minion-filtered-trimmed.fastq.gz")))
          # seq = FASTX.FASTQ.sequence(record)
          push!(vec_lengths, FASTX.FASTQ.seqlen(record))
     end
catch
     run(`gunzip minion-filtered-trimmed.fastq.gz`)
     for record in FASTX.FASTQ.Reader(open("minion-filtered-trimmed.fastq"))
          seq = FASTX.FASTQ.sequence(record)
          push!(vec_lengths, length(seq))
     end
     run(`gzip minion-filtered-trimmed.fastq`)
end
vec_lengths = convert(Array{Int64,1}, vec_lengths)
# using UnicodePlots
# UnicodePlots.histogram(vec_lengths)
hist = Plots.histogram(vec_lengths, legend=false);
savefig("minion-filtered-trimmed.readlen.hist.svg")
' > readlen_hist.jl
julia readlen_hist.jl

