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


./delta-filter -l 1000 -q Lope_on_Lori.delta > Lope_on_Lori_filter.delta
./show-coords -c -l -L 1000 -r -T Lope_on_Lori_filter.delta > Lope_on_Lori_filter_coords.txt
```

## Visualise output in R
```{R}
install.packages(c("devtools", "ggplot2"))

library(devtools)
install_github("timflutre/rutilstimflutre")

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(c("ggbio", "GenomicRanges"))



load_coords <- function(coords_file, perc.id) {
  library(rutilstimflutre)
  library(GenomicRanges)
  coords <- loadMummer(coords_file, algo = "nucmer")                    # read in nucmer results as a GRanges obj
  coords <- coords[(elementMetadata(coords)[ , "perc.id"] >= perc.id)]  # filter entries with perc.id != 100
  seqlevels(coords) <- seqlevelsInUse(coords)                           # drop seq levels that are no longer used
  coords$qry.name <- names(coords)                                      # set column with query contig name
  coords$qry <- rep("Query", nrow(values(coords)))
  return(coords)
}

faidx_to_GRanges <- function(faidx_file){
  faidx <- read.table(file = faidx_file, header = F, stringsAsFactors = F,
                      col.names = c("name", "contig.len", "offset", 
                                    "linebases", "linewidth"))
  # create a GRanges object for reference sequences
  grange <- GRanges(seqnames = Rle(faidx$name), 
                    ranges = IRanges(start = rep(1, nrow(faidx)), 
                                     end = faidx$contig.len))
  # add seqlengths to the reference GRanges object
  seqlengths(grange) <- faidx$contig.len
  genome(grange) <- "Reference"
  grange <- sortSeqlevels(grange)
  grange <- sort(grange)
  return(grange)
}

circular_plot_w_ref <- function(reference_GRange, NUCmer_coords){
  library(ggplot2)
  library(ggbio)
  # reference_GRange: reference sequence GRanges obj
  # NUCmer_coords: a GRanges object produced by reading in a show-coords processed NUCmer object.
  p <- ggbio() + 
    circle(NUCmer_coords, geom = "rect", 
           aes(color = "steelblue", fill = "steelblue")) +  # NUCmer obj
    circle(reference_GRange, geom = "ideo", 
           aes(color = "gray70", fill = "gray70")) +        # Ideogram of ref genome
    circle(reference_GRange, geom = "scale", 
           scale.type = "M", size = 1.5) +                  # Scale from seqlen of ref genome
    # circle(reference_GRange, geom = "text", 
    #       aes(label = seqnames), size = 2) +              # Uncomment for sequence label
    scale_color_manual(name = "Sequence Origin", 
                       labels = c("Reference", "Query"), 
                       values = c("gray70", "steelblue")) +
    scale_fill_manual(name = "Sequence Origin", 
                      labels = c("Reference", "Query"), 
                      values = c("gray70", "steelblue")) +
    ggtitle("Reference vs. Query")
  return(p)
}

load_and_plot_nucmer_w_ref <- function(NUCmer_coords_file, ref_faidx_file, perc.id) {
  # NUCmer_coords_file: string for file path of NUCmer output produced by show coords 
  # ref_faidx_file: reference sequence GRanges obj
  # perc.id: percent id cutoff to show on plot
  
  library(GenomicRanges)
  
  NUCmer_coords <- load_coords(NUCmer_coords_file, perc.id = perc.id)   # Make GRanges obj of nucmer output file
  referenceGR <- faidx_to_GRanges(faidx_file = ref_faidx_file)
  plot <- circular_plot_w_ref(reference_GRange = referenceGR, NUCmer_coords = NUCmer_coords)
  return(plot)
}


NUCmer_coords_file = "Lope_on_Lori_filter_coords.txt"
ref_faidx_file = "/data/Lolium/Quantitative_Genetics/02_FASTQ/REFERENCE/Reference.fasta.fai"
perc.id = 0.85

png("Lope_on_Lori_nucmer_map.png", width=2000, height=2000)
load_and_plot_nucmer_w_ref(NUCmer_coords_file=NUCmer_coords_file, 
                           ref_faidx_file=ref_faidx_file, 
                           perc.id=perc.id)
dev.off()
```

### Find unmapped scaffolds: 7,667 unmapped Lope scaffolds on the Lori assembly
```{sh}
cut -f13 Lope_on_Lori_filter_coords.txt | sort | uniq | wc -l ### 40,748
rg "^>" ${REF_LOPE} | wc -l                                   ### 48,415
mkdir Lope_on_Lori_LOPE_STATS/
cp ${REF_LOP} Lope_on_Lori_LOPE_STATS/
cd Lope_on_Lori_LOPE_STATS/
```

```{julia}
using ProgressMeter
function fun_fasta_lengths_GC_N(str_filename_fasta, n_int_keep_top_seq=0, bool_write_files=false)
    vec_str_sequence_names = String.([])
    vec_int_sequence_lengths = Int.([])
    vec_int_sequence_GC = Int.([])
    vec_int_sequence_N = Int.([])
    regex_GgCc = Regex("[GgCc]")
    regex_Nn = Regex("[Nn]")
    FILE = open(str_filename_fasta, "r")
    ProgressMeter.seekend(FILE) # find file size by navigating to the end of the file
    n_int_FILE_size = ProgressMeter.position(FILE) # find the size of the file which is not equal to the the number of lines in the file
    pb = ProgressMeter.Progress(n_int_FILE_size, 1) # initialise the progress bar
    ProgressMeter.seekstart(FILE) # reset to the begining of the file to initial the progress bar and the while loop
    while !eof(FILE)
        line = readline(FILE)
        if line[1]=='>'
            try
                close(file_seq)
                close(file_GC)
            catch
                nothing
            end
            str_sequence_name = line[2:end]
            push!(vec_str_sequence_names, str_sequence_name)
            append!(vec_int_sequence_lengths, 0)
            append!(vec_int_sequence_GC, 0)
            append!(vec_int_sequence_N, 0)
            if bool_write_files
                global file_seq = open(string(str_sequence_name, ".fasta"), "w")
                global file_GC  = open(string(str_sequence_name, "-GC_content.txt"), "w")
                write(file_seq, string(line, '\n'))
                # write(file_GC, string(line, '\n'))
            end
        else
            n_int_total = length(line)
            n_int_GC = length(collect(eachmatch(regex_GgCc, line)))
            n_int_N = sum(match.(regex_Nn, string.(collect(line))) .!= nothing)
            vec_int_sequence_lengths[end] = vec_int_sequence_lengths[end] + n_int_total
            vec_int_sequence_GC[end] = vec_int_sequence_GC[end] + n_int_GC
            vec_int_sequence_N[end] = vec_int_sequence_N[end] + n_int_N
            if bool_write_files
                write(file_seq, string(line, '\n'))
                write(file_GC, string(round(n_int_GC / (n_int_total-n_int_N+1e-100), digits=4), '\n'))
            end
        end
        ProgressMeter.update!(pb, ProgressMeter.position(FILE))
    end
    try
        close(file_seq)
        close(file_GC)
    catch
        nothing
    end        
    close(FILE)
    ### Calculate the assembly size including only the n largest sequences (if n=0, then use all the sequences)
    if n_int_keep_top_seq == 0
        n_int_keep_top_seq = length(vec_str_sequence_names)
    end
    vec_int_idx_sortperm = sortperm(vec_int_sequence_lengths, rev=true)
    vec_str_chromosome_names = vec_str_sequence_names[vec_int_idx_sortperm][1:n_int_keep_top_seq]
    vec_int_chromosome_lengths = vec_int_sequence_lengths[vec_int_idx_sortperm][1:n_int_keep_top_seq]
    vec_int_chromosome_GC = vec_int_sequence_GC[vec_int_idx_sortperm][1:n_int_keep_top_seq]
    vec_int_chromosome_N = vec_int_sequence_N[vec_int_idx_sortperm][1:n_int_keep_top_seq]
    vec_int_sort_by_name = sortperm(vec_str_chromosome_names)
    vec_str_chromosome_names = String.(vec_str_chromosome_names[vec_int_sort_by_name])
    vec_int_chromosome_lengths = Int.(vec_int_chromosome_lengths[vec_int_sort_by_name])
    vec_int_chromosome_GC = Int.(vec_int_chromosome_GC[vec_int_sort_by_name])
    vec_int_chromosome_N = Int.(vec_int_chromosome_N[vec_int_sort_by_name])
    n_int_assembly_size = sum(vec_int_chromosome_lengths)
    n_int_assembly_GC = sum(vec_int_chromosome_GC)
    n_int_assembly_N = sum(vec_int_chromosome_N)
    ### Remove the uncateogised contigs
    for name in vec_str_sequence_names[vec_int_idx_sortperm][(n_int_keep_top_seq+1):end]
        try
            rm(string(name, ".fasta"))
            rm(string(name, "-GC_content.txt"))
        catch
            nothing
        end
    end
    ### Return assembly size, chromosome names, and lengths
    return(n_int_assembly_size, n_int_assembly_GC, n_int_assembly_N, vec_str_chromosome_names, vec_int_chromosome_lengths, vec_int_chromosome_GC, vec_int_chromosome_N)
end

@time n_int_assembly_size, 
      n_int_assembly_GC,
      n_int_assembly_N,
      vec_str_chromosome_names,
      vec_int_chromosome_lengths,
      vec_int_chromosome_GC, 
      vec_int_chromosome_N = fun_fasta_lengths_GC_N("Reference.fasta")

FILE = open("../Lope_on_Lori_filter_coords.txt", "r")
for i in 1:4
  header = readline(FILE)
end
vec_str_scaffold_names = []
while !eof(FILE)
  line = readline(FILE)
  push!(vec_str_scaffold_names, split(line, '\t')[end])
end
close(FILE)
vec_str_scaffold_names = unique(vec_str_scaffold_names)

vec_idx = []
@showprogress for x in vec_str_chromosome_names
  append!(vec_idx, sum(x .== vec_str_scaffold_names))
end

missed_scaffolds = sum(.!Bool.(vec_idx))                        ###         7,667 scaffolds
captured = sum(vec_int_chromosome_lengths[Bool.(vec_idx)])      ### 1,112,454,952 bp
no_captured = sum(vec_int_chromosome_lengths[.!Bool.(vec_idx)]) ###    15,589,113 bp
### We captured almost the entire L perenne assembly!
### Then, why are we getting less SNPs?!?!?!
### Is the:
###@@@  BAM=03_BAM/ACC09.bam ### ACC09 mapped to our new Lolium rigidum genome assembly
###@@@  BAMBASE=$(echo ${BAM%.bam*})
###@@@  time samtools stats ${BAMBASE}.bam > ${BAMBASE}.stat
###@@@  plot-bamstats ${BAMBASE}.stat -p ${BAMBASE}
###@@@  eog ${BAMBASE}-coverage.png
### good indcator of actual coverage and downstrean SNP depths?
```

## Mapping with bowtie2 instead of bwa mem
```{sh}
DIR=/data/Lolium/Quantitative_Genetics/
REF_LORI=${DIR}/02_FASTQ/REFERENCE/Reference.fasta
REF_LOPE=${DIR}/02_FASTQ/REFERENCE/BK_Lperenne/Reference.fasta
cd $DIR
wget https://github.com/BenLangmead/bowtie2/releases/download/v2.4.5/bowtie2-2.4.5-linux-x86_64.zip
unzip bowtie2-2.4.5-linux-x86_64.zip

BOWTIE_BUILD_INDEX=bowtie2-2.4.5-linux-x86_64/bowtie2-build
${BOWTIE_BUILD_INDEX} ${REF_LORI} ${REF_LORI%.fasta*}
${BOWTIE_BUILD_INDEX} ${REF_LOPE} ${REF_LOPE%.fasta*}

BOWTIE=bowtie2-2.4.5-linux-x86_64/bowtie2
FNAME_READ1=${DIR}/02_FASTQ/ACC01_combined_R1.fastq.gz
FNAME_READ2=${DIR}/02_FASTQ/ACC01_combined_R2.fastq.gz
BAM=${DIR}/ACC01.bam
# BAM=${DIR}/ACC01_LOPE.bam
MAPQ=20
# ${BOWTIE} \
#     -x ${REF_LOPE%.fasta*} \
time \
${BOWTIE} \
    -x ${REF_LORI%.fasta*} \
    -1 ${FNAME_READ1} \
    -2 ${FNAME_READ2} | \
    samtools view -q ${MAPQ} -b | \
    samtools sort > ${BAM}

samtools stat ${BAM} > ${BAM%.bam*}.stat
plot-bamstats ${BAM%.bam*}.stat -p ${BAM%.bam*}
eog ${BAM%.bam*}-coverage.png
```


