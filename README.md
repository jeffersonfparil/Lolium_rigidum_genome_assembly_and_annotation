# *Lolium rigidum* genome assembly, annotation, and comparative genomics
Assembly of the *Lolium rigidum*, L. genome, genome annotations, and comparative genomics

|**Lab Website**|**CSIRO Website**|**License**|
|:-------:|:--------:|:--------:|
| <a href="https://adaptive-evolution.biosciences.unimelb.edu.au/"><img src="misc/Adaptive Evolution Logo mod.png" width="150"> | <a href="https://www.csiro.au/en/"><img src="misc/CSIRO_Logo.svg" width="70"> | [![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0) |

Our objective is to generate a better genome assembly and genome annotations than that of the closely-related forage crop species [*Lolium perenne*](https://onlinelibrary.wiley.com/doi/full/10.1111/tpj.13037).

## NCBI WGS Accession identifiers

| SUBID | BioProject | BioSample | Accession | Organism |
|:-----:|:----------:|:---------:|:---------:|:--------:|
| SUB10975047 | PRJNA799061 | SAMN25144995 | JAKKIG000000000 | Lolium rigidum FL_2022 |

This Whole Genome Shotgun project has been deposited at DDBJ/ENA/GenBank under the accession **JAKKIG000000000**. The version described in this paper is version **JAKKIG010000000**. Notice the **01** in JAKKIG**01**0000000 to denote the version - probably?

## Sequencing
- Illumina Hiseq X
- MinION and PromethION
- Hi-C (Illumina Hiseq X)

## Genome assembly

### CSIRO pipeline
Short- and long-read DNA sequence data were generated and scaffolded using Hi-C sequence information. Short-read sequencing libraries were constructed using NEBNext Ultra II DNA Library Prep kit for Illumina (NEB, USA) and sequenced using HiSeq X platform (Illumina, Inc., San Diego, USA) ran in 150-bp paired-end mode. Adapter sequences were removed from the resulting reads using TrimGalore (v 0.6.6). Long read sequencing was carried out on MinION and PromethION platforms. Basecalling was performed using guppy (v5.1; Wick et al, 2019) under the dna_r9.4.1_450bps_sup.cfg model. The long-read sequences were trimmed using Porechop (v0.2.4; Wick et al, 2017) and filtered using filtlong (v0.2.1) to obtain high quality reads. The long-reads were assembled using Flye (v2.9; Kolmogorov et al, 2020) with the minimum overlap parameter set to 6,000, kmer size of 17, genome size of 2.25 Gb, and with no scaffolding. Duplicate contigs were purged using purge_dups (v1.2.5; Guan et al, 2020) with the default settings. The long-reads were error-corrected and trimmed using Canu (v2.2; Koren et al, 2017) under default settings, and used in three rounds of contig polishing using Racon (v1.4.22; Vaser et al, 2017) under default settings. This was followed by three rounds of short-read-based polishing using Polca (MaSURCA v4.0.7; Zimin et al, 2013)  to obtain the final contig assembly using default settings. This assembly was assessed using BUSCO (v5; Simão et al, 2015) against the Viridiplantae and Poales lineages’ gene sets.
 
A Hi-C library was prepared using 20 mg of leaf tissue and the Arima HiC kit following the manufacturer’s instructions. The library was sequenced on Novaseq 6000 platform (Illumina, Inc., San Diego, USA) to generate 500 million reads. The final contig assembly was scaffolded based on the genomic topological information using ALLHiC (v1; Zhang et al, 2019) and manually curated using JuiceBox (v1.9.8; Dudchenko et al, 2017).


### Genome statistics and CIRCOS-like figure
1. Install Julia and clone this repository
```shell
DIR=/data/Lolium_rigidum_ASSEMBLY/GENOME_ASSEMBLY
cd $DIR
sudo apt install -y julia
git clone https://github.com/jeffersonfparil/Lolium_rigidum_genome_assembly_and_annotation.git
```

2. Install Julia packages
```julia
using Pkg
Pkg.add(["Plots", "DataFrames", "CSV", "ProgressMeter", "JLD2"])
```

3. Some statistics first
```shell
### Count the number of each feature
grep -v "^#" Lolium_rigidum.gff | cut -f 3 | sort | uniq -c > Lolium_rigidum_annotation_counts.txt
### Extract the coordinates of genes
awk '{if ($3=="gene") print}' Lolium_rigidum.gff | cut -f4-5 > Lolium_rigidum_annotation_gene_coordinates.txt
### Find the mean length of gene models
echo 'dat = read.table("Lolium_rigidum_annotation_gene_coordinates.txt", header=FALSE)
MEAN_GENE_MODEL_LENGTH = mean(apply(dat, MARGIN=1, FUN=function(x){abs(diff(x))}))
print(MEAN_GENE_MODEL_LENGTH)
' > find_mean_gene_model_length.R
Rscript find_mean_gene_model_length.R
```

4. Plot genome assembly circos diagram
```julia
include("Lolium_rigidum_genome_assembly_and_annotation/genome_statistics.jl")
function execute(; recompute=false, cleanup=false)
    str_filename_fasta = "Lolium_rigidum.fasta"
    str_fname_output_svg = "Lolium_rigidum.svg"
    str_filename_LTR_COPIA = "Lolium_rigidum-LTR_COPIA.csv"
    str_filename_LTR_GYPSY = "Lolium_rigidum-LTR_GYPSY.csv"
    str_filename_groupings_and_coordinates = "Lolium_rigidum-for_plotting.plg"
    vec_idx_groups_chr_pos = [4, 2, 3, 3]
    n = 7 ## haploid chromosome number
    n_int_tick_length_bp = 100*1e+6 # make this adjustable
    n_int_tick_label_size=12
    n_int_chrom_name_size=12
    # n_int_stats_label_size=10

    ### To (re)run the genome statistics computation or to use the saved genome statistics data
    str_filename_output_jld2 = string(join(split(str_filename_fasta, ".")[1:(end-1)], "."), "-statistics.jld2")
    if recompute
        ### Split into chromosomes, count assembly size, GC content, Ns, and GC fractiopn per line
        @time n_int_assembly_size,
            n_int_assembly_GC,
            n_int_assembly_N,
            vec_str_chromosome_names,
            vec_int_chromosome_lengths,
            vec_int_chromosome_GC,
            vec_int_chromosome_N = PlotGenome.fun_fasta_lengths_GC_N(str_filename_fasta)

        ### Some assembly stats including the pseudo-chromosomes and small contigs
        vec_int_idx_sort_by_decreasing_size = sortperm(vec_int_chromosome_lengths, rev=true)
        vec_str_chromosome_names = vec_str_chromosome_names[vec_int_idx_sort_by_decreasing_size]
        vec_int_chromosome_lengths = vec_int_chromosome_lengths[vec_int_idx_sort_by_decreasing_size]
        vec_int_cumsum_length = cumsum(vec_int_chromosome_lengths)
        vec_bool_idx_50th_percentile = vec_int_cumsum_length .>= (n_int_assembly_size*0.5)
        vec_bool_idx_90th_percentile = vec_int_cumsum_length .>= (n_int_assembly_size*0.9)
        n_int_size_of_largest_chromosome = maximum(vec_int_chromosome_lengths)
        str_largest_chromosome = vec_str_chromosome_names[vec_int_chromosome_lengths.==maximum(vec_int_chromosome_lengths)]
        int_n_chromosomes_whole_assembly = length(vec_str_chromosome_names)
        int_size_whole_assembly = sum(vec_int_chromosome_lengths)
        L50_whole_assembly = collect(1:length(vec_str_chromosome_names))[vec_bool_idx_50th_percentile][1]
        L90_whole_assembly = collect(1:length(vec_str_chromosome_names))[vec_bool_idx_90th_percentile][1]
        N50_whole_assembly = vec_int_chromosome_lengths[vec_bool_idx_50th_percentile][1]
        N90_whole_assembly = vec_int_chromosome_lengths[vec_bool_idx_90th_percentile][1]

        ### Using only the 7 pseudo-chromosomes, i.e. 7 largest sequences
        vec_str_chromosome_names = vec_str_chromosome_names[1:n]
        vec_int_chromosome_lengths = vec_int_chromosome_lengths[1:n]
        vec_int_chromosome_GC = vec_int_chromosome_GC[1:n]
        vec_int_chromosome_N = vec_int_chromosome_N[1:n]

        ### Some assembly stats after excluding the small contigs
        vec_int_idx_sort_by_decreasing_size = sortperm(vec_int_chromosome_lengths, rev=true)
        vec_str_chromosome_names = vec_str_chromosome_names[vec_int_idx_sort_by_decreasing_size]
        vec_int_chromosome_lengths = vec_int_chromosome_lengths[vec_int_idx_sort_by_decreasing_size]
        vec_int_cumsum_length = cumsum(vec_int_chromosome_lengths)
        vec_bool_idx_50th_percentile = vec_int_cumsum_length .>= (n_int_assembly_size*0.5)
        vec_bool_idx_90th_percentile = vec_int_cumsum_length .>= (n_int_assembly_size*0.9)
        n_int_size_of_largest_chromosome = maximum(vec_int_chromosome_lengths)
        str_largest_chromosome = vec_str_chromosome_names[vec_int_chromosome_lengths.==maximum(vec_int_chromosome_lengths)]
        int_n_chromosomes = n
        int_size = sum(vec_int_chromosome_lengths)
        L50 = collect(1:length(vec_str_chromosome_names))[vec_bool_idx_50th_percentile][1]
        L90 = collect(1:length(vec_str_chromosome_names))[vec_bool_idx_90th_percentile][1]
        N50 = vec_int_chromosome_lengths[vec_bool_idx_50th_percentile][1]
        N90 = vec_int_chromosome_lengths[vec_bool_idx_90th_percentile][1]

        ### Sort back according to chromosome names
        vec_int_idx_sort_by_name = sortperm(vec_str_chromosome_names)
        vec_str_chromosome_names = vec_str_chromosome_names[vec_int_idx_sort_by_name]
        vec_int_chromosome_lengths = vec_int_chromosome_lengths[vec_int_idx_sort_by_name]

        ### Save the genome statistics variables
        if save_temp_files_and_data
            vec_stats = [n_int_assembly_size,
                            n_int_assembly_GC,
                            n_int_assembly_N,
                            int_n_chromosomes_whole_assembly,
                            int_size_whole_assembly,
                            L50_whole_assembly,
                            L90_whole_assembly,
                            N50_whole_assembly,
                            N90_whole_assembly,
                            vec_str_chromosome_names,
                            vec_int_chromosome_lengths,
                            vec_int_chromosome_GC,
                            vec_int_chromosome_N,
                            str_largest_chromosome,
                            int_n_chromosomes,
                            int_size,
                            L50,
                            L90,
                            N50,
                            N90]
            jldsave(str_filename_output_jld2; vec_stats)
        end
    else
        ### If we're not rerunning the statistics computation, then just load the save statistics
        n_int_assembly_size,
        n_int_assembly_GC,
        n_int_assembly_N,
        int_n_chromosomes_whole_assembly,
        int_size_whole_assembly,
        L50_whole_assembly,
        L90_whole_assembly,
        N50_whole_assembly,
        N90_whole_assembly,
        vec_str_chromosome_names,
        vec_int_chromosome_lengths,
        vec_int_chromosome_GC,
        vec_int_chromosome_N,
        str_largest_chromosome,
        int_n_chromosomes,
        int_size,
        L50,
        L90,
        N50,
        N90 = load_object(str_filename_output_jld2)
    end

    ### Some assembly stats including the pseudo-chromosomes and small contigs
    vec_int_idx_sort_by_decreasing_size = sortperm(vec_int_chromosome_lengths, rev=true)
    vec_str_chromosome_names = vec_str_chromosome_names[vec_int_idx_sort_by_decreasing_size]
    vec_int_chromosome_lengths = vec_int_chromosome_lengths[vec_int_idx_sort_by_decreasing_size]
    vec_int_cumsum_length = cumsum(vec_int_chromosome_lengths)
    vec_bool_idx_50th_percentile = vec_int_cumsum_length .>= (n_int_assembly_size*0.5)
    vec_bool_idx_90th_percentile = vec_int_cumsum_length .>= (n_int_assembly_size*0.9)
    n_int_size_of_largest_chromosome = maximum(vec_int_chromosome_lengths)
    str_largest_chromosome = vec_str_chromosome_names[vec_int_chromosome_lengths.==maximum(vec_int_chromosome_lengths)]
    int_n_chromosomes_whole_assembly = length(vec_str_chromosome_names)
    int_size_whole_assembly = sum(vec_int_chromosome_lengths)
    L50_whole_assembly = collect(1:length(vec_str_chromosome_names))[vec_bool_idx_50th_percentile][1]
    L90_whole_assembly = collect(1:length(vec_str_chromosome_names))[vec_bool_idx_90th_percentile][1]
    N50_whole_assembly = vec_int_chromosome_lengths[vec_bool_idx_50th_percentile][1]
    N90_whole_assembly = vec_int_chromosome_lengths[vec_bool_idx_90th_percentile][1]

    ### Using only the 7 pseudo-chromosomes, i.e. 7 largest sequences
    vec_str_chromosome_names = vec_str_chromosome_names[1:n]
    vec_int_chromosome_lengths = vec_int_chromosome_lengths[1:n]
    vec_int_chromosome_GC = vec_int_chromosome_GC[1:n]
    vec_int_chromosome_N = vec_int_chromosome_N[1:n]

    ### Some assembly stats after excluding the small contigs
    vec_int_idx_sort_by_decreasing_size = sortperm(vec_int_chromosome_lengths, rev=true)
    vec_str_chromosome_names = vec_str_chromosome_names[vec_int_idx_sort_by_decreasing_size]
    vec_int_chromosome_lengths = vec_int_chromosome_lengths[vec_int_idx_sort_by_decreasing_size]
    vec_int_cumsum_length = cumsum(vec_int_chromosome_lengths)
    vec_bool_idx_50th_percentile = vec_int_cumsum_length .>= (n_int_assembly_size*0.5)
    vec_bool_idx_90th_percentile = vec_int_cumsum_length .>= (n_int_assembly_size*0.9)
    n_int_size_of_largest_chromosome = maximum(vec_int_chromosome_lengths)
    str_largest_chromosome = vec_str_chromosome_names[vec_int_chromosome_lengths.==maximum(vec_int_chromosome_lengths)]
    int_n_chromosomes = n
    int_size = sum(vec_int_chromosome_lengths)
    L50 = collect(1:length(vec_str_chromosome_names))[vec_bool_idx_50th_percentile][1]
    L90 = collect(1:length(vec_str_chromosome_names))[vec_bool_idx_90th_percentile][1]
    N50 = vec_int_chromosome_lengths[vec_bool_idx_50th_percentile][1]
    N90 = vec_int_chromosome_lengths[vec_bool_idx_90th_percentile][1]

    ### Sort back according to chromosome names
    vec_int_idx_sort_by_name = sortperm(vec_str_chromosome_names)
    vec_str_chromosome_names = vec_str_chromosome_names[vec_int_idx_sort_by_name]
    vec_int_chromosome_lengths = vec_int_chromosome_lengths[vec_int_idx_sort_by_name]

    ### Base plot
    l = @layout [a b]
    plt1 = plot(xlim=(-1.5, 1.5), ylim=(-1.5, 1.5), size=(700,700), axis=([], false), title="")
    ### Plot 1 for the tick alone, i.e. will be occluded by the GC layer
    r=1.00; w=0.15
    PlotGenome.fun_plot_chrom_length_layer!(plt1, vec_str_chromosome_names, vec_int_chromosome_lengths;
                                r=r, w=w,
                                n_int_tick_length_bp=n_int_tick_length_bp,
                                n_int_tick_label_size=n_int_tick_label_size,
                                n_int_chrom_name_size=n_int_chrom_name_size,
                                add_legend=false)
    ### Plot 1 layer 1: GC content
    r=1.00; w=0.15
    annotate!(plt1, 0.0, (r-w/2), ("a", 10, :gray, :center))
    PlotGenome.fun_plot_GC_layer!(plt1, vec_str_chromosome_names, vec_int_chromosome_lengths;
                    r=r, w=w,
                    n_int_total_chunks_across_genome=1000,
                    vec_colours_GC=palette(:thermometer, 7),
                    n_int_tick_label_size=n_int_tick_label_size,
                    lane_id="a")
    int_n_chromosomes_whole_assembly = PlotGenome.fun_add_comma_separator_on_large_positive_integers(int_n_chromosomes_whole_assembly)
    int_size_whole_assembly = PlotGenome.fun_add_comma_separator_on_large_positive_integers(int_size_whole_assembly)
    int_n_chromosomes = PlotGenome.fun_add_comma_separator_on_large_positive_integers(int_n_chromosomes)
    int_size = PlotGenome.fun_add_comma_separator_on_large_positive_integers(int_size)
    L50 = PlotGenome.fun_add_comma_separator_on_large_positive_integers(L50)
    N50 = PlotGenome.fun_add_comma_separator_on_large_positive_integers(N50)
    L90 = PlotGenome.fun_add_comma_separator_on_large_positive_integers(L90)
    N90 = PlotGenome.fun_add_comma_separator_on_large_positive_integers(N90)
    # annotate!(plt1, -1.5, 1.25,(string("Whole assembly:\n  n=", int_n_chromosomes_whole_assembly,
    #                                 "\n  size=", int_size_whole_assembly, " bp",
    #                                 "\n",
    #                                 "\nPseudo-chromosomes:\n  n=", int_n_chromosomes,
    #                                 "\n  size=", int_size, " bp",
    #                                 "\n  L50=", L50, " chromosomes",
    #                                 "\n  N50=", N50, " bp",
    #                                 "\n  L90=", L90, " chromosomes",
    #                                 "\n  N90=", N90, " bp"), n_int_stats_label_size, :gray, :left))
    ### Plot 1 layer 2: Ty1-Copia LTR histogram
    r=0.75; w=0.15
    annotate!(plt1, 0.0, (r-w/2), ("b", 10, :gray, :center))
    PlotGenome.fun_plot_hits_histogram_layer!(plt1, vec_str_chromosome_names, vec_int_chromosome_lengths,
                                str_filename_LTR_COPIA,
                                r=r, w=w,
                                n_int_window_size = 1e6,
                                col=:black,
                                col_background=:lightgray)
    ### Plot 1 layer 3: Ty1-Gypsy LTR histogram
    r=0.50; w=0.15
    annotate!(plt1, 0.0, (r-w/2), ("c", 10, :gray, :center))
    PlotGenome.fun_plot_hits_histogram_layer!(plt1, vec_str_chromosome_names, vec_int_chromosome_lengths,
                                str_filename_LTR_GYPSY,
                                r=r, w=w,
                                n_int_window_size = 1e6,
                                col=:gray,
                                col_background=:lightgray)
    ### Plot 2: Chord diagram
    plt2 = plot(xlim=(-1.5, 1.5), ylim=(-1.5, 1.5), size=(700,700), axis=([], false), title="")
    r = 1.0; w=0.10
    PlotGenome.fun_plot_chrom_length_layer!(plt2, vec_str_chromosome_names, vec_int_chromosome_lengths;
                            r=r, w=w,
                            n_int_tick_length_bp=n_int_tick_length_bp,
                            n_int_tick_label_size=n_int_tick_label_size,
                            n_int_chrom_name_size=n_int_chrom_name_size,
                            add_legend=false)
                            r = r-w; w=0.05
    PlotGenome.fun_add_chords!(plt2,
                    str_filename_groupings_and_coordinates,
                    vec_str_chromosome_names,
                    vec_int_chromosome_lengths,
                    delim='\t',
                    vec_idx_groups_chr_pos=vec_idx_groups_chr_pos,
                    colour_per_chrom=true,
                    linewidth=1,
                    r=r, w=w, header=false)
    plt3 = plot(plt1, plt2, layout=l, size=(1400,700))
    ### Clean-up
    if cleanup
        vec_files = readdir()
        idx_fasta = match.(Regex(".fasta"), vec_files) .!= nothing
        idx_GC = match.(Regex("-GC_content.txt"), vec_files) .!= nothing
        idx_not_input_fasta = match.(Regex(str_filename_fasta), vec_files) .== nothing
        for f in vec_files[(idx_fasta .& idx_not_input_fasta) .| idx_GC]
            rm(f)
        end
    end

    ### Save as svg
    savefig(plt3, str_fname_output_svg)
end
@time execute()
### MISC: C-value into bp
_2C_in_pg_ = 5.494
0.978e9 * (_2C_in_pg_/2)
```

### Genome size estimation with jellyfish and GenomeScope
```shell
sudo apt intall -y jellyfish
git clone https://github.com/schatzlab/genomescope.git

for f in $(find /data/Lolium_rigidum_ASSEMBLY/lolium_illumina/LOL-WGS2/N2009012_FA_30-434329113_SEQ/201110-X4A_L007 -name 'LOL-WGS2-*.fastq.gz')
do
f_new=$(basename $f)
f_new=${f_new%.gz*}
gunzip --keep $f --stdout > $f_new
done
readlen=$(head -n 10 $(ls LOL-WGS2-*.fastq | head -n1) | grep -A1 "^@" | head -n2 | tail -n1 | sed 's/N//g' | wc -c)
time \
for kmer in $(seq 15 25)
do
echo "###################################"
echo ${kmer}
# kmer=21
jellyfish count \
    --mer-len ${kmer} \
    --size 1000000000 \
    --threads 10 \
    --output kmer-${kmer}.jf \
    LOL-WGS2-*.fastq
jellyfish histo \
    kmer-${kmer}.jf \
    --threads 10 \
    --output kmer-${kmer}.hist
Rscript genomescope/genomescope.R \
    kmer-${kmer}.hist \
    ${kmer} \
    ${readlen} \
    GenomeScope_OUT-kmer-${kmer}
mv kmer-* GenomeScope_OUT-kmer-${kmer}
done
### Result: at kmer=16 the error was minimised:which results to 2,261,492,244 (2.26 Gb)
```

### Assess the contitguity of the assembly using the ratio of intact LTR-RT (LTR retrotransposons) and the total LTR-RT
1. Download and index the reference genome
```shell
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/022/539/505/GCF_022539505.1_APGP_CSIRO_Lrig_0.1/GCF_022539505.1_APGP_CSIRO_Lrig_0.1_genomic.fna.gz
gunzip -c GCF_022539505.1_APGP_CSIRO_Lrig_0.1_genomic.fna.gz > Lolium_rigidum.fasta

time \
gt suffixerator \
   -db Lolium_rigidum.fasta \
   -indexname Lolium_rigidum.fasta \
   -tis -suf -lcp -des -ssp -sds -dna
```

2. Identify the LTRs and compute the LAIs
```shell
time \
gt ltrharvest \
   -index Lolium_rigidum.fasta \
   -similar 90 -vic 10 -seed 20 -seqids yes \
   -minlenltr 100 -maxlenltr 10000 -mintsd 4 -maxtsd 6 \
   -motif TGCA -motifmis 1 > Lolium_rigidum_GENOME_HARVEST.scn

time \
/usr/local/LTR_retriever-2.9.0/LTR_retriever \
   -genome Lolium_rigidum.fasta \
   -inharvest Lolium_rigidum_GENOME_HARVEST.scn \
   -threads 15
```

3. Assess the genome-wide LAI distribution
```R
dat = read.table("Lolium_rigidum.fasta.mod.out.LAI", header=TRUE)
# str(dat)
# dat$len = dat$To - dat$From
X = droplevels(dat[grepl("NC", dat$Chr), ])
X$Chr = as.factor(X$Chr)
chr = unique(X$Chr)

svg("Lolium_rigidum-LAI_distribution.svg", width=10, height=7)
par(mfrow=c(2,2))
hist(X$raw_LAI, xlab="Raw LAI", main="Raw LAI")
legend("topleft", legend=paste0("µ=", round(dat$raw_LAI[1],2)), bty="n")
hist(X$LAI, xlab="Standardised LAI", main="Standardised LAI")
legend("topleft", legend=paste0("µ=", round(dat$LAI[1],2)), bty="n")
colours = c("#e41a1c", "#377eb8", "#4daf4a", "#984ea3", "#ff7f00", "#ffff33", "#a65628")[as.numeric(X$Chr)]
xticks = aggregate(c(1:nrow(X)) ~ X$Chr, FUN=median)
plot(X$raw_LAI, col=colours, xaxt="n", ylab="Raw LAI", xlab="Chromosome")
grid()
axis(side=1, at=xticks[,2], labels=c(1:7), tick=FALSE)
plot(X$LAI, col=colours, xaxt="n", ylab="Standardized LAI", xlab="Chromosome")
grid()
axis(side=1, at=xticks[,2], labels=c(1:7), tick=FALSE)
dev.off()
```


### Map the filtered RNAseq data into the assembly
```shell
DIR=/data/Lolium_rigidum_ASSEMBLY/TRANSCRIPTOME_ASSEMBLY
DIR_READS=/data/Lolium_rigidum_ASSEMBLY/TRANSCRIPTOME_ASSEMBLY/rRNAdepleted_reads
DIR_REF=/data/Lolium_rigidum_ASSEMBLY/GENOME_ASSEMBLY
cd $DIR
### Install Hisat2
git clone https://github.com/DaehwanKimLab/hisat2.git
cd hisat2
make
PATH=${PATH}:$(pwd)
sudo ln -s /usr/bin/python3 /usr/bin/python ### hisat2 needs python and python2 is now fully deprecated and we may have to specify that the default python is python3
hisat2 -h
cd $DIR
### Prepare the reference genome index
REF=${DIR_REF}/Lolium_rigidum.fasta
time \
hisat2-build \
    ${REF} \
    ${REF%.fasta*}
### Align
echo '#!/bin/bash
R1=$1
R2=$2
REF=$3
f=$(basename $R1)
SAMOUT=${f%.fq*}.sam
BAMOUT=${f%.fq*}.bam
hisat2 \
    -x ${REF} \
    -1 ${R1} \
    -2 ${R2} \
    -S ${SAMOUT}
samtools \
    sort ${SAMOUT} \
    -O BAM \
    > ${BAMOUT}
' > map_RNAseq.sh
chmod +x map_RNAseq.sh
time \
parallel \
--link \
--jobs 11 \
./map_RNAseq.sh \
    {1} \
    {2} \
    ${DIR_REF}/Lolium_rigidum \
    ::: $(find ${DIR_READS} -name '*.fq.1.gz' | sort) \
    ::: $(find ${DIR_READS} -name '*.fq.2.gz' | sort)
### Assess mapping rate
for f in $(ls *.bam)
do
samtools flagstat $f > ${f%.bam*}.stat
done
grep "mapped (" *.stat | cut -d'(' -f2 | cut -d':' -f1 | sort | sed 's/%//g' > mapping_rate.stat
mean_mapping_rate=$(echo "scale=2; (" $(cat mapping_rate.stat | sed -z 's/\n/+/g' | sed 's/+$//g') ") / " $(cat mapping_rate.stat | wc -l) | bc)
echo "MEAN MAPPING RATE = " ${mean_mapping_rate} "%"
```

## Transcriptome assembly
```shell
##########################################################################################################
#
#                                          Transcriptome assembly pipeline
#                              (based on Trinity best practice from Harvard FAS)
# https://informatics.fas.harvard.edu/best-practices-for-de-novo-transcriptome-assembly-with-trinity.html
##########################################################################################################

R1=`ls /data/Lolium_rigidum_ASSEMBLY/jeff_lolium/FASTQ/ILLUMINA/2020.9/30-409520043/N2006601_FA_30-409520043_SEQ/200908-X4A_L005/*R1* | grep -v determined`

for f in $R1; do 
    f_dest=`echo $f | cut -d'/' -f11 | sed 's/.gz//'`
    gunzip -c $f > ./${f_dest} &
    gunzip -c ${f/R1/R2} > ./${f_dest/R1/R2} 
done

R1=`ls *R1*fastq | grep -v WGS`
R1=`echo $R1 | sed 's/ /,/g'`
R2=`echo $R1 | sed 's/R1/R2/g'`

#How to count number of reads in fastq
for i in *_R1.fastq; do
    echo $(cat $i | wc -l)/4|bc
done

#Error correction and flagging.(https://github.com/mourisl/Rcorrector)
perl /data/Software/rcorrector/run_rcorrector.pl -t 48 -1 $R1 -2 $R2

#Filter out uncorrectable reads (https://github.com/harvardinformatics/TranscriptomeAssemblyTools)
R1=`ls *R1*cor.fq | grep -v WGS`
R2=`echo $R1 | sed 's/R1/R2/g'`

parallel -j 12 --link 'python2.7 /data/Software/TranscriptomeAssemblyTools/FilterUncorrectabledPEfastq.py -1 {1} -2 {2} -s {3}' ::: `echo $R1` ::: `echo $R2` ::: `echo $R1 | sed 's/_combined_R1.cor.fq/_filter.log/g'`

#Renaming/Overwritting for space keeping
for i in unfixrm_*_combined_R*.cor.fq; do
	mv $i ${i/unfixrm_/}
done

#Trim adapter and low quality bases from fastq files
R1=`ls *R1*cor.fq | grep -v WGS`
R2=`echo $R1 | sed 's/R1/R2/g'`

#The command given by the workflow
#./data/Software/TrimGalore-0.6.0/trim_galore --paired --retain_unpaired --phred33 --output_dir trimmed_reads --length 36 -q 5 --stringency 1 -e 0.1 $1 $2

parallel -j 12 --link '/data/Software/TrimGalore-0.6.0/trim_galore --paired --retain_unpaired --phred33 --output_dir trimmed_reads --length 36 -q 5 --stringency 1 -e 0.1 {1} {2}' ::: `echo $R1` ::: `echo $R2`



#remove rRNA using mapping mapping strategy to SILVA database
# rRNA fasta were sourced from SILVA (https://www.arb-silva.de/no_cache/download/archive/release_138_1/Exports/) SSU+SSUParc

#in RefGenomes directory
zcat SILVA_138.1_* | sed 's/U/T/g' > SILVA_138.1_LSU+SSU.fasta
#index the database
bowtie2-build SILVA_138.1_LSU+SSU_noIUPAC.fasta SILVA_138.1_LSU+SSU_noIUPAC

#from trimmed_reads dir
for i in *trimmed_R1.fq; do
	bowtie2 --quiet --very-sensitive-local --phred33  \
	        -x /data/RefGenomes/SILVA_rRNA/SILVA_138.1 \
	        -1 $i -2 ${i/R1/R2} \
	        --threads 48 \
	        --met-file ../rRNAdepleted_reads/${i/_trimmed_R1.fq/}_bowtie2_metrics.txt \
	        --al-conc-gz ../rRNAdepleted_reads/${i/_trimmed_R1.fq/}_paired_aligned.fq.gz --un-conc-gz ../rRNAdepleted_reads/${i/_trimmed_R1.fq/}_paired_unaligned.fq.gz  \
	        --al-gz ../rRNAdepleted_reads/${i/_trimmed_R1.fq/}_unpaired_aligned.fq.gz --un-gz ../rRNAdepleted_reads/${i/_trimmed_R1.fq/}_unpaired_unaligned.fq.gz \
	        > 
done


# from trinity_assembly/ dir
cat ../*1_paired_unaligned.fq.gz > ALL_TISSUES_R1.fq.gz
cat *2_paired_unaligned.fq.gz > ALL_TISSUES_R2.fq.gz

#Number of reads that were used for assembly
echo $(zcat ALL_TISSUES_R1.fq.gz | wc -l)/4|bc


for i in `ls ../rRNAdepleted_reads/*1_paired_unaligned.fq.1* | cut -d- -f1`; do 
	out_prefix=`echo $i | cut -d'/' -f3`
	#cat ${i}*paired_unaligned.fq.1.gz > ${out_prefix}_R1.fq.gz &
	#cat ${i}*paired_unaligned.fq.2.gz > ${out_prefix}_R2.fq.gz
	Trinity --seqType fq --max_memory 100G \
	        --left ${out_prefix}_R1.fq.gz  --right ${out_prefix}_R2.fq.gz \
	        --CPU 24 --output ${out_prefix}_trinity
done


#BUSCO FOR COMPLETENESS
source /data/miniconda3/bin/activate busco

busco -m transcriptome -i ../trinity_assembly/trinity_output.Trinity.fasta \
      -c 24 -o ALL_TISSUES -l poales_odb10

#Compacting an assembly with DRAP (http://www.sigenae.org/drap/index.html)
#Pull the docker image
sudo docker pull sigenae/drap
#create a container
sudo docker create --name drap --privileged -v /data/Lolium_rigidum_ASSEMBLY:/docker/scratch -i -t sigenae/drap:latest /bin/bash
#run the container
docker start drap
docker exec -i -t drap /bin/bash

#For the All tissues run
/usr/local/src/drap/runDrap  \
        --R1 ../trinity_assembly/ALL_TISSUES_R1.fq.gz  --R2 ../trinity_assembly/ALL_TISSUES_R2.fq.gz  \
        --ref Oryza_sativa.pep.fasta  --dbg trinity  --outdir /docker/scratch/drap_out/ALL_TISSUES_drap_wTrinity \
        --dbg-mem 24 --norm-mem 24

/usr/local/src/drap/runDrap  \
        --R1 ../trinity_assembly/ALL_TISSUES_R1.fq.gz  --R2 ../trinity_assembly/ALL_TISSUES_R2.fq.gz  \
        --ref Oryza_sativa.pep.fasta  --dbg oases  --outdir /docker/scratch/drap_out/ALL_TISSUES_drap_wOases \
        --dbg-mem 24 --norm-mem 24

#for the tissue specific run
for i in `ls ../trinity_assembly/*_R1.fq.gz | grep -v ALL_`; do
    j=`echo $i | cut -d"/" -f3 | cut -d_ -f1`
	/usr/local/src/drap/runDrap  \
	        --R1 $i  --R2 ${i/R1/R2}  \
	        --ref Oryza_sativa.pep.fasta  --dbg oases  --outdir /docker/scratch/drap_out/${j}_wOases \
	        --dbg-mem 24 --norm-mem 24
done


#from outside the container to copy
cp ../drap_out/*fa ../TRANSCRIPTOME_ASSEMBLY/docker_outpoint/drap_wTrinity/


# The meta-assembler to merge tissue specific runs of Trinity + Oases
for i in `ls ../trinity_assembly/*_R1.fq.gz | grep -v ALL_`; do
    j=`echo $i | cut -d"/" -f3 | cut -d_ -f1`
	/usr/local/src/drap/runMeta  \
	        --drap-dirs /docker/scratch/drap_out/${j}_wOases,/docker/scratch/drap_out/${j}_wTrinity \
	        --ref Oryza_sativa.pep.fasta \
	        --outdir /docker/scratch/drap_out/${j}_Meta
done

# The meta assembly across all tisues
j=`ls -d /docker/scratch/drap_out/*_Meta/ | sed 's/ /,/g'`
/usr/local/src/drap/runMeta  \
	        --drap-dirs `echo $j | sed 's/\/ /,/g'` \
	        --ref Oryza_sativa.pep.fasta \
	        --outdir /docker/scratch/drap_out/ALLTISSUES_Meta

#And the meta of the ALL tissues with Trinity+Oases
/usr/local/src/drap/runMeta  \
	        --drap-dirs /docker/scratch/drap_out/ALL_TISSUES_drap_wTrinity,/docker/scratch/drap_out/ALL_TISSUES_drap_wOases \
	        --ref Oryza_sativa.pep.fasta \
	        --outdir /docker/scratch/drap_out/ALLTISSUES_drap_Trinity+Oases

#Busco for completneess
#activate the conda environment to get salmon et dependencies loaded for Trinity
source /data/miniconda3/bin/activate busco

for i in ../docker_outpoint/*/*fa; do
    j=`echo $i | cut -d"/" -f3-4 | sed 's/\//_/g' | sed 's/.fa//g'`
    busco -m transcriptome -i $i -c 24 -o $j -l poales_odb10
done
```

## Annotation
The reads were demultiplexed and  error-corrected using Rcorrector (v1.0.4; Song et al, 2015). Adapters and low quality base pairs were trimmed using TrimGalore (v0.6.0). Ribosomal RNA sequences were discarded when one of the paired-end reads mapped to the sequences present in the SILVA database (v138.1; Quast et al, 2012) using Bowtie2 (v2.3; Langmead & Salzberg, 2012). After filtering, ~197 million reads were used for de novo transcriptome assembly using the De novo RNA-Seq Assembly Pipeline (Cabau et al, 2017) including the rice protein sequences (release 51) as guide and using both Trinity (v2.8.4; Haas et al, 2013) and Oases (v0.2.09; Schulz et al, 2012) as assemblers. The resulting two assemblies were merged into a single compacted meta-assembly. The filtered reads were re-mapped against the meta-assembly and transcripts with FPKM>1 were included in the transcriptome.
The genome was annotated using NCBI’s genome annotation pipeline using the de novo assembled transcriptome.
Transposable elements were identified using RepeatMasker and RepeatModeller (v4.1.2 and v2.0.3, respectively; Flynn et al, 2020).

## Comparative genomics

### Set working directories, executables, and output files
```shell
DIR=/data/Lolium_rigidum_ASSEMBLY/COMPARATIVE_GENOMICS
SRC=${DIR}/Lolium_rigidum_genome_assembly_and_annotation
PATH=${PATH}:${DIR}/OrthoFinder
PATH=${PATH}:${DIR}/CAFE5/bin
MACSE=${DIR}/MACSE/macse_v2.06.jar
PATH=${PATH}:${DIR}/iqtree-2.0.7-Linux/bin
PATH=${PATH}:${DIR}/paml4.9j/bin
PATH=${PATH}:${DIR}/kakscalculator2-2.0.1/src

DIR_ORTHOGROUPS=${DIR}/ORTHOGROUPS/OrthoFinder/Results_May13
DIR_ORTHOGROUP_SEQS=${DIR}/ORTHOGROUPS/OrthoFinder/Results_*/Orthogroup_Sequences
DIR_PANTHER=${DIR}/PantherHMM_17.0/famlib/rel/PANTHER17.0_altVersion/hmmscoring/PANTHER17.0/books
GOT_PATHER=${DIR}/PantherHMM_17.0/PANTHER17.0_HMM_classifications
MERGED_ORTHOGROUPS=${DIR}/ORTHOGROUPS/orthogroups.faa
ORTHOUT=${DIR}/ORTHOGROUPS/orthogroups_gene_counts_families_go.out
TREE=${DIR_ORTHOGROUPS}/Species_Tree/SpeciesTree_rooted.txt
DIR_ORTHOGROUP_SEQS=${DIR_ORTHOGROUPS}/Orthogroup_Sequences
DIR_GENES=${DIR}/TSR_NTSR_GENES
cd $DIR
```

### Download genomes, genome annotations, and predicted CDS sequences
**Note:** The first column of the annotation files may not correspond to the chromosome IDs.
```shell
### Arabidopsis thaliana
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/735/GCF_000001735.4_TAIR10.1/GCF_000001735.4_TAIR10.1_genomic.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/735/GCF_000001735.4_TAIR10.1/GCF_000001735.4_TAIR10.1_genomic.gff.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/735/GCF_000001735.4_TAIR10.1/GCF_000001735.4_TAIR10.1_cds_from_genomic.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/735/GCF_000001735.4_TAIR10.1/GCF_000001735.4_TAIR10.1_protein.faa.gz
gunzip -c GCF_000001735.4_TAIR10.1_genomic.fna.gz > Arabidopsis_thaliana.fasta
gunzip -c GCF_000001735.4_TAIR10.1_genomic.gff.gz > Arabidopsis_thaliana.gff
gunzip -c GCF_000001735.4_TAIR10.1_cds_from_genomic.fna.gz > Arabidopsis_thaliana.cds
gunzip -c GCF_000001735.4_TAIR10.1_protein.faa.gz > Arabidopsis_thaliana.faa
### Oryza sativa
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/001/433/935/GCF_001433935.1_IRGSP-1.0/GCF_001433935.1_IRGSP-1.0_genomic.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/001/433/935/GCF_001433935.1_IRGSP-1.0/GCF_001433935.1_IRGSP-1.0_genomic.gff.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/001/433/935/GCF_001433935.1_IRGSP-1.0/GCF_001433935.1_IRGSP-1.0_cds_from_genomic.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/001/433/935/GCF_001433935.1_IRGSP-1.0/GCF_001433935.1_IRGSP-1.0_protein.faa.gz
gunzip -c GCF_001433935.1_IRGSP-1.0_genomic.fna.gz > Oryza_sativa.fasta
gunzip -c GCF_001433935.1_IRGSP-1.0_genomic.gff.gz > Oryza_sativa.gff
gunzip -c GCF_001433935.1_IRGSP-1.0_cds_from_genomic.fna.gz > Oryza_sativa.cds
gunzip -c GCF_001433935.1_IRGSP-1.0_protein.faa.gz > Oryza_sativa.faa
### Zea mays
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/902/167/145/GCF_902167145.1_Zm-B73-REFERENCE-NAM-5.0/GCF_902167145.1_Zm-B73-REFERENCE-NAM-5.0_genomic.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/902/167/145/GCF_902167145.1_Zm-B73-REFERENCE-NAM-5.0/GCF_902167145.1_Zm-B73-REFERENCE-NAM-5.0_genomic.gff.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/902/167/145/GCF_902167145.1_Zm-B73-REFERENCE-NAM-5.0/GCF_902167145.1_Zm-B73-REFERENCE-NAM-5.0_cds_from_genomic.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/902/167/145/GCF_902167145.1_Zm-B73-REFERENCE-NAM-5.0/GCF_902167145.1_Zm-B73-REFERENCE-NAM-5.0_protein.faa.gz
gunzip -c GCF_902167145.1_Zm-B73-REFERENCE-NAM-5.0_genomic.fna.gz > Zea_mays.fasta
gunzip -c GCF_902167145.1_Zm-B73-REFERENCE-NAM-5.0_genomic.gff.gz > Zea_mays.gff
gunzip -c GCF_902167145.1_Zm-B73-REFERENCE-NAM-5.0_cds_from_genomic.fna.gz > Zea_mays.cds
gunzip -c GCF_902167145.1_Zm-B73-REFERENCE-NAM-5.0_protein.faa.gz > Zea_mays.faa
### Sorghum bicolor
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/003/195/GCF_000003195.3_Sorghum_bicolor_NCBIv3/GCF_000003195.3_Sorghum_bicolor_NCBIv3_genomic.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/003/195/GCF_000003195.3_Sorghum_bicolor_NCBIv3/GCF_000003195.3_Sorghum_bicolor_NCBIv3_genomic.gff.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/003/195/GCF_000003195.3_Sorghum_bicolor_NCBIv3/GCF_000003195.3_Sorghum_bicolor_NCBIv3_cds_from_genomic.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/003/195/GCF_000003195.3_Sorghum_bicolor_NCBIv3/GCF_000003195.3_Sorghum_bicolor_NCBIv3_protein.faa.gz
gunzip -c GCF_000003195.3_Sorghum_bicolor_NCBIv3_genomic.fna.gz > Sorghum_bicolor.fasta
gunzip -c GCF_000003195.3_Sorghum_bicolor_NCBIv3_genomic.gff.gz > Sorghum_bicolor.gff
gunzip -c GCF_000003195.3_Sorghum_bicolor_NCBIv3_cds_from_genomic.fna.gz > Sorghum_bicolor.cds
gunzip -c GCF_000003195.3_Sorghum_bicolor_NCBIv3_protein.faa.gz > Sorghum_bicolor.faa
### Lolium rigidum
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/022/539/505/GCF_022539505.1_APGP_CSIRO_Lrig_0.1/GCF_022539505.1_APGP_CSIRO_Lrig_0.1_genomic.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/022/539/505/GCF_022539505.1_APGP_CSIRO_Lrig_0.1/GCF_022539505.1_APGP_CSIRO_Lrig_0.1_genomic.gff.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/022/539/505/GCF_022539505.1_APGP_CSIRO_Lrig_0.1/GCF_022539505.1_APGP_CSIRO_Lrig_0.1_cds_from_genomic.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/022/539/505/GCF_022539505.1_APGP_CSIRO_Lrig_0.1/GCF_022539505.1_APGP_CSIRO_Lrig_0.1_protein.faa.gz
gunzip -c GCF_022539505.1_APGP_CSIRO_Lrig_0.1_genomic.fna.gz > Lolium_rigidum.fasta
gunzip -c GCF_022539505.1_APGP_CSIRO_Lrig_0.1_genomic.gff.gz > Lolium_rigidum.gff
gunzip -c GCF_022539505.1_APGP_CSIRO_Lrig_0.1_cds_from_genomic.fna.gz > Lolium_rigidum.cds
gunzip -c GCF_022539505.1_APGP_CSIRO_Lrig_0.1_protein.faa.gz > Lolium_rigidum.faa
### Lolium perenne
# wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/019/359/855/GCA_019359855.1_MPB_Lper_Kyuss_1697/GCA_019359855.1_MPB_Lper_Kyuss_1697_genomic.fna.gz
# gunzip -c GCA_019359855.1_MPB_Lper_Kyuss_1697_genomic.fna.gz > Lolium_perenne.fasta
wget https://de.cyverse.org/anon-files//iplant/home/shared/commons_repo/curated/Copetti_Kyuss_assembly_annotation_March_2021/Kyuss_1697_assembly.fa
wget https://de.cyverse.org/anon-files//iplant/home/shared/commons_repo/curated/Copetti_Kyuss_assembly_annotation_March_2021/Kyuss_1697_KYUS.gff
wget https://de.cyverse.org/anon-files//iplant/home/shared/commons_repo/curated/Copetti_Kyuss_assembly_annotation_March_2021/Kyuss_1697_KYUS_CDS.fa
wget https://de.cyverse.org/anon-files//iplant/home/shared/commons_repo/curated/Copetti_Kyuss_assembly_annotation_March_2021/Kyuss_1697_KYUS_proteins.fa
mv Kyuss_1697_assembly.fa Lolium_perenne.fasta
mv Kyuss_1697_KYUS.gff Lolium_perenne.gff
mv Kyuss_1697_KYUS_CDS.fa Lolium_perenne.cds
mv Kyuss_1697_KYUS_proteins.fa Lolium_perenne.faa
### Clean-up
rm GCF*.gz GCA*.gz
```

### Install R and Julia
```shell
sudo apt install -y r-base julia
```

### Clone this repository
```shell
git clone https://github.com/jeffersonfparil/Lolium_rigidum_genome_assembly_and_annotation.git
```

### Install OrthoFinder for classifying genes into orthologs, and paralogs, as well as to build a tree for the analysis of gene family evolution
```shell
wget https://github.com/davidemms/OrthoFinder/releases/download/2.5.4/OrthoFinder.tar.gz
tar -xvzf OrthoFinder.tar.gz
cd OrthoFinder/
./orthofinder -h
PATH=${PATH}:${DIR}/OrthoFinder
cd -
rm OrthoFinder.tar.gz
```

### Install HMMER for mapping CDS to PantherHMM gene family models
```shell
sudo apt install -y hmmer
```

### Install CAFE5 to analyse gene family evolution
```shell
wget https://github.com/hahnlab/CAFE5/releases/download/v5.0/CAFE5-5.0.0.tar.gz
tar -xvzf CAFE5-5.0.0.tar.gz
cd CAFE5/
./configure
make
PATH=${PATH}:${DIR}/CAFE5/bin
cafe5 -h
cd -
rm CAFE5-5.0.0.tar.gz
```

### Install MACSE: Multiple Alignment of Coding SEquences Accounting for Frameshifts and Stop Codons
```shell
mkdir MACSE/
cd MACSE/
wget https://bioweb.supagro.inra.fr/macse/releases/macse_v2.06.jar
MACSE=${DIR}/MACSE/macse_v2.06.jar
java -Xmx250G -jar ${MACSE} -help
cd -
```

### Install IQ-TREE for building trees with fossil root dates
```shell
sudo apt install libeigen3-dev
wget https://github.com/Cibiv/IQ-TREE/releases/download/v2.0.7/iqtree-2.0.7-Linux.tar.gz
tar -xvzf iqtree-2.0.7-Linux.tar.gz
PATH=${PATH}:${DIR}/iqtree-2.0.7-Linux/bin
iqtree2 -h
rm iqtree-2.0.7-Linux.tar.gz
```

### Install KaKs_Calculator2.0 to assess signatures of selection
```shell
wget https://github.com/kullrich/kakscalculator2/archive/refs/tags/v2.0.1.tar.gz
tar -xvzf v2.0.1.tar.gz
cd kakscalculator2-2.0.1/src
make
PATH=${PATH}:${DIR}/kakscalculator2-2.0.1/src
KaKs_Calculator -h
cd -
rm v2.0.1.tar.gz
```

### Install R::ape
1. Install BLAS and LAPACK libraries for R::ape package
```shell
sudo apt install -y libblas-dev liblapack-dev
```

2. Install R::ape library
```R
install.packages("ape")
```

### Install Julia packages: Plots, DataFrames, CSV, ProgressMeter, and JLD2
```julia
using Pkg
Pkg.add(["Plots", "DataFrames", "CSV", "ProgressMeter", "JLD2"])
```

### Download PantherHMM library including 15,619 protein family HMMs and their GO terms
```shell
wget http://data.pantherdb.org/ftp/panther_library/current_release/PANTHER17.0_hmmscoring.tgz
tar -xvzf PANTHER17.0_hmmscoring.tgz
mv target/ PantherHMM_17.0/
cd PantherHMM_17.0/
wget http://data.pantherdb.org/ftp/hmm_classifications/current_release/PANTHER17.0_HMM_classifications
### Isolate family codes and names, i.e. exclude subfamily info
grep -v ':SF' PANTHER17.0_HMM_classifications > Panther17.0_HMM_familyIDs.txt
rm PANTHER17.0_hmmscoring.tgz
cd -
```

### Use OrthoFinder to find orthologs and paralogs
```shell
mkdir ORTHOGROUPS/
cp *.faa ORTHOGROUPS/
for f in ORTHOGROUPS/*.faa
do
    # f=$(ls ORTHOGROUPS/*.faa | head -n1)
    fname=$(basename $f)
    species=${fname%.faa*}
    sed -i "s/^>/>$species|/g" $f
done
time \
orthofinder \
    -f ORTHOGROUPS/ \
    -t 31

DIR_ORTHOGROUPS=${DIR}/ORTHOGROUPS/OrthoFinder/Results_May11/
```

### Assign orthogroups into gene families
1. Define the location of the 15,619 protein family HMMs
```shell
DIR_PANTHER=${DIR}/PantherHMM_17.0/famlib/rel/PANTHER17.0_altVersion/hmmscoring/PANTHER17.0/books
GOT_PATHER=${DIR}/PantherHMM_17.0/PANTHER17.0_HMM_classifications
```

2. Iteratively, for each genome's predicted protein sequences run hmmsearch in paralel for each PantherHMM protein family
```shell
### Add Orthogroup name to each protein sequence name and merge so that we can be more efficient with hmmsearch
echo '#!/bin/bash
f=$1
ORTHOGROUP=$(basename ${f} | sed s/.fa$//g)
sed "s/^>/>${ORTHOGROUP}:/g" ${f} > ${ORTHOGROUP}.tmp
' > rename_sequences_with_orthogroup_ID.sh
chmod +x rename_sequences_with_orthogroup_ID.sh

### Split the large list of orthogroup filenames (too long for bash commands)
find ${DIR_ORTHOGROUPS}/Orthogroup_Sequences/ -name '*.fa' > \
    orthogroup_filenames.tmp
split -l 10000 orthogroup_filenames.tmp orthogroup_filenames-SPLIT-

### For each chunk of orthogroup filenames append the orthogroup names into the protein sequence names in parallel
time \
for F in $(ls orthogroup_filenames-SPLIT-*)
do
    echo $F
    parallel \
        ./rename_sequences_with_orthogroup_ID.sh \
        {} \
        ::: $(cat $F)
done

### Merge orthogroup sequences into a single fasta file
MERGED_ORTHOGROUPS=${DIR}/ORTHOGROUPS/orthogroups.faa
touch $MERGED_ORTHOGROUPS
time \
for f in $(ls | grep "^OG" | grep "tmp$")
do
    cat $f >> ${MERGED_ORTHOGROUPS}
    rm $f
done

### Clean-up
rm orthogroup_filenames*
```

3. Find PantherHMM protein families for each orthogroup
```shell
### Prepare parallelisable HMMER search script
echo '#!/bin/bash
PROTFA=$1
DIR_PANTHER=$2
d=$3
HMMODL=${DIR_PANTHER}/${d}/hmmer.hmm
OUTEMP=${PROTFA}-hhmer_gene_family_hits-${d}.tmp
hmmsearch -E 0.0001 --tblout ${OUTEMP} ${HMMODL} ${PROTFA}
sed "/^#/d" ${OUTEMP} | awk @{print $1,$3,$5}@ > ${OUTEMP}.tmp
if [ $(cat ${OUTEMP}.tmp | wc -l) -eq 0 ]
then
    rm ${OUTEMP} ${OUTEMP}.tmp
else
    mv ${OUTEMP}.tmp ${OUTEMP}
fi
' | sed "s/@/'/g" > hmmsearch_for_parallel_execution.sh
chmod +x hmmsearch_for_parallel_execution.sh

### Find PantherHMM protein families for each orthogroup
time \
parallel \
./hmmsearch_for_parallel_execution.sh \
    ${MERGED_ORTHOGROUPS} \
    ${DIR_PANTHER} \
    {} \
    ::: $(ls $DIR_PANTHER)

### Concatenate hmmsearch output for each protein family into a single output file
PANTHER_ORTHOGROUPS=$(echo ${MERGED_ORTHOGROUPS} | sed s/.faa$//g).pthr
cat ${MERGED_ORTHOGROUPS}-hhmer_gene_family_hits-* > ${PANTHER_ORTHOGROUPS}
rm ${MERGED_ORTHOGROUPS}-hhmer_gene_family_hits-*
```

4. Find the best fitting gene family to each unique sequence per orthogroup. This means that each orthogroup can have multiple gene families. Next, add family name and GO terms to each gene family.
```shell
grep "^>" ${MERGED_ORTHOGROUPS} | cut -d':' -f1 | sed 's/>//g' | sort | uniq > all_orthogroups.tmp
time \
julia orthogroup_classification_gene_family_GO_terms.jl \
        ORTHOGROUPS/orthogroups.pthr \
        PantherHMM_17.0/Panther17.0_HMM_familyIDs.txt \
        all_orthogroups.tmp \
        ${DIR_ORTHOGROUPS}/Orthogroups/Orthogroups.GeneCount.tsv \
        ${DIR_ORTHOGROUPS}/Orthogroups/Orthogroups_UnassignedGenes.tsv \
        ORTHOGROUPS/orthogroups_gene_counts_families_go.out
```

### Preliminary assessment of the distribution of the genes, orthogroups and gene family classifications.
```shell
time \
julia count_genes_per_ortholog_paralog_classes.jl \
        ORTHOGROUPS/orthogroups_gene_counts_families_go.out \
        ORTHOGROUPS/orthogroups_summarised_gene_counts.csv
```

### Infer gene family expansion and contraction in each species using CAFE5 (at alpha=1%)
```shell
ORTHOUT=${DIR}/ORTHOGROUPS/orthogroups_gene_counts_families_go.out
rev ${ORTHOUT} | cut -f5- | rev > col2_to_coln.tmp
awk -F'\t' '{print $(NF-1)}' ${ORTHOUT} > col1.tmp
paste -d'\t' col1.tmp col2_to_coln.tmp > counts.tmp
TREE=${DIR_ORTHOGROUPS}/Species_Tree/SpeciesTree_rooted.txt

time \
cafe5 \
    --infile counts.tmp \
    --tree ${TREE} \
    --n_gamma_cats 1000 \
    --cores 31 \
    --pvalue 0.01 \
    --output_prefix CAFE_Gamma1000_results
### Output using n_gamma_cats=1,000
echo -e "Species\tExpansion\tContraction" > CONTRACTION_EXPANSION.txt
grep -v "^#" ${DIR}/CAFE_Gamma1000_results/Gamma_clade_results.txt | \
    grep -v "^<" | \
    sed 's/<..>//g' | \
    sed 's/<.>//g' >> CONTRACTION_EXPANSION.txt
```

### GO term enrichment analysis of expanded gene families
1. Extract gene names of the expanded and contracted gene families in *Lolium rigidum*
```shell
ORTHOUT=${DIR}/ORTHOGROUPS/orthogroups_gene_counts_families_go.out
n=$(head -n1 ${DIR}/CAFE_Gamma1000_results/Gamma_change.tab | sed -z "s/\t/\n/g" | grep -n "Lolium_rigidum" | cut -d":" -f1)
cut -f1,${n} ${DIR}/CAFE_Gamma1000_results/Gamma_change.tab | grep -v "+0" | grep "+" | cut -f1 > expanded_orthogroups_for_grep.tmp
cut -f1,${n} ${DIR}/CAFE_Gamma1000_results/Gamma_change.tab | grep -v "+" | cut -f1 > contracted_orthogroups_for_grep.tmp

grep -wf expanded_orthogroups_for_grep.tmp $ORTHOUT | cut -f$(head -n1 $ORTHOUT | awk '{printf NF-2}') | grep -v "^$" > expanded_orthogroups.pthr.tmp
grep -wf contracted_orthogroups_for_grep.tmp $ORTHOUT | cut -f$(head -n1 $ORTHOUT | awk '{printf NF-2}') | grep -v "^$" > contracted_orthogroups.pthr.tmp
wget http://data.pantherdb.org/ftp/sequence_classifications/current_release/PANTHER_Sequence_Classification_files/PTHR17.0_rice
grep -wf expanded_orthogroups.pthr.tmp PTHR17.0_rice | cut -f3 > expanded_orthogroups.forgo
grep -wf contracted_orthogroups.pthr.tmp PTHR17.0_rice | cut -f3 > contracted_orthogroups.forgo
```

2. Run gene ontology enrichment analysis
- Open browser to: `http://geneontology.org/`.
- Paste the gene names in `expanded_orthogroups.forgo` into the text box.
- Set to "Biological process"
- Set the reference gene list to "Oryza sativa"
- Launch the analysis
- Set the correction to "Bonferroni correction" and relaunch analysis
- Sort by decreasing "Fold Enrichment"
- Export as table

3. `expanded_orthogroups.goea`
```shell
Analysis Type:	PANTHER Overrepresentation Test (Released 20220202)
Annotation Version and Release Date:	GO Ontology database DOI:  10.5281/zenodo.6399963 Released 2022-03-22
Analyzed List:	upload_1 (Oryza sativa)
Reference List:	Oryza sativa (all genes in database)
Test Type:	FISHER
Correction:	BONFERRONI
Bonferroni count:	2108
GO biological process complete	Oryza sativa - REFLIST (43658)	upload_1 (15118)	upload_1 (expected)	upload_1 (over/under)	upload_1 (fold Enrichment)	upload_1 (P-value)
xylan biosynthetic process (GO:0045492)	39	39	13.51	+	2.89	1.20E-02
xenobiotic detoxification by transmembrane export across the plasma membrane (GO:1990961)	39	38	13.51	+	2.81	2.08E-02
xenobiotic export from cell (GO:0046618)	39	38	13.51	+	2.81	2.08E-02
lignin metabolic process (GO:0009808)	55	53	19.05	+	2.78	4.28E-04
xyloglucan metabolic process (GO:0010411)	53	51	18.35	+	2.78	8.30E-04
cell wall polysaccharide biosynthetic process (GO:0070592)	74	71	25.62	+	2.77	5.36E-06
hydrogen peroxide catabolic process (GO:0042744)	144	137	49.86	+	2.75	2.55E-13
hydrogen peroxide metabolic process (GO:0042743)	144	137	49.86	+	2.75	2.55E-13
cellulose metabolic process (GO:0030243)	73	69	25.28	+	2.73	1.25E-05
xenobiotic transport (GO:0042908)	50	47	17.31	+	2.71	3.90E-03
cellular component macromolecule biosynthetic process (GO:0070589)	77	71	26.66	+	2.66	1.51E-05
cell wall macromolecule biosynthetic process (GO:0044038)	77	71	26.66	+	2.66	1.51E-05
reactive oxygen species metabolic process (GO:0072593)	159	146	55.06	+	2.65	1.74E-13
cellulose biosynthetic process (GO:0030244)	46	42	15.93	+	2.64	2.58E-02
cellular response to ethylene stimulus (GO:0071369)	64	58	22.16	+	2.62	4.55E-04
ethylene-activated signaling pathway (GO:0009873)	63	57	21.82	+	2.61	6.24E-04
cellular response to abscisic acid stimulus (GO:0071215)	83	75	28.74	+	2.61	9.17E-06
cellular response to alcohol (GO:0097306)	83	75	28.74	+	2.61	9.17E-06
plant-type cell wall organization (GO:0009664)	62	56	21.47	+	2.61	8.60E-04
oligopeptide transport (GO:0006857)	51	46	17.66	+	2.60	8.43E-03
peptide transport (GO:0015833)	51	46	17.66	+	2.60	8.43E-03
abscisic acid-activated signaling pathway (GO:0009738)	79	71	27.36	+	2.60	3.47E-05
diterpenoid metabolic process (GO:0016101)	66	59	22.85	+	2.58	6.00E-04
cell wall biogenesis (GO:0042546)	168	150	58.18	+	2.58	4.17E-13
plant-type secondary cell wall biogenesis (GO:0009834)	56	50	19.39	+	2.58	4.11E-03
protein refolding (GO:0042026)	46	41	15.93	+	2.57	4.20E-02
flavonoid biosynthetic process (GO:0009813)	46	41	15.93	+	2.57	4.20E-02
flavonoid metabolic process (GO:0009812)	46	41	15.93	+	2.57	4.20E-02
phosphorelay signal transduction system (GO:0000160)	117	104	40.52	+	2.57	1.66E-08
plant-type cell wall organization or biogenesis (GO:0071669)	156	138	54.02	+	2.55	8.89E-12
hemicellulose metabolic process (GO:0010410)	112	99	38.78	+	2.55	7.98E-08
plant-type cell wall biogenesis (GO:0009832)	108	95	37.40	+	2.54	1.97E-07
auxin-activated signaling pathway (GO:0009734)	124	109	42.94	+	2.54	1.05E-08
cell wall polysaccharide metabolic process (GO:0010383)	122	107	42.25	+	2.53	2.05E-08
cell wall macromolecule metabolic process (GO:0044036)	146	128	50.56	+	2.53	1.41E-10
polysaccharide biosynthetic process (GO:0000271)	177	155	61.29	+	2.53	3.90E-13
diterpenoid biosynthetic process (GO:0016102)	56	49	19.39	+	2.53	1.05E-02
cytokinin-activated signaling pathway (GO:0009736)	62	54	21.47	+	2.52	3.73E-03
phenylpropanoid biosynthetic process (GO:0009699)	54	47	18.70	+	2.51	1.32E-02
beta-glucan metabolic process (GO:0051273)	84	73	29.09	+	2.51	4.57E-05
cellular polysaccharide biosynthetic process (GO:0033692)	158	137	54.71	+	2.50	4.70E-11
phenylpropanoid metabolic process (GO:0009698)	112	97	38.78	+	2.50	2.23E-07
cellular response to auxin stimulus (GO:0071365)	126	109	43.63	+	2.50	2.28E-08
amide transport (GO:0042886)	66	57	22.85	+	2.49	1.67E-03
cell wall organization (GO:0071555)	314	271	108.73	+	2.49	9.77E-24
pectin catabolic process (GO:0045490)	65	56	22.51	+	2.49	2.28E-03
secondary metabolic process (GO:0019748)	176	151	60.95	+	2.48	3.51E-12
cellular response to cytokinin stimulus (GO:0071368)	63	54	21.82	+	2.48	4.29E-03
export across plasma membrane (GO:0140115)	62	53	21.47	+	2.47	5.94E-03
pectin metabolic process (GO:0045488)	89	76	30.82	+	2.47	3.72E-05
galacturonan metabolic process (GO:0010393)	89	76	30.82	+	2.47	3.72E-05
response to cytokinin (GO:0009735)	67	57	23.20	+	2.46	2.94E-03
SCF-dependent proteasomal ubiquitin-dependent protein catabolic process (GO:0031146)	91	77	31.51	+	2.44	4.81E-05
secondary metabolite biosynthetic process (GO:0044550)	97	82	33.59	+	2.44	1.73E-05
peptidyl-serine phosphorylation (GO:0018105)	84	71	29.09	+	2.44	1.87E-04
response to ethylene (GO:0009723)	84	71	29.09	+	2.44	1.87E-04
regulation of cyclin-dependent protein serine/threonine kinase activity (GO:0000079)	70	59	24.24	+	2.43	1.85E-03
protein autophosphorylation (GO:0046777)	70	59	24.24	+	2.43	1.85E-03
external encapsulating structure organization (GO:0045229)	340	286	117.74	+	2.43	5.61E-24
response to oomycetes (GO:0002239)	79	66	27.36	+	2.41	5.96E-04
defense response to oomycetes (GO:0002229)	79	66	27.36	+	2.41	5.96E-04
peptidyl-serine modification (GO:0018209)	85	71	29.43	+	2.41	2.13E-04
glucan biosynthetic process (GO:0009250)	115	96	39.82	+	2.41	1.27E-06
cell wall organization or biogenesis (GO:0071554)	415	346	143.71	+	2.41	6.25E-29
regulation of cyclin-dependent protein kinase activity (GO:1904029)	72	60	24.93	+	2.41	2.30E-03
transmembrane receptor protein serine/threonine kinase signaling pathway (GO:0007178)	77	64	26.66	+	2.40	1.13E-03
enzyme-linked receptor protein signaling pathway (GO:0007167)	77	64	26.66	+	2.40	1.13E-03
response to toxic substance (GO:0009636)	295	244	102.15	+	2.39	1.74E-19
xylan metabolic process (GO:0045491)	63	52	21.82	+	2.38	1.73E-02
detoxification (GO:0098754)	276	226	95.57	+	2.36	1.90E-17
protein targeting to membrane (GO:0006612)	77	63	26.66	+	2.36	1.80E-03
hormone-mediated signaling pathway (GO:0009755)	410	335	141.98	+	2.36	7.21E-27
cellular response to hormone stimulus (GO:0032870)	418	341	144.75	+	2.36	3.16E-27
response to water (GO:0009415)	70	57	24.24	+	2.35	6.93E-03
mitotic cell cycle phase transition (GO:0044772)	59	48	20.43	+	2.35	4.04E-02
cell cycle phase transition (GO:0044770)	59	48	20.43	+	2.35	4.04E-02
cellular oxidant detoxification (GO:0098869)	218	177	75.49	+	2.34	4.19E-13
glucan metabolic process (GO:0044042)	202	164	69.95	+	2.34	7.10E-12
cellular glucan metabolic process (GO:0006073)	196	159	67.87	+	2.34	1.96E-11
response to water deprivation (GO:0009414)	68	55	23.55	+	2.34	1.31E-02
cell surface receptor signaling pathway (GO:0007166)	261	211	90.38	+	2.33	8.91E-16
cellular polysaccharide metabolic process (GO:0044264)	259	209	89.69	+	2.33	1.67E-15
hormone metabolic process (GO:0042445)	93	75	32.20	+	2.33	2.43E-04
regulation of protein serine/threonine kinase activity (GO:0071900)	81	65	28.05	+	2.32	1.91E-03
cellular response to endogenous stimulus (GO:0071495)	425	341	147.17	+	2.32	2.71E-26
cellular response to toxic substance (GO:0097237)	234	187	81.03	+	2.31	2.04E-13
cellular detoxification (GO:1990748)	234	187	81.03	+	2.31	2.04E-13
carbohydrate derivative transport (GO:1901264)	79	63	27.36	+	2.30	3.55E-03
potassium ion transmembrane transport (GO:0071805)	64	51	22.16	+	2.30	3.14E-02
cellular carbohydrate biosynthetic process (GO:0034637)	205	163	70.99	+	2.30	2.41E-11
polysaccharide metabolic process (GO:0005976)	395	314	136.78	+	2.30	1.12E-23
cellular response to organic substance (GO:0071310)	483	381	167.25	+	2.28	1.35E-28
potassium ion transport (GO:0006813)	70	55	24.24	+	2.27	2.59E-02
protein deubiquitination (GO:0016579)	97	76	33.59	+	2.26	4.09E-04
response to abscisic acid (GO:0009737)	151	118	52.29	+	2.26	2.25E-07
response to auxin (GO:0009733)	215	168	74.45	+	2.26	3.16E-11
cellular response to chemical stimulus (GO:0070887)	758	592	262.48	+	2.26	6.57E-45
steroid metabolic process (GO:0008202)	102	79	35.32	+	2.24	4.43E-04
chloroplast organization (GO:0009658)	133	102	46.06	+	2.21	1.06E-05
response to cold (GO:0009409)	90	69	31.17	+	2.21	3.41E-03
regulation of protein catabolic process (GO:0042176)	77	59	26.66	+	2.21	1.61E-02
chaperone-mediated protein folding (GO:0061077)	77	59	26.66	+	2.21	1.61E-02
positive regulation of cellular catabolic process (GO:0031331)	111	85	38.44	+	2.21	2.21E-04
positive regulation of catabolic process (GO:0009896)	128	98	44.32	+	2.21	2.14E-05
cellular response to lipid (GO:0071396)	161	123	55.75	+	2.21	2.82E-07
response to acid chemical (GO:0001101)	76	58	26.32	+	2.20	2.17E-02
carbohydrate transmembrane transport (GO:0034219)	88	67	30.47	+	2.20	4.39E-03
response to alcohol (GO:0097305)	155	118	53.67	+	2.20	7.73E-07
regulation of hormone levels (GO:0010817)	134	102	46.40	+	2.20	1.18E-05
endocytosis (GO:0006897)	71	54	24.59	+	2.20	4.52E-02
amino acid transport (GO:0006865)	99	75	34.28	+	2.19	1.68E-03
response to hormone (GO:0009725)	613	463	212.27	+	2.18	3.96E-32
organic acid transport (GO:0015849)	138	104	47.79	+	2.18	1.19E-05
phyllome development (GO:0048827)	100	75	34.63	+	2.17	1.86E-03
response to oxidative stress (GO:0006979)	275	206	95.23	+	2.16	6.25E-13
regulation of protein kinase activity (GO:0045859)	95	71	32.90	+	2.16	3.84E-03
sterol metabolic process (GO:0016125)	87	65	30.13	+	2.16	1.28E-02
response to endogenous stimulus (GO:0009719)	620	463	214.70	+	2.16	2.61E-31
glutathione metabolic process (GO:0006749)	94	70	32.55	+	2.15	5.15E-03
carbohydrate transport (GO:0008643)	113	84	39.13	+	2.15	6.03E-04
cellular response to oxygen-containing compound (GO:1901701)	206	153	71.33	+	2.14	8.30E-09
regulation of protein phosphorylation (GO:0001932)	101	75	34.97	+	2.14	2.16E-03
signal transduction (GO:0007165)	1097	807	379.87	+	2.12	8.74E-55
polysaccharide catabolic process (GO:0000272)	188	138	65.10	+	2.12	1.65E-07
peptidyl-tyrosine modification (GO:0018212)	94	69	32.55	+	2.12	7.92E-03
peptidyl-tyrosine phosphorylation (GO:0018108)	94	69	32.55	+	2.12	7.92E-03
regulation of kinase activity (GO:0043549)	97	71	33.59	+	2.11	6.86E-03
response to chemical (GO:0042221)	1189	868	411.73	+	2.11	3.66E-58
regulation of phosphorylation (GO:0042325)	103	75	35.67	+	2.10	3.80E-03
protein dephosphorylation (GO:0006470)	187	136	64.75	+	2.10	3.26E-07
export from cell (GO:0140352)	154	112	53.33	+	2.10	1.17E-05
carbohydrate biosynthetic process (GO:0016051)	278	202	96.27	+	2.10	1.23E-11
signaling (GO:0023052)	1112	807	385.07	+	2.10	5.06E-53
response to organic substance (GO:0010033)	757	549	262.14	+	2.09	3.94E-35
response to lipid (GO:0033993)	250	180	86.57	+	2.08	6.58E-10
response to bacterium (GO:0009617)	159	114	55.06	+	2.07	1.76E-05
glycosylation (GO:0070085)	208	149	72.03	+	2.07	9.95E-08
response to temperature stimulus (GO:0009266)	194	138	67.18	+	2.05	6.16E-07
establishment of protein localization to membrane (GO:0090150)	117	83	40.52	+	2.05	2.76E-03
protein modification by small protein removal (GO:0070646)	120	85	41.55	+	2.05	1.75E-03
regulation of phosphate metabolic process (GO:0019220)	154	109	53.33	+	2.04	5.33E-05
regulation of phosphorus metabolic process (GO:0051174)	154	109	53.33	+	2.04	5.33E-05
positive regulation of transcription by RNA polymerase II (GO:0045944)	155	109	53.67	+	2.03	8.12E-05
cell division (GO:0051301)	196	137	67.87	+	2.02	1.55E-06
response to salt stress (GO:0009651)	96	67	33.24	+	2.02	4.22E-02
plastid organization (GO:0009657)	165	115	57.14	+	2.01	4.28E-05
regulation of catabolic process (GO:0009894)	155	108	53.67	+	2.01	1.19E-04
intracellular signal transduction (GO:0035556)	352	245	121.89	+	2.01	5.70E-13
response to osmotic stress (GO:0006970)	115	80	39.82	+	2.01	7.12E-03
cell communication (GO:0007154)	1274	885	441.16	+	2.01	1.31E-52
regulation of transferase activity (GO:0051338)	111	77	38.44	+	2.00	1.09E-02
carbohydrate catabolic process (GO:0016052)	293	202	101.46	+	1.99	4.05E-10
response to oxygen-containing compound (GO:1901700)	428	294	148.21	+	1.98	2.13E-15
fatty acid biosynthetic process (GO:0006633)	127	87	43.98	+	1.98	4.27E-03
regulation of protein modification process (GO:0031399)	168	115	58.18	+	1.98	1.10E-04
defense response to bacterium (GO:0042742)	136	93	47.09	+	1.97	2.02E-03
root system development (GO:0022622)	120	82	41.55	+	1.97	7.59E-03
root development (GO:0048364)	120	82	41.55	+	1.97	7.59E-03
peptidyl-amino acid modification (GO:0018193)	469	319	162.41	+	1.96	2.26E-16
proteasome-mediated ubiquitin-dependent protein catabolic process (GO:0043161)	246	167	85.19	+	1.96	1.16E-07
proteasomal protein catabolic process (GO:0010498)	275	186	95.23	+	1.95	1.41E-08
positive regulation of RNA metabolic process (GO:0051254)	386	261	133.67	+	1.95	8.75E-13
fatty acid metabolic process (GO:0006631)	194	131	67.18	+	1.95	2.22E-05
response to abiotic stimulus (GO:0009628)	597	403	206.73	+	1.95	7.88E-21
regulation of cellular catabolic process (GO:0031329)	132	89	45.71	+	1.95	4.47E-03
carbohydrate metabolic process (GO:0005975)	1039	697	359.79	+	1.94	3.86E-37
positive regulation of cellular biosynthetic process (GO:0031328)	348	233	120.51	+	1.93	5.86E-11
organelle localization (GO:0051640)	124	83	42.94	+	1.93	1.40E-02
positive regulation of biosynthetic process (GO:0009891)	350	234	121.20	+	1.93	6.62E-11
positive regulation of macromolecule biosynthetic process (GO:0010557)	340	227	117.74	+	1.93	1.82E-10
cellular carbohydrate metabolic process (GO:0044262)	427	285	147.86	+	1.93	1.27E-13
positive regulation of protein metabolic process (GO:0051247)	132	88	45.71	+	1.93	8.72E-03
response to inorganic substance (GO:0010035)	174	116	60.25	+	1.93	2.38E-04
positive regulation of nucleobase-containing compound metabolic process (GO:0045935)	399	266	138.17	+	1.93	1.53E-12
transmembrane transport (GO:0055085)	1365	908	472.68	+	1.92	2.34E-48
protein folding (GO:0006457)	229	152	79.30	+	1.92	2.97E-06
positive regulation of nitrogen compound metabolic process (GO:0051173)	522	346	180.76	+	1.91	1.06E-16
regulation of RNA biosynthetic process (GO:2001141)	2028	1344	702.26	+	1.91	4.14E-73
regulation of transcription, DNA-templated (GO:0006355)	2028	1344	702.26	+	1.91	4.14E-73
regulation of nucleic acid-templated transcription (GO:1903506)	2028	1344	702.26	+	1.91	4.14E-73
positive regulation of RNA biosynthetic process (GO:1902680)	317	210	109.77	+	1.91	2.25E-09
positive regulation of transcription, DNA-templated (GO:0045893)	317	210	109.77	+	1.91	2.25E-09
positive regulation of nucleic acid-templated transcription (GO:1903508)	317	210	109.77	+	1.91	2.25E-09
shoot system development (GO:0048367)	265	175	91.76	+	1.91	2.16E-07
cellular response to stimulus (GO:0051716)	1912	1261	662.09	+	1.90	2.46E-67
plant organ development (GO:0099402)	255	168	88.30	+	1.90	5.90E-07
regulation of RNA metabolic process (GO:0051252)	2135	1406	739.31	+	1.90	1.58E-75
response to stimulus (GO:0050896)	3698	2432	1280.55	+	1.90	2.13E-138
protein phosphorylation (GO:0006468)	1515	996	524.62	+	1.90	1.17E-51
defense response (GO:0006952)	1040	680	360.13	+	1.89	1.02E-33
response to biotic stimulus (GO:0009607)	346	226	119.81	+	1.89	1.03E-09
protein glycosylation (GO:0006486)	147	96	50.90	+	1.89	6.17E-03
macromolecule glycosylation (GO:0043413)	147	96	50.90	+	1.89	6.17E-03
positive regulation of cellular metabolic process (GO:0031325)	515	336	178.34	+	1.88	1.89E-15
anion transport (GO:0006820)	204	133	70.64	+	1.88	6.27E-05
metal ion transport (GO:0030001)	224	146	77.57	+	1.88	1.58E-05
regulation of macromolecule biosynthetic process (GO:0010556)	2237	1458	774.63	+	1.88	1.90E-76
positive regulation of macromolecule metabolic process (GO:0010604)	528	344	182.84	+	1.88	8.39E-16
regulation of nucleobase-containing compound metabolic process (GO:0019219)	2192	1427	759.05	+	1.88	2.14E-74
regulation of cellular biosynthetic process (GO:0031326)	2258	1467	781.91	+	1.88	3.47E-76
response to external biotic stimulus (GO:0043207)	319	207	110.46	+	1.87	1.33E-08
response to other organism (GO:0051707)	319	207	110.46	+	1.87	1.33E-08
glycoprotein biosynthetic process (GO:0009101)	148	96	51.25	+	1.87	6.61E-03
regulation of biosynthetic process (GO:0009889)	2265	1469	784.33	+	1.87	5.36E-76
positive regulation of metabolic process (GO:0009893)	544	352	188.38	+	1.87	8.38E-16
reproductive shoot system development (GO:0090567)	147	95	50.90	+	1.87	8.59E-03
terpenoid metabolic process (GO:0006721)	137	88	47.44	+	1.85	2.37E-02
positive regulation of cellular process (GO:0048522)	578	371	200.15	+	1.85	1.95E-16
regulation of cellular process (GO:0050794)	3867	2481	1339.07	+	1.85	1.92E-132
monocarboxylic acid metabolic process (GO:0032787)	391	250	135.40	+	1.85	3.39E-10
sulfur compound metabolic process (GO:0006790)	305	195	105.62	+	1.85	1.41E-07
regulation of nitrogen compound metabolic process (GO:0051171)	2644	1690	915.57	+	1.85	3.86E-85
flower development (GO:0009908)	141	90	48.83	+	1.84	2.19E-02
cellular macromolecule biosynthetic process (GO:0034645)	951	607	329.31	+	1.84	7.78E-28
regulation of primary metabolic process (GO:0080090)	2669	1701	924.23	+	1.84	5.87E-85
reproductive structure development (GO:0048608)	291	185	100.77	+	1.84	5.78E-07
response to external stimulus (GO:0009605)	455	289	157.56	+	1.83	7.22E-12
protein modification process (GO:0036211)	3176	2013	1099.79	+	1.83	5.92E-101
cellular protein modification process (GO:0006464)	3176	2013	1099.79	+	1.83	5.92E-101
reproductive system development (GO:0061458)	292	185	101.11	+	1.83	8.24E-07
chemical homeostasis (GO:0048878)	270	171	93.50	+	1.83	3.92E-06
macromolecule modification (GO:0043412)	3495	2212	1210.26	+	1.83	8.88E-112
monocarboxylic acid biosynthetic process (GO:0072330)	174	110	60.25	+	1.83	3.32E-03
regulation of cellular metabolic process (GO:0031323)	2690	1699	931.50	+	1.82	1.08E-82
regulation of gene expression (GO:0010468)	2484	1565	860.17	+	1.82	6.16E-75
regulation of biological process (GO:0050789)	4286	2700	1484.17	+	1.82	1.44E-138
biological regulation (GO:0065007)	4867	3065	1685.36	+	1.82	1.08E-160
glycoprotein metabolic process (GO:0009100)	159	100	55.06	+	1.82	1.28E-02
negative regulation of catalytic activity (GO:0043086)	229	144	79.30	+	1.82	1.06E-04
negative regulation of molecular function (GO:0044092)	229	144	79.30	+	1.82	1.06E-04
cellular catabolic process (GO:0044248)	1349	846	467.14	+	1.81	9.52E-38
homeostatic process (GO:0042592)	303	190	104.92	+	1.81	8.25E-07
phosphorylation (GO:0016310)	1841	1154	637.51	+	1.81	7.69E-53
macromolecule catabolic process (GO:0009057)	952	596	329.66	+	1.81	1.23E-25
regulation of transcription by RNA polymerase II (GO:0006357)	516	323	178.68	+	1.81	6.17E-13
RNA modification (GO:0009451)	302	189	104.58	+	1.81	1.06E-06
cellular protein metabolic process (GO:0044267)	4088	2555	1415.60	+	1.80	4.94E-127
regulation of macromolecule metabolic process (GO:0060255)	2862	1785	991.06	+	1.80	3.57E-84
dephosphorylation (GO:0016311)	374	232	129.51	+	1.79	1.89E-08
regulation of metabolic process (GO:0019222)	2943	1823	1019.11	+	1.79	1.98E-84
protein metabolic process (GO:0019538)	4589	2839	1589.09	+	1.79	2.38E-139
defense response to other organism (GO:0098542)	275	170	95.23	+	1.79	1.33E-05
system development (GO:0048731)	492	304	170.37	+	1.78	1.62E-11
catabolic process (GO:0009056)	1640	1013	567.90	+	1.78	7.28E-44
regulation of catalytic activity (GO:0050790)	518	319	179.37	+	1.78	5.66E-12
developmental process involved in reproduction (GO:0003006)	320	197	110.81	+	1.78	1.13E-06
regulation of molecular function (GO:0065009)	527	323	182.49	+	1.77	5.00E-12
positive regulation of biological process (GO:0048518)	687	421	237.90	+	1.77	2.96E-16
protein catabolic process (GO:0030163)	622	379	215.39	+	1.76	3.44E-14
response to stress (GO:0006950)	2142	1303	741.74	+	1.76	4.97E-55
proteolysis involved in cellular protein catabolic process (GO:0051603)	612	372	211.92	+	1.76	9.09E-14
cellular protein catabolic process (GO:0044257)	613	372	212.27	+	1.75	1.25E-13
protein modification by small protein conjugation or removal (GO:0070647)	915	555	316.85	+	1.75	2.36E-21
establishment of localization (GO:0051234)	2519	1527	872.29	+	1.75	9.06E-65
regulation of protein metabolic process (GO:0051246)	527	319	182.49	+	1.75	2.87E-11
cellular macromolecule catabolic process (GO:0044265)	793	480	274.60	+	1.75	4.50E-18
transport (GO:0006810)	2485	1503	860.51	+	1.75	3.22E-63
organic substance transport (GO:0071702)	1246	752	431.47	+	1.74	2.37E-29
peptide metabolic process (GO:0006518)	715	430	247.59	+	1.74	1.19E-15
localization (GO:0051179)	2618	1570	906.57	+	1.73	1.36E-64
multicellular organism development (GO:0007275)	651	390	225.43	+	1.73	1.01E-13
cellular macromolecule metabolic process (GO:0044260)	5139	3067	1779.55	+	1.72	1.40E-136
response to light stimulus (GO:0009416)	213	127	73.76	+	1.72	5.55E-03
organic substance catabolic process (GO:1901575)	1434	855	496.57	+	1.72	1.66E-32
nitrogen compound transport (GO:0071705)	1005	599	348.01	+	1.72	6.56E-22
macromolecule biosynthetic process (GO:0009059)	1321	787	457.44	+	1.72	1.32E-29
biological_process (GO:0008150)	17397	10359	6024.28	+	1.72	0.00E00
biological process involved in interspecies interaction between organisms (GO:0044419)	353	210	122.24	+	1.72	3.12E-06
regulation of biological quality (GO:0065008)	585	348	202.58	+	1.72	1.04E-11
regulation of cellular protein metabolic process (GO:0032268)	501	298	173.49	+	1.72	8.46E-10
organonitrogen compound metabolic process (GO:1901564)	5580	3319	1932.26	+	1.72	1.38E-148
organic hydroxy compound metabolic process (GO:1901615)	285	169	98.69	+	1.71	1.26E-04
methylation (GO:0032259)	334	198	115.66	+	1.71	9.54E-06
phosphorus metabolic process (GO:0006793)	2618	1551	906.57	+	1.71	3.21E-61
regulation of cell cycle (GO:0051726)	211	125	73.07	+	1.71	6.83E-03
protein modification by small protein conjugation (GO:0032446)	794	470	274.95	+	1.71	2.59E-16
modification-dependent macromolecule catabolic process (GO:0043632)	547	323	189.42	+	1.71	1.73E-10
phosphate-containing compound metabolic process (GO:0006796)	2581	1523	893.76	+	1.70	4.07E-59
proteolysis (GO:0006508)	1219	719	422.12	+	1.70	6.83E-26
ion transmembrane transport (GO:0034220)	529	312	183.18	+	1.70	5.43E-10
modification-dependent protein catabolic process (GO:0019941)	537	316	185.95	+	1.70	4.46E-10
protein ubiquitination (GO:0016567)	755	444	261.44	+	1.70	7.03E-15
metabolic process (GO:0008152)	10316	6038	3572.25	+	1.69	3.84E-308
organonitrogen compound catabolic process (GO:1901565)	786	460	272.18	+	1.69	2.93E-15
isoprenoid metabolic process (GO:0006720)	188	110	65.10	+	1.69	4.07E-02
cellular amide metabolic process (GO:0043603)	826	483	286.03	+	1.69	4.46E-16
macromolecule metabolic process (GO:0043170)	6945	4056	2404.93	+	1.69	2.76E-179
cellular process (GO:0009987)	12820	7486	4439.34	+	1.69	0.00E00
lipid metabolic process (GO:0006629)	961	559	332.78	+	1.68	1.29E-18
RNA splicing, via transesterification reactions (GO:0000375)	217	126	75.14	+	1.68	1.53E-02
organic substance metabolic process (GO:0071704)	9447	5475	3271.33	+	1.67	7.39E-260
cellular metabolic process (GO:0044237)	8971	5195	3106.50	+	1.67	4.39E-241
ion transport (GO:0006811)	631	365	218.50	+	1.67	4.28E-11
post-embryonic development (GO:0009791)	381	220	131.93	+	1.67	7.81E-06
response to radiation (GO:0009314)	220	127	76.18	+	1.67	1.73E-02
translation (GO:0006412)	584	337	202.23	+	1.67	4.72E-10
RNA splicing, via transesterification reactions with bulged adenosine as nucleophile (GO:0000377)	215	124	74.45	+	1.67	2.42E-02
primary metabolic process (GO:0044238)	8914	5129	3086.76	+	1.66	7.21E-232
ubiquitin-dependent protein catabolic process (GO:0006511)	513	295	177.64	+	1.66	2.13E-08
mitotic cell cycle (GO:0000278)	222	127	76.87	+	1.65	2.45E-02
peptide biosynthetic process (GO:0043043)	596	340	206.38	+	1.65	1.26E-09
post-transcriptional regulation of gene expression (GO:0010608)	270	154	93.50	+	1.65	3.02E-03
amide biosynthetic process (GO:0043604)	662	376	229.24	+	1.64	1.09E-10
anatomical structure development (GO:0048856)	768	436	265.94	+	1.64	1.11E-12
aromatic compound catabolic process (GO:0019439)	242	137	83.80	+	1.63	1.43E-02
vesicle-mediated transport (GO:0016192)	493	278	170.72	+	1.63	3.82E-07
biosynthetic process (GO:0009058)	2747	1549	951.24	+	1.63	2.74E-51
regulation of cellular macromolecule biosynthetic process (GO:2000112)	229	129	79.30	+	1.63	3.22E-02
protein transport (GO:0015031)	689	388	238.59	+	1.63	1.05E-10
nitrogen compound metabolic process (GO:0006807)	7413	4172	2566.99	+	1.63	1.36E-163
lipid biosynthetic process (GO:0008610)	548	308	189.76	+	1.62	4.82E-08
inorganic ion transmembrane transport (GO:0098660)	347	195	120.16	+	1.62	2.37E-04
organic substance biosynthetic process (GO:1901576)	2612	1463	904.49	+	1.62	7.02E-47
establishment of protein localization (GO:0045184)	697	390	241.36	+	1.62	1.80E-10
developmental process (GO:0032502)	899	502	311.31	+	1.61	5.84E-14
organic acid biosynthetic process (GO:0016053)	429	239	148.56	+	1.61	1.67E-05
organic cyclic compound catabolic process (GO:1901361)	247	137	85.53	+	1.60	3.66E-02
cellular biosynthetic process (GO:0044249)	2505	1389	867.44	+	1.60	1.42E-42
multicellular organismal process (GO:0032501)	830	459	287.41	+	1.60	4.59E-12
cellular component biogenesis (GO:0044085)	1116	616	386.45	+	1.59	7.93E-17
inorganic cation transmembrane transport (GO:0098662)	299	165	103.54	+	1.59	5.21E-03
regulation of response to stimulus (GO:0048583)	306	168	105.96	+	1.59	5.17E-03
ribonucleoprotein complex biogenesis (GO:0022613)	450	247	155.83	+	1.59	2.63E-05
RNA splicing (GO:0008380)	261	143	90.38	+	1.58	3.56E-02
cellular component organization or biogenesis (GO:0071840)	2322	1265	804.07	+	1.57	1.06E-35
carboxylic acid biosynthetic process (GO:0046394)	405	219	140.24	+	1.56	3.75E-04
protein localization (GO:0008104)	767	412	265.60	+	1.55	2.75E-09
cellular component organization (GO:0016043)	1998	1073	691.87	+	1.55	6.41E-28
small molecule biosynthetic process (GO:0044283)	578	310	200.15	+	1.55	2.26E-06
ribosome biogenesis (GO:0042254)	358	192	123.97	+	1.55	3.82E-03
cation transmembrane transport (GO:0098655)	362	194	125.35	+	1.55	3.40E-03
organic acid metabolic process (GO:0006082)	883	473	305.77	+	1.55	9.69E-11
cellular lipid metabolic process (GO:0044255)	684	366	236.86	+	1.55	7.77E-08
reproduction (GO:0000003)	576	308	199.46	+	1.54	3.38E-06
intracellular protein transport (GO:0006886)	477	255	165.18	+	1.54	7.68E-05
cation transport (GO:0006812)	387	206	134.01	+	1.54	2.21E-03
generation of precursor metabolites and energy (GO:0006091)	440	233	152.36	+	1.53	5.20E-04
macromolecule localization (GO:0033036)	954	498	330.35	+	1.51	4.17E-10
carboxylic acid metabolic process (GO:0019752)	840	438	290.88	+	1.51	1.55E-08
reproductive process (GO:0022414)	556	289	192.53	+	1.50	8.49E-05
oxoacid metabolic process (GO:0043436)	851	442	294.69	+	1.50	1.66E-08
cellular protein localization (GO:0034613)	553	285	191.49	+	1.49	1.84E-04
cellular macromolecule localization (GO:0070727)	555	286	192.19	+	1.49	1.57E-04
organonitrogen compound biosynthetic process (GO:1901566)	1372	706	475.10	+	1.49	3.64E-14
gene expression (GO:0010467)	1725	879	597.34	+	1.47	2.53E-17
carbohydrate derivative metabolic process (GO:1901135)	623	317	215.73	+	1.47	7.73E-05
mRNA metabolic process (GO:0016071)	474	241	164.14	+	1.47	3.54E-03
small molecule metabolic process (GO:0044281)	1436	727	497.26	+	1.46	1.85E-13
establishment of localization in cell (GO:0051649)	703	355	243.44	+	1.46	2.14E-05
cell cycle (GO:0007049)	399	201	138.17	+	1.45	4.77E-02
cellular localization (GO:0051641)	819	412	283.61	+	1.45	2.06E-06
RNA metabolic process (GO:0016070)	1458	729	504.88	+	1.44	1.42E-12
RNA processing (GO:0006396)	805	399	278.76	+	1.43	1.59E-05
intracellular transport (GO:0046907)	671	329	232.36	+	1.42	7.77E-04
organic cyclic compound metabolic process (GO:1901360)	2840	1389	983.44	+	1.41	1.60E-23
cellular nitrogen compound biosynthetic process (GO:0044271)	1298	633	449.47	+	1.41	4.30E-09
cellular nitrogen compound metabolic process (GO:0034641)	3203	1557	1109.14	+	1.40	5.82E-26
organic cyclic compound biosynthetic process (GO:1901362)	868	421	300.57	+	1.40	4.06E-05
cellular response to stress (GO:0033554)	639	309	221.27	+	1.40	3.96E-03
cellular aromatic compound metabolic process (GO:0006725)	2745	1326	950.55	+	1.39	1.00E-20
ncRNA metabolic process (GO:0034660)	518	250	179.37	+	1.39	4.56E-02
aromatic compound biosynthetic process (GO:0019438)	769	369	266.29	+	1.39	7.98E-04
cellular component assembly (GO:0022607)	670	314	232.01	+	1.35	2.47E-02
organelle organization (GO:0006996)	1279	598	442.90	+	1.35	3.08E-06
negative regulation of biological process (GO:0048519)	707	330	244.82	+	1.35	2.11E-02
nucleic acid metabolic process (GO:0090304)	1995	928	690.83	+	1.34	1.57E-10
heterocycle metabolic process (GO:0046483)	2613	1198	904.84	+	1.32	6.19E-13
nucleobase-containing compound metabolic process (GO:0006139)	2375	1082	822.42	+	1.32	7.73E-11
Unclassified (UNCLASSIFIED)	26261	4759	9093.72	-	.52	0.00E00
```

4. Findings: Stress response and xenobiotic metabolism genes are significantly expanded and significantly enriched
- (01) xylan biosynthetic process (most enriched) - integral to the integrity of cell walls and plays a role in defence against herbivory and mechanical stress
- (02) xenobiotic detoxification by transmembrane export across the plasma membrane
- (03) xenobiotic export from cell
- (04) lignin metabolic process
- (05) xyloglucan metabolic process
- (06) cell wall polysaccharide biosynthetic process
- (07) hydrogen peroxide catabolic process
- (08) hydrogen peroxide metabolic process
- (09) cellulose metabolic process
- (10) xenobiotic transport
- (13) reactive oxygen species metabolic process 

5. Clean-up
```shell
rm *.tmp
```

### Build tree using single-gene orthogroups
1. Identify single-copy orthogroups and their respective gene names across all **6 species**:
```shell
awk '($2 == 1) && ($3 == 1) && ($4 == 1) && ($5 == 1) && ($6 == 1) && ($7 == 1)' $ORTHOUT | cut -f1 > single_gene_list.grep
grep -f single_gene_list.grep ${DIR_ORTHOGROUPS}/Orthogroups/Orthogroups.tsv > single_gene_list.geneNames
```

2. Extract the CDS of these genes (Outputs: ${ORTHONAME}.fasta [includes sequences from each species]):
```shell
echo '#!/bin/bash
i=$1
line=$(head -n${i} single_gene_list.geneNames | tail -n1)
ORTHONAME=$(echo $line | cut -d" " -f1)
for name in $(echo $line | cut -d" " -f2-)
do
    SPECIES=$(echo $name | cut -d"|" -f1)
    GENE_NAME=$(echo $name | cut -d"|" -f2)
    julia extract_sequence_using_name_query.jl \
        ${SPECIES}.cds \
        ${GENE_NAME} \
        ${ORTHONAME}-${SPECIES}.fasta \
        ${SPECIES} \
        false
done
cat ${ORTHONAME}-*.fasta > ${ORTHONAME}.fasta
rm ${ORTHONAME}-*.fasta
' > parallel_extract_single_gene_orthogroups.sh
chmod +x parallel_extract_single_gene_orthogroups.sh
time \
parallel \
./parallel_extract_single_gene_orthogroups.sh {} \
::: $(seq 1 $(cat single_gene_list.geneNames | wc -l))
```

3. Align CDS (Outputs: ${ORTHOLOG}.NT.cds [nucleotide alignments] and ${ORTHOLOG}.AA.prot [amino acid alignments])
```shell
echo '#!/bin/bash
f=$1
MACSE=$2
ORTHOLOG=${f%.fasta*}
# Align the CDS across species
java -Xmx8G \
    -jar ${MACSE} \
    -prog alignSequences \
    -seq ${f} \
    -out_NT ${ORTHOLOG}.aligned.unsorted.cds.tmp \
    -out_AA ${ORTHOLOG}.aligned.unsorted.prot.tmp
# Convert stop codons and frameshifts as "---" for compatibility with downstream tools
java -Xmx8G \
    -jar ${MACSE} \
    -prog exportAlignment \
    -align ${ORTHOLOG}.aligned.unsorted.cds.tmp \
    -codonForFinalStop --- \
    -codonForInternalStop NNN \
    -codonForExternalFS --- \
    -codonForInternalFS --- \
    -out_NT ${ORTHOLOG}.NT.cds \
    -out_AA ${ORTHOLOG}.AA.prot
# Clean-up
rm ${ORTHOLOG}*.tmp ${ORTHOLOG}.AA.prot # we are not using the amino acid sequences
' > parallel_align_cds.sh
chmod +x parallel_align_cds.sh
time \
parallel \
./parallel_align_cds.sh {} ${MACSE} \
::: $(ls OG*.fasta)
```

4. Build the tree (Output: ORTHOGROUPS_SINGLE_GENE.NT.timetree.nex)
```shell
TYPE=NT.cds
### Extract sequences per species (Outputs: ${ORTHONAME}-${SPECIES}.fasta)
parallel \
julia extract_sequence_using_name_query.jl \
    {1}.${TYPE} \
    {2} \
    {1}-{2}.fasta \
    {1}-{2} \
    false \
::: $(ls *.${TYPE} | sed "s/.$TYPE//g") \
::: $(grep "^>" $(ls *.${TYPE} | head -n1) | sed 's/^>//g')

### Concatenate alignments per species (Outputs: ${SPECIES}.aln)
for SPECIES in $(grep "^>" $(ls *.${TYPE} | head -n1) | sed 's/^>//g')
do
    echo $SPECIES
    # Concatenate sequences
    echo ">${SPECIES}" > ${SPECIES}.aln.tmp
    for f in $(ls *-${SPECIES}.fasta)
    do
        sed '/^>/d' $f | sed -z 's/\n//g' >> ${SPECIES}.aln.tmp
    done
    echo "" >> ${SPECIES}.aln.tmp
    julia reformat_fasta_sequence.jl ${SPECIES}.aln.tmp 50 ${SPECIES}-${TYPE%.*}.aln
    rm ${SPECIES}.aln.tmp
done

### Extract sequence lengths to build the sequence partitioning nexus file (Output: alignment_parition.${TYPE%.*}.nex)
SPECIES=$(grep "^>" $(ls *.${TYPE} | head -n1) | sed 's/^>//g' | head -n1)
echo '#nexus
begin sets;' > alignment_parition.${TYPE%.*}.nex
N0=0
for f in $(ls *-${SPECIES}.fasta)
do
    # f=$(ls *-${SPECIES}.fasta | head -n1)
    NAME=$(head -n1 $f | sed 's/>//g' | cut -d'-' -f1)
    N1=$(cat $f | sed '/^>/d' | sed -z 's/\n//g' | wc -c)
    START=$(echo "$N0 + 1" | bc)
    END=$(echo "$N0 + $N1" | bc)
    echo "charset $NAME = $START-$END;" >> alignment_parition.${TYPE%.*}.nex
    N0=$END
done
echo 'end;' >> alignment_parition.${TYPE%.*}.nex

### Concatenate species alignments (Output: ORTHOGROUPS_SINGLE_GENE.${TYPE%.*}.aln)
cat *-${TYPE%.*}.aln > ORTHOGROUPS_SINGLE_GENE.${TYPE%.*}.aln.tmp
mv ORTHOGROUPS_SINGLE_GENE.${TYPE%.*}.aln.tmp ORTHOGROUPS_SINGLE_GENE.${TYPE%.*}.aln

### Lookup divergence times between species (Output: dates.txt)
echo 'Arabidopsis_thaliana,Oryza_sativa     -160.00
Sorghum_bicolor,Lolium_perenne               -62.00
Lolium_perenne,Lolium_rigidum                 -2.74' > dates.txt

### Build tree
BOOTSTRAP_REPS=1000
THREADS=20
TIP_DATE=0
time \
iqtree2 \
    -s ORTHOGROUPS_SINGLE_GENE.${TYPE%.*}.aln \
    -p alignment_parition.${TYPE%.*}.nex \
    -B ${BOOTSTRAP_REPS} \
    -T ${THREADS} \
    --date dates.txt \
    --date-tip ${TIP_DATE} \
    --prefix ORTHOGROUPS_SINGLE_GENE.${TYPE%.*} \
    --redo
```

5. **Additional**: Compute pairwise 4DTv (Output: ORTHOGROUPS_SINGLE_GENE.NT.4DTv)
```shell
### Compute the transversion rate among 4-fold degenerate sites (Output: ${ORTHOLOG}.NT.cds.4DTv.tmp)
time \
parallel \
julia calculate_4DTv.jl {1} {1}.4DTv.tmp \
    ::: $(ls *.NT.cds)
### Concatenate pairwise 4DTv among species across single-copy orthogroups
echo -e "ORTHOGROUP\tSPECIES_1\tSPECIES_2\tn4D_sites\tnTv4D_sites\t4DTv" > ORTHOGROUPS_SINGLE_GENE.${TYPE%.*}.4DTv
for f in $(ls *.cds.4DTv.tmp)
do
    # f=$(ls *.cds.4DTv.tmp | head -n10 | tail -n1)
    n=$(cat $f | wc -l)
    printf "${f%.${TYPE}*}\n%.0s" $(seq 1 $n) > col1.tmp
    sed -z "s/ /\t/g" $f | sed -z "s/:/\t/g"  > col2n.tmp
    paste col1.tmp col2n.tmp >> ORTHOGROUPS_SINGLE_GENE.${TYPE%.*}.4DTv
done
```

6. Clean-up
```shell
rm OG*.fasta
rm OG*.NT.cds
rm single_gene_list.*
rm dates.txt
rm *-${TYPE%.*}.aln
rm *.tmp
```

### Assess whole genome duplication (WGD) events 
We will use the distribution of four-fold degenerate sites (4DTv) across multi-copy paralogs within genomes and across sing-copy gene orthologs between pairs of species
1. Prepare script to extract CDS, align, calculate 4DTv (pairwise), and divergence time (pairwise) in parallel
```shell
echo '#!/bin/bash
### NOTE: The file: multi_gene_list.geneNames contains the gene names of one species. The gene names are in the second column (space-delimited), and each gene is comma-space-delimited
j=$1
MACSE=$2
# j=1017
line=$(head -n${j} multi_gene_list.geneNames | tail -n1)
ORTHONAME=$(echo $line | cut -d" " -f1)
# Extract CDS of these multi-copy genes
for name in $(echo $line | cut -d" " -f2- | sed -z "s/, /\n/g")
do
    # name=$(echo $line | cut -d" " -f2- | sed -z "s/, /\n/g" | head -n2 | tail -n1)
    SPECIES=$(echo $name | cut -d"|" -f1)
    GENE_NAME=$(echo $name | cut -d"|" -f2 | sed -z "s/\r//g")
    julia extract_sequence_using_name_query.jl \
        ${SPECIES}.cds \
        ${GENE_NAME} \
        ${ORTHONAME}-${GENE_NAME}.fasta \
        ${GENE_NAME} \
        false
done
# Concatenate the sequences
if [ $(ls ${ORTHONAME}-*.fasta | wc -l) -gt 1 ]
then
    cat ${ORTHONAME}-*.fasta > ${ORTHONAME}.fasta
fi
# Align the CDS across species
java -Xmx8G \
    -jar ${MACSE} \
    -prog alignSequences \
    -seq ${ORTHONAME}.fasta \
    -out_NT ${ORTHONAME}.aligned.unsorted.cds.tmp \
    -out_AA ${ORTHONAME}.aligned.unsorted.prot.tmp
# Convert stop codons and frameshifts as "---" for compatibility with downstream tools
java -Xmx8G \
    -jar ${MACSE} \
    -prog exportAlignment \
    -align ${ORTHONAME}.aligned.unsorted.cds.tmp \
    -codonForFinalStop --- \
    -codonForInternalStop NNN \
    -codonForExternalFS --- \
    -codonForInternalFS --- \
    -out_NT ${ORTHONAME}.NT.cds \
    -out_AA ${ORTHONAME}.AA.prot
# Calculate 4DTv (i.e. the ratio of the number of 4-fold degenerate codons with transversion and the total number of 4-fold degenerate codons)
julia calculate_4DTv.jl ${ORTHONAME}.NT.cds ${ORTHONAME}.4DTv.tmp
##### # Calculate divergence time with KaKs_Calculator
##### grep "^>" ${ORTHONAME}.NT.cds | \
#####     sed "s/^>//g" | \
#####     sed -z "s/\n/:/g" | \
#####     sed -z "s/:$/\n/g" > ${ORTHONAME}.axt.tmp
##### grep -v "^>" ${ORTHONAME}.NT.cds >> ${ORTHONAME}.axt.tmp
##### KaKs_Calculator -i ${ORTHONAME}.axt.tmp \
#####                 -m MA \
#####                 -o ${ORTHONAME}.kaks.tmp
##### divergence_time=$(tail -n1 ${ORTHONAME}.kaks.tmp | cut -f16)
##### sed -i -z "s/\n/\t${divergence_time}\n/g" ${ORTHONAME}.4DTv.tmp
# Clean-up
rm ${ORTHONAME}*.fasta
rm ${ORTHONAME}.aligned.unsorted*.tmp ${ORTHONAME}.NT.cds ${ORTHONAME}.AA.prot
##### rm ${ORTHONAME}.axt.tmp ${ORTHONAME}.kaks.tmp
' > parallel_extract_multi_gene_orthogroups.sh
chmod +x parallel_extract_multi_gene_orthogroups.sh
```

2. Identify multi-copy paralogs (2 to 5 copies) per species, align, and estimate 4DTv
```shell
head -n1 ${ORTHOUT} | rev | cut -f5- | rev | cut -f2- | sed -z "s/\t/\n/g" > species_names.tmp
time \
for i in $(seq 1 $(cat species_names.tmp | wc -l))
do
    # i=1
    ### Extract species name
    SPECIES=$(head -n${i} species_names.tmp | tail -n1)
    idx=$(echo $i + 1 | bc)
    ### Exract names of orthogroups with 2 copies in the current species
    awk -v col="$idx" '($col >= 2) && ($col <= 5)' $ORTHOUT | cut -f1 > multi_gene_list.grep
    ### Extract names of the genes of these multi-copy orthogroups
    grep -f multi_gene_list.grep ${DIR_ORTHOGROUPS}/Orthogroups/Orthogroups.tsv | cut -f1,${idx} > multi_gene_list.geneNames
    ### Extract CDS, and align in parallel
    parallel \
    ./parallel_extract_multi_gene_orthogroups.sh {} ${MACSE} \
    ::: $(seq 1 $(cat multi_gene_list.geneNames | wc -l))
    ### Concatenate 4DTv estimates
    cat *.4DTv.tmp > ${SPECIES}.4DTv
    ### Clean-up
    rm *.4DTv.tmp
    rm multi_gene_list.grep
    rm multi_gene_list.geneNames
done
### Clean-up
rm species_names.tmp
```

### Identify EPSPS and detoxification genes
1. Download protein sequences of genes from UniProt (https://www.uniprot.org) (Outputs: ${GENE}.faa)
```shell
###############################
### SET WORKING DIRECTORIES ###
###############################
DIR_ORTHOGROUP_SEQS=${DIR}/ORTHOGROUPS/OrthoFinder/Results_*/Orthogroup_Sequences
DIR_GENES=${DIR}/TSR_NTSR_GENES
mkdir $DIR_GENES
cd $DIR_GENES

##################
### GLYPHOSATE ###
##################
### TARGET: 5-enolpyruvylshikimate-3-phosphate synthase
echo 'enolpyruvylshikimate phosphate synthase AND organism:"Arabidopsis thaliana (Mouse-ear cress) [3702]"'
wget https://www.uniprot.org/uniprot/P05466.fasta
wget https://www.uniprot.org/uniprot/P57720.fasta
echo 'enolpyruvylshikimate phosphate synthase AND organism:"Oryza sativa (Rice) [4530]"'
echo "NONE"
echo 'enolpyruvylshikimate phosphate synthase AND organism:"Zea mays (Maize) [4577]"'
echo "NONE"
echo 'enolpyruvylshikimate phosphate synthase AND organism:"Lolium rigidum (Annual ryegrass) [89674]"'
wget https://www.uniprot.org/uniprot/F4YBG5.fasta
wget https://www.uniprot.org/uniprot/F4YBG4.fasta
wget https://www.uniprot.org/uniprot/B1P936.fasta
wget https://www.uniprot.org/uniprot/B1P932.fasta
wget https://www.uniprot.org/uniprot/Q2PT59.fasta
wget https://www.uniprot.org/uniprot/Q2PT49.fasta
echo 'enolpyruvylshikimate phosphate synthase AND organism:"Lolium multiflorum (Italian ryegrass) (Lolium perenne subsp. multiflorum) [4521]"'
wget https://www.uniprot.org/uniprot/A0A2D2AP30.fasta
wget https://www.uniprot.org/uniprot/A0A2D2AP54.fasta
cat *.fasta > EPSPS.faa
rm *.fasta

# ################
# ### PARAQUAT ###
# ################
############################
### DETOXIFICATION GENES ###
############################

##################################
### 1.) superoxide dismutase (SOD)
echo 'superoxide dismutase AND organism:"Arabidopsis thaliana (Mouse-ear cress) [3702]"'
wget https://www.uniprot.org/uniprot/P24704.fasta
wget https://www.uniprot.org/uniprot/O78310.fasta
echo 'superoxide dismutase AND organism:"Oryza sativa (Rice) [4530]"'
wget https://www.uniprot.org/uniprot/Q43803.fasta
wget https://www.uniprot.org/uniprot/Q01JW6.fasta
echo 'superoxide dismutase AND organism:"Zea mays (Maize) [4577]"'
wget https://www.uniprot.org/uniprot/P09233.fasta
wget https://www.uniprot.org/uniprot/P23346.fasta
echo 'superoxide dismutase AND organism:"Lolium rigidum (Annual ryegrass) [89674]"'
echo "NONE!"
echo 'superoxide dismutase AND organism:"Lolium multiflorum (Italian ryegrass) (Lolium perenne subsp. multiflorum) [4521]"'
echo "NONE!"
cat *.fasta > SOD.faa
rm *.fasta

##################################
### 2.) ascorbate peroxidase (APX)
echo 'ascorbate peroxidase AND organism:"Arabidopsis thaliana (Mouse-ear cress) [3702]"'
wget https://www.uniprot.org/uniprot/Q42593.fasta
wget https://www.uniprot.org/uniprot/Q42564.fasta
echo 'ascorbate peroxidase AND organism:"Oryza sativa (Rice) [4530]"'
wget https://www.uniprot.org/uniprot/Q9ST81.fasta
wget https://www.uniprot.org/uniprot/Q25AM6.fasta
echo 'ascorbate peroxidase AND organism:"Zea mays (Maize) [4577]"'
wget https://www.uniprot.org/uniprot/B4G031.fasta
wget https://www.uniprot.org/uniprot/C0HEE6.fasta
echo 'ascorbate peroxidase AND organism:"Lolium rigidum (Annual ryegrass) [89674]"'
echo "NONE!"
echo 'ascorbate peroxidase AND organism:"Lolium multiflorum (Italian ryegrass) (Lolium perenne subsp. multiflorum) [4521]"'
echo "NONE!"
cat *.fasta > APX.faa
rm *.fasta

#######################################
### 3.) glutathione S-transferase (GST)
echo 'glutathione S-transferase AND organism:"Arabidopsis thaliana (Mouse-ear cress) [3702]"'
wget https://www.uniprot.org/uniprot/Q96266.fasta
wget https://www.uniprot.org/uniprot/Q9FRL8.fasta
echo 'glutathione S-transferase AND organism:"Oryza sativa (Rice) [4530]"'
echo "NONE!"
echo 'glutathione S-transferase AND organism:"Zea mays (Maize) [4577]"'
wget https://www.uniprot.org/uniprot/P12653.fasta
wget https://www.uniprot.org/uniprot/P46420.fasta
wget https://www.uniprot.org/uniprot/P04907.fasta
wget https://www.uniprot.org/uniprot/Q9ZP62.fasta
wget https://www.uniprot.org/uniprot/C0PG07.fasta
wget https://www.uniprot.org/uniprot/Q9FQC6.fasta
wget https://www.uniprot.org/uniprot/A0A3L6E9E2.fasta
wget https://www.uniprot.org/uniprot/A0A1D6LVU5.fasta
echo 'glutathione S-transferase AND organism:"Lolium rigidum (Annual ryegrass) [89674]"'
echo "NONE!"
echo 'glutathione S-transferase AND organism:"Lolium multiflorum (Italian ryegrass) (Lolium perenne subsp. multiflorum) [4521]"'
echo "NONE!"
cat *.fasta > GST.faa
rm *.fasta

#############################################
### 4.) monodehydroascorbate reductase (MDAR)
echo 'monodehydroascorbate reductase AND organism:"Arabidopsis thaliana (Mouse-ear cress) [3702]"'
wget https://www.uniprot.org/uniprot/P92947.fasta
wget https://www.uniprot.org/uniprot/Q9LK94.fasta
wget https://www.uniprot.org/uniprot/Q9LFA3.fasta
wget https://www.uniprot.org/uniprot/F4J849.fasta
echo 'monodehydroascorbate reductase AND organism:"Oryza sativa (Rice) [4530]"'
echo "NONE!"
echo 'monodehydroascorbate reductase AND organism:"Zea mays (Maize) [4577]"'
wget https://www.uniprot.org/uniprot/C4J4E4.fasta
wget https://www.uniprot.org/uniprot/B8A028.fasta
wget https://www.uniprot.org/uniprot/A0A3L6G1K4.fasta
wget https://www.uniprot.org/uniprot/B4FQK0.fasta
wget https://www.uniprot.org/uniprot/A0A3L6ESY8.fasta
wget https://www.uniprot.org/uniprot/C0P4M0.fasta
echo 'monodehydroascorbate reductase AND organism:"Lolium rigidum (Annual ryegrass) [89674]"'
echo "NONE!"
echo 'monodehydroascorbate reductase AND organism:"Lolium multiflorum (Italian ryegrass) (Lolium perenne subsp. multiflorum) [4521]"'
echo "NONE!"
cat *.fasta > MDAR.faa
rm *.fasta

####################################
### 5.) glutathione peroxidase (GPX)
echo 'glutathione peroxidase AND organism:"Arabidopsis thaliana (Mouse-ear cress) [3702]"'
wget https://www.uniprot.org/uniprot/O22850.fasta
wget https://www.uniprot.org/uniprot/Q9LYB4.fasta
echo 'glutathione peroxidase AND organism:"Oryza sativa (Rice) [4530]"'
wget https://www.uniprot.org/uniprot/Q259Q9.fasta
wget https://www.uniprot.org/uniprot/Q9FEV2.fasta
wget https://www.uniprot.org/uniprot/Q8L8G3.fasta
echo 'glutathione peroxidase AND organism:"Zea mays (Maize) [4577]"'
wget https://www.uniprot.org/uniprot/Q6JAH6.fasta
wget https://www.uniprot.org/uniprot/B4FRF0.fasta
wget https://www.uniprot.org/uniprot/B6T5N2.fasta
wget https://www.uniprot.org/uniprot/A0A1D6K2F9.fasta
wget https://www.uniprot.org/uniprot/A0A3L6DIL3.fasta
echo 'glutathione peroxidase AND organism:"Lolium rigidum (Annual ryegrass) [89674]"'
echo "NONE!"
echo 'glutathione peroxidase AND organism:"Lolium multiflorum (Italian ryegrass) (Lolium perenne subsp. multiflorum) [4521]"'
echo "NONE!"
cat *.fasta > GPX.faa
rm *.fasta

################################
### 6.) cytochrome P450 (CYP450)
echo 'cytochrome P450 AND organism:"Arabidopsis thaliana (Mouse-ear cress) [3702]"'
wget https://www.uniprot.org/uniprot/O65782.fasta
wget https://www.uniprot.org/uniprot/Q9SMP5.fasta
echo 'cytochrome P450 AND organism:"Oryza sativa (Rice) [4530]"'
wget https://www.uniprot.org/uniprot/Q712G4.fasta
wget https://www.uniprot.org/uniprot/Q8S3F7.fasta
echo 'cytochrome P450 AND organism:"Zea mays (Maize) [4577]"'
wget https://www.uniprot.org/uniprot/Q43246.fasta
wget https://www.uniprot.org/uniprot/A0A1D6IZ40.fasta
echo 'cytochrome P450 AND organism:"Lolium rigidum (Annual ryegrass) [89674]"'
wget https://www.uniprot.org/uniprot/A5LFV4.fasta
wget https://www.uniprot.org/uniprot/Q6F4F2.fasta
wget https://www.uniprot.org/uniprot/Q6F4F4.fasta
wget https://www.uniprot.org/uniprot/A0A858NLU3.fasta
echo 'cytochrome P450 AND organism:"Lolium multiflorum (Italian ryegrass) (Lolium perenne subsp. multiflorum) [4521]"'
echo "NONE!"
cat *.fasta > CYP450.faa
rm *.fasta

##############################################
### 7.) ATP-binding cassette transporter (ABC)
echo 'ABC transporter AND organism:"Arabidopsis thaliana (Mouse-ear cress) [3702]"'
wget https://www.uniprot.org/uniprot/Q9LVM1.fasta
wget https://www.uniprot.org/uniprot/Q94FB9.fasta
echo 'ABC transporter AND organism:"Oryza sativa (Rice) [4530]"'
wget https://www.uniprot.org/uniprot/A0A650FN06.fasta
wget https://www.uniprot.org/uniprot/Q94HL6.fasta
wget https://www.uniprot.org/uniprot/Q01HZ4.fasta
wget https://www.uniprot.org/uniprot/Q25AJ5.fasta
echo 'ABC transporter AND organism:"Zea mays (Maize) [4577]"'
wget https://www.uniprot.org/uniprot/A7KVC2.fasta
wget https://www.uniprot.org/uniprot/C0P683.fasta
wget https://www.uniprot.org/uniprot/A0A317YHC9.fasta
wget https://www.uniprot.org/uniprot/A0A1D6E8S2.fasta
echo 'ABC transporter AND organism:"Lolium rigidum (Annual ryegrass) [89674]"'
echo "NONE!"
echo 'ABC transporter AND organism:"Lolium multiflorum (Italian ryegrass) (Lolium perenne subsp. multiflorum) [4521]"'
echo "NONE!"
cat *.fasta > ABC.faa
rm *.fasta

# #################
# ### CLETHODIM ###
# #################
# ### TARGET: Acetyl Co-A carboxylase
# echo 'acetyl-coenzyme a carboxylase AND organism:"Arabidopsis thaliana (Mouse-ear cress) [3702]" '
# wget https://www.uniprot.org/uniprot/Q38970.fasta
# wget https://www.uniprot.org/uniprot/Q9LD43.fasta
# echo 'acetyl-coenzyme a carboxylase AND organism:"Oryza sativa (Rice) [4530]"'
# wget https://www.uniprot.org/uniprot/P0C2Y2.fasta
# echo 'acetyl-coenzyme a carboxylase AND organism:"Zea mays (Maize) [4577]"'
# wget https://www.uniprot.org/uniprot/Q41743.fasta
# wget https://www.uniprot.org/uniprot/A0A1D6J3Q3.fasta
# echo 'acetyl-coenzyme a carboxylase AND organism:"Lolium rigidum (Annual ryegrass) [89674]"'
# wget https://www.uniprot.org/uniprot/Q3HWZ8.fasta
# wget https://www.uniprot.org/uniprot/A5JJU5.fasta
# wget https://www.uniprot.org/uniprot/Q7XHK2.fasta
# echo 'acetyl-coenzyme a carboxylase AND organism:"Lolium multiflorum (Italian ryegrass) (Lolium perenne subsp. multiflorum) [4521]"'
# wget https://www.uniprot.org/uniprot/A0A4D6X5K8.fasta
# wget https://www.uniprot.org/uniprot/Q8VWG1.fasta
# cat *.fasta > ACCase.faa
# rm *.fasta

# #################
# ### INTERCEPT ### Imazamox, and
# ################# Imazapyr
# ### TARGET: Acetolactate synthase a.k.a. acetohydroxy acid synthase
# echo 'acetolactate synthase AND organism:"Arabidopsis thaliana (Mouse-ear cress) [3702]"'
# wget https://www.uniprot.org/uniprot/A0A178VL64.fasta
# wget https://www.uniprot.org/uniprot/P17597.fasta
# echo 'acetolactate synthase AND organism:"Oryza sativa (Rice) [4530]"'
# wget https://www.uniprot.org/uniprot/A0A5J6D4R6.fasta
# wget https://www.uniprot.org/uniprot/Q01LD9.fasta
# echo 'acetolactate synthase AND organism:"Zea mays (Maize) [4577]"'
# wget https://www.uniprot.org/uniprot/Q41769.fasta
# wget https://www.uniprot.org/uniprot/K7TWQ8.fasta
# echo 'acetolactate synthase AND organism:"Lolium rigidum (Annual ryegrass) [89674]"'
# wget https://www.uniprot.org/uniprot/A0A5B9T5W5.fasta
# wget https://www.uniprot.org/uniprot/A7XBQ0.fasta
# echo 'acetolactate synthase AND organism:"Lolium multiflorum (Italian ryegrass) (Lolium perenne subsp. multiflorum) [4521]"'
# wget https://www.uniprot.org/uniprot/Q9FUD0.fasta
# wget https://www.uniprot.org/uniprot/A0A2R4NC54.fasta
# cat *.fasta > ALS.faa
# rm *.fasta

# ###############
# ### LUXIMAX ###
# ############### Cinmethylin
# ### TARGET: acyl-acp thioesterase
# echo 'acyl-acp thioesterase AND organism:"Arabidopsis thaliana (Mouse-ear cress) [3702]"'
# wget https://www.uniprot.org/uniprot/Q42561.fasta
# wget https://www.uniprot.org/uniprot/Q9SJE2.fasta
# wget https://www.uniprot.org/uniprot/Q9SV64.fasta
# wget https://www.uniprot.org/uniprot/Q8W583.fasta
# wget https://www.uniprot.org/uniprot/F4HX80.fasta
# echo 'acyl-acp thioesterase AND organism:"Oryza sativa (Rice) [4530]"'
# echo "NONE!"
# echo 'acyl-acp thioesterase AND organism:"Zea mays (Maize) [4577]"'
# wget https://www.uniprot.org/uniprot/A0A3L6FUI9.fasta
# wget https://www.uniprot.org/uniprot/A0A077D597.fasta
# wget https://www.uniprot.org/uniprot/A0A3L6FVN3.fasta
# wget https://www.uniprot.org/uniprot/A0A3L6E7J9.fasta
# wget https://www.uniprot.org/uniprot/K7V747.fasta
# echo 'acyl-acp thioesterase AND organism:"Lolium rigidum (Annual ryegrass) [89674]"'
# echo "NONE!"
# echo 'acyl-acp thioesterase AND organism:"Lolium multiflorum (Italian ryegrass) (Lolium perenne subsp. multiflorum) [4521]"'
# echo "NONE!"
# cat *.fasta > AACPT.faa
# rm *.fasta

# #################
# ### OVERWATCH ###
# ################# Bixlozone
# ### TARGET: deoxyxylulose 5-phosphate synthase (DXS)
# echo 'deoxyxylulose 5-phosphate synthase AND organism:"Arabidopsis thaliana (Mouse-ear cress) [3702]"'
# wget https://www.uniprot.org/uniprot/Q38854.fasta
# wget https://www.uniprot.org/uniprot/Q8LAD0.fasta
# wget https://www.uniprot.org/uniprot/Q9XFS9.fasta
# wget https://www.uniprot.org/uniprot/Q0WUB4.fasta
# echo 'deoxyxylulose 5-phosphate synthase AND organism:"Oryza sativa (Rice) [4530]"'
# echo "NONE!"
# echo 'deoxyxylulose 5-phosphate synthase AND organism:"Zea mays (Maize) [4577]"'
# echo "NONE!"
# echo 'deoxyxylulose 5-phosphate synthase AND organism:"Lolium rigidum (Annual ryegrass) [89674]"'
# echo "NONE!"
# echo 'deoxyxylulose 5-phosphate synthase AND organism:"Lolium multiflorum (Italian ryegrass) (Lolium perenne subsp. multiflorum) [4521]"'
# echo "NONE!"
# cat *.fasta > DXS.faa
# rm *.fasta

# #################
# ### BOXERGOLD ### S-metolachlor, and
# ################# Prosulfocarb
# ##############
# ### SAKURA ###
# ############## Pyroxasulfone
# #################
# ### TRIALLATE ###
# #################
# ### TARGET: fatty acid elongase (FAE)
# echo 'fatty acid elongase AND organism:"Arabidopsis thaliana (Mouse-ear cress) [3702]"'
# wget https://www.uniprot.org/uniprot/Q38860.fasta
# wget https://www.uniprot.org/uniprot/V9Z5W3.fasta
# wget https://www.uniprot.org/uniprot/Q5XEP9.fasta
# wget https://www.uniprot.org/uniprot/Q9MAM3.fasta
# wget https://www.uniprot.org/uniprot/Q570B4.fasta
# echo 'fatty acid elongase AND organism:"Oryza sativa (Rice) [4530]"'
# echo "NONE!"
# echo 'fatty acid elongase AND organism:"Zea mays (Maize) [4577]"'
# wget https://www.uniprot.org/uniprot/A0A1D6MV48.fasta
# wget https://www.uniprot.org/uniprot/A0A1D6L0Y4.fasta
# wget https://www.uniprot.org/uniprot/A0A3L6FEW8.fasta
# wget https://www.uniprot.org/uniprot/Q6A4M2.fasta
# wget https://www.uniprot.org/uniprot/A0A1D6P2P8.fasta
# echo 'fatty acid elongase AND organism:"Lolium rigidum (Annual ryegrass) [89674]"'
# echo "NONE!"
# echo 'fatty acid elongase AND organism:"Lolium multiflorum (Italian ryegrass) (Lolium perenne subsp. multiflorum) [4521]"'
# echo "NONE!"
# cat *.fasta > FAE.faa
# rm *.fasta

# ################
# ### ATRAZINE ###
# ################
# ### TARGET: Photosystem II protein D1
# echo 'Photosystem II protein D1 AND organism:"Arabidopsis thaliana (Mouse-ear cress) [3702]"'
# wget https://www.uniprot.org/uniprot/P83755.fasta
# wget https://www.uniprot.org/uniprot/A0A1B1W4S7.fasta
# echo 'Photosystem II protein D1 AND organism:"Oryza sativa (Rice) [4530]"'
# wget https://www.uniprot.org/uniprot/P0C432.fasta
# wget https://www.uniprot.org/uniprot/A0A0K0LK42.fasta
# echo 'Photosystem II protein D1 AND organism:"Zea mays (Maize) [4577]"'
# wget https://www.uniprot.org/uniprot/P48183.fasta
# wget https://www.uniprot.org/uniprot/A0A3L6DET6.fasta
# echo 'Photosystem II protein D1 AND organism:"Lolium rigidum (Annual ryegrass) [89674]"'
# echo "NONE!"
# echo 'Photosystem II protein D1 AND organism:"Lolium multiflorum (Italian ryegrass) (Lolium perenne subsp. multiflorum) [4521]"'
# wget https://www.uniprot.org/uniprot/B8Q684.fasta
# cat *.fasta > psbA.faa
# rm *.fasta
```

2. Generate BLAST database for each orthogroup (Outputs: ${ORTHOGROUP}.*)
```shell
echo '#!/bin/bash
f=$1
# f=$(find $DIR_ORTHOGROUP_SEQS -name "*fa" | head -n1)
makeblastdb \
    -in $f \
    -dbtype prot
' > prepare_blastdb_per_orthogroup.sh
chmod +x prepare_blastdb_per_orthogroup.sh

find $DIR_ORTHOGROUP_SEQS -name "*.fa" > orthogroup_faa_list.tmp
split --additional-suffix=".tmp" --lines 10000 orthogroup_faa_list.tmp

time \
for f in $(ls x*.tmp)
do
parallel ./prepare_blastdb_per_orthogroup.sh {} ::: $(cat $f)
done
```

3. Blastp (Outputs: ${GENE}-${ORTHOGROUP}.blastout)
```shell
echo '#!/bin/bash
GENE=$1
f=$2
# GENE=EPSPS
# f=$(find $DIR_ORTHOGROUP_SEQS -name "*.fa" | head -n10 | tail -n1)
ortname=$(basename $f)
blastp \
    -db $f \
    -query ${GENE}.faa \
    -out ${GENE}-${ortname%.fa*}.blastout \
    -evalue 1e-10 \
    -max_hsps 1 \
    -outfmt "6 qseqid sseqid slen qlen length evalue bitscore qcovs qcovhsp"
if [ $(cat ${GENE}-${ortname%.fa*}.blastout | wc -l) -eq 0 ]
then
rm ${GENE}-${ortname%.fa*}.blastout
fi
' > blastp_and_remove_no_hits.sh
chmod +x blastp_and_remove_no_hits.sh

time \
for f in $(ls x*.tmp)
do
parallel ./blastp_and_remove_no_hits.sh \
    {1} \
    {2} \
    ::: $(ls *.faa | sed 's/.faa$//g') \
    ::: $(cat $f)
done

rm *.tmp
```

### Enrichment of stress-related genes: Do we have more ortholog members for herbicide and stress-related genes in Lolium rigidum compared with the other species?
1. Extract orthologs per gene (Outputs: ${GENE}.ortho)
```shell
echo 'args = commandArgs(trailingOnly=TRUE)
# args = c("EPSPS")
gene = args[1]
fnames = list.files()
fnames = fnames[grep(paste0(gene, "-"), fnames)]
fnames = fnames[grep("blastout$", fnames)]
orthoout = c()
for (f in fnames){
    # f = fnames[1]
    orthogroup = unlist(strsplit(unlist(strsplit(f, "-"))[2], "[.]"))[1]
    df = read.delim(f, header=FALSE)
    if (mean(df[,ncol(df)]) >= 50){
        orthoout = c(orthoout, orthogroup)
    }
}
out = data.frame(gene=rep(gene, times=length(orthoout)), orthogroup=orthoout)
write.table(out, file=paste0(gene, ".ortho"), col.names=FALSE, row.names=FALSE, quote=FALSE, sep="\t")
' > extract_orthogroups_per_gene.R

time \
parallel Rscript extract_orthogroups_per_gene.R \
    {} ::: $(ls *.blastout | cut -d"-" -f1 | sort | uniq)

mkdir BLASTOUT/
mv *.blastout BLASTOUT/
```

2. Infer gene family expansion and contraction (Outpus: ${GENE}.conex)
```shell
### Input files
ORTHOUT=${DIR}/ORTHOGROUPS/orthogroups_gene_counts_families_go.out
TREE=${DIR_ORTHOGROUPS}/Species_Tree/SpeciesTree_rooted.txt

### Extract TSR and NTSR gene families (Output: ${GENE}.orthocounts)
rev ${ORTHOUT} | cut -f5- | rev > coln.tmp
awk -F'\t' '{print $(NF-1)}' ${ORTHOUT} > col1.tmp
paste -d'\t' col1.tmp coln.tmp > counts.tmp

### Infer expansion and contraction of gene families
time \
for GENE in $(ls *.ortho | sed 's/.ortho//g')
do
    # GENE=$(ls *.ortho | sed 's/.ortho//g' | head -n1)
    ### Extract gene family counts
    cat ${GENE}.ortho | cut -f2 > ${GENE}.ortho.tmp
    head -n1 counts.tmp > ${GENE}.orthocounts
    grep -f ${GENE}.ortho.tmp counts.tmp >> ${GENE}.orthocounts
    ### Run with lambda_i ~ Gamma(alpha), for each of the i_th gene family category (Output: ${GENE}_CAFE_Gamma100_results/)
    cafe5 \
        --infile ${GENE}.orthocounts \
        --tree ${TREE} \
        --n_gamma_cats 100 \
        --cores 31 \
        --pvalue 0.01 \
        --output_prefix ${GENE}_CAFE_Gamma100_results
    ### Output(Output: ${GENE}.conex)
    echo -e "Species\tExpansion\tContraction" > ${GENE}.conex
    grep -v "^#" ${GENE}_CAFE_Gamma100_results/Gamma_clade_results.txt | \
        grep -v "^<" | \
        sed 's/<..>//g' | \
        sed 's/<.>//g' >> ${GENE}.conex
done

### Clean-up
rm *.tmp
```

### dN/dS assessment of EPSPS gene
Note: we use "gene" to refer to TSR and NTSR genes, as well as orthogroups within TSR/NTSR genes per species. Apologies for any misunderstanding.

1. Extract EPSPS CDS (i.e. all orthologs and paralogs within blast-hit orthologs) (Outputs: ${species}-${gene}-${ortho}.cds)
```shell
### Extract species names and number of species
head -n1 ${DIR_ORTHOGROUPS}/Orthogroups/Orthogroups.tsv | cut -f2- > species_names.tmp
NSPECIES=$(awk -F"\t" "{print NF}" species_names.tmp)

### Extract gene names
echo '#!/bin/bash
f=$1
DIR_ORTHOGROUPS=$2
NSPECIES=$3
# f=$(ls *.ortho | head -n13 | tail -n1)
gene=$(echo ${f%.ortho*})
for ortho in $(cut -f2 $f)
do
    # ortho=$(cut -f2 $f | head -n1 | tail -n1)
    grep "$ortho" ${DIR_ORTHOGROUPS}/Orthogroups/Orthogroups.tsv > \
        ${gene}-${ortho}-list_gene_names.tmp
    for i in $(seq 1 $NSPECIES)
    do
        # i=1
        species=$(cut -f$i species_names.tmp | sed -z "s/\r//g" | sed -z "s/\n//g")
        col=$(echo $i + 1 | bc)
        cut -f${col} ${gene}-${ortho}-list_gene_names.tmp | \
            sed -z "s/, /\n/g" | \
            sed "s/$species|//g" | \
            sed "/^$/d" | \
            sed "/\r/d" > \
            ${species}-${gene}-${ortho}-list_gene_names.tmp
    if [ $(cat ${species}-${gene}-${ortho}-list_gene_names.tmp | wc -l) -lt 1 ]
    then
        rm ${species}-${gene}-${ortho}-list_gene_names.tmp
    fi
    done
    rm  ${gene}-${ortho}-list_gene_names.tmp
done
' > extract_gene_names.sh
chmod +x extract_gene_names.sh
# time \
# parallel ./extract_gene_names.sh \
#     {} \
#     ${DIR_ORTHOGROUPS} \
#     ${NSPECIES} \
#     ::: $(ls *.ortho)
time ./extract_gene_names.sh \
    EPSPS.ortho \
    ${DIR_ORTHOGROUPS} \
    ${NSPECIES}

### Extract gene sequences
echo '#!/bin/bash
f=$1
DIR=$2
SRC=$3
species=$(echo $f | cut -d"-" -f1)
gene=$(echo $f | cut -d"-" -f2)
ortho=$(echo $f | cut -d"-" -f3)
for query in $(cat $f)
do
    julia ${SRC}/extract_sequence_using_name_query.jl \
                    ${DIR}/${species}.cds \
                    ${query} \
                    ${species}-${gene}-${ortho}-${query}.cds.tmp \
                    ${species}-${gene}-${ortho}-${query} \
                    false
done
cat ${species}-${gene}-${ortho}-*.cds.tmp > ${species}-${gene}-${ortho}.cds
rm $f
' > extract_sequences_in_parallel.sh
chmod +x extract_sequences_in_parallel.sh
# time \
# parallel ./extract_sequences_in_parallel.sh \
#     {} \
#     ${DIR} \
#     ${SRC} \
#     ::: $(ls *-list_gene_names.tmp)
time \
parallel ./extract_sequences_in_parallel.sh \
    {} \
    ${DIR} \
    ${SRC} \
    ::: $(ls *-EPSPS-*-list_gene_names.tmp)
```

2. Merge per orthogropup prior to alignment (Outputs: ${gene}-${ortho}.cds)
```shell
for ortho in $(ls *.cds | cut -d"-" -f3 | cut -d"." -f1 | sort | uniq)
do
    # ortho=$(ls *.cds | cut -d"-" -f3 | cut -d"." -f1 | sort | uniq | head -n13 | tail -n1)
    fname_output=$(ls *-${ortho}.cds | cut -d"-" -f2 | sort | uniq | sed -z "s/\n/-/g")${ortho}.cds.tmp
    cat *-${ortho}.cds > ${fname_output}
    rm *-${ortho}.cds
    mv ${fname_output} ${fname_output%.tmp*}
done
rm *.tmp
```

3. Align CDS per orthogroup per gene (Outputs: ${gene}-${ortho}.aln)
```shell
echo '#!/bin/bash
f=$1
ext=$2
DIR=$3
MACSE=$4
# f=Zea_mays-SOD-OG0026358.cds; ext=cds
java -Xmx8G \
    -jar ${MACSE} \
    -prog alignSequences \
    -seq ${f} \
    -out_NT ${f%.${ext}*}.aligned.unsorted.cds.tmp \
    -out_AA ${f%.${ext}*}.aligned.unsorted.prot.tmp
# Convert stop codons and frameshifts as "---" for compatibility with downstream tools
java -Xmx8G \
    -jar ${MACSE} \
    -prog exportAlignment \
    -align ${f%.${ext}*}.aligned.unsorted.cds.tmp \
    -codonForFinalStop --- \
    -codonForInternalStop NNN \
    -codonForExternalFS --- \
    -codonForInternalFS --- \
    -out_NT ${f%.${ext}*}.aln \
    -out_AA ${f%.${ext}*}.AA.prot
' > align_in_parallel.sh
chmod +x align_in_parallel.sh
time \
parallel ./align_in_parallel.sh {} cds ${DIR} ${MACSE} \
    ::: $(ls *.cds)

rm *.tmp *.AA.prot
```

4. Remove alignments (also the corresponding cds) without Lolium rigidum genes (Outputs: ${gene}-${ortho}.aln)
```shell
for f in $(ls *.aln)
do
    # f=$(ls *.aln | head -n1 | tail -n1)
    x=$(grep "^>Lolium_rigidum" $f | wc -l)
    if [ $x -eq 0 ]
    then
        rm $f
        rm ${f%.aln*}.cds
    fi
done
```

5. Create pairwise cds alignments using the first Lolium rigidum alignment as the focal alignment per orthogroup per gene (Outputs: ${gene}-${ortho}.aln.pw)
```shell
echo '#!/bin/bash
f=$1
# f=$(ls *.aln | head -n1 | tail -n1)

for i in $(seq 1 $(grep "^>Lolium_rigidum" $f | wc -l))
do
    focal_aln=$(grep "^>Lolium_rigidum" $f | head -n${i} | tail -n1 | sed "s/^>//g")
    grep -A1 "^>Lolium_rigidum" $f | head -n$(echo "$i * 2" | bc) | tail -n1 > ${f}.FOCAL_SEQ.tmp
    touch ${f}-${i}.pw.tmp
    ### temp file without the focal alignment
    sed -e "/$focal_aln/,+1d" $f > ${f}.tmp
    for line in $(seq 2 2 $(cat ${f}.tmp | wc -l))
    do
        # line=2
        curr_aln_name=$(head -n$(echo $line -1 | bc) ${f}.tmp | tail -n1 | sed "s/>//g")
        if [ ${focal_aln} != ${curr_aln_name} ]
        then
            name=${focal_aln}--:--${curr_aln_name}
            echo $name >> ${f}-${i}.pw.tmp                           ### alignment pair name
            cat ${f}.FOCAL_SEQ.tmp >> ${f}-${i}.pw.tmp               ### focal alignment
            head -n${line} ${f}.tmp | tail -n1 >> ${f}-${i}.pw.tmp   ### current alignment
            echo "" >> ${f}-${i}.pw.tmp
        fi
    done
    ### Clean-up
    mv ${f}-${i}.pw.tmp ${f}-${i}.pw
    rm ${f}.FOCAL_SEQ.tmp ${f}.tmp
    ### Remove single alignments
    if [ $(cat ${f}-${i}.pw | wc -l) -eq 0 ]
    then
        rm ${f}-${i}.pw
    fi
done
' > prepare_pairwise_alignments_with_Lolium_rigidum_as_focus_in_parallel.sh
chmod +x prepare_pairwise_alignments_with_Lolium_rigidum_as_focus_in_parallel.sh
time \
parallel ./prepare_pairwise_alignments_with_Lolium_rigidum_as_focus_in_parallel.sh \
    {} ::: $(ls *.aln)
```

6. KaKs_calculator2 for pairwise orthogroup gene comparisons with sliding 15-bp windows
```shell
echo '#!/bin/bash
f=$1
SRC=$2
julia ${SRC}/split_alignment_pairs.jl \
    ${f} \
    15 \
    15 \
    ${f}.windows.tmp

KaKs_Calculator \
    -m MS \
    -i ${f}.windows.tmp \
    -o ${f%.tmp*}.kaks.tmp

Rscript ${SRC}/plot_KaKs_across_windows.R \
    ${f%.tmp*}.kaks.tmp\
    0.001
' > KaKs_per_window_and_plot_in_parallel.sh
chmod +x KaKs_per_window_and_plot_in_parallel.sh
time \
parallel ./KaKs_per_window_and_plot_in_parallel.sh \
    {} ${SRC} ::: $(ls EPSPS-*.aln-*.pw)
```

### Locate paralogs in the Lolium rigidum genome for the Circos-like figure
```shell
ORT=${DIR_ORTHOGROUPS}/Orthogroups/Orthogroups.tsv
GFF=${DIR}/Lolium_rigidum.gff

echo '
using ProgressMeter
fname_paralogs = ARGS[1]
fname_annotations = ARGS[2]
fname_output = ARGS[3]
# fname_paralogs = "ORTHOGROUPS/OrthoFinder/Results_May13/Orthogroups/Orthogroups.tsv"
# fname_annotations = "Lolium_rigidum.gff"
# fname_output = "Lolium_rigidum.plg"

# Load gene names and start position into temproray vectors
function load_gene_names_and_coordinates(fname_annotations)
    temp_geneIDs = []
    temp_genes = []
    temp_chromosomes = []
    temp_positions = []
    file = open(fname_annotations, "r")
    seekend(file); n=position(file); seekstart(file)
    pb = Progress(n)
    while !eof(file)
        line = split(readline(file), "\t"[1])
        update!(pb, position(file))
        if line[1][1] .!= "#"[1]
            if line[3] == "CDS"
                if match(Regex("Name="), line[end]) != nothing
                    desc = split(line[end], ";"[1])
                    id = try
                            split(desc[match.(Regex("Dbxref=GeneID:"), desc) .!= nothing][1], ","[1])[1]
                        catch
                            split(desc[match.(Regex("Dbxref=Genbank:"), desc) .!= nothing][1], ","[1])[2]
                        end
                    geneID = replace(id, "Dbxref=GeneID:" => "")
                    geneID = replace(id, "GeneID:" => "")
                    name = split(desc[match.(Regex("Name="), desc) .!= nothing][1], ","[1])[1]
                    gene = replace(name, "Name=" => "")
                    push!(temp_geneIDs, geneID)
                    push!(temp_genes, gene)
                    push!(temp_chromosomes, line[1])
                    push!(temp_positions, parse(Int, line[4]))
                end
            else
                continue
            end
        else
            continue
        end
    end
    close(file)
    @time idx = sortperm(temp_genes)
    temp_geneIDs = temp_geneIDs[idx]
    temp_genes = temp_genes[idx]
    temp_chromosomes = temp_chromosomes[idx]
    temp_positions = temp_positions[idx]

    # Use only one CDS/gene starting position
    genes = []
    chromosomes = []
    positions = []
    n = length(temp_genes)
    i = 1
    pb = Progress(n)
    while i < n
        id = temp_geneIDs[i]
        g = temp_genes[i]
        c = temp_chromosomes[i]
        p = temp_positions[i]
        push!(genes, g)
        push!(chromosomes, c)
        push!(positions, p)
        while (i < n) & ( (g == temp_genes[i]) | (id == temp_geneIDs[i]) )
            i += 1
        end
        update!(pb, i)
    end
    if genes[end] != temp_genes[n]
        push!(genes, temp_genes[n])
        push!(chromosomes, temp_chromosomes[n])
        push!(positions, temp_positions[n])
    end

    return(genes, chromosomes, positions)
end
genes, chromosomes, positions = load_gene_names_and_coordinates(fname_annotations)

# Add orthogroup labels, i.e. paralog classification
function add_paralogs(fname_paralogs, genes)
    paralogs = repeat(["None"], length(genes))
    file = open(fname_paralogs, "r")
    seekend(file); n=position(file); seekstart(file)
    pb = Progress(n)
    header = split(readline(file), "\t"[1])
    idx = header .== "Lolium_rigidum"
    while !eof(file)
        line = split(readline(file), "\t"[1])
        orthogroup = line[1]
        gene_names = replace.(split(line[idx][1], ", "), "Lolium_rigidum|" => "")
        for g in gene_names
            # g = gene_names[1]
            paralogs[genes .== g] .= orthogroup
        end
        update!(pb, position(file))
    end
    close(file)
    return(paralogs)
end
paralogs = add_paralogs(fname_paralogs, genes)

# Remove remove unclassified genes
idx = paralogs .!= "None"
paralogs = paralogs[idx]
genes = genes[idx]
chromosomes = chromosomes[idx]
positions = positions[idx]

# Save the list of gene names, chromosome, position, and ortholog info into a file
file = open(fname_output, "a")
for i in 1:length(genes)
    line = string(join([genes[i], chromosomes[i], positions[i], paralogs[i]], "\t"), "\n")
    write(file, line)
end
close(file)

# Save only the top 5 paralogs with the most genes
vec_p = unique(paralogs)
vec_p_counts = []
@showprogress for p in vec_p
    push!(vec_p_counts, sum(p .== paralogs))
end

idx = sortperm(vec_p_counts)
vec_p = vec_p[idx][(end-4):end]
vec_p_counts = vec_p_counts[idx][(end-3):end]

file = open(replace(fname_output, ".plg"=>"-for_plotting.plg"), "a")
for i in 1:length(paralogs)
    if sum(paralogs[i] .== vec_p) > 0
        line = string(join([genes[i], chromosomes[i], positions[i], paralogs[i]], "\t"), "\n")
        write(file, line)
    end
end
close(file)
' > locate_paralogs.jl

time \
julia locate_paralogs.jl \
    ${ORT} \
    ${GFF} \
    ${GFF%.gff*}.plg

### These top 5 most numureous orthogroups are:
# OG0000001	47	PTHR44586	F-BOX DOMAIN CONTAINING PROTEIN, EXPRESSED
# OG0000003	44	PTHR33070	OS06G0725500 PROTEIN
# OG0000004	54	PTHR32141	FAMILY NOT NAMED
# OG0000006	43	PTHR35546	F-BOX PROTEIN INTERACTION DOMAIN PROTEIN-RELATED
# OG0000013	45	PTHR31175	AUXIN-RESPONSIVE FAMILY PROTEIN
```

### Synteny maps between Lolium rigidum and Oryza sativa
1. Extract orthologs with the most paralogs in each of the 2 genomes
```shell
ORT=${DIR_ORTHOGROUPS}/Orthogroups/Orthogroups.tsv
GFF_Lr_=${DIR}/Lolium_rigidum.gff
GFF_Os=${DIR}/Oryza-sativa.gff

echo '
using ProgressMeter
fname_paralogs = ARGS[1]
fname_annotations_Lr = ARGS[2]
fname_annotations_Os = ARGS[3]
fname_output = ARGS[4]
# fname_paralogs = "ORTHOGROUPS/OrthoFinder/Results_May13/Orthogroups/Orthogroups.tsv"
# fname_annotations_Lr = "Lolium_rigidum.gff"
# fname_annotations_Os = "Oryza_sativa.gff"
# fname_output = "Synteny_Lolium_rigidum_and_Oryza_sativa.plg"

# Load gene names and start position into temproray vectors
function load_gene_names_and_coordinates(fname_annotations)
    temp_geneIDs = []
    temp_genes = []
    temp_chromosomes = []
    temp_positions = []
    file = open(fname_annotations, "r")
    seekend(file); n=position(file); seekstart(file)
    pb = Progress(n)
    while !eof(file)
        line = split(readline(file), "\t"[1])
        update!(pb, position(file))
        if line[1][1] .!= "#"[1]
            if line[3] == "CDS"
                if match(Regex("Name="), line[end]) != nothing
                    desc = split(line[end], ";"[1])
                    idx = match.(Regex("GeneID:"), desc) .!= nothing
                    if sum(idx) > 0
                        id = split(desc[idx][1], ","[1])[end]
                        geneID = split(id, ":")[end]
                        name = split(desc[match.(Regex("Name="), desc) .!= nothing][1], ","[1])[1]
                        gene = split(name, "Name=")[end]
                        push!(temp_geneIDs, geneID)
                        push!(temp_genes, gene)
                        push!(temp_chromosomes, line[1])
                        push!(temp_positions, parse(Int, line[4]))
                    end
                end
            else
                continue
            end
        else
            continue
        end
    end
    close(file)
    @time idx = sortperm(temp_genes)
    temp_geneIDs = temp_geneIDs[idx]
    temp_genes = temp_genes[idx]
    temp_chromosomes = temp_chromosomes[idx]
    temp_positions = temp_positions[idx]

    # Use only one CDS/gene starting position
    genes = []
    chromosomes = []
    positions = []
    n = length(temp_genes)
    i = 1
    pb = Progress(n)
    while i < n
        id = temp_geneIDs[i]
        g = temp_genes[i]
        c = temp_chromosomes[i]
        p = temp_positions[i]
        push!(genes, g)
        push!(chromosomes, c)
        push!(positions, p)
        while (i < n) & ( (g == temp_genes[i]) | (id == temp_geneIDs[i]) )
            i += 1
        end
        update!(pb, i)
    end
    if genes[end] != temp_genes[n]
        push!(genes, temp_genes[n])
        push!(chromosomes, temp_chromosomes[n])
        push!(positions, temp_positions[n])
    end

    return(genes, chromosomes, positions)
end
genes_Lr, chromosomes_Lr, positions_Lr = load_gene_names_and_coordinates(fname_annotations_Lr)
genes_Os, chromosomes_Os, positions_Os = load_gene_names_and_coordinates(fname_annotations_Os)

# Add orthogroup labels, i.e. paralog classification
function add_paralogs(fname_paralogs, genes; species="Lolium_rigidum")
    paralogs = repeat(["None"], length(genes))
    file = open(fname_paralogs, "r")
    seekend(file); n=position(file); seekstart(file)
    pb = Progress(n)
    header = split(readline(file), "\t"[1])
    idx = header .== species
    while !eof(file)
        line = split(readline(file), "\t"[1])
        orthogroup = line[1]
        gene_names = replace.(split(line[idx][1], ", "), string(species, "|") => "")
        for g in gene_names
            # g = gene_names[1]
            paralogs[genes .== g] .= orthogroup
        end
        update!(pb, position(file))
    end
    close(file)
    return(paralogs)
end
paralogs_Lr = add_paralogs(fname_paralogs, genes_Lr)
paralogs_Os = add_paralogs(fname_paralogs, genes_Os, species="Oryza_sativa")

# Remove remove unclassified genes
idx = paralogs_Lr .!= "None"
paralogs_Lr = paralogs_Lr[idx]
genes_Lr = genes_Lr[idx]
chromosomes_Lr = chromosomes_Lr[idx]
positions_Lr = positions_Lr[idx]

idx = paralogs_Os .!= "None"
paralogs_Os = paralogs_Os[idx]
genes_Os = genes_Os[idx]
chromosomes_Os = chromosomes_Os[idx]
positions_Os = positions_Os[idx]

# Find common paralogs
idx = [x ∈ paralogs_Os for x in paralogs_Lr]
paralogs_Lr = paralogs_Lr[idx]
genes_Lr = genes_Lr[idx]
chromosomes_Lr = chromosomes_Lr[idx]
positions_Lr = positions_Lr[idx]

idx = [x ∈ paralogs_Lr for x in paralogs_Os]
paralogs_Os = paralogs_Os[idx]
genes_Os = genes_Os[idx]
chromosomes_Os = chromosomes_Os[idx]
positions_Os = positions_Os[idx]

### Sort
idx = sortperm(paralogs_Lr)
paralogs_Lr = paralogs_Lr[idx]
genes_Lr = genes_Lr[idx]
chromosomes_Lr = chromosomes_Lr[idx]
positions_Lr = positions_Lr[idx]

idx = sortperm(paralogs_Os)
paralogs_Os = paralogs_Os[idx]
genes_Os = genes_Os[idx]
chromosomes_Os = chromosomes_Os[idx]
positions_Os = positions_Os[idx]

# Save the list of gene names, chromosome, position, and ortholog info into a file
file = open(fname_output, "a")
for i in 1:length(genes_Lr)
    line = string(join([paralogs_Lr[i], genes_Lr[i], chromosomes_Lr[i], positions_Lr[i]], "\t"), "\n")
    write(file, line)
    end
for i in 1:length(genes_Os)
    line = string(join([paralogs_Os[i], genes_Os[i], chromosomes_Os[i], positions_Os[i]], "\t"), "\n")
    write(file, line)
end
close(file)

# Save only the top 5 paralogs with the most genes
vec_p_Lr = unique(paralogs_Lr)
vec_p_counts_Lr = []
@showprogress for p in vec_p_Lr
    push!(vec_p_counts_Lr, sum(p .== paralogs_Lr))
end
vec_p_Os = unique(paralogs_Os)
vec_p_counts_Os = []
@showprogress for p in vec_p_Os
    push!(vec_p_counts_Os, sum(p .== paralogs_Os))
end

idx = sortperm(vec_p_counts_Lr)
vec_p_Lr = vec_p_Lr[idx][(end-4):end]
vec_p_counts_Lr = vec_p_counts_Lr[idx][(end-3):end]
idx = sortperm(vec_p_counts_Os)
vec_p_Os = vec_p_Os[idx][(end-4):end]
vec_p_counts_Os = vec_p_counts_Os[idx][(end-3):end]

vec_p = copy(vec_p_Lr)
append!(vec_p, vec_p_Os)
vec_p = unique(vec_p)

file = open(replace(fname_output, ".plg"=>"-for_plotting.plg"), "a")
for i in 1:length(paralogs_Lr)
    if sum(paralogs_Lr[i] .== vec_p) > 0
        line = string(join([genes_Lr[i], chromosomes_Lr[i], positions_Lr[i], paralogs_Lr[i]], "\t"), "\n")
        write(file, line)
    end
end
for i in 1:length(paralogs_Os)
    if sum(paralogs_Os[i] .== vec_p) > 0
        line = string(join([genes_Os[i], chromosomes_Os[i], positions_Os[i], paralogs_Os[i]], "\t"), "\n")
        write(file, line)
    end
end
close(file)
' > locate_paralogs.jl

time \
julia locate_paralogs.jl \
    ${ORT} \
    ${GFF_Lr} \
    ${GFF_Os} \
    Synteny_Lolium_rigidum_and_Oryza_sativa.plg

```

2. Merge the two genomes' chromosome sequences
```shell
cat Lolium_rigidum.fasta Oryza_sativa.fasta > Synteny_Lolium_rigidum_and_Oryza_sativa.fasta
```

3. Plot genome synteny circos plot
```julia
include("Lolium_rigidum_genome_assembly_and_annotation/genome_statistics.jl")
function execute(; recompute=false, cleanup=false)
    str_filename_fasta = "Synteny_Lolium_rigidum_and_Oryza_sativa.fasta"
    str_fname_output_svg = "Synteny_Lolium_rigidum_and_Oryza_sativa.svg"
    str_filename_groupings_and_coordinates = "Synteny_Lolium_rigidum_and_Oryza_sativa-for_plotting.plg"
    vec_idx_groups_chr_pos = [4, 2, 3, 3] ### for drawing the chords, i.e. paralog positions: specify the columns of orthogroup ID, chromosome ID and start or end position
    n = 7 + 12 ## haploid chromosome number
    n_int_tick_length_bp = 100*1e+6 # make this adjustable
    n_int_tick_label_size=12
    n_int_chrom_name_size=12
    # n_int_stats_label_size=10

    ### To (re)run the genome statistics computation or to use the saved genome statistics data
    str_filename_output_jld2 = string(join(split(str_filename_fasta, ".")[1:(end-1)], "."), "-statistics.jld2")
    if recompute
        ### Split into chromosomes, count assembly size, GC content, Ns, and GC fractiopn per line
        @time n_int_assembly_size,
            n_int_assembly_GC,
            n_int_assembly_N,
            vec_str_chromosome_names,
            vec_int_chromosome_lengths,
            vec_int_chromosome_GC,
            vec_int_chromosome_N = PlotGenome.fun_fasta_lengths_GC_N(str_filename_fasta)

        ### Some assembly stats including the pseudo-chromosomes and small contigs
        vec_int_idx_sort_by_decreasing_size = sortperm(vec_int_chromosome_lengths, rev=true)
        vec_str_chromosome_names = vec_str_chromosome_names[vec_int_idx_sort_by_decreasing_size]
        vec_int_chromosome_lengths = vec_int_chromosome_lengths[vec_int_idx_sort_by_decreasing_size]
        vec_int_cumsum_length = cumsum(vec_int_chromosome_lengths)
        vec_bool_idx_50th_percentile = vec_int_cumsum_length .>= (n_int_assembly_size*0.5)
        vec_bool_idx_90th_percentile = vec_int_cumsum_length .>= (n_int_assembly_size*0.9)
        n_int_size_of_largest_chromosome = maximum(vec_int_chromosome_lengths)
        str_largest_chromosome = vec_str_chromosome_names[vec_int_chromosome_lengths.==maximum(vec_int_chromosome_lengths)]
        int_n_chromosomes_whole_assembly = length(vec_str_chromosome_names)
        int_size_whole_assembly = sum(vec_int_chromosome_lengths)
        L50_whole_assembly = collect(1:length(vec_str_chromosome_names))[vec_bool_idx_50th_percentile][1]
        L90_whole_assembly = collect(1:length(vec_str_chromosome_names))[vec_bool_idx_90th_percentile][1]
        N50_whole_assembly = vec_int_chromosome_lengths[vec_bool_idx_50th_percentile][1]
        N90_whole_assembly = vec_int_chromosome_lengths[vec_bool_idx_90th_percentile][1]

        ### Using only the 7 pseudo-chromosomes, i.e. 7 largest sequences
        vec_str_chromosome_names = vec_str_chromosome_names[1:n]
        vec_int_chromosome_lengths = vec_int_chromosome_lengths[1:n]
        vec_int_chromosome_GC = vec_int_chromosome_GC[1:n]
        vec_int_chromosome_N = vec_int_chromosome_N[1:n]

        ### Some assembly stats after excluding the small contigs
        vec_int_idx_sort_by_decreasing_size = sortperm(vec_int_chromosome_lengths, rev=true)
        vec_str_chromosome_names = vec_str_chromosome_names[vec_int_idx_sort_by_decreasing_size]
        vec_int_chromosome_lengths = vec_int_chromosome_lengths[vec_int_idx_sort_by_decreasing_size]
        vec_int_cumsum_length = cumsum(vec_int_chromosome_lengths)
        vec_bool_idx_50th_percentile = vec_int_cumsum_length .>= (n_int_assembly_size*0.5)
        vec_bool_idx_90th_percentile = vec_int_cumsum_length .>= (n_int_assembly_size*0.9)
        n_int_size_of_largest_chromosome = maximum(vec_int_chromosome_lengths)
        str_largest_chromosome = vec_str_chromosome_names[vec_int_chromosome_lengths.==maximum(vec_int_chromosome_lengths)]
        int_n_chromosomes = n
        int_size = sum(vec_int_chromosome_lengths)
        L50 = collect(1:length(vec_str_chromosome_names))[vec_bool_idx_50th_percentile][1]
        L90 = collect(1:length(vec_str_chromosome_names))[vec_bool_idx_90th_percentile][1]
        N50 = vec_int_chromosome_lengths[vec_bool_idx_50th_percentile][1]
        N90 = vec_int_chromosome_lengths[vec_bool_idx_90th_percentile][1]

        ### Sort back according to chromosome names
        vec_int_idx_sort_by_name = sortperm(vec_str_chromosome_names)
        vec_str_chromosome_names = vec_str_chromosome_names[vec_int_idx_sort_by_name]
        vec_int_chromosome_lengths = vec_int_chromosome_lengths[vec_int_idx_sort_by_name]

        ### Save the genome statistics variables
        if save_temp_files_and_data
            vec_stats = [n_int_assembly_size,
                            n_int_assembly_GC,
                            n_int_assembly_N,
                            int_n_chromosomes_whole_assembly,
                            int_size_whole_assembly,
                            L50_whole_assembly,
                            L90_whole_assembly,
                            N50_whole_assembly,
                            N90_whole_assembly,
                            vec_str_chromosome_names,
                            vec_int_chromosome_lengths,
                            vec_int_chromosome_GC,
                            vec_int_chromosome_N,
                            str_largest_chromosome,
                            int_n_chromosomes,
                            int_size,
                            L50,
                            L90,
                            N50,
                            N90]
            jldsave(str_filename_output_jld2; vec_stats)
        end
    else
        ### If we're not rerunning the statistics computation, then just load the save statistics
        n_int_assembly_size,
        n_int_assembly_GC,
        n_int_assembly_N,
        int_n_chromosomes_whole_assembly,
        int_size_whole_assembly,
        L50_whole_assembly,
        L90_whole_assembly,
        N50_whole_assembly,
        N90_whole_assembly,
        vec_str_chromosome_names,
        vec_int_chromosome_lengths,
        vec_int_chromosome_GC,
        vec_int_chromosome_N,
        str_largest_chromosome,
        int_n_chromosomes,
        int_size,
        L50,
        L90,
        N50,
        N90 = load_object(str_filename_output_jld2)
    end

    ### Some assembly stats including the pseudo-chromosomes and small contigs
    vec_int_idx_sort_by_decreasing_size = sortperm(vec_int_chromosome_lengths, rev=true)
    vec_str_chromosome_names = vec_str_chromosome_names[vec_int_idx_sort_by_decreasing_size]
    vec_int_chromosome_lengths = vec_int_chromosome_lengths[vec_int_idx_sort_by_decreasing_size]
    vec_int_cumsum_length = cumsum(vec_int_chromosome_lengths)
    vec_bool_idx_50th_percentile = vec_int_cumsum_length .>= (n_int_assembly_size*0.5)
    vec_bool_idx_90th_percentile = vec_int_cumsum_length .>= (n_int_assembly_size*0.9)
    n_int_size_of_largest_chromosome = maximum(vec_int_chromosome_lengths)
    str_largest_chromosome = vec_str_chromosome_names[vec_int_chromosome_lengths.==maximum(vec_int_chromosome_lengths)]
    int_n_chromosomes_whole_assembly = length(vec_str_chromosome_names)
    int_size_whole_assembly = sum(vec_int_chromosome_lengths)
    L50_whole_assembly = collect(1:length(vec_str_chromosome_names))[vec_bool_idx_50th_percentile][1]
    L90_whole_assembly = collect(1:length(vec_str_chromosome_names))[vec_bool_idx_90th_percentile][1]
    N50_whole_assembly = vec_int_chromosome_lengths[vec_bool_idx_50th_percentile][1]
    N90_whole_assembly = vec_int_chromosome_lengths[vec_bool_idx_90th_percentile][1]

    ### Using only the 7 pseudo-chromosomes, i.e. 7 largest sequences
    vec_str_chromosome_names = vec_str_chromosome_names[1:n]
    vec_int_chromosome_lengths = vec_int_chromosome_lengths[1:n]
    vec_int_chromosome_GC = vec_int_chromosome_GC[1:n]
    vec_int_chromosome_N = vec_int_chromosome_N[1:n]

    ### Some assembly stats after excluding the small contigs
    vec_int_idx_sort_by_decreasing_size = sortperm(vec_int_chromosome_lengths, rev=true)
    vec_str_chromosome_names = vec_str_chromosome_names[vec_int_idx_sort_by_decreasing_size]
    vec_int_chromosome_lengths = vec_int_chromosome_lengths[vec_int_idx_sort_by_decreasing_size]
    vec_int_cumsum_length = cumsum(vec_int_chromosome_lengths)
    vec_bool_idx_50th_percentile = vec_int_cumsum_length .>= (n_int_assembly_size*0.5)
    vec_bool_idx_90th_percentile = vec_int_cumsum_length .>= (n_int_assembly_size*0.9)
    n_int_size_of_largest_chromosome = maximum(vec_int_chromosome_lengths)
    str_largest_chromosome = vec_str_chromosome_names[vec_int_chromosome_lengths.==maximum(vec_int_chromosome_lengths)]
    int_n_chromosomes = n
    int_size = sum(vec_int_chromosome_lengths)
    L50 = collect(1:length(vec_str_chromosome_names))[vec_bool_idx_50th_percentile][1]
    L90 = collect(1:length(vec_str_chromosome_names))[vec_bool_idx_90th_percentile][1]
    N50 = vec_int_chromosome_lengths[vec_bool_idx_50th_percentile][1]
    N90 = vec_int_chromosome_lengths[vec_bool_idx_90th_percentile][1]

    ### Sort back according to chromosome names
    vec_int_idx_sort_by_name = sortperm(vec_str_chromosome_names)
    vec_str_chromosome_names = vec_str_chromosome_names[vec_int_idx_sort_by_name]
    vec_int_chromosome_lengths = vec_int_chromosome_lengths[vec_int_idx_sort_by_name]

    ### Standardise the chromosome length of the rice chromosomes
    vec_int_chromosome_lengths[1:12] = 10 .* vec_int_chromosome_lengths[1:12]
    str_filename_groupings_and_coordinates_STD = string(str_filename_groupings_and_coordinates, ".Osativa-chrom-lengths-standardised.tmp")
    file = open(str_filename_groupings_and_coordinates, "r")
    file_out = open(str_filename_groupings_and_coordinates_STD, "w")
    while !eof(file)
        line = split(readline(file), "\t")
        if isnothing(match(Regex("NC_02"), line[2])) == false
            line[3] = string(10 * parse(Int, line[3]))
        end
        write(file_out, string(join(line, "\t"), "\n"))
    end
    close(file)
    close(file_out)

    ### Set two colours for the two genomes
    vec_colours_chrom = vcat([repeat([palette(:default, n)[1]], 12), repeat([palette(:default, n)[2]], 7)]...)

    ### Revise the chromosome names
    vec_str_chromosome_names_NEW = vcat([string.("Os", collect(1:12)), string.("Lr", collect(1:7))]...)

    ### Plot
    plt = plot(xlim=(-1.5, 1.5), ylim=(-1.5, 1.5), size=(700,700), axis=([], false), title="")
    r=1.00; w=0.15
    PlotGenome.fun_plot_chrom_length_layer!(plt, vec_str_chromosome_names_NEW, vec_int_chromosome_lengths;
                                r=r, w=w,
                                vec_colours_chrom=vec_colours_chrom,
                                n_int_tick_length_bp=n_int_tick_length_bp,
                                n_int_tick_label_size=n_int_tick_label_size,
                                n_int_chrom_name_size=n_int_chrom_name_size,
                                add_legend=false)
    PlotGenome.fun_add_chords!(plt,
                    str_filename_groupings_and_coordinates_STD,
                    vec_str_chromosome_names,
                    vec_int_chromosome_lengths,
                    delim='\t',
                    vec_idx_groups_chr_pos=vec_idx_groups_chr_pos,
                    vec_colours=vec_colours_chrom,
                    add_legend=false,
                    colour_per_chrom=true,
                    linewidth=1,
                    r=r-w, w=w, header=false)
    ### Save as svg
    savefig(plt, str_fname_output_svg)
end
@time execute()
```

## Extract and assemble plastid genomes
```shell
python3 -m pip install
git clone https://github.com/Kinggerm/GetOrganelle.git

curl -L https://github.com/Kinggerm/GetOrganelleDB/releases/download/0.0.1/v0.0.1.tar.gz | tar zx
python3 GetOrganelle/Utilities/get_organelle_config.py \
    -a embplant_pt,embplant_mt \
    --use-local ./0.0.1

time \
python3 GetOrganelle/get_organelle_from_reads.py \
    -F embplant_pt,embplant_mt \
    -1 /data/Lolium_rigidum_ASSEMBLY/lolium_illumina/LOL-WGS2/N2009012_FA_30-434329113_SEQ/201110-X4A_L007/LOL-WGS2-1_combined_R1.fastq.gz \
    -2 /data/Lolium_rigidum_ASSEMBLY/lolium_illumina/LOL-WGS2/N2009012_FA_30-434329113_SEQ/201110-X4A_L007/LOL-WGS2-1_combined_R2.fastq.gz \
    -t 20 \
    -o Lolium_rigidum_PLASTIDS
```
