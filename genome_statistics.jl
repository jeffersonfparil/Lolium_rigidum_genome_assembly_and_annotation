println("#################################################")
println("###                                           ###")
println("### Lolium rigidum genome assembly statistics ###")
println("###                                           ###")
println("#################################################")

### Input
str_filename_fasta = "APGP_CSIRO_Lrig_flye-racon-polca-allhic-juicebox_v0.1n.fasta"
n = 7 ## haploid chromosome number

### Split the genome asembly file into reads, i.e. psuedochromosomes and uncategorised tiny contigs
### Also, measure the size of each pseudochromosome and contig
n_int_number_of_sequences = 0
vec_str_sequence_names = []
vec_int_sequence_lengths = []
FILE = open(str_filename_fasta, "r")
line = readline(FILE) ### read the first line assumes we start with the very first sequence in the first line, i.e. no need to test for "^>"
@time while !eof(FILE)
    n_int_number_of_sequences += 1
    str_sequence_name = line[2:end]
    push!(vec_str_sequence_names, str_sequence_name)
    append!(vec_int_sequence_lengths, 0)
    file = open(string(str_sequence_name, ".fasta"), "w")
    write(file, string(line, '\n'))
    line = readline(FILE)
    while (!eof(FILE)) & (line[1] != '>')
        vec_int_sequence_lengths[end] = vec_int_sequence_lengths[end] + length(line)
        write(file, string(line, '\n'))
        line = readline(FILE)
    end
    close(file)
end
close(FILE)

### Calculate the assembly size including only the 7 largest sequences (n=7 expectation)
vec_int_idx_sortperm = sortperm(vec_int_sequence_lengths, rev=true)
vec_str_chromosome_names = vec_str_sequence_names[vec_int_idx_sortperm][1:n]
vec_int_chromosome_lengths = vec_int_sequence_lengths[vec_int_idx_sortperm][1:n]
vec_int_sort_by_name = sortperm(vec_str_chromosome_names)
@show vec_str_chromosome_names = vec_str_chromosome_names[vec_int_sort_by_name]
@show vec_int_chromosome_lengths = vec_int_chromosome_lengths[vec_int_sort_by_name]
@show n_int_assembly_size = sum(vec_int_chromosome_lengths)

### GC content per line per chromosome and across the entire genome
### Ouput a file per chromosome where each line correspond to the GC count per line of the fasta file
n_int_GC_count_genomewide = 0
n_flt_min_GC_perc_genomewide = 1.0
n_flt_max_GC_perc_genomewide = 0.0
@time for str_chrom in vec_str_chromosome_names
    # str_chrom = vec_str_chromosome_names[1]
    FILE = open(string(str_chrom, ".fasta"), "r")
    file = open(string(str_chrom, "-GC_content.txt"), "w")
    line = readline(FILE)
    str_previous_base = 'N'
    while !eof(FILE)
        line = readline(FILE)
        n_int_GC_count = length(collect(eachmatch(r"[GC]", line)))
        n_int_total_bases = length(line)
        write(file, string(round(n_int_GC_count/n_int_total_bases, digits=4), '\n'))
        n_int_GC_count_genomewide = n_int_GC_count_genomewide + n_int_GC_count
    end
    close(FILE)
    close(file)
end
@show n_flt_percent_GC_content = n_int_GC_count_genomewide / n_int_assembly_size

### Circular plot

using Plots
gr()

function fun_arcshape(θ1, θ2; r=1.0, w=0.2, n_int_points=10)
    ##################
    ### TEST
    # θ1 = 0.00
    # θ2 = 0.88
    # r = 1.0
    # w = 0.2
    # n_int_points = 10
    ##################
    vec_flt_coor_top_arc =    Plots.partialcircle(θ1, θ2, n_int_points, r)
    vec_flt_coor_bottom_arc = reverse(Plots.partialcircle(θ1, θ2, n_int_points, r-w))
    Shape(vcat(vec_flt_coor_top_arc, vec_flt_coor_bottom_arc))
end

function fun_plot_circos_layer!(plt, vec_str_chromosome_names, vec_int_chromosome_lengths;
                                r=1.00, w=0.20,
                                bool_chrom=true, vec_colours_chrom = palette(:default), n_int_tick_length_bp=100*1e+6,
                                bool_GC_heatmap=false)
    ########################
    ### TEST
    # r=0.75
    # w=0.20
    # bool_chrom=true; vec_colours_chrom = palette(:default); n_int_tick_length_bp=100*1e+6
    # bool_GC_heatmap=false; vec_colours_GC=palette(:thermal)
    ########################
    ### chromosome plotting
    vec_colours_chrom = palette(:default)
    vec_colours_GC = palette(:thermal)
    n_int_tick_length_bp = 100 * 1e+6
    ### chromosome slice properties
    n_int_assembly_size = sum(vec_int_chromosome_lengths)
    n_flt_empty_slice = 0.05(2π)
    n_flt_available_circumference = 2*π - n_flt_empty_slice
    n_flt_fraction_filled = 0.9*n_flt_available_circumference
    n_flt_fraction_spacer = 0.1*n_flt_available_circumference
    vec_flt_arc_lengths_filled = n_flt_fraction_filled * (vec_int_chromosome_lengths ./ n_int_assembly_size)
    n_flt_arc_lengths_spacer = n_flt_fraction_spacer / (n-1)
    n_flt_arch_length_per_tick = n_flt_fraction_filled * n_int_tick_length_bp / n_int_assembly_size
    ### Find minimum and maximum GC content per line in each chromosome sequence file 
    if bool_GC_heatmap
        n_flt_min_GC_perc = 1.0
        n_flt_max_GC_perc = 0.0
        n_int_bp_per_chunk = Int(round(n_int_assembly_size / (length(vec_str_chromosome_names)*1e+6)))
        for i in 1:n
            # i = 1
            vec_flt_GC_content = parse.(Float64, readlines(string(vec_str_chromosome_names[i], "-GC_content.txt")))
            n_int_lines_per_chunk = Int(round(length(vec_flt_GC_content) / n_int_bp_per_chunk))
            for j in 1:n_int_bp_per_chunk
                # j = 1
                n_int_idx_start = ((j-1)*n_int_lines_per_chunk)+1
                n_int_idx_end = j*n_int_lines_per_chunk
                n_int_idx_end > length(vec_flt_GC_content) ? n_int_idx_end = length(vec_flt_GC_content) : nothing
                n_flt_mean_GC_content = sum(vec_flt_GC_content[n_int_idx_start:n_int_idx_end]) / n_int_lines_per_chunk
                n_flt_min_GC_perc > n_flt_mean_GC_content ? n_flt_min_GC_perc = n_flt_mean_GC_content : nothing
                n_flt_max_GC_perc < n_flt_mean_GC_content ? n_flt_max_GC_perc = n_flt_mean_GC_content : nothing
              end
        end
        @show n_flt_min_GC_perc
        @show n_flt_max_GC_perc
    end
    ### Intial position (unrotated)
    θ1 = 2π - (n_flt_empty_slice/2)
    ### plot per chromosome
    for i in 1:n
        # i = 1
        ### update track tail coordinate
        θ2 = θ1 - vec_flt_arc_lengths_filled[i]
        ### rotate to proper orientation
        (θ1 + π/2) > 2π ? θ1 = (θ1+π/2)-2π : nothing
        (θ2 + π/2) > 2π ? θ2 = (θ2+π/2)-2π : nothing
        ### plot outer chromosome ring
        if bool_chrom
            ### plot each chromosome
            plot!(plt, fun_arcshape(θ1,
                                    θ2,
                                    r=r,
                                    w=w),
                legend=false,
                color=vec_colours_chrom[i],
                linecolor=vec_colours_chrom[i])
            ### print ticks and tick labels
            for j in 0:floor(vec_flt_arc_lengths_filled[i] / n_flt_arch_length_per_tick)
                # j = 1
                x = θ1-(n_flt_arch_length_per_tick*j)
                plot!(plt, fun_arcshape(x,
                                        x,
                                        r=r+(w/4),
                                        w=w/4,
                                        n_int_points=1),
                      legend=false,
                      linecolor=:gray)
                vec_flt_coor_text = fun_arcshape(x, x, r=r+(w/2), w=0.00, n_int_points=1)
                annotate!(plt,
                          vec_flt_coor_text.x[1],
                          vec_flt_coor_text.y[1],
                          (string(Int(j)), 7, :gray, :center))
            end
            ### print chromosome names
            vec_flt_coor_text = fun_arcshape(θ1-(vec_flt_arc_lengths_filled[i]/2),
                                             θ1-(vec_flt_arc_lengths_filled[i]/2),
                                             r=r+(r/4), w=0.00, n_int_points=1)
            str_chromosome_name = replace(vec_str_chromosome_names[i], "Chromosome"=>"Chr")
            annotate!(plt, vec_flt_coor_text.x[1],
                      vec_flt_coor_text.y[1],
                      (str_chromosome_name, 7, :gray, :center))
        end
        ### plot GC content heatmap
        if bool_GC_heatmap
            vec_colours_GC = palette(:thermometer, 256) ## 256 colours
            vec_flt_GC_content = parse.(Float64, readlines(string(vec_str_chromosome_names[i], "-GC_content.txt")))

            n_int_bp_per_chunk = Int(round(n_int_assembly_size / (length(vec_str_chromosome_names)*1e+6)))
            n_int_lines_per_chunk = Int(round(length(vec_flt_GC_content) / n_int_bp_per_chunk))

            for j in 1:n_int_bp_per_chunk
                # j = 1
                n_int_idx_start = ((j-1)*n_int_lines_per_chunk)+1
                n_int_idx_end = j*n_int_lines_per_chunk
                n_int_idx_end > length(vec_flt_GC_content) ? n_int_idx_end = length(vec_flt_GC_content) : nothing
                n_flt_mean_GC_content = sum(vec_flt_GC_content[n_int_idx_start:n_int_idx_end]) / n_int_lines_per_chunk
                n_flt_remapped_GC = (n_flt_mean_GC_content - n_flt_min_GC_perc) / (n_flt_max_GC_perc - n_flt_min_GC_perc)
                col = vec_colours_GC[Int(round(n_flt_remapped_GC*255)+1)]
                x = θ1 - (abs(θ1-θ2)*j/n_int_bp_per_chunk)
                y = θ1 - (abs(θ1-θ2)*(j-1)/n_int_bp_per_chunk)
                vec_flt_coor_heat_slice = fun_arcshape(x, y, r=r, w=w, n_int_points=10)
                plot!(plt,
                      vec_flt_coor_heat_slice,
                      legend=false,
                      color=col,
                      linecolor=col)
            end
        end
        ### update track head coordinate
        θ1 = θ2 - n_flt_arc_lengths_spacer
    end
    return(plt)
end


plt = plot(xlim=(-1.5, 1.5), ylim=(-1.5, 1.5), size=(700,700), axis=([], false), title="Lolium rigidum genome")
### Layer 1: chromosome lengths
r=1.00; w=0.20
annotate!(plt, 0.0, (r-w/2), ("a", 10, :gray, :center))
fun_plot_circos_layer!(plt, vec_str_chromosome_names, vec_int_chromosome_lengths,
                       r=r, w=w,
                       bool_chrom=true,
                       bool_GC_heatmap=false)
### Layer 2: GC content
r=0.75; w=0.15
annotate!(plt, 0.0, (r-w/2), ("b", 10, :gray, :center))
fun_plot_circos_layer!(plt, vec_str_chromosome_names, vec_int_chromosome_lengths,
                       r=r, w=w,
                       bool_chrom=false,
                       bool_GC_heatmap=true)

savefig(plt, "test.svg")

