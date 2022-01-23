println("#################################################")
println("###                                           ###")
println("### Lolium rigidum genome assembly statistics ###")
println("###                                           ###")
println("#################################################")

### Load plotting library and set backend to gr(): GKS and openGL
using ProgressMeter
using Plots
gr()

### Split the genome asembly file into reads, i.e. psuedochromosomes and uncategorised tiny contigs
### Measure the size of each pseudochromosome and contigs
### Retains only the 1:nth largest sequences and removes contigs split file
### Returns assembly size, chromosome names, and chromosome lengths
function fun_split_fasta_count_lengths(str_filename_fasta, n)
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
    vec_str_chromosome_names = vec_str_chromosome_names[vec_int_sort_by_name]
    vec_int_chromosome_lengths = vec_int_chromosome_lengths[vec_int_sort_by_name]
    n_int_assembly_size = sum(vec_int_chromosome_lengths)
    ### Remove the uncateogised contigs
    for name in vec_str_sequence_names[vec_int_idx_sortperm][n+1:end]
        rm(string(name, ".fasta"))
    end
    ### Return assembly size, chromosome names, and lengths
    return(n_int_assembly_size, vec_str_chromosome_names, vec_int_chromosome_lengths)
end
### precompile with tiny test file
temp_N, temp_M, temp_L = fun_split_fasta_count_lengths("test.fasta", 3)

### GC content per line per chromosome and across the entire genome
### Ouputs a file per chromosome where each line correspond to the GC count per line of the fasta file
### Returns the genomewide GC content, minimum and maximum GC contents per line
function fun_estimate_GC_content(vec_str_chromosome_names, n_int_assembly_size)
    n_int_GC_count_genomewide = 0
    ProgressMeter.@showprogress for str_chrom in vec_str_chromosome_names
        # str_chrom = vec_str_chromosome_names[1]
        FILE = open(string(str_chrom, ".fasta"), "r")
        file = open(string(str_chrom, "-GC_content.txt"), "w")
        line = readline(FILE)
        str_previous_base = 'N'
        while !eof(FILE)
            line = readline(FILE)
            n_int_GC_count = length(collect(eachmatch(r"[GC]", line)))
            n_int_total_bases = length(line)
            n_flt_GC_perc = round(n_int_GC_count/n_int_total_bases, digits=4)
            write(file, string(n_flt_GC_perc, '\n'))
            n_int_GC_count_genomewide = n_int_GC_count_genomewide + n_int_GC_count
        end
        close(FILE)
        close(file)
    end
    n_flt_percent_GC_content = n_int_GC_count_genomewide / n_int_assembly_size
    return(n_flt_percent_GC_content)
end
### precompile
temp_GC = fun_estimate_GC_content(temp_M, temp_N)

### Compute the coordinates for an arc polygon
function fun_arcshape(θ1, θ2; r=1.0, w=0.2, n_int_points=50)
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
### precompile
fun_arcshape(0.0, π)

### Plot chromosome lengths layer (outer-most layer)
function fun_plot_chrom_length_layer!(plt, vec_str_chromosome_names, vec_int_chromosome_lengths;
                                r=1.00, w=0.10,
                                vec_colours_chrom=palette(:default),
                                n_int_tick_length_bp=100*1e+6,
                                n_int_tick_label_size=7,
                                n_int_chrom_name_size=7)
    ########################
    ### TEST
    # r=0.75
    # w=0.20
    # vec_colours_chrom=palette(:default)
    # n_int_tick_length_bp=100*1e+6
    # n_int_tick_label_size=7
    # n_int_chrom_name_size=7
    ########################
    ### number of chromosomes
    n = length(vec_str_chromosome_names)
    ### chromosome slice properties
    n_int_assembly_size = sum(vec_int_chromosome_lengths)
    n_flt_empty_slice = 0.05(2π)
    n_flt_available_circumference = 2*π - n_flt_empty_slice
    n_flt_fraction_filled = 0.9*n_flt_available_circumference
    n_flt_fraction_spacer = 0.1*n_flt_available_circumference
    vec_flt_arc_lengths_filled = n_flt_fraction_filled * (vec_int_chromosome_lengths ./ n_int_assembly_size)
    n_flt_arc_lengths_spacer = n_flt_fraction_spacer / (n-1)
    ### Tick arc length
    n_flt_arc_length_per_tick = n_flt_fraction_filled * n_int_tick_length_bp / n_int_assembly_size
    ### Intial position (unrotated)
    θ1 = 2π - (n_flt_empty_slice/2)
    ### plot per chromosome
    ProgressMeter.@showprogress for i in 1:n
        # i = 1
        ### update track tail coordinate
        θ2 = θ1 - vec_flt_arc_lengths_filled[i]
        ### rotate to proper orientation
        if ((θ1 + π/2) > 2π) | ((θ2 + π/2) > 2π)
            θ1 = (θ1+π/2)-2π
            θ2 = (θ2+π/2)-2π
        end
        ### plot each chromosome
        plot!(plt, fun_arcshape(θ1,
                                θ2,
                                r=r,
                                w=w),
            legend=false,
            color=vec_colours_chrom[i],
            linecolor=vec_colours_chrom[i])
        ### print ticks and tick labels
        for j in 0:floor(vec_flt_arc_lengths_filled[i] / n_flt_arc_length_per_tick)
            # j = 1
            x = θ1-(n_flt_arc_length_per_tick*j)
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
                        (string(Int(j)), n_int_tick_label_size, :gray, :center))
        end
        ### print chromosome names
        vec_flt_coor_text = fun_arcshape(θ1-(vec_flt_arc_lengths_filled[i]/2),
                                            θ1-(vec_flt_arc_lengths_filled[i]/2),
                                            r=r+(r/4), w=0.00, n_int_points=1)
        str_chromosome_name = replace(vec_str_chromosome_names[i], "Chromosome"=>"Chr")
        annotate!(plt, vec_flt_coor_text.x[1],
                    vec_flt_coor_text.y[1],
                    (str_chromosome_name, n_int_chrom_name_size, :gray, :center))
        ### update track head coordinate
        θ1 = θ2 - n_flt_arc_lengths_spacer
    end
    ### add legend
    annotate!(plt, -1.5, -1.5, (string("×", n_int_tick_length_bp, " bases"), n_int_chrom_name_size, :gray, :left))
    ### return plot
    return(plt)
end
### precompile
plt = plot(xlim=(-1.5, 1.5), ylim=(-1.5, 1.5), size=(700,700), axis=([], false))
fun_plot_chrom_length_layer!(plt, temp_M, temp_L, n_int_tick_length_bp=150)

### Plot GC content heatmap layer
function fun_plot_GC_layer!(plt, vec_str_chromosome_names, vec_int_chromosome_lengths;
                                r=0.80, w=0.10,
                                n_int_total_chunks_across_genome=1000,
                                vec_colours_GC=palette(:thermometer, 256),
                                n_int_tick_label_size=7)
    ########################
    ### TEST
    # r=0.75
    # w=0.20
    # n_int_total_chunks_across_genome=1000
    # vec_colours_GC=palette(:thermal)
    # n_int_tick_label_size=7
    ########################
    ### number of chromosomes
    n = length(vec_str_chromosome_names)
    ### chromosome slice properties
    n_int_assembly_size = sum(vec_int_chromosome_lengths)
    n_flt_empty_slice = 0.05(2π)
    n_flt_available_circumference = 2*π - n_flt_empty_slice
    n_flt_fraction_filled = 0.9*n_flt_available_circumference
    n_flt_fraction_spacer = 0.1*n_flt_available_circumference
    vec_flt_arc_lengths_filled = n_flt_fraction_filled * (vec_int_chromosome_lengths ./ n_int_assembly_size)
    n_flt_arc_lengths_spacer = n_flt_fraction_spacer / (n-1)
    #################################################################################
    ### Find minimum and maximum GC per chunk and plot legend
    ### Intial position (unrotated)
    θ1 = 2π - (n_flt_empty_slice/2)
    ### Initialise minimum and maximum GC content per chunk
    n_flt_min_GC_perc_genomewide = 1.0
    n_flt_max_GC_perc_genomewide = 0.0
    ### Iterate across chromosomes to find the minumum and maximum GC contents per chunk
    for i in 1:n
        # i = 1
        ### update track tail coordinate
        θ2 = θ1 - vec_flt_arc_lengths_filled[i]
        ### rotate to proper orientation
        if ((θ1 + π/2) > 2π) | ((θ2 + π/2) > 2π)
            θ1 = (θ1+π/2)-2π
            θ2 = (θ2+π/2)-2π
        end
        ### plot GC heatmap with 1 coloured slice per GC content chucnk
        vec_flt_GC_content = parse.(Float64, readlines(string(vec_str_chromosome_names[i], "-GC_content.txt")))
        n_int_bp_per_chunk = Int(ceil(n_int_assembly_size / n_int_total_chunks_across_genome))
        n_int_chunks_in_this_chrom = Int(ceil(vec_int_chromosome_lengths[i] / n_int_bp_per_chunk))
        n_int_lines_per_chunk = Int(ceil(length(vec_flt_GC_content) / n_int_chunks_in_this_chrom))
        ### correct the number of chunks in this chromosome to that we have a maximum of the number of lines == number of chunks
        if (n_int_lines_per_chunk*n_int_chunks_in_this_chrom) > length(vec_flt_GC_content)
            # n_int_chunks_in_this_chrom = length(vec_flt_GC_content)
            n_int_chunks_in_this_chrom = Int(ceil(length(vec_flt_GC_content) / n_int_lines_per_chunk))
        end
        for j in 1:n_int_chunks_in_this_chrom
            # j = 1
            n_int_idx_start = ((j-1)*n_int_lines_per_chunk)+1
            n_int_idx_end = j*n_int_lines_per_chunk
            n_int_idx_start > length(vec_flt_GC_content) ? n_int_idx_start = length(vec_flt_GC_content) : nothing
            n_int_idx_end > length(vec_flt_GC_content) ? n_int_idx_end = length(vec_flt_GC_content) : nothing
            n_flt_mean_GC_content = sum(vec_flt_GC_content[n_int_idx_start:n_int_idx_end]) / length(collect(n_int_idx_start:n_int_idx_end))
            ### Update minimum and maximum GC content per chunk
            n_flt_min_GC_perc_genomewide > n_flt_mean_GC_content ? n_flt_min_GC_perc_genomewide=n_flt_mean_GC_content : nothing
            n_flt_max_GC_perc_genomewide < n_flt_mean_GC_content ? n_flt_max_GC_perc_genomewide=n_flt_mean_GC_content : nothing
        end
        ### update track head coordinate
        θ1 = θ2 - n_flt_arc_lengths_spacer
    end
    n_int_remainder_min = round(n_flt_min_GC_perc_genomewide*100) % 5
    n_flt_min_GC_perc_genomewide = round(n_flt_min_GC_perc_genomewide + ((0 - n_int_remainder_min)/100), digits=2)
    n_int_remainder_max = round(n_flt_max_GC_perc_genomewide*100) % 5
    n_flt_max_GC_perc_genomewide = round(n_flt_max_GC_perc_genomewide + ((5 - n_int_remainder_max)/100), digits=2)
    ### plot legend
    x0 = 0.75; y0 = -1.25
    x1 = 1.50; y1 = -1.35
    n_int_colours = length(vec_colours_GC)
    ### heatmap
    for i in 1:n_int_colours
        dx0 = x0+(x1-x0)*((i-1)/n_int_colours)
        dx1 = x1+(x1-x0)*((i-0)/n_int_colours)
        shp_rectangle = Shape(vcat((dx0, y0),
                                   (dx1, y0),
                                   (dx1, y1),
                                   (dx0, y1)))
        plot!(plt, shp_rectangle, color=vec_colours_GC[i], linecolor=vec_colours_GC[i])
    end
    ### ticks
    vec_flt_legend_ticks = round.(collect(range(n_flt_min_GC_perc_genomewide, n_flt_max_GC_perc_genomewide, length=5)), digits=2)
    for i in 1:5
        dx0 = x0+(x1-x0)*((i-1)/4)
        shp_tick = Shape(vcat((dx0, y1),
                              (dx0, y1+(y1*0.03))))
        plot!(plt, shp_tick, color=:gray, linecolor=:gray)
        annotate!(plt,
                  dx0,
                  y1+(y1*0.05),
                  (vec_flt_legend_ticks[i], n_int_tick_label_size, :gray, :center))
    end
    ### GC content label
    annotate!(plt,
              x0+((x1-x0)/2),
              -1.5,
              ("GC content", n_int_tick_label_size, :gray, :bottom))
    #################################################################################
    ### Intial position (unrotated)
    θ1 = 2π - (n_flt_empty_slice/2)
    ### plot per chromosome
    ProgressMeter.@showprogress for i in 1:n
        # i = 1
        ### update track tail coordinate
        θ2 = θ1 - vec_flt_arc_lengths_filled[i]
        ### rotate to proper orientation
        if ((θ1 + π/2) > 2π) | ((θ2 + π/2) > 2π)
            θ1 = (θ1+π/2)-2π
            θ2 = (θ2+π/2)-2π
        end
        ### plot GC heatmap with 1 coloured slice per GC content chucnk
        vec_flt_GC_content = parse.(Float64, readlines(string(vec_str_chromosome_names[i], "-GC_content.txt")))
        n_int_bp_per_chunk = Int(ceil(n_int_assembly_size / n_int_total_chunks_across_genome))
        n_int_chunks_in_this_chrom = Int(ceil(vec_int_chromosome_lengths[i] / n_int_bp_per_chunk))
        n_int_lines_per_chunk = Int(ceil(length(vec_flt_GC_content) / n_int_chunks_in_this_chrom))
        ### correct the number of chunks in this chromosome to that we have a maximum of the number of lines == number of chunks
        if (n_int_lines_per_chunk*n_int_chunks_in_this_chrom) > length(vec_flt_GC_content)
            # n_int_chunks_in_this_chrom = length(vec_flt_GC_content)
            n_int_chunks_in_this_chrom = Int(ceil(length(vec_flt_GC_content) / n_int_lines_per_chunk))
        end
        for j in 1:n_int_chunks_in_this_chrom
            # j = 1
            n_int_idx_start = ((j-1)*n_int_lines_per_chunk)+1
            n_int_idx_end = j*n_int_lines_per_chunk
            n_int_idx_start > length(vec_flt_GC_content) ? n_int_idx_start = length(vec_flt_GC_content) : nothing
            n_int_idx_end > length(vec_flt_GC_content) ? n_int_idx_end = length(vec_flt_GC_content) : nothing
            n_flt_mean_GC_content = sum(vec_flt_GC_content[n_int_idx_start:n_int_idx_end]) / length(collect(n_int_idx_start:n_int_idx_end))
            n_flt_remapped_GC = (n_flt_mean_GC_content - n_flt_min_GC_perc_genomewide) / (n_flt_max_GC_perc_genomewide - n_flt_min_GC_perc_genomewide)
            col = vec_colours_GC[Int(ceil(n_flt_remapped_GC*(length(vec_colours_GC)-1))+1)] ### less 1 to the multiplier, and added 1 because n_flt_remapped_GC ranges from 0 to 1
            x = θ1 - (abs(θ1-θ2)*j/n_int_chunks_in_this_chrom)
            y = θ1 - (abs(θ1-θ2)*(j-1)/n_int_chunks_in_this_chrom)
            vec_flt_coor_heat_slice = fun_arcshape(x, y, r=r, w=w, n_int_points=10)
            plot!(plt,
                  vec_flt_coor_heat_slice,
                  legend=false,
                  color=col,
                  linecolor=col)
        end
        ### update track head coordinate
        θ1 = θ2 - n_flt_arc_lengths_spacer
    end
    return(plt)
end
### precompile
fun_plot_GC_layer!(plt, temp_M, temp_L, n_int_total_chunks_across_genome=20, vec_colours_GC=palette(:thermometer, 256))

### Plot histogram layer
function fun_plot_histogram_layer!(plt, vec_str_chromosome_names, vec_int_chromosome_lengths, vec_bool_hits;
                                      r=0.60, w=0.10,
                                      n_int_window_size = 10,
                                      col=:black,
                                      col_bckground=:lightgray)
    ########################
    ### TEST
    # vec_bool_hits = Bool.(round.(rand(n_int_chrom_size)))
    # r=0.50
    # w=0.15
    # n_int_window_size = 10
    # col=:black
    # col_bckground=:lightgray
    ########################
    ### number of chromosomes
    n = length(vec_str_chromosome_names)
    ### chromosome slice properties
    n_int_assembly_size = sum(vec_int_chromosome_lengths)
    n_flt_empty_slice = 0.05(2π)
    n_flt_available_circumference = 2*π - n_flt_empty_slice
    n_flt_fraction_filled = 0.9*n_flt_available_circumference
    n_flt_fraction_spacer = 0.1*n_flt_available_circumference
    vec_flt_arc_lengths_filled = n_flt_fraction_filled * (vec_int_chromosome_lengths ./ n_int_assembly_size)
    n_flt_arc_lengths_spacer = n_flt_fraction_spacer / (n-1)
    ### Intial position (unrotated)
    θ1 = 2π - (n_flt_empty_slice/2)
    ### plot per chromosome
    ProgressMeter.@showprogress for i in 1:n
        # i = 1
        ### update track tail coordinate
        θ2 = θ1 - vec_flt_arc_lengths_filled[i]
        ### rotate to proper orientation
        if ((θ1 + π/2) > 2π) | ((θ2 + π/2) > 2π)
            θ1 = (θ1+π/2)-2π
            θ2 = (θ2+π/2)-2π
        end
        ### PENDING: load "vec_bool_hits" per chromosome
        ### Currently we're using just one randomly simulated "vec_bool_hits" 
        vec_bool_hits = vec_bool_hits
        ### measure frequencies per window
        n_int_chrom_size = temp_L[i]
        n_int_windows_count = Int(round(n_int_chrom_size/n_int_window_size))
        vec_flt_frequency = []
        for i in 1:n_int_windows_count
            n_int_start = ((i-1)*n_int_window_size)+1
            n_int_end = ((i-0)*n_int_window_size)+0
            n_int_end > n_int_chrom_size ? n_int_end = n_int_chrom_size : nothing
            append!(vec_flt_frequency, sum(vec_bool_hits[n_int_start:n_int_end])/(n_int_end-n_int_start+1))
        end
        ### draw background as ink to the main histogram
        vec_flt_coor_background_slice = fun_arcshape(θ1, θ2, r=r, w=w)
        plot!(plt, vec_flt_coor_background_slice, legend=false, color=col, linecolor=col)
        ### remove upper part to draw the histograms since we are drawing the arc polygons from top to bottom
        for i in 1:n_int_windows_count
            # i = 1
            x = θ1 - (abs(θ1-θ2)*i/n_int_windows_count)
            y = θ1 - (abs(θ1-θ2)*(i-1)/n_int_windows_count)
            vec_flt_coor_hist_inverse_slice = fun_arcshape(x,
                                                y,
                                                r=r,
                                                w=w*(1-vec_flt_frequency[i]),
                                                n_int_points=10)
            plot!(plt, vec_flt_coor_hist_inverse_slice, legend=false, color=col_bckground, linecolor=col_bckground)
        end
        ### update track head coordinate
        θ1 = θ2 - n_flt_arc_lengths_spacer
    end
    return(plt)
end
### precompile
vec_bool_hits = Bool.(round.(rand(1000)))
fun_plot_histogram_layer!(plt, temp_M, temp_L, vec_bool_hits)





### Precompilation clean-up
for f in readdir()[match.(r"FOR_TESTING-", readdir()) .!= nothing]
    rm(f)
end



###############
### EXECUTE ###
###############
str_filename_fasta = "APGP_CSIRO_Lrig_flye-racon-polca-allhic-juicebox_v0.1n.fasta"
n = 7 ## haploid chromosome number
n_int_tick_length_bp = 100*1e+6
n_int_tick_label_size=7
n_int_chrom_name_size=7

### Split into chromosomes
@time n_int_assembly_size,
      vec_str_chromosome_names,
      vec_int_chromosome_lengths = fun_split_fasta_count_lengths(str_filename_fasta, n)
### Measure GC content
@time n_flt_percent_GC_content = fun_estimate_GC_content(vec_str_chromosome_names, n_int_assembly_size)
### Base plot
plt = plot(xlim=(-1.5, 1.5), ylim=(-1.5, 1.5), size=(700,700), axis=([], false), title="Lolium rigidum genome")
### Layer 1: chromosome lengths
r=1.00; w=0.10
annotate!(plt, 0.0, (r-w/2), ("a", 10, :gray, :center))
fun_plot_chrom_length_layer!(plt, vec_str_chromosome_names, vec_int_chromosome_lengths;
                             r=r, w=w,
                             vec_colours_chrom=palette(:default),
                             n_int_tick_length_bp=n_int_tick_length_bp,
                             n_int_tick_label_size=n_int_tick_label_size,
                             n_int_chrom_name_size=n_int_chrom_name_size)
### Layer 2: GC content
r=0.80; w=0.10
annotate!(plt, 0.0, (r-w/2), ("b", 10, :gray, :center))
fun_plot_GC_layer!(plt, vec_str_chromosome_names, vec_int_chromosome_lengths;
                   r=r, w=w,
                   n_int_total_chunks_across_genome=1000,
                   vec_colours_GC=palette(:thermometer, 25))

savefig(plt, "Lolium_rigidum_genome.svg")

