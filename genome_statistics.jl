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
### Creates fasta files for wach chromosome where the filenames are the Chromosome or contig names with .fasta extension
### Measure the size of each pseudochromosome and contigs
### Retains only the 1:nth largest sequences and removes contigs split file
### Returns assembly size, chromosome names, and chromosome lengths
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
### precompile with tiny test file
temp_N, temp_GC, temp_X, temp_M, temp_L, temp_gc, temp_n = fun_fasta_lengths_GC_N("test.fasta", 3, true)

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
fun_plot_chrom_length_layer!(plt, temp_M, temp_L, n_int_tick_length_bp=10000)

### Plot GC content heatmap layer using the per chromosom GC count files generated by "fun_estimate_GC_content()"
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
fun_plot_GC_layer!(plt, temp_M, temp_L, n_int_total_chunks_across_genome=1000, vec_colours_GC=palette(:thermometer, 256))

### Plot histogram layer
function fun_plot_hits_histogram_layer!(plt, vec_str_chromosome_names, vec_int_chromosome_lengths, str_filename_coor;
                                      r=0.60, w=0.10,
                                      n_int_window_size = 10,
                                      col=:black,
                                      col_background=:lightgray)
    ########################
    ### TEST
    # str_filename_coor = "APGP_CSIRO_Lrig_flye-racon-polca-allhic-juicebox_v0.1n_clean1.fasta.out-LTR_coordinates-COPIA.csv"
    # r=0.50
    # w=0.15
    # n_int_window_size = 10 # 1e6
    # col=:black
    # col_background=:lightgray
    ########################
    ### number of chromosomes
    n = length(vec_str_chromosome_names)
    ### find the maximum hit count across the genome
    n_int_max_hit_count = 0
    for i in 1:n
        ### extraction positions per chromosome
        FILE = open(str_filename_coor, "r")
        vec_int_position = []
        while !eof(FILE)
            str_chrom, str_pos = split(readline(FILE), ',')[[1,2]]
            if str_chrom == vec_str_chromosome_names[i]
                append!(vec_int_position, parse(Int, str_pos))
            end
        end
        close(FILE)
        ### measure frequencies per window per chromosome
        n_int_chrom_size = vec_int_chromosome_lengths[i]
        n_int_windows_count = Int(ceil(n_int_chrom_size/n_int_window_size))
        for j in 1:n_int_windows_count
            # j = 1
            n_int_start = ((j-1)*n_int_window_size)+1
            n_int_end = ((j-0)*n_int_window_size)+0 
            n_int_hit_count = sum((vec_int_position .>= n_int_start) .& (vec_int_position .<= n_int_end))
            n_int_max_hit_count<n_int_hit_count ? n_int_max_hit_count=n_int_hit_count : nothing
        end
    end
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
        ### Extract the end-positions for the current chromosome
        FILE = open(str_filename_coor, "r")
        vec_int_position = []
        while !eof(FILE)
            str_chrom, str_pos = split(readline(FILE), ',')[[1,2]]
            if str_chrom == vec_str_chromosome_names[i]
                append!(vec_int_position, parse(Int, str_pos))
            end
        end
        close(FILE)
        ### measure frequencies per window
        n_int_chrom_size = vec_int_chromosome_lengths[i]
        n_int_windows_count = Int(ceil(n_int_chrom_size/n_int_window_size))
        vec_flt_hits_freq = []
        for j in 1:n_int_windows_count
            # j = 1
            n_int_start = ((j-1)*n_int_window_size)+1
            n_int_end = ((j-0)*n_int_window_size)+0
            append!(vec_flt_hits_freq, sum((vec_int_position .>= n_int_start) .& (vec_int_position .<= n_int_end)) / n_int_max_hit_count)
        end
        ### draw background
        vec_flt_coor_background_slice = fun_arcshape(θ1, θ2, r=r, w=w)
        plot!(plt, vec_flt_coor_background_slice, legend=false, color=col_background, linecolor=col_background)
        ### draw histogram
        for i in 1:n_int_windows_count
            # i = 1
            x = θ1 - (abs(θ1-θ2)*i/n_int_windows_count)
            y = θ1 - (abs(θ1-θ2)*(i-1)/n_int_windows_count)
            vec_flt_coor_hist_inverse_slice = fun_arcshape(x, y,
                                                r=r-(w*(1-vec_flt_hits_freq[i])),
                                                w=w*vec_flt_hits_freq[i],
                                                n_int_points=10)
            plot!(plt, vec_flt_coor_hist_inverse_slice, legend=false, color=col, linecolor=col)
        end
        ### update track head coordinate
        θ1 = θ2 - n_flt_arc_lengths_spacer
    end
    return(plt)
end
### precompile
str_filename_coor = "test-coor.temp"
FILE = open(str_filename_coor, "w")
### random segment setting up sa features
vec_int_n_features = [10, 20, 100]
for i in 1:length(temp_M)
    # i = 1
    str_chr = temp_M[i]
    int_len = temp_L[i]
    int_features = vec_int_n_features[i]
    vec_int_end1 = sort(rand(collect(1:int_len), int_features))
    vec_int_end2 = sort(rand(collect(1:int_len), int_features))
    for j in 1:length(vec_int_end1)
        if vec_int_end1[j] < vec_int_end2[j]
            write(FILE, string(join([str_chr, string(vec_int_end1[j]), string(vec_int_end2[j]), "DUMMY"], ','), '\n'))
        else
        end
    end
end
close(FILE)
fun_plot_hits_histogram_layer!(plt, temp_M, temp_L, str_filename_coor)
### Precompilation clean-up
for f in readdir()[match.(r"FOR_TESTING", readdir()) .!= nothing]
    rm(f)
end

### add comma separators on large positive integers
function fun_add_comma_separator_on_large_positive_integers(int_positive_number)
    if match(Regex(","), string(int_positive_number)) != nothing
        str_number = int_positive_number
    else
        str_number = ""
        vec_char_number = reverse(collect(string(int_positive_number)))
        i = 0
        while i < length(vec_char_number)
            i += 1
            char_digit = vec_char_number[i]
            str_number = string(char_digit, str_number)
            int_counter = 1
            while (int_counter < 3) & (i < length(vec_char_number))
                int_counter += 1
                i += 1
                char_digit = vec_char_number[i]
                str_number = string(char_digit, str_number)
            end
            str_number = string(",", str_number)
        end
        str_number = str_number[2:end]
    end
    return(str_number)
end









###############
### EXECUTE ###
###############
str_filename_fasta = "APGP_CSIRO_Lrig_flye-racon-polca-allhic-juicebox_v0.1n_clean1.fasta"
str_filename_LTR_COPIA = "APGP_CSIRO_Lrig_flye-racon-polca-allhic-juicebox_v0.1n_clean1.fasta.out-LTR_coordinates-COPIA.csv"
str_filename_LTR_GYPSY = "APGP_CSIRO_Lrig_flye-racon-polca-allhic-juicebox_v0.1n_clean1.fasta.out-LTR_coordinates-GYPSY.csv"
n = 7 ## haploid chromosome number
n_int_tick_length_bp = 100*1e+6
n_int_tick_label_size=7
n_int_chrom_name_size=7

### Split into chromosomes, count assembly size, GC content, Ns, and GC fractiopn per line
@time n_int_assembly_size,
      n_int_assembly_GC,
      n_int_assembly_N,
      vec_str_chromosome_names,
      vec_int_chromosome_lengths,
      vec_int_chromosome_GC,
      vec_int_chromosome_N = fun_fasta_lengths_GC_N(str_filename_fasta, 0, true)

### Some assembly stats including the pseudo-chromosomes and small contigs
vec_int_idx_sort_by_decreasing_size = sortperm(vec_int_chromosome_lengths, rev=true)
vec_str_chromosome_names = vec_str_chromosome_names[vec_int_idx_sort_by_decreasing_size]
vec_int_chromosome_lengths = vec_int_chromosome_lengths[vec_int_idx_sort_by_decreasing_size]
vec_int_cumsum_length = reverse(cumsum(reverse(vec_int_chromosome_lengths)))
vec_int_cumsum_less_median = abs.(vec_int_cumsum_length .- n_int_assembly_size/2)
n_int_size_of_largest_chromosome = maximum(vec_int_chromosome_lengths)
str_largest_chromosome = vec_str_chromosome_names[vec_int_chromosome_lengths.==maximum(vec_int_chromosome_lengths)]
int_n_chromosomes_whole_assembly = length(vec_str_chromosome_names)
int_size_whole_assembly = sum(vec_int_chromosome_lengths)
L50_whole_assembly = collect(1:length(vec_str_chromosome_names))[vec_int_cumsum_less_median .== minimum(vec_int_cumsum_less_median)][1]
N50_whole_assembly = vec_int_chromosome_lengths[vec_int_cumsum_less_median .== minimum(vec_int_cumsum_less_median)][1]

### Using only the 7 pseudo-chromosomes, i.e. 7 largest sequences
vec_str_chromosome_names = vec_str_chromosome_names[1:n]
vec_int_chromosome_lengths = vec_int_chromosome_lengths[1:n]
vec_int_chromosome_GC = vec_int_chromosome_GC[1:n]
vec_int_chromosome_N = vec_int_chromosome_N[1:n]

### Some assembly stats after excluding the small contigs
vec_int_idx_sort_by_decreasing_size = sortperm(vec_int_chromosome_lengths, rev=true)
vec_str_chromosome_names = vec_str_chromosome_names[vec_int_idx_sort_by_decreasing_size]
vec_int_chromosome_lengths = vec_int_chromosome_lengths[vec_int_idx_sort_by_decreasing_size]
vec_int_cumsum_length = reverse(cumsum(reverse(vec_int_chromosome_lengths)))
vec_int_cumsum_less_median = abs.(vec_int_cumsum_length .- n_int_assembly_size/2)
n_int_size_of_largest_chromosome = maximum(vec_int_chromosome_lengths)
str_largest_chromosome = vec_str_chromosome_names[vec_int_chromosome_lengths.==maximum(vec_int_chromosome_lengths)]
int_n_chromosomes = n
int_size = sum(vec_int_chromosome_lengths)
L50 = collect(1:length(vec_str_chromosome_names))[vec_int_cumsum_less_median .== minimum(vec_int_cumsum_less_median)][1]
N50 = vec_int_chromosome_lengths[vec_int_cumsum_less_median .== minimum(vec_int_cumsum_less_median)][1]

### Sort back according to chromosome names
vec_int_idx_sort_by_name = sortperm(vec_str_chromosome_names)
vec_str_chromosome_names = vec_str_chromosome_names[vec_int_idx_sort_by_name]
vec_int_chromosome_lengths = vec_int_chromosome_lengths[vec_int_idx_sort_by_name]

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
int_n_chromosomes_whole_assembly = fun_add_comma_separator_on_large_positive_integers(int_n_chromosomes_whole_assembly)
int_size_whole_assembly = fun_add_comma_separator_on_large_positive_integers(int_size_whole_assembly)
int_n_chromosomes = fun_add_comma_separator_on_large_positive_integers(int_n_chromosomes)
int_size = fun_add_comma_separator_on_large_positive_integers(int_size)
L50 = fun_add_comma_separator_on_large_positive_integers(L50)
N50 = fun_add_comma_separator_on_large_positive_integers(N50)
annotate!(plt, -1.5, 1.25,(string("Whole assembly:\n  n=", int_n_chromosomes_whole_assembly,
                                  "\n  size=", int_size_whole_assembly, " bp",
                                  "\n",
                                  "\nPseudo-chromosomes:\n  n=", int_n_chromosomes,
                                  "\n  size=", int_size, " bp",
                                  "\n  L50=", L50, " chromosomes",
                                  "\n  N50=", N50, " bp"), 8, :gray, :left))
### Layer 2: GC content
r=0.80; w=0.10
annotate!(plt, 0.0, (r-w/2), ("b", 10, :gray, :center))
fun_plot_GC_layer!(plt, vec_str_chromosome_names, vec_int_chromosome_lengths;
                   r=r, w=w,
                   n_int_total_chunks_across_genome=1000,
                   vec_colours_GC=palette(:thermometer, 7))
### Layer 3: Ty1-Copia LTR histogram
r=0.60; w=0.10
annotate!(plt, 0.0, (r-w/2), ("c", 10, :gray, :center))
fun_plot_hits_histogram_layer!(plt, vec_str_chromosome_names, vec_int_chromosome_lengths,
                               str_filename_LTR_COPIA,
                               r=r, w=w,
                               n_int_window_size = 1e6,
                               col=:black,
                               col_background=:lightgray)
### Layer 3: Ty1-Gypsy LTR histogram
r=0.40; w=0.10
annotate!(plt, 0.0, (r-w/2), ("d", 10, :gray, :center))
fun_plot_hits_histogram_layer!(plt, vec_str_chromosome_names, vec_int_chromosome_lengths,
                               str_filename_LTR_GYPSY,
                               r=r, w=w,
                               n_int_window_size = 1e6,
                               col=:gray,
                               col_background=:lightgray)
### Save as svg
savefig(plt, "Lolium_rigidum_genome.svg")

