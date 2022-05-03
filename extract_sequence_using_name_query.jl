using ProgressMeter

### BE SURE TO REPLACE "|" in the sequence_name_query input

fasta_input = ARGS[1]
sequence_name_query = ARGS[2]
fasta_output = try 
                    ARGS[3]
                catch
                    ""
                end
new_sequence_name = try 
                    ARGS[4]
                catch
                    ""
                end
add_gene_coordinates = try 
                    parse(Bool, ARGS[5])
                catch
                    true
                end

if fasta_output == ""
    fasta_output = string(join(split(fasta_input, ".")[1:(end-1)], "."), "-", sequence_name_query, ".out")
end

### Add escape characters in front of "."
if match(Regex("\\."), sequence_name_query) != nothing
    sequence_name_query = replace(sequence_name_query, "."=> "\\.")
end
### Remove return character "\r"
if match(Regex("\\r"), sequence_name_query) != nothing
    sequence_name_query = replace(sequence_name_query, "\r"=> "")
end

file_input = open(fasta_input, "r")
seekend(file_input); n = position(file_input)
seekstart(file_input)
pb = Progress(n)
while !eof(file_input)
    line = readline(file_input)
    if line[1] == '>'
        while match(Regex(sequence_name_query), replace(line, "|"=>":")) != nothing
            file_output = open(fasta_output, "a")
            vec_line = split(line, " ")
            if (new_sequence_name != "")
                line = string(">", new_sequence_name)
            end
            if add_gene_coordinates
                coordinates = try
                        vec_line[(match.(Regex("interval="), vec_line) .!= nothing) .| (match.(Regex("location="), vec_line) .!= nothing)][1]
                    catch
                        ""
                    end
                line = string(line, "(", coordinates, ")")
            end
            write(file_output, string(line, '\n'))
            line = readline(file_input)
            bool_test = line[1] != '>'
            while bool_test
                write(file_output, line)
                line = readline(file_input)
                bool_test = try
                    line[1] != '>'
                catch
                    false
                end
                update!(pb, position(file_input))
            end
            write(file_output, '\n')
            close(file_output)
        end
    end
    update!(pb, position(file_input))
end
close(file_input)
