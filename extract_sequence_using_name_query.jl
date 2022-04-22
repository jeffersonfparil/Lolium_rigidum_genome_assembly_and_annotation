using ProgressMeter
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
@show fasta_output
if match(Regex("\\."), sequence_name_query) != nothing
    sequence_name_query = replace(sequence_name_query, "."=> "\\.")
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
            if (new_sequence_name != "")
                vec_line = split(line, " ")
                line = string(">", new_sequence_name)
            end
            write(file_output, string(line, '\n'))
            line = readline(file_input)
            while line[1] != '>'
                write(file_output, line)
                line = readline(file_input)
                update!(pb, position(file_input))
            end
            write(file_output, '\n')
            close(file_output)
        end
    end
    update!(pb, position(file_input))
end
close(file_input)
