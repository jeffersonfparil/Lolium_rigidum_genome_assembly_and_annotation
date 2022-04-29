filename_input = ARGS[1]
filename_output = string(join(split(filename_input, '.')[1:(end-1)], '.'), ".phylip")
# count number of sequences
function COUNT_SEQUENCES(filename_input)
    FILE = open(filename_input, "r")
    count_sequences = 0
    while !eof(FILE)
        if readline(FILE)[1] == '>'
        count_sequences += 1
        end
    end
    close(FILE)
    return(count_sequences)
end
count_sequences = COUNT_SEQUENCES(filename_input)
# count sequence length
function SEQUENCE_LENGTH(filename_input)
    FILE = open(filename_input, "r")
    _ = readline(FILE)
    line = readline(FILE)
    sequence_length = length(collect(line))
    while line[1] != '>'
        line = readline(FILE)
        sequence_length = sequence_length + length(collect(line))
    end
    sequence_length = sequence_length - length(collect(line))
    close(FILE)
    return(sequence_length)
end
sequence_length = SEQUENCE_LENGTH(filename_input)
# output phylip format
function CONVERT_FASTA_TO_PHYLIP(filename_input, filename_output, count_sequences, sequence_length)
    FILE_INPUT = open(filename_input, "r")
    FILE_OUTPUT = open(filename_output, "a")
    write(FILE_OUTPUT, string(count_sequences, " ", sequence_length, '\n'))
    while !eof(FILE_INPUT)
        line = readline(FILE_INPUT)
        if line[1] == '>'
            line = line[2:end]
        end
        write(FILE_OUTPUT, string(line, '\n'))
    end
    close(FILE_INPUT)
    close(FILE_OUTPUT)
    return(0)
end
CONVERT_FASTA_TO_PHYLIP(filename_input, filename_output, count_sequences, sequence_length)