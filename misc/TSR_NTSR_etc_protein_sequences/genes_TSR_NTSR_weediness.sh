DIR=/data/Lolium_rigidum_ASSEMBLY/COMPARATIVE_GENOMICS
mkdir ${DIR}/TSR_NTSR_GENES
cd ${DIR}/TSR_NTSR_GENES




wget https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.dat.gz
wget https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz
gunzip -c uniprot_sprot.dat.gz   > uniprot.txt
gunzip -c uniprot_sprot.fasta.gz > uniprot.faa
rm uniprot_sprot.dat.gz
rm uniprot_sprot.fasta.gz

grep -i "acetyl-coenzyme A carboxylase" uniprot.txt
Mesangiospermae [1437183] acetyl-coenzyme A carboxylase

julia

using ProgressMeter

function SPLIT(x::String)::Vector{String}
    y = split(x, " ")
    y[y .!= ""]
end

function EXTRACT_SEQ(filename, query)
    file = open(filename, "r")
    seekend(file); n=position(file); seekstart(file)
    pb = Progress(n)
    SEQUENCES = []
    while !eof(file)
        line = SPLIT(readline(file))
        update!(pb, position(file))
        if (line[1] == "ID")
            ID = try 
                    line[2]
                catch
                    "Unknown"
                end
            LEN = try
                    line[end-1]
                catch
                    "Unknown"
                end
            if !isa(query, Vector)
                query = [query]
            end
            bool = zeros(Bool, length(query))
            while line[1] != "SQ"
                line = SPLIT(readline(file))
                for i in 1:length(query)
                    bool[i] = bool[i] | (sum(match.(Regex(query[i]), line) .!= nothing) > 0)
                end
            end
            if prod(bool)
                sequence = string(">", ID, ":", LEN, "\n")
                while line[1] != "//"
                    line = SPLIT(readline(file))
                    sequence = string(sequence, join(line))
                end
                SEQUENCES = push!(SEQUENCES, sequence[1:(end-2)])
            end
        end
    end
    close(file)
    return(SEQUENCES)
end

x = EXTRACT_SEQ("test.tmp", "virus")

filename = "uniprot.txt"
query = ["acetyl", "coenzyme A", "carboxylase"]
seq = EXTRACT_SEQ(filename, query)





#################
### CLETHODIM ###
#################
### TARGET: Acetyle Co-A carboxylase
echo 'taxonomy:"Mesangiospermae [1437183]" acetyl-coenzyme A carboxylase'
echo 'LINK: https://www.uniprot.org/uniprot/?query=taxonomy%3A%22Mesangiospermae%20%5B1437183%5D%22%20acetyl-coenzyme%20A%20carboxylase&columns=id%2Centry%20name%2Creviewed%2Cprotein%20names%2Cgenes%2Corganism%2Clength&sort=score'
wget https://www.uniprot.org/uniprot/?query=taxonomy%3A%22Mesangiospermae+%5B1437183%5D%22+acetyl-coenzyme+a+carboxylase+database%3A%28type%3Aensemblplants%29&sort=score#
gunzip -c uniprot_sprot.fasta.gz > ACCase.faa





##################
### GLYPHOSATE ###
##################
### TARGET: 5-enolpyruvylshikimate-3-phosphate synthase
###         Synonyms: 
###             - EPSPS (in rice)
###             - 3-phosphoshikimate 1-carboxyvinyltransferase (in rice and Arabidopsis)
###             - epsp-s (in rice)
### UNIPROT QUERY (https://www.uniprot.org/uniprot/?query=taxonomy%3A%22Mesangiospermae+%5B1437183%5D%22+epsp&sort=score):
echo 'taxonomy:"Mesangiospermae [1437183]" epsp'
### Downloaded all 923 matches, and rename as:
echo 'Glyphosate-target_EPSPS_UniProt_Mesangiospermae.fasta'

wget https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz

#################
### INTERCEPT ### Imazamox, and
################# Imazapyr
### TARGET: Acetolactate synthase a.k.a. acetohydroxy acid synthase
###         Synonyms: 
###             - ALS
###             - AHAS
###             - 
### UNIPROT QUERY (https://www.uniprot.org/uniprot/?query=taxonomy%3A%22Mesangiospermae+%5B1437183%5D%22+als+ahas&sort=score):
echo 'taxonomy:"Mesangiospermae [1437183]" als ahas'
### Downloaded all 983 matches, and rename as:
echo 'Intercept-target_ALS_UniProt_Mesangiospermae.fasta'

###############
### LUXIMAX ###
############### Cinmethylin
### TARGET: acyl-acp thioesterase
###         Synonyms: 
###             - ALT
###             - 
### UNIPROT QUERY (https://www.uniprot.org/uniprot/?query=taxonomy%3A%22Mesangiospermae+%5B1437183%5D%22+acyl-ACP+thioesterase&sort=score):
echo 'taxonomy:"Mesangiospermae [1437183]" acyl-acp thioesterase'
### Downloaded all 1,978 matches, and rename as:
echo 'Luximax-target_ALT_UniProt_Mesangiospermae.fasta'

#################
### OVERWATCH ###
################# Bixlozone
### TARGET: deoxyxylulose 5-phosphate synthase
###         Synonyms: 
###             - DXS
###             - 
### UNIPROT QUERY (https://www.uniprot.org/uniprot/?query=taxonomy%3A%22Mesangiospermae+%5B1437183%5D%22+deoxyxylulose+5-phosphate+synthase&sort=score):
echo 'taxonomy:"Mesangiospermae [1437183]" deoxyxylulose 5-phosphate synthase'
### Downloaded all 295 matches, and rename as:
echo 'Overwatch-target_DXS_UniProt_Mesangiospermae.fasta'

##############
### SAKURA ###
############## Pyroxasulfone
### TARGET: very-long chain fatty acid elongase
###         Synonyms: 
###             - Very-long-chain enoyl-CoA reductase
###             - Elongation of fatty acids protein
###             - 3-ketoacyl-CoA synthase
###             - 
### UNIPROT QUERY (https://www.uniprot.org/uniprot/?query=taxonomy%3A%22Mesangiospermae+%5B1437183%5D%22+very-long+chain+fatty+acid+elongase&sort=score):
echo 'taxonomy:"Mesangiospermae [1437183]" very-long chain fatty acid elongase'
### Downloaded all 295 matches, and rename as:
echo 'Sakura-target_VLCFAE_UniProt_Mesangiospermae.fasta'

#################
### TRIALLATE ###
#################
### TARGET: fatty acid elongase
###         Synonyms: 
###             - 3-ketoacyl-CoA synthase (same proteins as Sakura's: pyroxasulfone)
###             - Elongation of fatty acids protein (same proteins as Sakura's: pyroxasulfone)
###             - 
### UNIPROT QUERY (https://www.uniprot.org/uniprot/?query=taxonomy%3A%22Mesangiospermae+%5B1437183%5D%22+fatty+acid+elongase&sort=score):
echo 'taxonomy:"Mesangiospermae [1437183]" fatty acid elongase'
### Downloaded all 311 matches, and rename as:
echo 'Triallate-target_FAE_UniProt_Mesangiospermae.fasta'

#################
### BOXERGOLD ### S-metolachlor, and
################# Prosulfocarb
### TARGETS: very-long chain fatty acid elongase
###         Synonyms: 
###             - Very-long-chain enoyl-CoA reductase
###             - 3-ketoacyl-CoA synthase (same proteins as Sakura's: pyroxasulfone)
###             - Elongation of fatty acids protein (same proteins as Sakura's: pyroxasulfone)
###             - 
### UNIPROT QUERY (https://www.uniprot.org/uniprot/?query=taxonomy%3A%22Mesangiospermae+%5B1437183%5D%22+very-long+chain+fatty+acid+elongase&sort=score):
echo 'taxonomy:"Mesangiospermae [1437183]" very-long chain fatty acid elongase'
### Downloaded all 295 matches, and rename as:
echo 'BoxerGold-target_VLCFAE_UniProt_Mesangiospermae.fasta'

################
### ATRAZINE ###
################
### TARGET: Photosystem II protein D1
###         Synonyms: 
###             - psbA
###             - 
### UNIPROT QUERY (https://www.uniprot.org/uniprot/?query=taxonomy%3A%22Mesangiospermae+%5B1437183%5D%22+photosystem+2+protein+d1+%22psbA%22&sort=score):
echo 'taxonomy:"Mesangiospermae [1437183]" photosystem 2 protein d1 psba'
### Downloaded all 54,044 matches, and rename as:
echo 'Atrazine-target_psbA_UniProt_Mesangiospermae.fasta'

################
### PARAQUAT ###
################ methylviologen
### TARGET: Photosytem I
###         Notes: 
###             - Paraquat steals electrons from photosynthesis, and diverts them to generate highly damaging reactive oxygen species.
###             - Paraquat generates superoxide from water via the photosynthetic electron transport chain.
###             - Chloroplastic superoxide dismutase (SOD) reduces superoxide into peroxide and oxygen.
###             - Peroxide reacts with other superoxides to generate hydroxyl radicals.
###             - These hydrogen radicals are highly reactive ad damages membrane lipids, DNA, proteins and even carbohydrates.
###             - References: Funderburk and Lawrence, 1964; and Li et al, 2013
###             - Since their is no specific target enzyme it directly inhibits then we will serach for NTSR,
###             - specifically detoxification enzymes:
###                 1.) superoxide dismutase (SOD) | UNIPROT QUERY: 'taxonomy:"Mesangiospermae [1437183]" superoxide dismutase' | 3,915 matches
###                 2.) ascorbate peroxidase (APX) | UNIPROT QUERY: 'taxonomy:"Mesangiospermae [1437183]" ascorbate peroxidase' | 2,823 matches
###                 3.) dehydroascorbate reductase (GST) | UNIPROT QUERY: 'taxonomy:"Mesangiospermae [1437183]" dehydroascorbate reductase' | 224 matches
###                 4.) monodehydroascorbate reductase (MDAR) | UNIPROT QUERY: 'taxonomy:"Mesangiospermae [1437183]" monodehydroascorbate reductase' | 820 matches
###                 5.) glutathione peroxidase | UNIPROT QUERY: 'taxonomy:"Mesangiospermae [1437183]" glutathione peroxidase' | 2,174 matches
###                 6.) L-type amino acid transporter (LAT) | UNIPROT QUERY: 'taxonomy:"Mesangiospermae [1437183]" polyamine transporter rmv1' | 1,658 matches
###                 7.) AtPDR11 (ATP-BINDING CASSETTE G39, ATPDR11, PDR11, PLEIOTROPIC DRUG RESISTANCE 11) | UNIPROT QUERY: 'taxonomy:"Mesangiospermae [1437183]" pdr abc' | 7,460 matches
###                 # tar.xz 8.) Cytochrome P450 | UNIPROT QUERY: 'taxonomy:"Mesangiospermae [1437183]" cytochrome p450 monooxygenase' | 74,868 matches | filename: 'uniprot-taxonomy Mesangiospermae+[1437183] +cytochrome+p450+monooxyge--.tar.xz'
###                 9.) 
echo 'Paraquat-target_NTSR_UniProt_Mesangiospermae.fasta'


