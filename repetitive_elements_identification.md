# Repeat element identification with:
- **RepeatModeler** to generate de novo TE family identification, and
- **RepeatMasker** to identify the coordinates of the identified TEs.

## Install software denpendencies
1. Perl (should be installed by default on linux)
2. Python 3 `sudo apt install python3 python3-pip;`
3. h5py: `pip3 install h5py`
4. TRF (randem repeat finder)
   ```{sh}
   wget https://github.com/Benson-Genomics-Lab/TRF/releases/download/v4.09.1/trf409.linux64
   chmod +x trf409.linux64
   ./trf409.linux64 -h
   sudo mv trf409.linux64 /usr/local/
   ```
5. RMBlast
   ```{sh}
   # NOTE: for some reason libgomp.so.1 is missing so I needed to install inskscape ¯\_(ツ)_/¯ on ubuntu 21.04: 'sudo apt install inkscape'
   wget http://www.repeatmasker.org/rmblast-2.11.0+-x64-linux.tar.gz
   tar -xzvf rmblast-2.11.0+-x64-linux.tar.gz
   sudo mv rmblast-2.11.0/ /usr/local/
   cd /usr/local/rmblast-2.11.0/bin
   ./windowmasker -h
   cd -
   rm rmblast-2.11.0+-x64-linux.tar.gz
   ```
6. RepeatMasker (Note: Do not move to /usr/local/ because we need to place the large Dfam.h5 library which will not fit in the main storage so we're keeing it on the vdb)

   ```{sh}
   wget http://www.repeatmasker.org/RepeatMasker/RepeatMasker-4.1.2-p1.tar.gz
   tar -xvzf RepeatMasker-4.1.2-p1.tar.gz
   # sudo mv RepeatMasker/ /usr/local/
   # cd /usr/local/RepeatMasker/
   cd RepeatMasker/
   ./RepeatMasker -h
   cd -
   ```
7. Replace the RepeatMasker library with Dfam families version 3.5
   ```{sh}
   wget https://www.dfam.org/releases/Dfam_3.5/families/Dfam.h5.gz ### Will take ~3 hours at ~1MB/s
   gunzip Dfam.h5.gz
   # sudo cp /usr/local/RepeatMasker/Libraries/Dfam.h5 /usr/local/RepeatMasker/Libraries/Dfam-3.2.h5 ### backup v3.2
   cp RepeatMasker/Libraries/Dfam.h5 RepeatMasker/Libraries/Dfam-3.2.h5 ### backup v3.2
   cp Dfam.h5 RepeatMasker/Libraries/
   ```
8. Configure RepeatMasker
   ```{sh}
   # cd /usr/local/RepeatMasker
   cd RepeatMasker/
   perl ./configure
   ### Set the proper directories:
   ### /usr/local/trf409.linux64
   ### /usr/local/rmblast-2.11.0/bin
   ### Need to explicitly exit with '5. Done' afterwards.
   cd -
   ```
9. RECON
   ```{sh}
   wget http://www.repeatmasker.org/RepeatModeler/RECON-1.08.tar.gz
   tar -xvzf RECON-1.08.tar.gz
   cd RECON-1.08/src
   make
   make install
   cd ../bin
   ./edgeredef -h
   cd ../..
   sudo mv RECON-1.08/ /usr/local/
   ```
10. RepeatScout
   ```{sh}
   wget http://www.repeatmasker.org/RepeatScout-1.0.6.tar.gz
   tar -xvzf RepeatScout-1.0.6.tar.gz
   cd RepeatScout-1.0.6/
   make
   ./RepeatScout -h
   cd ..
   sudo mv RepeatScout-1.0.6/ /usr/local/
   ```
11. CD-HIT
   ```{sh}
   wget https://github.com/weizhongli/cdhit/releases/download/V4.8.1/cd-hit-v4.8.1-2019-0228.tar.gz
   tar -xzvf cd-hit-v4.8.1-2019-0228.tar.gz
   cd cd-hit-v4.8.1-2019-0228/
   make
   ./cd-hit -h
   cd ..
   sudo mv cd-hit-v4.8.1-2019-0228/ /usr/local/
   ```
12. UCSC Tools
   ```{sh}
   mkdir UCSC_Tools/
   rsync -aP rsync://hgdownload.soe.ucsc.edu/genome/admin/exe/linux.x86_64/ UCSC_Tools/
   sudo mv UCSC_Tools/ /usr/local/
   ```
13. Genome Tools, HMMER, and MAFFT
   ```{sh}
   sudo apt install genometools hmmer mafft
   ```
14. LTR_Retriever
   ```{sh}
   wget https://github.com/oushujun/LTR_retriever/archive/refs/tags/v2.9.0.tar.gz
   tar -xzvf v2.9.0.tar.gz
   cd LTR_retriever-2.9.0/
   sed -i "s/^BLAST/# BLAST/g" paths
   sed -i "s/^Repeat/# Repeat/g" paths
   sed -i "s/^HMMER/# HMMER/g" paths
   sed -i "s/^CDHIT/# CDHIT/g" paths
   echo 'BLAST+=/usr/local/rmblast-2.11.0/bin' >> paths
   echo 'RepeatMasker=/usr/local/RepeatMasker' >> paths
   echo 'HMMER=/usr/bin' >> paths
   echo 'CDHIT=/usr/local/cd-hit-v4.8.1-2019-0228' >> paths
   cd ..
   sudo mv LTR_retriever-2.9.0/ /usr/local/
   ```
15. Ninja Phylogenetic Analysis
   ```{sh}
   sudo apt install default-jre
   wget http://wheelerlab.org/software/ninja/files/ninja.tgz
   tar -xvzf ninja.tgz
   cd ninja_1.2.2/
   ./ninja -h
   cd ..
   sudo mv ninja_1.2.2/ /usr/local/
   sudo mv /usr/local/ninja_1.2.2/ninja /usr/local/ninja_1.2.2/Ninja ### RepeatModeler looks for "Ninja" not "ninja"
   ```
13. RepeatModeler
   ```{sh}
   wget http://www.repeatmasker.org/RepeatModeler/RepeatModeler-2.0.2a.tar.gz
   tar -xvzf RepeatModeler-2.0.2a.tar.gz
   cd RepeatModeler-2.0.2a/
   ```
14. Configure RepeatModeler
   ```{sh}
   cpan JSON File::Which URI LWP::UserAgent Devel::Size
   perl ./configure
   ### /data-weedomics-3/RepeatMasker
   ### /usr/local/RECON-1.08/bin
   ### /usr/local/RepeatScout-1.0.6
   ### /usr/local/trf409.linux64
   ### /usr/local/cd-hit-v4.8.1-2019-0228
   ### /usr/local/UCSC_Tools
   ### /usr/local/rmblast-2.11.0/bin
   ### /usr/local/LTR_retriever-2.9.0
   ### /usr/local/ninja_1.2.2/
   ./RepeatModeler -h
   cd ..
   sudo mv RepeatModeler-2.0.2a/ /usr/local/
   ```

## Build the BLAST database using our genome sequence
```{sh}
DIR=/data-weedomics-3
REF=${DIR}/APGP_CSIRO_Lrig_flye-racon-polca-allhic-juicebox_v0.1n_clean1.fasta
cd $DIR
time \
/usr/local/RepeatModeler-2.0.2a/BuildDatabase \
   -name Lolium_rigidum \
   ${REF}
```

## Run RepeatModeler to identify and classify the repeats
```{sh}
time \
/usr/local/RepeatModeler-2.0.2a/RepeatModeler \
   -database Lolium_rigidum \
   -pa 14 \
   -LTRStruct
```
Ran for about 2216 minutes (~1.5 days) in a 15-threaded machine.

## Find the repeats by generating a masked genome
```{sh}
time \
RepeatMasker/RepeatMasker \
   -lib Lolium_rigidum-families.fa \
   -pa 14 \
   ${REF}
```
- **Summary output**: ${REF}.tbl
- **List of repeats and their locations**: ${REF}.out

## Find Copia and Gypsy LTR (2 out of 3 known LTR families, where the 3rd one BEL/pa family have only been found in animals)
- Copia (order of protein coding domains: *protease*, *integrase*, *reverse transcriptase*, and *ribonuclease H*)

- Gypsy (order of protein coding domains: *protease*, *reverse transcriptase*, *ribonuclease H*, and *integrase*)

- Extract the coordinates of the Copia and Gypsy long termina retrotransposons
   ```{julia}
   using ProgressMeter
   str_filename_stk = "APGP_CSIRO_Lrig_flye-racon-polca-allhic-juicebox_v0.1n_clean1.fasta.out"
   str_filename_output_COPIA = string(str_filename_stk, "-LTR_coordinates-COPIA.csv")
   str_filename_output_GYPSY = string(str_filename_stk, "-LTR_coordinates-GYPSY.csv")
   FILE = open(str_filename_stk, "r")
   file_output_COPIA = open(str_filename_output_COPIA, "w")
   file_output_GYPSY = open(str_filename_output_GYPSY, "w")
   ### Find file size by navigating to the end of the file
   ProgressMeter.seekend(FILE)
   n_int_FILE_size = ProgressMeter.position(FILE)
   ### Reset to the begining of the file to initial the progress bar and the while loop
   ProgressMeter.seekstart(FILE)
   pb = ProgressMeter.Progress(n_int_FILE_size, 1)
   while !eof(FILE)
      line = readline(FILE)
      if (match(Regex("Chromosome"), line) != nothing) & ( (match(Regex("Copia"), line) != nothing) | (match(Regex("Gypsy"), line) != nothing) )
         match(Regex("Copia"), line) != nothing ? str_LTR="Copia" : str_LTR="Gypsy"
         str_chrom, str_start, str_end = split(line)[5:7]
         if str_LTR=="Copia"
            write(file_output_COPIA, string(join([str_chrom, str_start, str_end, str_LTR], ','), '\n'))
         else
            write(file_output_GYPSY, string(join([str_chrom, str_start, str_end, str_LTR], ','), '\n'))
         end
         line = readline(FILE)
      end
      ProgressMeter.update!(pb, ProgressMeter.position(FILE))
   end
   close(FILE)
   close(file_output_COPIA)
   close(file_output_GYPSY)
   ```
