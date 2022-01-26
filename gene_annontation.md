# Gene annontation with:
- **ProtHint** to generate gff3 file
- Transciptome
- ???

## Install software dependencies

0. Set our working directory, and define our genome assembly, as well as the transcriptome assembly
```{sh}
DIR=/data-weedomics-3
REF=${DIR}/APGP_CSIRO_Lrig_flye-racon-polca-allhic-juicebox_v0.1n.fasta
TRA=${DIR}/LolRig_transcripts_fpkm_1.fa ### Pooled all tissues
# TRA=${DIR}/transcripts_fpkm_1.fa ### Tissue-specific
cd $DIR
```

1. Perl modules for GeneMark-EX and ProtHint
    ```{sh}
    ### NOTE: may require sudo
    cpan Hash::Merge MCE::Mutex Math::Utils Parallel:ForkManager ### GeneMark-EX dependencies
    cpan threads YAML Thread::Queue ### ProtHint dependencies
    ```
2. GeneMark-EX
Download **GeneMark-ES/ET/EP** manually from (http://exon.gatech.edu/GeneMark/license_download.cgi)[http://exon.gatech.edu/GeneMark/license_download.cgi]. Enter the credentials being required. You will need to download the software and its corresponding key.
    ```{sh}
    tar -xvzf gmes_linux_64.tar.gz ### decompress the software
    gunzip -d gm_key_64.gz; mv gm_key_64 gmes_linux_64/.gm_key ### decompress, rename, set as hidden, and move to the GeneMark-EX directory
    cd gmes_linux_64/
    ./check_install.bash ### check installation of GeneMark-EX
    cd -
    ```
3. Download, install, and configure ProtHint
    ```{sh}
    git clone https://github.com/gatech-genemark/ProtHint.git
    cd ProtHint/
    echo "export PROTHINT_PATH=${DIR}/ProtHint/bin/"  >> ~/.bashrc ### add to path
    source ~/.bashrc
    bin/prothint.py -h
    cd -
    cp -R gmes_linux_64/* ProtHint/dependencies/GeneMarkES/
    cp gmes_linux_64/.gm_key ProtHint/dependencies/GeneMarkES/
    cd -
    ```
4. Install Star transcriptome aligner
    ```{sh}
    sudo apt install -y rna-star
    ```

## Download OrthoDB protein sequencies and gene list
```{sh}
wget https://v101.orthodb.org/download/odb10v1_all_fasta.tab.gz ### ~22 minutes at ~7MB/s
wget https://v101.orthodb.org/download/odb10v1_genes.tab.gz ### <1 minute
gunzip -d odb10v1_genes.tab.gz
gunzip -d odb10v1_all_fasta.tab.gz
mv odb10v1_all_fasta.tab odb10v1_all.fasta
```

## *Ab initio* gene annotation
```{sh}
time \
ProtHint/bin/prothint.py \
    ${REF} \
    odb10v1_all.fasta
```


## Transcript-supported gene annontation
1. Align the transcriptome assembly to the genome assembly
    ```{sh}
    time \
    STAR --runMode genomeGenerate \
        --genomeDir $(dirname ${REF}) \
        --genomeFastaFiles ${REF} \
        --genomeSAindexNbases 13 \
        --runThreadN 31
    ```
2. Align
    ```{sh}
    DIR_RAW_RNASEQ=/data/Lolium_rigidum_ASSEMBLY/TRANSCRIPTOME_ASSEMBLY/raw_reads
    time \
    for tissue in INFLO LEAF MERI ROOT SEEDL STEM
    do
        for rep in 1 2
        do
            echo ${tissue}-${rep}
            STAR --genomeDir $(dirname ${REF}) \
                 --readFilesIn \
                    ${DIR_RAW_RNASEQ}/${tissue}-${rep}_combined_R1.fastq \
                    ${DIR_RAW_RNASEQ}/${tissue}-${rep}_combined_R2.fastq \
                 --runThreadN 27 \
                 --outFileNamePrefix Lolium_rigidum-transcriptome-${tissue}-${organ}-UNSORTED
        done
    done

time \
${STAR} --genomeDir /data/Lolium_rigidum_ASSEMBLY/assembly_annotation_pipeline_tests_20210104/ASSEMBLY/ \
        --readFilesIn \
            /data/Lolium_rigidum_ASSEMBLY/assembly_annotation_pipeline_tests_20210104/FASTQ/ILLUMINA/RNAseq/INFLO-1_combined_R1.fastq.gz \
            /data/Lolium_rigidum_ASSEMBLY/assembly_annotation_pipeline_tests_20210104/FASTQ/ILLUMINA/RNAseq/INFLO-1_combined_R2.fastq.gz \
        --readFilesCommand zcat \
        --runThreadN 12 \
        --outFileNamePrefix /data/Lolium_rigidum_ASSEMBLY/assembly_annotation_pipeline_tests_20210104/ASSEMBLY/Lori_hh_RNAseq


    ```

3. Sort, compress, and index
    ```{sh}
    time \
    samtools view \
        -q 20 \
        -b Lori_hh_RNAseqAligned.out.sam | \
    samtools sort > Lolium_rigidum_transcriptome_all_tissues.bam
    samtools index Lolium_rigidum_transcriptome_all_tissues.bam
    ```

4. bam to gff
    ```{sh}
    wget https://metacpan.org/raw/TJPARNELL/Bio-ToolBox-1.17/scripts/bam2gff_bed.pl?download=1
    mv 'bam2gff_bed.pl?download=1' bam2gff_bed.pl
    sudo cpanm Bio::ToolBox
    time \
    perl bam2gff_bed.pl \
        --in Lolium_rigidum_transcriptome_all_tissues.bam \
        --pe \
        --gff \
        --source RNAseq \
        --out Lolium_rigidum_transcriptome_all_tissues
    ```

5. GeneMark-EP+ (generate gtf annotations)
    ```{sh}
    GENEMARK_EPP=/data/Lolium_rigidum_ASSEMBLY/assembly_annotation_pipeline_tests_20210104/gmes_linux_64/gmes_petap.pl
    time \
    ${GENEMARK_EPP} \
        --EP prothint.gff \
        --evidence Lolium_rigidum_transcriptome_all_tissues.gff \
        --seq ${REF} \
        --soft_mask 1000 \
        --cores 12 \
        --verbose
    ```

