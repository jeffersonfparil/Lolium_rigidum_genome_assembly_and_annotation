# Reference genomes of other species

## Set working directory
```{sh}
DIR=/data/Lolium_rigidum_ASSEMBLY/COMPARATIVE_GENOMICS
# DIR=/data-weedomics-3
```

## Lolium rigidum (our genome assembly)
```{sh}
mkdir ${DIR}/Lolium_rigidum
cd ${DIR}/Lolium_rigidum
ln -s ${DIR}/../GENOME_ASSEMBLY/APGP_CSIRO_Lrig_flye-racon-polca-allhic-juicebox_v0.1n_clean1.fasta APGP_CSIRO_Lrig_flye-racon-polca-allhic-juicebox_v0.1n_clean1.fasta
cd -
```

## Lolium perenne
Frei, Daniel, Elisabeth Veekman, Daniel Grogg, Ingrid Stoffel-Studer, Aki Morishima, Rie Shimizu-Inatsugi, Steven Yates, et al. “Ultralong Oxford Nanopore Reads Enable the Development of a Reference-Grade Perennial Ryegrass Genome Assembly.” Genome Biology and Evolution 13, no. 8 (July 10, 2021): evab159. https://doi.org/10.1093/gbe/evab159.

```{sh}
mkdir ${DIR}/Lolium_perenne
cd ${DIR}/Lolium_perenne    
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/019/359/855/GCA_019359855.1_MPB_Lper_Kyuss_1697/GCA_019359855.1_MPB_Lper_Kyuss_1697_genomic.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/019/359/855/GCA_019359855.1_MPB_Lper_Kyuss_1697/GCA_019359855.1_MPB_Lper_Kyuss_1697_genomic.gbff.gz

gunzip -c GCA_019359855.1_MPB_Lper_Kyuss_1697_genomic.fna.gz > GCA_019359855.1_MPB_Lper_Kyuss_1697_genomic.fasta
gunzip -c GCA_019359855.1_MPB_Lper_Kyuss_1697_genomic.gbff.gz > GCA_019359855.1_MPB_Lper_Kyuss_1697_genomic.gbff
cd -
```

## Secale cereale
Li, Guangwei, Lijian Wang, Jianping Yang, Hang He, Huaibing Jin, Xuming Li, Tianheng Ren, et al. “A High-Quality Genome Assembly Highlights Rye Genomic Characteristics and Agronomically Important Genes.” Nature Genetics 53, no. 4 (April 2021): 574–84. https://doi.org/10.1038/s41588-021-00808-z.

```{sh}
mkdir ${DIR}/Secale_cereale
cd ${DIR}/Secale_cereale    
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/016/097/815/GCA_016097815.1_HAU_Weining_v1.0/GCA_016097815.1_HAU_Weining_v1.0_genomic.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/016/097/815/GCA_016097815.1_HAU_Weining_v1.0/GCA_016097815.1_HAU_Weining_v1.0_genomic.gbff.gz

gunzip -c GCA_016097815.1_HAU_Weining_v1.0_genomic.fna.gz > GCA_016097815.1_HAU_Weining_v1.0_genomic.fasta
gunzip -c GCA_016097815.1_HAU_Weining_v1.0_genomic.gbff.gz > GCA_016097815.1_HAU_Weining_v1.0_genomic.gbff
cd -
```

## Zea mays
Link: https://www.ncbi.nlm.nih.gov/assembly/GCA_016097815.1

```{sh}
mkdir ${DIR}/Zea_mays
cd ${DIR}/Zea_mays  
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/902/167/145/GCA_902167145.1_Zm-B73-REFERENCE-NAM-5.0/GCA_902167145.1_Zm-B73-REFERENCE-NAM-5.0_genomic.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/902/167/145/GCA_902167145.1_Zm-B73-REFERENCE-NAM-5.0/GCA_902167145.1_Zm-B73-REFERENCE-NAM-5.0_genomic.gbff.gz

gunzip -c GCA_902167145.1_Zm-B73-REFERENCE-NAM-5.0_genomic.fna.gz > GCA_902167145.1_Zm-B73-REFERENCE-NAM-5.0_genomic.fasta
gunzip -c GCA_902167145.1_Zm-B73-REFERENCE-NAM-5.0_genomic.gbff.gz > GCA_902167145.1_Zm-B73-REFERENCE-NAM-5.0_genomic.gbff
cd -
```

## Oryza sativa
Link: https://www.ncbi.nlm.nih.gov/assembly/GCF_001433935.1/

```{sh}
mkdir ${DIR}/Oryza_sativa
cd ${DIR}/Oryza_sativa  
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/001/433/935/GCA_001433935.1_IRGSP-1.0/GCA_001433935.1_IRGSP-1.0_genomic.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/001/433/935/GCA_001433935.1_IRGSP-1.0/GCA_001433935.1_IRGSP-1.0_genomic.gbff.gz

gunzip -c GCA_001433935.1_IRGSP-1.0_genomic.fna.gz > GCA_001433935.1_IRGSP-1.0_genomic.fasta
gunzip -c GCA_001433935.1_IRGSP-1.0_genomic.gbff.gz > GCA_001433935.1_IRGSP-1.0_genomic.gbff
cd -
```

## Arabidopsis thaliana
Link: https://www.ncbi.nlm.nih.gov/assembly/GCF_000001735.4

```{sh}
mkdir ${DIR}/Arabidopsis_thaliana
cd ${DIR}/Arabidopsis_thaliana  

wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/735/GCA_000001735.2_TAIR10.1/GCA_000001735.2_TAIR10.1_genomic.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/735/GCA_000001735.2_TAIR10.1/GCA_000001735.2_TAIR10.1_genomic.gbff.gz

gunzip -c GCA_000001735.2_TAIR10.1_genomic.fna.gz > GCA_000001735.2_TAIR10.1_genomic.fasta
gunzip -c GCA_000001735.2_TAIR10.1_genomic.gbff.gz > GCA_000001735.2_TAIR10.1_genomic.gbff
cd -
```

## Marchantia polymorpha
Link: https://www.ncbi.nlm.nih.gov/assembly/GCA_003032435.1#/def

```{sh}
mkdir ${DIR}/Marchantia_polymorpha
cd ${DIR}/Marchantia_polymorpha  

wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/003/032/435/GCA_003032435.1_Marchanta_polymorpha_v1/GCA_003032435.1_Marchanta_polymorpha_v1_genomic.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/003/032/435/GCA_003032435.1_Marchanta_polymorpha_v1/GCA_003032435.1_Marchanta_polymorpha_v1_genomic.gbff.gz

gunzip -c GCA_003032435.1_Marchanta_polymorpha_v1_genomic.fna.gz > GCA_003032435.1_Marchanta_polymorpha_v1_genomic.fasta
gunzip -c GCA_003032435.1_Marchanta_polymorpha_v1_genomic.gbff.gz > GCA_003032435.1_Marchanta_polymorpha_v1_genomic.gbff
cd -
```

