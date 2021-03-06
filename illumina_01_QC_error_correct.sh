#!/bin/bash
#####################################################################
### Quality check and Bayesian error correction of Illumina reads ###
#####################################################################

### Inputs:
### (1) Illumina reads in compressed fastq format (*.fastq.gz)
### Notes: - all reads from 2020.11 were renamed from LOL-WGS2-{2..5}.fastq.gz into simply LOL-WGS-{2..5}.fastq.gz
###        - LOL-WGS-*.fastq.gz from 2020.9 sequencing was split into 2 and renamed to LOL-WGS-{0..1}-*.fastq.gz
###        - LOL-WGS-1* was further dived into 2 because 280GB of RAM is insufficient, i.e. LOL-WGS-1.{0..1}*

### Outputs:
### (1) Quality check html and zip files (*_fastqc.html and *_fastqc.zip)
### (2) Error corrected reads in compressed fastq format (BayesHammer_output/*.cor.fastq.gz)
### (3) Number of bases sequenced after BayesHammer error correction (illumina-corrected.base.count)

### Parameters:
DIR=/data/Lolium_rigidum_ASSEMBLY/assembly_annotation_pipeline_tests_20210104/FASTQ/ILLUMINA
FASTQC=/data/Lolium_rigidum_ASSEMBLY/assembly_annotation_pipeline_tests_20210104/FastQC/fastqc
SPADES=/data/Lolium_rigidum_ASSEMBLY/assembly_annotation_pipeline_tests_20210104/SPAdes-3.14.1-Linux/bin/spades.py

### Navigate to working directory
cd $DIR

### Quality check
time parallel ${FASTQC} {} ::: $(ls *.fastq.gz)
mkdir QC/
mv *.html QC/
mv *.zip QC/

### Bayesian error correction (may need to link python3 via: sudo ln -s /usr/bin/python3 /usr/bin/python)
### (Computation time: 3 days, 6 hours, 14 minutes, and 47.296 seconds on 32 threads and 280 GB of RAM)
time \
for i in 0 1.0 1.1 2 3 4 5
do
# i=1
mkdir BayesHammer_output_WGS-${i}/
time \
$SPADES \
    --only-error-correction \
    --threads 32 \
    --memory 280 \
    -1 LOL-WGS-${i}_combined_R1.fastq.gz \
    -2 LOL-WGS-${i}_combined_R2.fastq.gz \
    -o BayesHammer_output_WGS-${i}/
done

### Clean-up and another QC (21 minutes on 32 threads and 280GB RAM)
mkdir BayesHammer_output/
mv BayesHammer_output_WGS-*/ BayesHammer_output/
cd BayesHammer_output/
mv BayesHammer_output_WGS-*/corrected/*R*.fastq.*.cor.fastq.gz . ### exclude the unpaired fastq.gz file from each read-pair
time parallel ${FASTQC} {} ::: $(ls *.cor.fastq.gz)
mkdir QC/
mv *.html QC/
mv *.zip QC/
cd -

### Count the total number of bases sequenced
echo '#!/bin/bash
f=$1
zcat ${f} | paste - - - - | cut -f 2 | tr -d "\n" | wc -c > ${f}.base.count
' > parallel_count_bases.sh
chmod +x parallel_count_bases.sh
time \
parallel ./parallel_count_bases.sh {} ::: $(ls BayesHammer_output/*.cor.fastq.gz)
cat BayesHammer_output/*.base.count | paste -sd+ | bc > illumina-corrected.base.count
### As of 2021-01-11 we have ~200.6 billion bases (200,617,065,600) sequenced which means ~100X depth of coverage

### Clean-up
rm parallel_count_bases.sh
rm BayesHammer_output/*.base.count
