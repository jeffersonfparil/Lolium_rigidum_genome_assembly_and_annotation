#!/bin/bash
################################################
### Assembly using Illumina reads via SPAdes ###
################################################

### Inputs:
### (1) Error-corrected Illumina reads in compressed fastq format (BayesHammer_output/*cor.fastq.gz)

### Outputs:
### (1) scaffolds.fasta - contains resulting scaffolds (recommended for use as resulting sequences)
### (2) contigs.fasta - contains resulting contigs
### (3) assembly_graph.gfa - contains SPAdes assembly graph and scaffolds paths in GFA 1.0 format
### (4) assembly_graph.fastg - contains SPAdes assembly graph in FASTG format


### Parameters:
INPUT_DIR=/data/Lolium_rigidum_ASSEMBLY/assembly_annotation_pipeline_tests_20210104/FASTQ/ILLUMINA/BayesHammer_output/
SPADES=/data/Lolium_rigidum_ASSEMBLY/assembly_annotation_pipeline_tests_20210104/SPAdes-3.14.1-Linux/bin/spades.py
OUTPUT_DIR=/data/Lolium_rigidum_ASSEMBLY/assembly_annotation_pipeline_tests_20210104/ASSEMBLY/

### Intialise and spades output folder
cd ${OUTPUT_DIR}
mkdir SPADES_Lori_i1/

### ERROR: 32 threads and 280 GB RAM
### ERROR STILL: Assemble with 16 threads and 280 RAM (number of threads reduced to avoid insufficient memory allocation per thread)
### ERROR STILL: 32 threads and 200 GB RAM and with the --isolate flag
### ERROR STILL: 32 threads, 280 GB RAM, and specifying k-mer length to 33 since no errors were found at this point
### ERROR STILL: 32 threads, 280 GB RAM, specifying k-mer length to 33, and removing --isolate flag
### ERROR STILL: 32 threads, 280 GB RAM, specifying k-mer length to 33, --isolate flag, but only 1 library
### ERROR STILL: 32 threads, 280 GB RAM, specifying k-mer length to 33, --isolate flag, and concatenated all libraries into a single file
### ERROR STILL: 32 threads, 280 GB RAM, --isolate flag, concatenated all libraries into a single file, and reduced k-mer length from 33 to 21
### TRYING: 32 threads, 250 GB RAM, --isolate flag, concatenated all libraries into a single file, and using default k range
# time \
# $SPADES \
#     --only-assembler \
#     --isolate \
#     --threads 32 \
#     --memory 280 \
#     -k 33 \
#     --pe1-1 ${INPUT_DIR}/LOL-WGS-0_combined_R1.fastq.00.0_0.cor.fastq.gz --pe1-2 ${INPUT_DIR}/LOL-WGS-0_combined_R2.fastq.00.0_0.cor.fastq.gz \
#     --pe2-1 ${INPUT_DIR}/LOL-WGS-1.0_combined_R1.fastq.00.0_0.cor.fastq.gz --pe2-2 ${INPUT_DIR}/LOL-WGS-1.0_combined_R2.fastq.00.0_0.cor.fastq.gz \
#     --pe2-1 ${INPUT_DIR}/LOL-WGS-1.1_combined_R1.fastq.00.0_0.cor.fastq.gz --pe2-2 ${INPUT_DIR}/LOL-WGS-1.1_combined_R2.fastq.00.0_0.cor.fastq.gz \
#     --pe3-1 ${INPUT_DIR}/LOL-WGS-2_combined_R1.fastq.00.0_0.cor.fastq.gz --pe3-2 ${INPUT_DIR}/LOL-WGS-2_combined_R2.fastq.00.0_0.cor.fastq.gz \
#     --pe4-1 ${INPUT_DIR}/LOL-WGS-3_combined_R1.fastq.00.0_0.cor.fastq.gz --pe4-2 ${INPUT_DIR}/LOL-WGS-3_combined_R2.fastq.00.0_0.cor.fastq.gz \
#     --pe5-1 ${INPUT_DIR}/LOL-WGS-4_combined_R1.fastq.00.0_0.cor.fastq.gz --pe5-2 ${INPUT_DIR}/LOL-WGS-4_combined_R2.fastq.00.0_0.cor.fastq.gz \
#     --pe6-1 ${INPUT_DIR}/LOL-WGS-5_combined_R1.fastq.00.0_0.cor.fastq.gz --pe6-2 ${INPUT_DIR}/LOL-WGS-5_combined_R2.fastq.00.0_0.cor.fastq.gz \
#     -o ${OUTPUT_DIR}/SPADES_Lori_i1/
cat ${INPUT_DIR}/LOL-WGS-*_R1.fastq.*.cor.fastq.gz > ${INPUT_DIR}/Lrigidum_illumina_150bp_R1.fastq.gz
cat ${INPUT_DIR}/LOL-WGS-*_R2.fastq.*.cor.fastq.gz > ${INPUT_DIR}/Lrigidum_illumina_150bp_R2.fastq.gz
time \
$SPADES \
    --only-assembler \
    --isolate \
    --threads 32 \
    --memory 250 \
    --pe1-1 ${INPUT_DIR}/Lrigidum_illumina_150bp_R1.fastq.gz \
    --pe1-2 ${INPUT_DIR}/Lrigidum_illumina_150bp_R2.fastq.gz \
    -o ${OUTPUT_DIR}/SPADES_Lori_i1/

### TEST: Run on spartan physical partition with 500 GB RAM
### concatenate the BayesHammer-corrected reads into a single fastq.gz read pair
cat *_R1.fastq.*.cor.fastq.gz > Lrigidum_illumina_150bp_R1.fastq.gz
cat *_R2.fastq.*.cor.fastq.gz > Lrigidum_illumina_150bp_R2.fastq.gz
### prepare output directory
mkdir OUTPUT/
### write slurm script
echo '#!/bin/bash
# Partition for the job:
#SBATCH --partition=physical
# Multithreaded (SMP) job: must run on one node and the cloud partition
#SBATCH --nodes=1
# The name of the job:
#SBATCH --job-name=Lolium_rigidum_genome_assembly
# The project ID which this job should run under:
#SBATCH --account=punim0543
# Maximum number of tasks/CPU cores used by the job:
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
# The amount of memory in megabytes per process in the job:
#SBATCH --mem=550GB
# The maximum running time of the job in days-hours:mins:sec
#SBATCH --time=5-0:0:0
# Send yourself an email when the job:
# aborts abnormally (fails)
#SBATCH --mail-type=FAIL
# begins
#SBATCH --mail-type=BEGIN
# ends successfully
#SBATCH --mail-type=END
# Use this email address:
#SBATCH --mail-user=jparil@student.unimelb.edu.au
### load modules
module load spades/3.14.0-python-3.7.4
### Execute
cd /data/gpfs/projects/punim0543/jparil/GENOME_ASSEMBLY_Lolium_rigidum/
time \
SPAdes-3.14.1-Linux/bin/spades.py \
    --only-assembler \
    --isolate \
    --threads 32 \
    --memory 500 \
    -k 21,33
    --pe1-1 Lrigidum_illumina_150bp_R1.fastq.gz --pe1-2 Lrigidum_illumina_150bp_R2.fastq.gz \
    -o OUTPUT/
' > Lrigidum_gassembly.slurm
### submit the job and execute
sbatch Lrigidum_gassembly.slurm
### monitoring
cd /data/gpfs/projects/punim0543/jparil/GENOME_ASSEMBLY_Lolium_rigidum/
squeue -u jparil
check_project_usage
