#!/bin/bash

### Assembly error correction with pilon or quiver/arrow

ASSEMBLY_FASTA=$1

### INPUT:
### (1) Genome assembly in fasta format

### OUTPUTS:
### (1) 

### Pilon
wget https://github.com/broadinstitute/pilon/releases/download/v1.24/pilon-1.24.jar

### Quiver/Arrow
git clone https://github.com/PacificBiosciences/GenomicConsensus.git
