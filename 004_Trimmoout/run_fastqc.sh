#!/bin/bash

set -e # Instructs a shell to exit if a command fails
set -x # Prints out command arguments during execution.

# save this script in $HOME/hw/002_bamToFastq.out_fastq if you process raw reads, 
# and in $HOME/hw/004_trimmomatic.out_trimmedFastq.gz if you process trimmed reads

#READS=BSF_examplefile.fastq
#OUTDIR=$HOME/hw/003_fastqc.out_html.zip # outdir for raw reads
#OUTDIR=$HOME/hw/005_fastqc.out_html.zip # outdir for trimmed reads
#FASTQC_PATH=/apps/anaconda3/envs/omics2023/bin/fastqc

READS=$1
OUTDIR=$2
FASTQC_PATH=$3

# fastqc
mkdir -p ${OUTDIR}
${FASTQC_PATH} -o $OUTDIR/ $READS
