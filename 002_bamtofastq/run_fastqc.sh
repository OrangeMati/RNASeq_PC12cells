#!/bin/bash

set -e # Instructs a shell to exit if a command fails
set -x # Prints out command arguments during execution.


#READS=BSF_examplefile.fastq
#OUTDIR=outdir for raw reads
#OUTDIR=outdir for trimmed reads
#FASTQC_PATH=/apps/anaconda3/envs/omics2023/bin/fastqc

READS=$1
OUTDIR=$2
FASTQC_PATH=$3

# fastqc
mkdir -p ${OUTDIR}
${FASTQC_PATH} -o $OUTDIR/ $READS
