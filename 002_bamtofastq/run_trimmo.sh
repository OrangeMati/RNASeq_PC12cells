#!/bin/bash

set -e # Instructs a shell to exit if a command fails
set -x # Prints out command arguments during execution.

#PATHTO_TRIMJAR=/apps/anaconda3/envs/omics2023/share/trimmomatic-0.39-2/trimmomatic.jar
#ADAPTORPATH=${PATHTO_TRIMJAR}/../adapters
#ADAPTER="TruSeq2-SE.fa"
#THREADS="16"
#READS1=BSF_examplefile.fastq
#BN_READS1=$( basename ${READS1} .fastq )
#OUTDIR="$( pwd )/../005_trimmomatic.out_trimmedFastq.gz"
#BASE_QUAL=20
#WINDOW_LEN=6
#WINDOW_QUAL=20
#LENGTH_MIN=25

PATHTO_TRIMJAR=$1
ADAPTORPATH=$2
ADAPTER=$3
THREADS=$4
READS1=$5 #BSF_examplefile.fastq
BN_READS1=$( basename ${READS1} .fastq )
OUTDIR=$6
BASE_QUAL=$7
WINDOW_LEN=$8
WINDOW_QUAL=$9
LENGTH_MIN=${10}

java -jar ${PATHTO_TRIMJAR} SE \
    -threads $THREADS -phred33 ${READS1} \
    ${OUTDIR}/${BN_READS1}_SE.fastq.gz \
    ILLUMINACLIP:${ADAPTORPATH}/${ADAPTER}:2:30:10:6:true \
    TRAILING:$BASE_QUAL \
    SLIDINGWINDOW:$WINDOW_LEN:$WINDOW_QUAL \
    MINLEN:$LENGTH_MIN
