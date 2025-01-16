#!/bin/bash
set -x
# this script has to be started in $HOME/RNASeq_LW/002_bamtofastq/
OUTDIR="$( pwd )/../004_Trimmoout"
PATHTO_TRIMJAR=/apps/anaconda3/envs/omics2023/share/trimmomatic-0.39-2/trimmomatic.jar
ADAPTORPATH=${PATHTO_TRIMJAR}/../adapters
ADAPTER="TruSeq3-SE.fa"
THREADS="3"
BASE_QUAL=20
WINDOW_LEN=6
WINDOW_QUAL=20
LENGTH_MIN=25

for READS1 in $( ls -1 *fastq )
do
BN_READS1=$( basename ${READS1} .fastq )
bash ./run_trimmo.sh ${PATHTO_TRIMJAR} ${ADAPTORPATH} ${ADAPTER} ${THREADS} ${READS1} ${OUTDIR} ${BASE_QUAL} ${WINDOW_LEN} ${WINDOW_QUAL} ${LENGTH_MIN};
done
