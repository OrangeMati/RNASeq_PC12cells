#!/bin/bash
set -x
# this script has to be started in $HOME/RNASeq_LW/004_Trimmoout
OUTDIR="$( pwd )/../005_FASTQCout2/"
FASTQC_PATH=/apps/anaconda3/envs/omics2023/bin/fastqc
for FASTQ_IN in $( ls -1 *fastq.gz )
do
bash ./run_fastqc.sh ${FASTQ_IN} ${OUTDIR} ${FASTQC_PATH}
done
