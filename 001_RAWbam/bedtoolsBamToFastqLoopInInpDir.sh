#!/bin/bash
set -x
# this script has to be started in $HOME/RNASeq_LW/001_RAWbam
OUTDIR="$( pwd )/../002_bamtofastq/"
BAMTOFASTQ=/apps/anaconda3/envs/omics2023/bin/bamToFastq
for BAM_IN in $( ls -1 *bam )
do
bash ./run_bam2fastq.sh ${BAM_IN} ${OUTDIR} ${BAMTOFASTQ}
done
