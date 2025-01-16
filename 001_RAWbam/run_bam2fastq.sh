#!/bin/bash

set -e # Instructs a shell to exit if a command fails
set -x # Prints out command arguments during execution.

# if this script is not called by bedtoolsBamToFastqLoopInInpDir.sh, but used for manually process single samples
# this script has to be started in 001_RAWbam. 
# The next 4 lines should be uncommented, while lines where parameters are passed should be commented

#BAM_IN=BSF_filename
#OUTDIR="$( pwd )/../002_bamtofastq/"
#BAMTOFASTQ=/apps/anaconda3/envs/omics2023/bin/bamToFastq
#PE=SE
#BAM_BASENAME=$( basename $BAM_IN .bam )

# this is uncommented when the script is called from bedtoolsBamToFastqLoopInInpDir.sh
BAM_IN=$1
OUTDIR=$2
BAMTOFASTQ=$3
PE=$4
BAM_BASENAME=$( basename $BAM_IN .bam )

mkdir -p ${OUTDIR}

$BAMTOFASTQ -i $BAM_IN -fq ${OUTDIR}/${BAM_BASENAME}_1.fastq

echo "BAM2FASTQ finished." >>$OUTDIR/log.txt
