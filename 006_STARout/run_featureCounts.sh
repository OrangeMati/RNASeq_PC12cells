#!/bin/bash

set -e # Instructs a shell to exit if a command fails
set -x # Prints out command arguments during execution.

FEATURECOUNTS_PATH=$1
BAM=$2
GTF=$3
OUTDIR=$4
BN_BAM=$( basename $BAM Aligned.sortedByCoord.out.bam )
OUTNAME=${OUTDIR}/${BN_BAM}
THREADS=8


# outdirectory must be created, otherwise error: cant create file
mkdir -p $OUTDIR

# single end
# S0
${FEATURECOUNTS_PATH} -s 0 -T $THREADS -t exon -g gene_id -a $GTF -o ${OUTNAME}_S0 $BAM
# S1
${FEATURECOUNTS_PATH} -s 1 -T $THREADS -t exon -g gene_id -a $GTF -o ${OUTNAME}_S1 $BAM
# S2
${FEATURECOUNTS_PATH} -s 2 -T $THREADS -t exon -g gene_id -a $GTF -o ${OUTNAME}_S2 $BAM
