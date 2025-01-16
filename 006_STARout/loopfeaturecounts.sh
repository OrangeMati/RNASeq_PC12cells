#!/bin/bash

# you need to be in 007_star_out

FEATURECOUNTS_PATH="/apps/anaconda3/envs/omics2023/bin/featureCounts"
GTF=/home/fhwn.ac.at/201138/RNASeq_LW/000_referenceseq/ncbi-genomes-2023-06-03/GCF_015227675.2_mRatBN7.2_genomic.gtf
OUTDIR=$( pwd )/../008_FeatureCountsout
THREADS=8

for STAR_IN in $( ls -1 *out.bam )
do
#BN_BAM=$( basename $STAR_IN Aligned.sortedByCoord.out.bam )
#OUTNAME=${OUTDIR}/${BN_BAM}
bash ./run_featureCounts.sh ${FEATURECOUNTS_PATH} ${STAR_IN} ${GTF} ${OUTDIR}
done