#!/bin/bash

PATH_QUALIMAP=/apps/anaconda3/envs/omics2023/bin/qualimap
#BAM=BSF_examplefile_sorted_1_SE_Aligned.sortedByCoord.out.bam
GTF=/home/fhwn.ac.at/201138/RNASeq_LW/000_referenceseq/ncbi-genomes-2023-06-03/GCF_015227675.2_mRatBN7.2_genomic.gtf
OUTDIR=$( pwd )/../007_QUALImapout
OUTFORMAT=HTML
MEMSIZE=32G
#BN_BAM=$( basename $BAM Aligned.sortedByCoord.out.bam )

for BAM_IN in $( ls -1 *bam )
do
bash ./run_qualimap.sh ${PATH_QUALIMAP} ${BAM_IN} ${GTF} ${OUTDIR} ${OUTFORMAT} ${MEMSIZE}
done
