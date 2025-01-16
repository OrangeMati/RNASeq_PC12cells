#!/bin/bash

set -e # Instructs a shell to exit if a command fails
set -x # Prints out command arguments during execution.

# star works with the unzipped fasta, gtf, and fastq files -> gzip before running star 


FASTA=/home/fhwn.ac.at/201138/RNASeq_LW/000_referenceseq/ncbi-genomes-2023-06-03/GCF_015227675.2_mRatBN7.2_genomic.fna
ANNO=/home/fhwn.ac.at/201138/RNASeq_LW/000_referenceseq/ncbi-genomes-2023-06-03/GCF_015227675.2_mRatBN7.2_genomic.gtf
INDEXDIR=/home/fhwn.ac.at/201138/RNASeq_LW/000_referenceseq/star_indices
INDEXNAME="" #it depends how the star index is done, whether this parameter is known; on highmem1, index is in primary_index subfolder, 
# on tubdsnode01, the files of the index are directly unter M32
OUTDIR=$( pwd )/../006_STARout
PATH_STAR=/apps/anaconda3/envs/omics2023/bin/STAR
PATH_SAMTOOLS=/apps/anaconda3/envs/omics2023/bin/samtools
#READS1=BSF_examplefile.fastq.gz
#READS1_BASE=$( basename $READS1 .fastq.gz )

# input
#FASTA=$1
#ANNO=$2
#INDEXDIR=$3
#INDEXNAME=$4
#OUTDIR=$4
#PATH_STAR=$5
#PATH_SAMTOOLS=$6
#READS1=$7

for FASTQGZ_IN in $( ls -1 *fastq.gz )
do
bash ./run_star.sh ${FASTA} ${ANNO} ${INDEXDIR} ${OUTDIR} ${PATH_STAR} ${PATH_SAMTOOLS} ${FASTQGZ_IN}
done
