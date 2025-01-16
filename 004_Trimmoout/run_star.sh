#!/bin/bash

set -e # Instructs a shell to exit if a command fails
set -x # Prints out command arguments during execution.

# star works with the unzipped fasta, gtf, and fastq files -> gzip before running star 

# input
#OUTDIR=$( pwd )/../006_star.out_bam
#PATH_STAR=/apps/anaconda3/envs/omics2023/bin/STAR
#PATH_SAMTOOLS=/apps/anaconda3/envs/omics2023/bin/samtools
#READS1=BSF_examplefile.fastq.gz
#READS1_BASE=$( basename $READS1 .fastq.gz )

# input
FASTA=$1
ANNO=$2
INDEXDIR=$3
#INDEXNAME=$4
OUTDIR=$4
PATH_STAR=$5
PATH_SAMTOOLS=$6
READS1=$7
READS1_BASE=$( basename $READS1 .fastq.gz )

mkdir -p $OUTDIR

# settings
THREADS=16

# star
# run star mapping gz files
if [[ $READS1 == *gz ]]; then
	UNZIP="--readFilesCommand zcat"
	R1="$READS1 $UNZIP"
else
	R1="$READS1"
fi

${PATH_STAR} \
	--runThreadN $THREADS \
	--outSAMtype BAM SortedByCoordinate \
	--genomeDir ${INDEXDIR} \
	--readFilesIn $R1 \
	--outFileNamePrefix ${OUTDIR}/${READS1_BASE}_

# index bam file
$PATH_SAMTOOLS index ${OUTDIR}/${READS1_BASE}_Aligned.sortedByCoord.out.bam


# this is an example call to star: 
#${PATH_STAR} --runThreadN $THREADS --outSAMtype BAM SortedByCoordinate --genomeDir ${INDEXDIR} \
#--readFilesIn $R1 --readFilesCommand zcat --outFileNamePrefix ${OUTDIR}/${READS1_BASE}_
#+ /apps/anaconda3/envs/omics2023/bin/STAR --runThreadN 16 --outSAMtype BAM SortedByCoordinate \
#--genomeDir /proj/courses/transcriptomics/SS23/localmirror/index/gencode/... \
#--readFilesIn BSF_examplefile_SE.fastq.gz --readFilesCommand zcat\
# --outFileNamePrefix \
#/proj/courses/transcriptomics/SS23/data/subsample/trimmedFastq/../mapping_star/BSF_examplefile_sorted_1_SE_
#Apr 16 13:36:01 ..... started STAR run
#Apr 16 13:36:01 ..... loading genome
#Apr 16 13:36:20 ..... started mapping
#Apr 16 13:36:28 ..... finished mapping
#Apr 16 13:36:30 ..... started sorting BAM
#Apr 16 13:36:31 ..... finished successfully
#
## this is an example call to samtools for indexing
#$PATH_SAMTOOLS index ${OUTDIR}/${READS1_BASE}_Aligned.sortedByCoord.out.bam
#+ /apps/anaconda3/envs/omics2023/bin/samtools index \
#/proj/courses/transcriptomics/SS23/data/subsample/trimmedFastq/../mapping_star/BSF_examplefile_\
#SE_Aligned.sortedByCoord.out.bam
