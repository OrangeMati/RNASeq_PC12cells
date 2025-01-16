#!/bin/bash

set -e # Instructs a shell to exit if a command fails
set -x # Prints out command arguments during execution.

# ACHTUNG: fuer qualimap braucht man java-11-openjdk.x86_64 (/usr/lib/jvm/java-11-openjdk-11.0.17.0.8-1.fc37.x86_64/bin/java)
# am highmem war urspruenglich nur java-17-openjdk.x86_64 (/usr/lib/jvm/java-17-openjdk-17.0.5.0.8-1.fc37.x86_64/bin/java)
# da gab es dann diesen Fehler
# Unrecognized VM option 'MaxPermSize=1024m' Error: Could not create the Java Virtual Machine. Error: A fatal exception has occurred. Program will exit.

# ich habe dann java 11 nach installiert als su, und dann als default gesetzt
# [root@i118highmem1 mapping_star]# dnf install java-11-openjdk
# [root@i118highmem1 mapping_star]# update-alternatives --config java

#There are 3 programs which provide 'java'.

#  Selection    Command
#-----------------------------------------------
#*+ 1           java-17-openjdk.x86_64 (/usr/lib/jvm/java-17-openjdk-17.0.5.0.8-1.fc37.x86_64/bin/java)
#   2           java-11-openjdk.x86_64 (/usr/lib/jvm/java-11-openjdk-11.0.17.0.8-1.fc37.x86_64/bin/java)
#   3           java-1.8.0-openjdk.x86_64 (/usr/lib/jvm/java-1.8.0-openjdk-1.8.0.352.b08-2.fc37.x86_64/jre/bin/java)

#Enter to keep the current selection[+], or type selection number: 2

#PATH_QUALIMAP=/apps/anaconda3/envs/omics2023/bin/qualimap
#BAM=BSF_examplefile_sorted_1_SE_Aligned.sortedByCoord.out.bam
#GTF=/proj/courses/transcriptomics/SS23/localmirror/genomes/gencode/ratannotation.gtf
#OUTDIR=$( pwd )/../007_QUALImap.out_html
#OUTFORMAT=HTML
#MEMSIZE=32G
#BN_BAM=$( basename $BAM Aligned.sortedByCoord.out.bam )
#


# input
PATH_QUALIMAP=$1
BAM=$2
GTF=$3
OUTDIR=$4
OUTFORMAT=$5 # PDF or HTML
MEMSIZE=$6
BN_BAM=$( basename $BAM Aligned.sortedByCoord.out.bam )

mkdir -p ${OUTDIR}/${BN_BAM}_out

${PATH_QUALIMAP} rnaseq -bam $BAM -gtf $GTF -outdir $OUTDIR/${BN_BAM}_out -outformat $OUTFORMAT --java-mem-size=$MEMSIZE


