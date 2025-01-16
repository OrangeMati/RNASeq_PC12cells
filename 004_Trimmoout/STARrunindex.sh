#!/bin/bash

set -e # Instructs a shell to exit if a command fails
set -x # Prints out command arguments during execution.

# Star index creation for Rattus norvegicus

/apps/anaconda3/envs/omics2023/bin/STAR \
 --runThreadN 16 \
 --runMode genomeGenerate \
 --genomeDir /home/fhwn.ac.at/201138/RNASeq_LW/000_referenceseq/star_indices \
 --genomeFastaFiles /home/fhwn.ac.at/201138/RNASeq_LW/000_referenceseq/ncbi-genomes-2023-06-03/GCF_015227675.2_mRatBN7.2_genomic.fna \
 --sjdbGTFfile /home/fhwn.ac.at/201138/RNASeq_LW/000_referenceseq/ncbi-genomes-2023-06-03/GCF_015227675.2_mRatBN7.2_genomic.gtf
 