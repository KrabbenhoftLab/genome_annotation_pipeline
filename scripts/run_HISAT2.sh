#!/bin/bash
# v0.2.4
# Pipeline to perform gene prediction and annotation
# author: Dan MacGuigan

# RUN HISAT2

# input variables passed from AISO_annotation_pipeline.sh
SPECIES=$1 # short name for your species
GENOME_DIR=$2 # directory containing your genome assembly
MASKED_GENOME_FILE=$3 # your genome assembly
RNA_R1=$4 # RNA read 1 files
RNA_R2=$5 # RNA read 2 files
RNA_UNPAIRED=$6 # unpaired RNA read files
HISAT_THREADS=$7

# create directory for HISAT2
mkdir ${SPECIES}_HISAT2


# first generate index files for reference genome
echo "starting to index reference genome"
cd ${GENOME_DIR}
hisat2-build -f -p ${HISAT_THREADS} ${MASKED_GENOME_FILE} ${MASKED_GENOME_FILE%.fasta}
cd ..

echo ""
echo "finished indexing, now mapping RNA-seq reads"
# now map RNA seq reads
cd ${SPECIES}_HISAT2
# if you DO NOT have unpaired RNA-seq reads
if [ "${RNA_UNPAIRED}" == "NONE" ]; then
	hisat2 -p ${HISAT_THREADS} -q -1 ${RNA_R1} -2 ${RNA_R2} -x ../${GENOME_DIR}/${MASKED_GENOME_FILE%.fasta} -S ${SPECIES}.rna.sam
elif [ "${RNA_1}" == "NONE" ]; then
	hisat2 -p ${HISAT_THREADS} -q -U ${RNA_UNPAIRED} -x ../${GENOME_DIR}/${MASKED_GENOME_FILE%.fasta} -S ${SPECIES}.rna.sam
else
	hisat2 -p ${HISAT_THREADS} -q -1 ${RNA_R1} -2 ${RNA_R2} -U ${RNA_UNPAIRED} -x ../${GENOME_DIR}/${MASKED_GENOME_FILE%.fasta} -S ${SPECIES}.rna.sam
fi
echo "finished RNA-seq read mapping"