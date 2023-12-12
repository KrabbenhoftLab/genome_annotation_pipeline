#!/bin/bash
# v0.3.2
# Pipeline to perform gene prediction and annotation
# author: Dan MacGuigan

# RUN HISAT2

# input variables passed from AISO_annotation_pipeline.sh
SPECIES=$1 # short name for your species
GENOME_DIR=$2 # directory containing your genome assembly
MASKED_GENOME_FILE=$3 # your genome assembly
RNA_DIR=$4 # RNA read 1 files
RNA_FILES=$5 # RNA read 2 files
HISAT_THREADS=$6

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

while read line; do
	line_arr=($line)
	if [ ${#line_arr[@]} -eq 1 ]; then # unpaired RNA-seq reads
		hisat2 -p ${HISAT_THREADS} -q -U ${RNA_DIR}/${line_arr[0]} -x ../${GENOME_DIR}/${MASKED_GENOME_FILE%.fasta} -S ${SPECIES}.${line_arr[0]}.rna.sam
	elif [ ${#line_arr[@]} -eq 2 ]; then # paired RNA-seq reads
		hisat2 -p ${HISAT_THREADS} -q -1 ${RNA_DIR}/${line_arr[0]} -2 ${RNA_DIR}/${line_arr[1]} -x ../${GENOME_DIR}/${MASKED_GENOME_FILE%.fasta} -S ${SPECIES}.${line_arr[0]}.rna.sam
	else
		echo "${line} --> contains more than two files"
	fi
	
	samtools view -S -b ${SPECIES}.${line_arr[0]}.rna.sam --threads ${HISAT_THREADS} > ${SPECIES}.${line_arr[0]}.rna.bam 
	samtools sort  ${SPECIES}.${line_arr[0]}.rna.bam -o ${SPECIES}.${line_arr[0]}.sorted.rna.bam --threads ${HISAT_THREADS}
	rm ${SPECIES}.${line_arr[0]}.rna.sam
	rm ${SPECIES}.${line_arr[0]}.rna.bam
	echo "finished mapping ${line}"

done < ../${RNA_FILES}

echo "finished mapping all files"
echo "merging bam files"
samtools merge ${SPECIES}.rna.bam ${SPECIES}.*.sorted.rna.bam --threads ${HISAT_THREADS}
samtools sort ${SPECIES}.rna.bam -o ${SPECIES}.sorted.rna.bam --threads ${HISAT_THREADS} -m 2G
rm ${SPECIES}.rna.bam

echo "merged bam file: ${SPECIES}.sorted.rna.bam"
echo "step complete"
