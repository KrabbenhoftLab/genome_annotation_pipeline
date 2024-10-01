#!/bin/bash -l
# v0.4.7
#SBATCH --qos=general-compute
#SBATCH --partition=general-compute
#SBATCH --account=tkrabben
#SBATCH --time=1:00:00
#SBATCH --mem=100G
#SBATCH --ntasks-per-node=1
#SBATCH --nodes=1
#SBATCH --export=NONE

ANNOTATION_DIR_CLUSTER=$1
SPECIES=$2 # short name for your species
GENOME_DIR=$3 # directory containing your genome assembly
GENOME_FILE=$4 # your genome assembly
MASKED_GENOME_FILE=$5 # name for soft masked genome assembly file

module load gcc/11.2.0 bedtools/2.30.0

mkdir -p "${ANNOTATION_DIR_CLUSTER}/logFiles" # create logFiles directory

bedtools maskfasta -bed "${ANNOTATION_DIR_CLUSTER}/${SPECIES}_RepeatMasker/final_repeat_mask/${SPECIES}.final_repeat_mask.gff3" -fi "${ANNOTATION_DIR_CLUSTER}/${GENOME_DIR}/${GENOME_FILE}" -fo "${ANNOTATION_DIR_CLUSTER}/${GENOME_DIR}/${MASKED_GENOME_FILE}" -soft
		
echo ""
echo "Step 3 COMPLETE"
echo "please check to see that a masked genome FASTA was created:"
echo "${ANNOTATION_DIR_CLUSTER}/${GENOME_DIR}/${MASKED_GENOME_FILE}"
echo ""
echo "if you have RNA-seq data, you may now proceed to Step 4: HISAT2 read mapping"
echo "if you do not have RNA-seq data, proceed to Step 5: downloading reference protein datasets from NCBI"
echo "if you have already downloaded reference protein datasets, proceed to Steps 6 and 7: BRAKER and GeMoMa"
echo "Steps 6 and 7 may be run simultaneously"