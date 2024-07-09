#!/bin/bash -l
# v0.4.3

#SBATCH --qos=general-compute
#SBATCH --partition=general-compute
#SBATCH --account=tkrabben
#SBATCH --time=72:00:00
#SBATCH --nodes=1


SPECIES=$1
ANNOTATION_DIR_CLUSTER=$2
BRAKER_WEIGHT_EGG=$3
GEMOMA_WEIGHT_EGG=$4
EGG_THREADS=$5
GEMOMA_SCORE_AA_FILTER=$6
EGGNOG_OPTIONS=$7

module load gcc/11.2.0 openmpi/4.1.1 eggnog-mapper/2.1.12
export EGGNOG_DATA_DIR="/projects/academic/tkrabben/software/eggNOG_db/eggnog-mapper-data"
# need to download new eggNOG DBs?
# load the modules, export the EGGNOG_DATA_DIR, and run download_eggnog_data.py

emapper.py --cpu ${EGG_THREADS} -i ${SPECIES}_EVM/B${BRAKER_WEIGHT_EGG}_G${GEMOMA_WEIGHT_EGG}-SCORE-${GEMOMA_SCORE_AA_FILTER}/${SPECIES}.EVM.pep --output_dir ${SPECIES}_eggNOG/B${BRAKER_WEIGHT_EGG}_G${GEMOMA_WEIGHT_EGG}-SCORE-${GEMOMA_SCORE_AA_FILTER} -o ${SPECIES}.EVM.B${BRAKER_WEIGHT_EGG}.G${GEMOMA_WEIGHT_EGG} --decorate_gff ${SPECIES}_EVM/B${BRAKER_WEIGHT_EGG}_G${GEMOMA_WEIGHT_EGG}-SCORE-${GEMOMA_SCORE_AA_FILTER}/${SPECIES}.EVM.mod.gff3 --excel ${EGGNOG_OPTIONS}

echo ""
echo "Step 9 COMPLETE"
echo "EggNOG-mapper should have produced several files here:"
echo "${ANNOTATION_DIR_CLUSTER}/${SPECIES}_eggNOG/B${BRAKER_WEIGHT_EGG}_G${GEMOMA_WEIGHT_EGG}-SCORE-${GEMOMA_SCORE_AA_FILTER}"
echo ""
echo "You can manually examine the functional annotations in Excel:"
echo "${SPECIES}.EVM.B${BRAKER_WEIGHT_EGG}.G${GEMOMA_WEIGHT_EGG}.emapper.annotations.xlsx"
echo ""
echo "There should also be an annotated GFF file:"
echo "${SPECIES}.EVM.B${BRAKER_WEIGHT_EGG}.G${GEMOMA_WEIGHT_EGG}.emapper.decorated.gff" 
