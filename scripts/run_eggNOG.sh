#!/bin/bash
# v0.3.0
#SBATCH --qos=general-compute
#SBATCH --partition=general-compute
#SBATCH --account=tkrabben
#SBATCH --time=72:00:00
#SBATCH --nodes=1
#SBATCH --reservation=ubhpc-future

SPECIES=$1
ANNOTATION_DIR_CLUSTER=$2
BRAKER_WEIGHT=$3
GEMOMA_WEIGHT=$4
EGG_THREADS=$5
EGGNOG_OPTIONS=$6

module load miniconda3/22.11.1-1
source activate eggnog-mapper

# use * after pep in case user has compressed the peptide fasta
emapper.py --cpu ${EGG_THREADS} -i ${SPECIES}_EVM_B${BRAKER_WEIGHT}_G${GEMOMA_WEIGHT}/${SPECIES}.EVM.pep* --output_dir {SPECIES}_EVM_B${BRAKER_WEIGHT}_G${GEMOMA_WEIGHT}_eggNOG -o Engr.EVM.B${BRAKER_WEIGHT}.G${GEMOMA_WEIGHT} --decorate_gff ${SPECIES}_EVM_B${BRAKER_WEIGHT}_G${GEMOMA_WEIGHT}/${SPECIES}.EVM.mod.gff3 --excel ${EGGNOG_OPTIONS}