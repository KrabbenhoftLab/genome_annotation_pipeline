#!/bin/bash
# v0.3.3
#SBATCH --qos=general-compute
#SBATCH --partition=general-compute
#SBATCH --account=tkrabben
#SBATCH --time=72:00:00
#SBATCH --nodes=1
#SBATCH --reservation=ubhpc-future

SPECIES=$1
ANNOTATION_DIR_CLUSTER=$2
BRAKER_WEIGHT_EGG=$3
GEMOMA_WEIGHT_EGG=$4
EGG_THREADS=$5
EGGNOG_OPTIONS=$6

module load gcc/11.2.0 openmpi/4.1.1 eggnog-mapper/2.1.12
export EGGNOG_DATA_DIR="/projects/academic/tkrabben/software/eggNOG_db/eggnog-mapper-data"
# need to download new eggNOG DBs?
# load the modules, export the EGGNOG_DATA_DIR, and run download_eggnog_data.py

emapper.py --cpu ${EGG_THREADS} -i ${SPECIES}_EVM_B${BRAKER_WEIGHT_EGG}_G${GEMOMA_WEIGHT_EGG}/${SPECIES}.EVM.pep --output_dir ${SPECIES}_EVM_B${BRAKER_WEIGHT_EGG}_G${GEMOMA_WEIGHT_EGG}_eggNOG -o Engr.EVM.B${BRAKER_WEIGHT_EGG}.G${GEMOMA_WEIGHT_EGG} --decorate_gff ${SPECIES}_EVM_B${BRAKER_WEIGHT_EGG}_G${GEMOMA_WEIGHT_EGG}/${SPECIES}.EVM.mod.gff3 --excel ${EGGNOG_OPTIONS}
