#!/bin/bash
# v0.3.1
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

module load gcc/11.2.0 openmpi/4.1.1 eggnog-mapper/2.1.12
export EGGNOG_DATA_DIR="/projects/academic/tkrabben/software/eggNOG_db/eggnog-mapper-data"
# need to download new eggNOG DBs?
# load the modules, export the EGGNOG_DATA_DIR, and run download_eggnog_data.py

emapper.py --cpu ${EGG_THREADS} -i ${SPECIES}_EVM_B${BRAKER_WEIGHT}_G${GEMOMA_WEIGHT}/${SPECIES}.EVM.pep --output_dir ${SPECIES}_EVM_B${BRAKER_WEIGHT}_G${GEMOMA_WEIGHT}_eggNOG -o Engr.EVM.B${BRAKER_WEIGHT}.G${GEMOMA_WEIGHT} --decorate_gff ${SPECIES}_EVM_B${BRAKER_WEIGHT}_G${GEMOMA_WEIGHT}/${SPECIES}.EVM.mod.gff3 --excel ${EGGNOG_OPTIONS}
