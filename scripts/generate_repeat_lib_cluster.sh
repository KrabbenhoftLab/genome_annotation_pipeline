#!/bin/bash
# v0.4.0
#SBATCH --qos=general-compute
#SBATCH --partition=general-compute
#SBATCH --account=tkrabben
#SBATCH --time=72:00:00
#SBATCH --nodes=1
#SBATCH --constraint=AVX512
#SBATCH --export=NONE
#SBATCH --reservation=ubhpc-future

# Pipeline to perform gene prediction and annotation
# author: Dan MacGuigan

# GENERATE CUSTOM REPEAT LIBRARY

# load modules
module load gcc/11.2.0 openmpi/4.1.1 repeatmodeler/2.0.4

# input variables passed from AISO_annotation_pipeline.sh
SPECIES=$1 # short name for your species
ANNOTATION_DIR_CLUSTER=$2
GENOME_DIR=$3 # directory containing your genome assembly
GENOME_FILE=$4 # your genome assembly
REPEAT_LIBRARY_NAME=$5 # name for your repeat library
RM_THREADS=$6 # Specify the number of parallel search jobs to run. RMBlast jobs will
              # use 4 cores each and ABBlast jobs will use a single core each. i.e.
              # on a machine with 12 cores and running with RMBlast you would use
RUN_LTR=$7
CONTINUE_RMODEL=$8
CONTINUE_RMODEL_DIR=$9

# move to correct cluster directory
cd ${ANNOTATION_DIR_CLUSTER}

# create directory for RepeatModeler
mkdir ${SPECIES}_RepeatModeler
cd ${SPECIES}_RepeatModeler

# load conda module
#source activate RepeatModeler-2.0.4

# make a blast database from the genome
makeblastdb -in ${ANNOTATION_DIR_CLUSTER}/${GENOME_DIR}/${GENOME_FILE} -dbtype nucl -parse_seqids
BuildDatabase -name ${REPEAT_LIBRARY_NAME} ${ANNOTATION_DIR_CLUSTER}/${GENOME_DIR}/${GENOME_FILE}

# run RepeatModeler
if [ ${CONTINUE_RMODEL} = "yes" ]
then
	if [ ${RUN_LTR} = "yes" ]
	then
		RepeatModeler -database ${REPEAT_LIBRARY_NAME} -threads ${RM_THREADS} -LTRStruct -recoverDir ${CONTINUE_RMODEL_DIR}
	else
		RepeatModeler -database ${REPEAT_LIBRARY_NAME} -threads ${RM_THREADS} -recoverDir ${CONTINUE_RMODEL_DIR}
	fi
else
	if [ ${RUN_LTR} = "yes" ]
	then
		RepeatModeler -database ${REPEAT_LIBRARY_NAME} -threads ${RM_THREADS} -LTRStruct
	else
		RepeatModeler -database ${REPEAT_LIBRARY_NAME} -threads ${RM_THREADS} 
	fi
fi

# continue RepeatModeler run if necessary, you will need to modify "-recoverDir"
#RepeatModeler -database ${REPEAT_LIBRARY_NAME} -threads ${RM_THREADS} -LTRStruct -recoverDir RM_18866.FriApr101507532020

# run RepeatClassifier if it fails, you will need to modify the line below
#RepeatClassifier -consensi ./RM_18866.FriApr101507532020/families.fa -stockholm ./RM_18866.FriApr101507532020/families.stk
