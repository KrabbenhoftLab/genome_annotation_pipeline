#!/bin/bash
# v0.4.7

# Pipeline to perform gene prediction and annotation 
# author: Dan MacGuigan

# GENERATE CUSTOM REPEAT LIBRARY

# input variables passed from AISO_annotation_pipeline.sh
SPECIES=$1 # short name for your species
GENOME_DIR=$2 # directory containing your genome assembly
GENOME_FILE=$3 # your genome assembly
REPEAT_LIBRARY_NAME=$4 # name for your repeat library
RM_THREADS=$5 # Specify the number of parallel search jobs to run. RMBlast jobs will
              # use 4 cores each and ABBlast jobs will use a single core each. i.e.
              # on a machine with 12 cores and running with RMBlast you would use

# create directory for RepeatModeler
mkdir ${SPECIES}_RepeatModeler
cd ${SPECIES}_RepeatModeler

# load conda module
source activate RepeatModeler-2.0.4

# make a blast database from the genome
makeblastdb -in ../${GENOME_DIR}/${GENOME_FILE} -dbtype nucl -parse_seqids
BuildDatabase -name ${REPEAT_LIBRARY_NAME} ../${GENOME_DIR}/${GENOME_FILE}

# run RepeatModeler
RepeatModeler -database ${REPEAT_LIBRARY_NAME} -threads ${RM_THREADS} -LTRStruct

# continue RepeatModeler run if necessary, you will need to modify "-recoverDir"
#RepeatModeler -database ${REPEAT_LIBRARY_NAME} -threads ${RM_THREADS} -LTRStruct -recoverDir RM_18866.FriApr101507532020

# run RepeatClassifier if it fails, you will need to modify the line below
#RepeatClassifier -consensi ./RM_18866.FriApr101507532020/families.fa -stockholm ./RM_18866.FriApr101507532020/families.stk
