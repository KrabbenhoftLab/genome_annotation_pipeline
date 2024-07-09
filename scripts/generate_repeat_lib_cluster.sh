#!/bin/bash -l
# v0.4.3

#SBATCH --qos=general-compute
#SBATCH --partition=general-compute
#SBATCH --account=tkrabben
#SBATCH --time=72:00:00
#SBATCH --nodes=1
##SBATCH --constraint=AVX512
#SBATCH --export=NONE

# Pipeline to perform gene prediction and annotation
# author: Dan MacGuigan

# GENERATE CUSTOM REPEAT LIBRARY

# load modules
module load gcc/11.2.0 openmpi/4.1.1 repeatmodeler/2.0.4.KRAB

whereis RepeatMasker
whereis RepeatModeler

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

# run RepeatClassifier if it fails, you will need to modify the line below
#RepeatClassifier -consensi ./RM_18866.FriApr101507532020/families.fa -stockholm ./RM_18866.FriApr101507532020/families.stk

echo ""
echo "Step 1 COMPLETE"
echo "Please check to make sure that your custom repeat family FASTA file was created:"
echo "${ANNOTATION_DIR_CLUSTER}/${SPECIES}_RepeatModeler/${REPEAT_LIBRARY_NAME}-families.fa"
echo ""
echo "If this file exists, you may proceed to Step 2: repeat masking"
echo ""
echo "If the repeat families file does not exist, try resuming Step 1 from a checkpoint."
echo "You will need to set CONTINUE_RMODEL=\"yes\""
echo "and specify the CONTINUE_RMODEL_DIR option"
echo "in your config file."
echo ""
echo "Occasionally, all 5 rounds of RepeatModeler will finish, but the LTR Structural Analysis"
echo "will not. Usually this is due to hitting a walltime limit, producing an error message like:"
echo 'slurmstepd: error: *** JOB 14904061 ON cpn-i16-24 CANCELLED AT 2024-02-11T11:09:23 ***'
echo ""
echo "If this happens, you will not see a ${REPEAT_LIBRARY_NAME}-families.fa file in your"
echo "RepeatModeler directory. Sadly, RepeatModeler is not smart enough to recognize this as"
echo "a failed run. If you restart Step 1, you will get a message like:"
echo "\"This directory appears to contain a successful run of RepeatModeler\""
echo ""
echo "In this case, in order to resume your failed RepeatModeler run, take the following steps:"
echo " - navigate to ${SPECIES}_RepeatModeler/RM_SOMENUMBER.DATE/round-5"
echo " - delete the file \"consensi.fa\""
echo " - rerun Step 1 of this pipeline"
echo ""
echo "Unfortunately, this workaround requires RepeatModeler to restart round 5 from scratch."
echo "After round 5 is complete, it will then run the LTR Structural Analysis, and eventually"
echo "produce a ${REPEAT_LIBRARY_NAME}-families.fa file."
echo ""