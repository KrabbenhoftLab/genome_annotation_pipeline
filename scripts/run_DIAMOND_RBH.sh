#!/bin/bash -l
# v0.4.0
#SBATCH --qos=general-compute
#SBATCH --partition=general-compute
#SBATCH --account=tkrabben
#SBATCH --time=72:00:00
#SBATCH --nodes=1
#SBATCH --export=NONE

### ABOUT ###
#Script to perform a reciprocal DIAMOND search
#modified from https://morphoscape.wordpress.com/2020/08/18/reciprocal-best-hits-DIAMOND-rbhb/

#if running on the cluster, uncomment the line below to load the DIAMOND conda environment
module load miniconda3/22.11.1-1
source activate DIAMOND

### INPUTS ###
ANNOTATION_DIR_CLUSTER=$1
SPECIES=$2
BRAKER_WEIGHT_DIAMOND=$3 # EVM weight for BRAKER predictions
GEMOMA_WEIGHT_DIAMOND=$4 # EVM weight for GeMoMa predictions
GEMOMA_SCORE_AA_FILTER=$5
DATABASE=$6 # path to protein FASTA file to use a search database, must be located within ANNOTATION_DIR_CLUSTER
E_VAL=$7 # 
DIAMOND_THREADS=$8

# set up directories
mkdir -p ${ANNOTATION_DIR_CLUSTER}/${SPECIES}_DIAMOND_RBH
mkdir -p ${ANNOTATION_DIR_CLUSTER}/${SPECIES}_DIAMOND_RBH/B${BRAKER_WEIGHT_DIAMOND}_G${GEMOMA_WEIGHT_DIAMOND}-SCORE-${GEMOMA_SCORE_AA_FILTER}

cd ${ANNOTATION_DIR_CLUSTER}/${SPECIES}_DIAMOND_RBH/B${BRAKER_WEIGHT_DIAMOND}_G${GEMOMA_WEIGHT_DIAMOND}-SCORE-${GEMOMA_SCORE_AA_FILTER}
 
#full path to query file, should be protein FASTA
inputQuery="${ANNOTATION_DIR_CLUSTER}/${SPECIES}_EVM/B${BRAKER_WEIGHT_DIAMOND}_G${GEMOMA_WEIGHT_DIAMOND}-SCORE-${GEMOMA_SCORE_AA_FILTER}/${SPECIES}.EVM.pep"
#full path to DB reciprocal file, should be protein FASTA
inputDB="${ANNOTATION_DIR_CLUSTER}/${DATABASE}"
#full path to output results
outputPath="${ANNOTATION_DIR_CLUSTER}/${SPECIES}_DIAMOND_RBH/B${BRAKER_WEIGHT_DIAMOND}_G${GEMOMA_WEIGHT_DIAMOND}-SCORE-${GEMOMA_SCORE_AA_FILTER}"
#DIAMOND threads
threads=${DIAMOND_THREADS}
#DIAMOND E value, see https://blast.ncbi.nlm.nih.gov/doc/blast-help/FAQ.html
e_val=${E_VAL}

### SCRIPT ###

#Move to query directory
queryPath=$(dirname "$inputQuery")
cd "${queryPath}"
#Make DIAMOND protein DB of the input query
diamond makedb --in "${inputQuery}" --db "${inputQuery}"

#Move to DB directory
dbPath=$(dirname "$inputDB")
cd "${dbPath}"
#Make DIAMOND protein DB of the input DB
diamond makedb --in "${inputDB}" --db "${inputDB}"

#Output start status message
echo "Beginning reciprocal DIAMOND..."

#Move to outputs folder
cd "${outputPath}"

#Use DIAMONDp to search a database
diamond blastp --ultra-sensitive --query "${inputQuery}" --db "${inputDB}" --outfmt 6 --max-target-seqs 1 --evalue ${e_val} --threads ${threads} > DIAMOND.outfmt6

#Switch query and search paths for reciprocal search
diamond blastp --ultra-sensitive --query "${inputDB}" --db "${inputQuery}" --outfmt 6 --max-target-seqs 1 --evalue "${e_val}" --threads "${threads}" > DIAMOND_reciprocal.outfmt6

#Output end status message
echo "Finished reciprocal DIAMOND!"

#now to filter reciprocal DIAMOND results for best hits

#Input query DIAMOND results file
queryPath="DIAMOND.outfmt6"
#Input DB reciprocal DIAMOND results file
dbPath="DIAMOND_reciprocal.outfmt6"

#Final output files
outFileRBH="DIAMOND_RBH.txt"
summaryFile="DIAMOND_RBH_summary.txt"

#Add headers to output RBH files
echo "queryHit,dbHit" > ${outFileRBH}

echo "queryHits,dbHits,bestHits" > ${summaryFile}

echo ${summaryFile}

#Output start status message
echo "Recording RBH..."

#Loop over query DIAMOND results
while IFS=$'\t' read -r f1 f2 f3 f4 f5 f6 f7 f8 f9 f10 f11 f12
do
	#Determine RBH to DB DIAMOND results
	if grep -q "$f2"$'\t'"$f1"$'\t' ${dbPath}; then #RBH
		echo "$f1,$f2" >> ${outFileRBH}
	fi
done < ${queryPath}

#Output summary of RBH
queryHits=$(wc -l "$queryPath" | cut -d ' ' -f 1)
dbHits=$(wc -l "$dbPath" | cut -d ' ' -f 1)
bestHits=$(($(wc -l "$outFileRBH" | cut -d ' ' -f 1)-1))
echo "${queryHits}","${dbHits}","${bestHits}" >> "${summaryFile}"

#Output end status message
echo "Finished recording RBH!"
