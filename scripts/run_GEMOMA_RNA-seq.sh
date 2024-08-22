#!/bin/bash -l
# v0.4.5

#SBATCH --qos=general-compute
#SBATCH --partition=general-compute
#SBATCH --account=tkrabben
#SBATCH --time=72:00:00
#SBATCH --nodes=1
##SBATCH --constraint=AVX512
#SBATCH --export=NONE

# Pipeline to perform gene prediction and annotation
# author: Dan MacGuigan

# RUN GEMOMA
## Script to run GeMoMa on UB CCR Vortex cluster
## Authors: Christopher Osborne, Dan MacGuigan

# input variables passed from AISO_annotation_pipeline.sh
SPECIES=$1 # short name for your species
GENOME_DIR=$2 # directory containing your genome assembly
MASKED_GENOME_FILE=$3 # your (soft) masked genome assembly
GEMOMA_THREADS=$4
GEMOMA_RAM=$5
ANNOTATION_DIR=$6
GEMOMA_REFS=$7
USE_HISAT2_RNA_SEQ=$8
HISAT2_RNA_LIB_TYPE=$9
USE_CUSTOM_RNA_SEQ=${10}
CUSTOM_RNA_BAM=${11}
CUSTOM_RNA_LIB_TYPE=${12}
GEMOMA_SCORE_AA_FILTER=${13} # filter value for score/aa, see https://www.jstacs.de/index.php/GeMoMa-Docs
GEMOMA_FILTER_ONLY=${14} # filter previous GeMoMa run (yes or no)

# following lines fix a problem with dconf, solution from https://groups.google.com/g/slurm-users/c/CDCICHG7yLo
unset XDG_RUNTIME_DIR
unset XDG_SESSION_ID
unset XDG_DATA_DIRS

module load miniconda3/22.11.1-1
source activate GeMoMa_1.9
module load gcc/11.2.0 samtools/1.16.1

if [[ "${USE_CUSTOM_RNA_SEQ}" == "YES" ]]; then
	# make sure custom RNA-seq BAM file is sorted
	samtools sort ${CUSTOM_RNA_BAM} -o ${CUSTOM_RNA_BAM}.sorted.bam --threads ${GEMOMA_THREADS} -m ${GEMOMA_RAM}
fi

# make GeMoMa output directory
cd ${ANNOTATION_DIR}
mkdir -p ${SPECIES}_GeMoMa
cd ${SPECIES}_GeMoMa

# output directory for combined predictions
outDir_combined="GeMoMa_combined"
mkdir -p ${outDir_combined}
mkdir -p ${outDir_combined}/filter_score_${GEMOMA_SCORE_AA_FILTER}

if [ ${GEMOMA_FILTER_ONLY} = "yes" ]
then
	cd ${outDir_combined}/filter_score_${GEMOMA_SCORE_AA_FILTER}

	# create .sh to run GeMoMa Annotation Filter
	echo "GeMoMa -Xmx${GEMOMA_RAM} GAF \\" > GeMoMa.filter.${GEMOMA_SCORE_AA_FILTER}.sh

	i=0
	for GFF in ${GEMOMA_REFS}/*.gff
	do
		echo "g=\"../unfiltered_predictions_from_species_${i}.gff\" \\" >> GeMoMa.filter.${GEMOMA_SCORE_AA_FILTER}.sh
		i=$((i+1))
	done
	
	echo "f=\"start=='M' and stop=='*' and score/aa>=${GEMOMA_SCORE_AA_FILTER}\";" >> GeMoMa.filter.${GEMOMA_SCORE_AA_FILTER}.sh

	# run GeMoMa Annotation Filter
	conda activate GeMoMa_1.9
	bash GeMoMa.filter.${GEMOMA_SCORE_AA_FILTER}.sh

	# create .sh to run GeMoMa Annotation Finalizer
	echo "GeMoMa -Xmx${GEMOMA_RAM} AnnotationFinalizer \\" > GeMoMa.finalizer.${GEMOMA_SCORE_AA_FILTER}.sh
	echo "a=filtered_predictions.gff \\" >> GeMoMa.finalizer.${GEMOMA_SCORE_AA_FILTER}.sh
	echo "g=${ANNOTATION_DIR}/${GENOME_DIR}/${MASKED_GENOME_FILE} \\" >> GeMoMa.finalizer.${GEMOMA_SCORE_AA_FILTER}.sh
	echo "rename=NO;" >> GeMoMa.finalizer.${GEMOMA_SCORE_AA_FILTER}.sh

	# run GeMoMa Annotation Finalizer
	bash GeMoMa.finalizer.${GEMOMA_SCORE_AA_FILTER}.sh

	# get longest isoform of each gene
	AGAT_SIF="/projects/academic/tkrabben/software/agat/agat_1.0.0--pl5321hdfd78af_0.sif"
	singularity exec -H ${PWD} ${AGAT_SIF} agat_sp_keep_longest_isoform.pl --gff final_annotation.gff -o final_annotation.longest_isoform.gff

	# generate protein predictions for longest isoform
	GeMoMa -Xmx${GEMOMA_RAM} Extractor \
	 Ambiguity=AMBIGUOUS \
	 p=true \
	 g=${ANNOTATION_DIR}/${GENOME_DIR}/${MASKED_GENOME_FILE} \
	 a=final_annotation.longest_isoform.gff;

	mv proteins.fasta proteins.longest_isoform.fasta

else
	# loop through each file in GEMOMA_REFS
	conda activate genometools-1.6.2
	for GFF in ${GEMOMA_REFS}/*.gff
	do 
		# get gff file prefix
		temp1=${GFF%.gff} # remove file suffix
		sp=${temp1##*/} # remove path
		echo ${sp}
		gt gff3 -tidy -o ${sp}.clean.gff ${GFF} # clean up GFF file, make sure it's GFF3 format
	done

	cd ${outDir_combined}/filter_score_${GEMOMA_SCORE_AA_FILTER}
	
	# make .sh file to run combined GEMOMA
	echo "GeMoMa -Xmx${GEMOMA_RAM} GeMoMaPipeline \\" > run_GeMoMa.combined.sh
	echo "r=MAPPED \\" >> run_GeMoMa.combined.sh
	if [[ "${USE_HISAT2_RNA_SEQ}" == "YES" ]]; then
		echo "ERE.m=${ANNOTATION_DIR}/${SPECIES}_HISAT2/${SPECIES}.sorted.rna.bam \\" >> run_GeMoMa.combined.sh
		echo "ERE.s=${HISAT2_RNA_LIB_TYPE} \\" >> run_GeMoMa.combined.sh
	fi
	if [[ "${USE_CUSTOM_RNA_SEQ}" == "YES" ]]; then
		echo "ERE.m=${CUSTOM_RNA_BAM}.sorted.bam \\" >> run_GeMoMa.combined.sh
		echo "ERE.s=${CUSTOM_RNA_LIB_TYPE} \\" >> run_GeMoMa.combined.sh
	fi
	echo "threads=${GEMOMA_THREADS} \\" >> run_GeMoMa.combined.sh
	echo "AnnotationFinalizer.r=NO \\" >> run_GeMoMa.combined.sh
	echo "p=false \\" >> run_GeMoMa.combined.sh
	echo "o=true \\" >> run_GeMoMa.combined.sh
	echo "outdir=\".\" \\" >> run_GeMoMa.combined.sh
	echo "GAF.f=\"start=='M' and stop=='*' and score/aa>=${GEMOMA_SCORE_AA_FILTER}\" \\" >> run_GeMoMa.combined.sh

	# loop through each file in GEMOMA_REFS
	for GFF in ${GEMOMA_REFS}/*.gff
	do 
		# get gff file prefix
		temp1=${GFF%.gff} # remove file suffix
		sp=${temp1##*/} # remove path
		echo ${sp}
		echo "s=own \\" >> run_GeMoMa.combined.sh
		echo "i=${sp} \\" >> run_GeMoMa.combined.sh
		echo "a=${ANNOTATION_DIR}/${SPECIES}_GeMoMa/${sp}.clean.gff \\" >> run_GeMoMa.combined.sh
		echo "g=${GEMOMA_REFS}/${sp}.fasta \\" >> run_GeMoMa.combined.sh
	done

	echo "t=${ANNOTATION_DIR}/${GENOME_DIR}/${MASKED_GENOME_FILE};" >> run_GeMoMa.combined.sh

	# run GeMoMa
	conda activate GeMoMa_1.9
	bash run_GeMoMa.combined.sh

	# generate protien predictions
	GeMoMa -Xmx${GEMOMA_RAM} Extractor \
	 Ambiguity=AMBIGUOUS \
	 p=true \
	 outdir="." \
	 g=${ANNOTATION_DIR}/${GENOME_DIR}/${MASKED_GENOME_FILE} \
	 a=final_annotation.gff;

	mv proteins.fasta proteins.all.fasta

	# get longest isoform of each gene
	AGAT_SIF="/projects/academic/tkrabben/software/agat/agat_1.0.0--pl5321hdfd78af_0.sif"
	singularity exec -H ${PWD} ${AGAT_SIF} agat_sp_keep_longest_isoform.pl --gff final_annotation.gff -o final_annotation.longest_isoform.gff

	# generate protein predictions for longest isoform
	GeMoMa -Xmx${GEMOMA_RAM} Extractor \
	 Ambiguity=AMBIGUOUS \
	 p=true \
	 outdir="." \
	 g=${ANNOTATION_DIR}/${GENOME_DIR}/${MASKED_GENOME_FILE} \
	 a=final_annotation.longest_isoform.gff;

	mv proteins_1.fasta proteins.longest_isoform.fasta
	
	cp unfiltered_predictions_from_species_*.gff ${ANNOTATION_DIR}/${SPECIES}_GeMoMa/GeMoMa_combined # copy unfiltered predictions in case user wants to filter post hoc later
fi


echo ""
echo "Step 7 COMPLETE"
echo "Please check to make sure that a combined GFF was produced by GeMoMa:"
echo "${ANNOTATION_DIR}/${SPECIES}_GeMoMa/GeMoMa_combined/filter_score_${GEMOMA_SCORE_AA_FILTER}/final_annotation.longest_isoform.gff"
echo ""
echo "If Step 6 (BRAKER) has also finished, you may proceed to Step 8: combining gene predictions with EVM."