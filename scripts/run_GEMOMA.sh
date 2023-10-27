#!/bin/bash
# v0.2.4
#SBATCH --qos=general-compute
#SBATCH --partition=general-compute
#SBATCH --account=tkrabben
#SBATCH --time=72:00:00
#SBATCH --nodes=1
#SBATCH --constraint=AVX512
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

# following lines fix a problem with dconf, solution from https://groups.google.com/g/slurm-users/c/CDCICHG7yLo
unset XDG_RUNTIME_DIR
unset XDG_SESSION_ID
unset XDG_DATA_DIRS

module load miniconda3/22.11.1-1
source activate GeMoMa_1.9

# make GeMoMa output directory
cd ${ANNOTATION_DIR}
mkdir ${SPECIES}_GeMoMa
cd ${SPECIES}_GeMoMa

# output directory for combined predictions
outDir_combined="GeMoMa_combined"
mkdir ${outDir_combined}

# make .sh file to run combined GEMOMA
echo "GeMoMa GeMoMaPipeline -XX${GEMOMA_RAM} \\" > run_GeMoMa.combined.sh
echo "threads=${GEMOMA_THREADS} \\" >> run_GeMoMa.combined.sh
echo "AnnotationFinalizer.r=NO \\" >> run_GeMoMa.combined.sh
echo "p=false \\" >> run_GeMoMa.combined.sh
echo "o=true \\" >> run_GeMoMa.combined.sh
echo "outdir=${outDir_combined} \\" >> run_GeMoMa.combined.sh

# loop through each file in GEMOMA_REFS
conda activate genometools-1.6.2
for GFF in ${GEMOMA_REFS}/*.gff
do 
	# get gff file prefix
	temp1=${GFF%.gff} # remove file suffix
	sp=${temp1##*/} # remove path
	echo ${sp}
	GENOME=${GEMOMA_REFS}/${sp}.fasta
	gt gff3 -tidy -o ${sp}.clean.gff ${GFF} # clean up GFF file, make sure it's GFF3 format
	echo "s=own \\" >> run_GeMoMa.combined.sh
	echo "i=${sp} \\" >> run_GeMoMa.combined.sh
	echo "a=${sp}.clean.gff \\" >> run_GeMoMa.combined.sh
	echo "g=${GENOME} \\" >> run_GeMoMa.combined.sh
done

echo "t=${ANNOTATION_DIR}/${GENOME_DIR}/${MASKED_GENOME_FILE};" >> run_GeMoMa.combined.sh

# run GeMoMa
conda activate GeMoMa_1.9
bash run_GeMoMa.combined.sh

# generate protien predictions
GeMoMa Extractor -XX${GEMOMA_RAM} \
 Ambiguity=AMBIGUOUS \
 p=true \
 outdir=${outDir_combined} \
 g=${ANNOTATION_DIR}/${GENOME_DIR}/${MASKED_GENOME_FILE} \
 a=${outDir_combined}/final_annotation.gff;

mv ${outDir_combined}/proteins_1.fasta ${outDir_combined}/proteins.all.fasta

# get longest isoform of each gene
AGAT_SIF="/projects/academic/tkrabben/software/agat/agat_1.0.0--pl5321hdfd78af_0.sif"
singularity run -H ${PWD} ${AGAT_SIF} agat_sp_keep_longest_isoform.pl --gff ${outDir_combined}/final_annotation.gff -o ${outDir_combined}/final_annotation.longest_isoform.gff

# generate protein predictions for longest isoform
GeMoMa Extractor -XX${GEMOMA_RAM} \
 Ambiguity=AMBIGUOUS \
 p=true \
 outdir=${outDir_combined} \
 g=${ANNOTATION_DIR}/${GENOME_DIR}/${MASKED_GENOME_FILE} \
 a=${outDir_combined}/final_annotation.longest_isoform.gff;

mv ${outDir_combined}/proteins_2.fasta ${outDir_combined}/proteins.longest_isoform.fasta
