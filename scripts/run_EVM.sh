#!/bin/bash
# v0.4.0
#SBATCH --qos=general-compute
#SBATCH --partition=general-compute
#SBATCH --account=tkrabben
#SBATCH --time=72:00:00
#SBATCH --nodes=1
#SBATCH --reservation=ubhpc-future

# Pipeline to perform gene prediction and annotation
# author: Dan MacGuigan

# RUN EVM 
# to combine annotations from Braker and GeMoMa

# input variables passed from AISO_annotation_pipeline.sh
SPECIES=$1 # short name for your species
GENOME_DIR=$2 # directory containing your genome assembly
MASKED_GENOME_FILE=$3 # your (soft) masked genome assembly
ANNOTATION_DIR=$4
BRAKER_WEIGHT=$5 # EVM weight for BRAKER predictions
GEMOMA_WEIGHT=$6 # EVM weight for GeMoMa predictions
EVM_THREADS=$7 # threads for EVM
GEMOMA_SCORE_AA_FILTER=$8 # which GeMoMa filtered dataset to use?

EVM_SIF="/projects/academic/tkrabben/software/EvidenceModeler/EVidenceModeler/Docker/EVidenceModeler.latest.simg"

cd ${ANNOTATION_DIR}

# combine predictions into one file
cat ${SPECIES}_EVM/B${BRAKER_WEIGHT}_G${GEMOMA_WEIGHT}-SCORE-${GEMOMA_SCORE_AA_FILTER}/evidence/*EVM.gff3 > ${SPECIES}_EVM/B${BRAKER_WEIGHT}_G${GEMOMA_WEIGHT}-SCORE-${GEMOMA_SCORE_AA_FILTER}/evidence/combined.gff3

# create evidence weights file
cd ${SPECIES}_EVM/B${BRAKER_WEIGHT}_G${GEMOMA_WEIGHT}-SCORE-${GEMOMA_SCORE_AA_FILTER}
echo "ABINITIO_PREDICTION	gmst	${BRAKER_WEIGHT}" > weights.txt
echo "ABINITIO_PREDICTION	AUGUSTUS	${BRAKER_WEIGHT}" >> weights.txt
echo "ABINITIO_PREDICTION	GeneMark.hmm3	${BRAKER_WEIGHT}" >> weights.txt
echo "OTHER_PREDICTION	GeMoMa	${GEMOMA_WEIGHT}" >> weights.txt

# run evidence modeler
echo "starting EVM run..."
singularity exec -H ${ANNOTATION_DIR} ${EVM_SIF} EVidenceModeler \
                   --sample_id ${SPECIES} \
                   --genome ${ANNOTATION_DIR}/${GENOME_DIR}/${MASKED_GENOME_FILE} \
                   --weights ./weights.txt \
                   --gene_predictions ./evidence/combined.gff3 \
                   --segmentSize 1000000 \
                   --overlapSize 100000 \
                   --CPU ${EVM_THREADS}

# get summary info for the transcripts
# from https://github.com/darencard/GenomeAnnotation
echo "annotation summary stats written to ${SPECIES}.EVM.summaryStats.txt"
echo "see ${SPECIES}.EVM.genestats.txt for list of stats by gene"
module load gcc/11.2.0
module load samtools/1.16.1
module load bedtools/2.30.0
module load bcftools/1.14
sed '/^$/d' ${SPECIES}.EVM.gff3 > ${SPECIES}.EVM.mod.gff3 # need to remove empty lines for tabix
echo -e 'transcript_ID\ttranscript_length\texons\ttotal_exon_length\tintrons\ttotal_intron_length\tCDS\ttotal_CDS_length' > ${SPECIES}.EVM.genestats.txt
/projects/academic/tkrabben/software/genestats/genestats ${SPECIES}.EVM.mod.gff3 >> ${SPECIES}.EVM.genestats.txt
echo "test"
cat ${SPECIES}.EVM.pep | awk '$0 ~ ">" {if (NR > 1) {print c;} c=0;printf substr($0,2,100) "\t"; } $0 !~ ">" {c+=length($0);} END { print c; }' > ${SPECIES}.EVM.pep.lens.txt
echo "test2"
rm ${SPECIES}.EVM.summaryStats.txt
module load miniconda3/22.11.1-1;source activate r_env
Rscript --vanilla ${ANNOTATION_DIR}/scripts/annotation_summary_stats.R ${SPECIES}.EVM.mod.gff3 introns.txt ${SPECIES}.EVM.genestats.txt ${SPECIES}.EVM.pep.lens.txt ${SPECIES}

