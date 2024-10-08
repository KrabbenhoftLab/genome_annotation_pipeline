#!/bin/bash -l
# v0.4.7

#SBATCH --qos=general-compute
#SBATCH --partition=general-compute
#SBATCH --account=tkrabben
#SBATCH --time=72:00:00
#SBATCH --nodes=1


# Pipeline to perform gene prediction and annotation
# author: Dan MacGuigan

# RUN BRAKER

# input variables passed from AISO_annotation_pipeline.sh
SPECIES=$1 # short name for your species
GENOME_DIR=$2 # directory containing your genome assembly
MASKED_GENOME_FILE=$3 # your (soft) masked genome assembly
PROT_FASTA=$4
BRAKER_THREADS=$5
ANNOTATION_DIR=$6
AUGUSTUS_SPECIES_NAME=$7

cd ${ANNOTATION_DIR}

BRAKER_SIF="/projects/academic/tkrabben/software/BRAKER3/braker3.0.6.sif" # location of BRAKER singularity image

# create directory for BRAKER
mkdir ${SPECIES}_BRAKER

# workaround to make Singularity function properly
mkdir .augustus
echo "cp -r /usr/share/augustus/config /mnt/.augustus" > temp.sh
singularity exec -B $PWD:/mnt ${BRAKER_SIF} /bin/bash /mnt/temp.sh
rm temp.sh

#from https://github.com/nf-core/chipseq/issues/123
export TMPDIR=/tmp/

if [ -f "${ANNOTATION_DIR}/${SPECIES}_HISAT2/${SPECIES}.sorted.rna.bam" ]; then  # if there is an RNA-seq BAM, supply RNA-seq evidence to BRAKER

    echo "starting BRAKER with RNA-seq evidence and protein evidence, sit tight"
    singularity exec -H ${PWD} -B ${PWD}:/mnt ${BRAKER_SIF} braker.pl --genome=${GENOME_DIR}/${MASKED_GENOME_FILE} --bam=${SPECIES}_HISAT2/${SPECIES}.sorted.rna.bam --prot_seq=${PROT_FASTA} --species=${AUGUSTUS_SPECIES_NAME} --workingdir=${SPECIES}_BRAKER --threads ${BRAKER_THREADS} --gff3 --AUGUSTUS_CONFIG_PATH=/mnt/.augustus/config

    # create annotation keeping all genes predicted by Augustus

    # need to make a TSEBRA config file, taken from their GitHub
    # could try tinkering with this sometime
    echo "# Weight for each hint source
# Values have to be >= 0
P 1
E 20
C 1
M 1
# Required fraction of supported introns or supported start/stop-codons for a transcript
# Values have to be in [0,1]
intron_support 1.0
stasto_support 2
# Allowed difference for each feature 
# Values have to be in [0,1]
e_1 0.1
e_2 0.5
e_3 0.05
e_4 0.18" > ${SPECIES}_BRAKER/default.cfg

    singularity exec -H ${PWD} ${BRAKER_SIF} tsebra.py -g ${SPECIES}_BRAKER/braker.gtf -k ${SPECIES}_BRAKER/Augustus/augustus.hints.gtf -e ${SPECIES}_BRAKER/hintsfile.gff -c ${SPECIES}_BRAKER/default.cfg -o ${SPECIES}_BRAKER/braker.allAugustus.gtf

else # if there is no RNA-seq BAM, run BRAKER only with protein evidence

    echo "starting BRAKER with only protein evidence, sit tight"
    singularity exec -H ${PWD} -B ${PWD}:/mnt ${BRAKER_SIF} braker.pl --genome=${GENOME_DIR}/${MASKED_GENOME_FILE} --prot_seq=${PROT_FASTA} --species=${AUGUSTUS_SPECIES_NAME} --workingdir=${SPECIES}_BRAKER --threads ${BRAKER_THREADS} --gff3 --AUGUSTUS_CONFIG_PATH=/mnt/.augustus/config
    
fi


# GET LONGEST ISOFORM OF EACH GENE FOR BUSCO ANALYSES WITH GVOLANTE
AGAT_SIF="/projects/academic/tkrabben/software/agat/agat_1.0.0--pl5321hdfd78af_0.sif"
singularity exec -H ${PWD} ${AGAT_SIF} agat_sp_keep_longest_isoform.pl --gff ${SPECIES}_BRAKER/braker.gff3 -o ${SPECIES}_BRAKER/braker.longest_isoform.gff3
singularity exec -H ${PWD} ${AGAT_SIF} agat_sp_extract_sequences.pl -g ${SPECIES}_BRAKER/braker.longest_isoform.gff3 -f ${GENOME_DIR}/${MASKED_GENOME_FILE} -t cds -p -o ${SPECIES}_BRAKER/braker.longest_isoform.aa 


echo ""
echo "Step 6 COMPLETE"
echo "Please check to make sure that a combined GFF was produced by BRAKER:"
echo "${ANNOTATION_DIR_CLUSTER}/${SPECIES}_BRAKER/braker.longest_isoform.gff3"
echo ""
echo "If Step 7 (GeMoMa) has also finished, you may proceed to Step 8: combining gene predictions with EVM."
