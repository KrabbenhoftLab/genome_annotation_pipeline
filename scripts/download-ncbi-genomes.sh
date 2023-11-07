#!/bin/bash
# v0.2.5
#SBATCH --qos=general-compute
#SBATCH --partition=general-compute
#SBATCH --account=tkrabben
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --mem=100G
#SBATCH --constraint=AVX512
#SBATCH --export=NONE

# Pipeline to perform gene prediction and annotation
# author: Dan MacGuigan

# DOWNLOAD REFERENCE DATA FROM NCBI

# input variables passed from AISO_annotation_pipeline.sh
NCBI_TARGET_GENOMES=$1
NCBI_DOWNLOAD_DIR=$2
ANNOTATION_DIR=$3

# create download dir
mkdir ${ANNOTATION_DIR}/${NCBI_DOWNLOAD_DIR}
cd ${ANNOTATION_DIR}/${NCBI_DOWNLOAD_DIR}

module load miniconda3
source activate ncbi_datasets

while read p
do
  IFS=$'\t'
  tmp=($p)
  mkdir -p ${tmp[0]}
  cd ${tmp[0]}

  # download from NCBI and unzip
  datasets download genome accession ${tmp[0]} --filename ${tmp[0]}.zip --include genome,gff3,protein
  unzip ${tmp[0]}.zip

  # make copies of fastas and gffs with new names
  cp ncbi_dataset/data/${tmp[0]}/*genomic.fna ../${tmp[1]}.fasta
  cp ncbi_dataset/data/${tmp[0]}/*.gff ../${tmp[1]}.gff
  cp ncbi_dataset/data/${tmp[0]}/protein.faa ../${tmp[1]}.protein.fasta
  
  cd ..
  rm -rf ${tmp[0]}
done <${NCBI_TARGET_GENOMES}
