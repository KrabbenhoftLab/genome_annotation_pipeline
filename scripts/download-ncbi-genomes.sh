#!/bin/bash -l
# v0.4.3

#SBATCH --qos=general-compute
#SBATCH --partition=general-compute
#SBATCH --account=tkrabben
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --mem=100G
##SBATCH --constraint=AVX512
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

module load miniconda3/22.11.1-1
source activate ncbi_datasets

# make sure the input file is correctly formatted, two columns and tab delimited
bad=$(awk -F "\t" 'NF != 2' ${NCBI_TARGET_GENOMES})

if [ -z "$bad" ]
then
    echo "${NCBI_TARGET_GENOMES} is properly formatted, downloading now"
else
    echo "${NCBI_TARGET_GENOMES} is not formatted correctly"
    echo "check to make sure the file contains two tab-delimited columns"
    echo "your columns might be space-delimited"
    exit 1
fi


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
