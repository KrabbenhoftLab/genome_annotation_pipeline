#!/bin/bash
# v0.3.1
# Pipeline to perform gene prediction and annotation
# author: Dan MacGuigan

# RUN REPEAT MASKER

# input variables passed from AISO_annotation_pipeline.sh
SPECIES=$1 # short name for your species
GENOME_DIR=$2 # directory containing your genome assembly
GENOME_FILE=$3 # your genome assembly
REPEAT_LIBRARY_NAME=$4 # name for your repeat library
RM_THREADS=$5 # number of threads to parallelize RepeatMasker
RM_SPECIES=$6 # species or clade to use for repeat masking
              # to find a species or clade of interest, run the following three commands
			  # 'source activate maker-3.01.03'
			  # 'cd /home/krablab/Documents/apps/RepeatMasker/Libraries"
			  # 'famdb.py -i Dfam.h5 names YOUR_SPECIES_OR_CLADE'
			  # this should return a list of potential matches for you to choose from
BLAST_CPUS=$7 # for use with steps 1B and 12

# create directory for RepeatMasker
mkdir ${SPECIES}_RepeatMasker
cd ${SPECIES}_RepeatMasker

# check for coding sequences
# download the uniprot database
wget https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz
gunzip uniprot_sprot.fasta.gz
# filter out transposable element from uniprot database
/home/krablab/Documents/apps/bbmap/filterbyname.sh in=uniprot_sprot.fasta out=uniprot_sprot.TEfiltered.fasta names=transposable,transposase substring=name ignorejunk
# blast against uniprot
# need to use an older version of BLAST for compatibility with ProtExcluder
/home/krablab/Documents/apps/ncbi-blast-2.4.0+/bin/makeblastdb -in uniprot_sprot.TEfiltered.fasta -dbtype prot
/home/krablab/Documents/apps/ncbi-blast-2.4.0+/bin/blastx -db uniprot_sprot.TEfiltered.fasta -query ../${SPECIES}_RepeatModeler/${REPEAT_LIBRARY_NAME}-families.fa -num_threads ${BLAST_CPUS} -evalue 10 -out ${SPECIES}_blastx_uniprot_TEfiltered_e10.out
# filter out sequences with blast hit, excluding the 50 bp flanking the matched region on each side
# note that ProtExcluder requires "esl-sfetch" from the Easel library
# on Bottlerocket, located at "/home/krablab/miniconda2/envs/hmmer/bin/esl-sfetch"
/home/krablab/Documents/apps/ProtExcluder/ProtExcluder.pl ${SPECIES}_blastx_uniprot_TEfiltered_e10.out ../${SPECIES}_RepeatModeler/${REPEAT_LIBRARY_NAME}-families.fa

# Run RepeatMasker to get GFF
# use coding gene filtered denovo predicted repeats
mkdir RMask_denovoPrediction_protFiltered
RepeatMasker -pa ${RM_THREADS} -gff -lib ../${SPECIES}_RepeatModeler/${REPEAT_LIBRARY_NAME}-families.fanoProtFinal -dir ./RMask_denovoPrediction_protFiltered ../${GENOME_DIR}/${GENOME_FILE}

# run RepeatMasker again using a clade or species from the Dfam database
mkdir RMask_denovoPlusDfam
RepeatMasker -pa ${RM_THREADS} -gff -species ${RM_SPECIES} -dir ./RMask_denovoPlusDfam ./RMask_denovoPrediction_protFiltered/${GENOME_FILE}.masked

# Combine two rounds of RepeatMasker
mkdir final_repeat_mask
gunzip ./RMask_denovoPrediction_protFiltered/*.cat.gz ./RMask_denovoPlusDfam/*.cat.gz
cat  ./RMask_denovoPrediction_protFiltered/*.cat  ./RMask_denovoPlusDfam/*.cat > final_repeat_mask/${SPECIES}.final_repeat_mask.cat
gzip ./RMask_denovoPrediction_protFiltered/*.cat
gzip ./RMask_denovoPlusDfam/*.cat
cd final_repeat_mask
ProcessRepeats -species ${RM_SPECIES} ${SPECIES}.final_repeat_mask.cat

# create GFF
/home/krablab/Documents/apps/RepeatMasker/util/rmOutToGFF3.pl ${SPECIES}.final_repeat_mask.out > ${SPECIES}.final_repeat_mask.gff3
# isolate complex repeats
grep -v -e "Satellite" -e ")n" -e "-rich" ${SPECIES}.final_repeat_mask.gff3 > ${SPECIES}.final_repeat_mask.complex.gff3
# reformat to work with MAKER
cat ${SPECIES}.final_repeat_mask.complex.gff3 | \
  perl -ane '$id; if(!/^\#/){@F = split(/\t/, $_); chomp $F[-1];$id++; $F[-1] .= "\;ID=$id"; $_ = join("\t", @F)."\n"} print $_' \
  > ${SPECIES}.final_repeat_mask.complex.reformat.gff3
