#!/bin/bash -l
# v0.4.2
#SBATCH --qos=general-compute
#SBATCH --partition=general-compute
#SBATCH --account=tkrabben
#SBATCH --time=72:00:00
#SBATCH --nodes=1
##SBATCH --constraint=AVX512
#SBATCH --export=NONE


# Pipeline to perform gene prediction and annotation
# author: Dan MacGuigan

# RUN REPEAT MASKER

# input variables passed from AISO_annotation_pipeline.sh
SPECIES=$1 # short name for your species
ANNOTATION_DIR_CLUSTER=$2
GENOME_DIR=$3 # directory containing your genome assembly
GENOME_FILE=$4 # your genome assembly
REPEAT_LIBRARY_NAME=$5 # name for your repeat library
RMASK_THREADS=$6 # number of threads to parallelize RepeatMasker
RM_SPECIES=$7 # species or clade to use for repeat masking
              # to find a species or clade of interest, run the following three commands
			  # 'source activate maker-3.01.03'
			  # 'cd /home/krablab/Documents/apps/RepeatMasker/Libraries"
			  # 'famdb.py -i Dfam.h5 names YOUR_SPECIES_OR_CLADE'
			  # this should return a list of potential matches for you to choose from

GENOME_FILE_PREFIX=${GENOME_FILE%".fasta"}


BLAST_CPUS=${RMASK_THREADS}

## RMASK_THREADS specifies the number of parallel search jobs to run for RepeatMasker.
## RMBlast jobs will use 4 cores each thread, so we need to do some math
RM_THREADS=$(( RMASK_THREADS / 4 )) 

echo "RM_THREADS = ${RM_THREADS}"

# create directory for RepeatMasker
cd ${ANNOTATION_DIR_CLUSTER}
mkdir ${SPECIES}_RepeatMasker
cd ${SPECIES}_RepeatMasker


if ! [ -f ./uniprot_sprot.TEfiltered.fasta ]; then
	# check for coding sequences
	# download the uniprot database
	wget https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz
	gunzip uniprot_sprot.fasta.gz
	# filter out transposable element from uniprot database
	# load bbmap module
	module load gcc/11.2.0 openmpi/4.1.1 bbmap/38.98
	filterbyname.sh in=uniprot_sprot.fasta out=uniprot_sprot.TEfiltered.fasta names=transposable,transposase substring=name ignorejunk
fi

if ! [ -f ./${SPECIES}_blastx_uniprot_TEfiltered_e10.out ]; then
	# blast against uniprot
	# need to use an older version of BLAST for compatibility with ProtExcluder
	/projects/academic/tkrabben/software/ncbi-blast-2.4.0+/bin/makeblastdb -in uniprot_sprot.TEfiltered.fasta -dbtype prot
	/projects/academic/tkrabben/software/ncbi-blast-2.4.0+/bin/blastx -db uniprot_sprot.TEfiltered.fasta -query ${ANNOTATION_DIR_CLUSTER}/${SPECIES}_RepeatModeler/${REPEAT_LIBRARY_NAME}-families.fa -num_threads ${BLAST_CPUS} -evalue 10 -out ${SPECIES}_blastx_uniprot_TEfiltered_e10.out
fi

# load repeatmasker module
module purge
module load gcc/11.2.0 openmpi/4.1.1 repeatmodeler/2.0.4.KRAB

if ! [ -f ${ANNOTATION_DIR_CLUSTER}/${SPECIES}_RepeatModeler/${REPEAT_LIBRARY_NAME}-families.fanoProtFinal ]; then
	# filter out sequences with blast hit, excluding the 50 bp flanking the matched region on each side
	# note that ProtExcluder requires "esl-sfetch" from the Easel library
	# on Bottlerocket, located at "/home/krablab/miniconda2/envs/hmmer/bin/esl-sfetch"
	/projects/academic/tkrabben/software/ProtExcluder/ProtExcluder.pl ${SPECIES}_blastx_uniprot_TEfiltered_e10.out ${ANNOTATION_DIR_CLUSTER}/${SPECIES}_RepeatModeler/${REPEAT_LIBRARY_NAME}-families.fa
fi

if ! [ -f ./RMask_denovoPrediction_protFiltered/*.rmalign.gz ]; then
	# Run RepeatMasker to get GFF
	# use coding gene filtered denovo predicted repeats
	mkdir RMask_denovoPrediction_protFiltered
	
	module load nextflow/21.10.6
	
	nextflow run /projects/academic/tkrabben/software/RepeatMasker_Nextflow/RepeatMasker_Nextflow.nf \
		--inputSequence ${ANNOTATION_DIR_CLUSTER}/${GENOME_DIR}/${GENOME_FILE} \
		--inputLibrary ${ANNOTATION_DIR_CLUSTER}/${SPECIES}_RepeatModeler/${REPEAT_LIBRARY_NAME}-families.fanoProtFinal \
		--cluster UB \
		--outputDir ${ANNOTATION_DIR_CLUSTER}/${SPECIES}_RepeatMasker/RMask_denovoPrediction_protFiltered \
		-resume;
fi

if ! [ -f ./RMask_denovoPrediction_protFiltered/*.masked ]; then
	# generate masked fasta
	cd RMask_denovoPrediction_protFiltered
	cp ${ANNOTATION_DIR_CLUSTER}/${GENOME_DIR}/${GENOME_FILE} .
	gunzip ${GENOME_FILE_PREFIX}.rmout.gz
	echo "masking genome assembly ${GENOME_FILE_PREFIX}.fasta"
	/projects/academic/tkrabben/modules_KrabLab/easybuild/2023.01/software/avx512/MPI/gcc/11.2.0/openmpi/4.1.1/repeatmasker/4.1.5/util/maskFile.pl -fasta ${GENOME_FILE} -annotations ${GENOME_FILE_PREFIX}.rmout
	echo "done masking genome"
	gzip ${GENOME_FILE_PREFIX}.rmout
	cd ..
	#RepeatMasker -pa ${RM_THREADS} -gff -lib ${ANNOTATION_DIR_CLUSTER}/${SPECIES}_RepeatModeler/${REPEAT_LIBRARY_NAME}-families.fanoProtFinal -dir ./RMask_denovoPrediction_protFiltered ${ANNOTATION_DIR_CLUSTER}/${GENOME_DIR}/${GENOME_FILE}
fi

if ! [ -f ./RMask_denovoPlusDfam/*.rmalign.gz ]; then
	# run RepeatMasker again using a clade or species from the Dfam database
	mkdir RMask_denovoPlusDfam
	
	module load nextflow/21.10.6
	
	# mask using Dfam, use masked FASTA generated using custom repeat library
	nextflow run /projects/academic/tkrabben/software/RepeatMasker_Nextflow/RepeatMasker_Nextflow.nf \
		--inputSequence ${ANNOTATION_DIR_CLUSTER}/${SPECIES}_RepeatMasker/RMask_denovoPrediction_protFiltered/${GENOME_FILE}.masked \
		--cluster UB \
		--species "${RM_SPECIES}" \
		--outputDir ${ANNOTATION_DIR_CLUSTER}/${SPECIES}_RepeatMasker/RMask_denovoPlusDfam \
		-resume;
	
	#RepeatMasker -pa ${RM_THREADS} -gff -species ${RM_SPECIES} -dir ./RMask_denovoPlusDfam ./RMask_denovoPrediction_protFiltered/${GENOME_FILE}.masked
fi

if ! [ -f ./final_repeat_mask/*.tbl ]; then
	cd ${ANNOTATION_DIR_CLUSTER}/${SPECIES}_RepeatMasker

	# Combine two rounds of RepeatMasker
	mkdir final_repeat_mask
	
	# prepare files
	gunzip ./RMask_denovoPlusDfam/${GENOME_FILE}.rmalign.gz 
	mv ./RMask_denovoPlusDfam/${GENOME_FILE}.rmalign ./RMask_denovoPlusDfam/${GENOME_FILE}.align
	
	gunzip ./RMask_denovoPlusDfam/${GENOME_FILE}.rmout.gz
	mv ./RMask_denovoPlusDfam/${GENOME_FILE}.rmout ./RMask_denovoPlusDfam/${GENOME_FILE}.out
	
	gunzip ./RMask_denovoPrediction_protFiltered/${GENOME_FILE_PREFIX}.rmalign.gz
	mv ./RMask_denovoPrediction_protFiltered/${GENOME_FILE_PREFIX}.rmalign ./RMask_denovoPrediction_protFiltered/${GENOME_FILE_PREFIX}.align 
	
	gunzip ./RMask_denovoPrediction_protFiltered/${GENOME_FILE_PREFIX}.rmout.gz
	mv ./RMask_denovoPrediction_protFiltered/${GENOME_FILE_PREFIX}.rmout ./RMask_denovoPrediction_protFiltered/${GENOME_FILE_PREFIX}.out
	
	# combine results of the two runs
	/projects/academic/tkrabben/modules_KrabLab/easybuild/2023.01/software/avx512/MPI/gcc/11.2.0/openmpi/4.1.1/repeatmasker/4.1.5/util/combineRMFiles.pl ./RMask_denovoPlusDfam/${GENOME_FILE} ./RMask_denovoPrediction_protFiltered/${GENOME_FILE_PREFIX} ./final_repeat_mask/${SPECIES}.final_repeat_mask
	
	#cat  ./RMask_denovoPrediction_protFiltered/*.rmalign  ./RMask_denovoPlusDfam/*.rmalign > final_repeat_mask/${SPECIES}.final_repeat_mask.rmalign
	# compress and clean up
	gzip ./RMask_denovoPlusDfam/*.align
	mv ./RMask_denovoPlusDfam/${GENOME_FILE}.align.gz ./RMask_denovoPlusDfam/${GENOME_FILE}.rmalign.gz
	
	gzip ./RMask_denovoPlusDfam/*.out
	mv ./RMask_denovoPlusDfam/${GENOME_FILE}.out.gz ./RMask_denovoPlusDfam/${GENOME_FILE}.rmout.gz
	
	gzip ./RMask_denovoPrediction_protFiltered/*.align
	mv ./RMask_denovoPrediction_protFiltered/${GENOME_FILE_PREFIX}.align.gz ./RMask_denovoPrediction_protFiltered/${GENOME_FILE_PREFIX}.rmalign.gz
	
	gzip ./RMask_denovoPrediction_protFiltered/*.out
	mv ./RMask_denovoPrediction_protFiltered/${GENOME_FILE_PREFIX}.out.gz ./RMask_denovoPrediction_protFiltered/${GENOME_FILE_PREFIX}.rmout.gz
	
	# generate repeat GFF and summary table
	cd final_repeat_mask
	#ProcessRepeats -species ${RM_SPECIES} ${SPECIES}.final_repeat_mask.align
	gzip ${SPECIES}.final_repeat_mask.align
	
	# get lengths of sequences in genome fasta
	cat ${ANNOTATION_DIR_CLUSTER}/${GENOME_DIR}/${GENOME_FILE} | awk '$0 ~ ">" {if (NR > 1) {print c;} c=0;printf substr($0,2,100) "\t"; } $0 !~ ">" {c+=length($0);} END { print c; }' > ${ANNOTATION_DIR_CLUSTER}/${GENOME_DIR}/${GENOME_FILE}.lens.tsv

	# generate detailed summary table
	/projects/academic/tkrabben/modules_KrabLab/easybuild/2023.01/software/avx512/MPI/gcc/11.2.0/openmpi/4.1.1/repeatmasker/4.1.5/util/buildSummary.pl -useAbsoluteGenomeSize -genome ${ANNOTATION_DIR_CLUSTER}/${GENOME_DIR}/${GENOME_FILE}.lens.tsv ${SPECIES}.final_repeat_mask.out > ${SPECIES}.final_repeat_mask.summary 
	
	# create GFF
	/projects/academic/tkrabben/modules_KrabLab/easybuild/2023.01/software/avx512/MPI/gcc/11.2.0/openmpi/4.1.1/repeatmasker/4.1.5/util/rmOutToGFF3.pl ${SPECIES}.final_repeat_mask.out > ${SPECIES}.final_repeat_mask.gff3
	# isolate complex repeats
	grep -v -e "Satellite" -e ")n" -e "-rich" ${SPECIES}.final_repeat_mask.gff3 > ${SPECIES}.final_repeat_mask.complex.gff3
	# reformat to work with MAKER
	cat ${SPECIES}.final_repeat_mask.complex.gff3 | \
	  perl -ane '$id; if(!/^\#/){@F = split(/\t/, $_); chomp $F[-1];$id++; $F[-1] .= "\;ID=$id"; $_ = join("\t", @F)."\n"} print $_' \
	  > ${SPECIES}.final_repeat_mask.complex.reformat.gff3
fi

echo ""
echo "Step 2 COMPLETE"
echo "please check to make sure that your final repeat GFF was generated"
echo "${ANNOTATION_DIR_CLUSTER}/${SPECIES}_RepeatMasker/final_repeat_mask/${SPECIES}.final_repeat_mask.gff3"
echo "you can inspect the repeat table to see a summary of the repeat masking:"
echo "${ANNOTATION_DIR_CLUSTER}/${SPECIES}_RepeatMasker/final_repeat_mask/${SPECIES}.final_repeat_mask.tbl"
echo ""
echo "if the repeat GFF does not exist, you can rerun Step 2"
echo "the script will attempt to resume the repeat masking run"
echo ""
echo "you may now proceed to Step 3: soft masking your genome FASTA"