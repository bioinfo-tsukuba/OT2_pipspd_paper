#!/bin/bash

# Path of OT2_yeast
path_ot2_yeast="/mnt/d/200_GitHub_Repository/OT2_yeast"

# Path of large data
path_large="/mnt/k/20240325_pipspd2nd"

# Path of RSEM reference for W303 CDS
path_ref="$path_large/RSEM/reference_W303_SGD_2015_JRIU00000000"

# Path of STAR index
path_index = "$path_large/STAR/index"

# Path of metadata TSV
path_metadata="$path_ot2_yeast/code/RNAseq_analysis/analysis_2nd/metadata.tsv"

# Path of GTF file
path_genemap="$path_large/reference_genome_W303//transcripts.gtf"

# Path of fasta file
path_fasta = "$path_large/reference_genome_W303/W303_SGD_2015_JRIU00000000.fsa"

# Path of input for star (dir S01, S02, ...)
path_indir="$path_large/STAR"

# Path of trimmomatic output directory
path_outdir="$path_large/RSEM"



###############################################
# Check whether metadata TSV exists
if [ ! -e ${path_metadata} ]; then
	echo "$path_metadata is missing!"
	exit
fi

# Check whether indir exists
if [ ! -d ${path_indir} ]; then
	echo "${path_indir} is missing!"
	exit
fi


# Make output directory
if [ ! -d ${path_outdir} ]; then
	mkdir -p ${path_outdir}
fi


# Make log file
datetime=$(date '+%Y%m%d-%H%M%S')
path_logfile=${path_large}/log_rsem_${datetime}.log
touch $path_logfile


###############################################


# Activate salmon virtual env
conda activate labopt


cd ${path_indir}

COUNT=0
while read LINE; do
	if [ $COUNT -gt 0 ]; then
		list=(`echo "$LINE"`)

		# Get Fastq_prefix and Sample_ID
        Fastq_prefix=${list[0]}
		Sample_ID=${list[1]}

		# Log start time
		datetime=$(date '+%Y%m%d-%H%M%S')
		message="[${datetime}] Start rsem ${Sample_ID}..."
		echo $message >> ${path_logfile}
		echo $message

		# # Log command
        echo "rsem-calculate-expression --alignments --paired-end --num-threads 16 \
        --strandedness reverse --append-names --estimate-rspd \
        ${path_large}/STAR/${Sample_ID}/Aligned.toTranscriptome.out.bam \
		${path_ref} \
		${path_outdir}/${Sample_ID} \
		--seed 12345" >> ${path_logfile}

 		  echo "rsem-calculate-expression --alignments --paired-end --num-threads 16 \
        --strandedness reverse --append-names --estimate-rspd \
        ${path_large}/STAR/${Sample_ID}/Aligned.toTranscriptome.out.bam \
		${path_ref} \
		${path_outdir}/${Sample_ID} \
		--seed 12345"
		
		# RSEM
        rsem-calculate-expression --alignments --paired-end --num-threads 16 \
        --strandedness reverse --append-names --estimate-rspd \
        ${path_large}/STAR/${Sample_ID}/Aligned.toTranscriptome.out.bam \
		${path_ref} \
		${path_outdir}/${Sample_ID} \
		--seed 12345
		


		# Log end time
		datetime=$(date '+%Y%m%d-%H%M%S')
		message="[${datetime}] Finished rsem ${Sample_ID}!"
		echo $message >> ${path_logfile}
		echo $message
	fi
	COUNT=$(( COUNT + 1 ))
done < ${path_metadata}
