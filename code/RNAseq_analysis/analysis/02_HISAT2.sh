#!/bin/bash

# Path of OT2_yeast
path_ot2_yeast="/mnt/d/200_GitHub_Repository/OT2_yeast"

# Path of large data
path_large="/mnt/k/20240325_pipspd2nd"

# Path of indexed for W303 exon
path_index="$path_large/HISAT2/index/index"

# Path of metadata TSV
path_metadata="$path_ot2_yeast/code/RNAseq_analysis/analysis_2nd/metadata.tsv"

# Path of input for STAR (separated by *R1.fastq.gz and *R2.fastq.gz)
path_indir="$path_large/fastp"

# Path of GTF file
path_genemap="$path_large/reference_genome_W303/W303_JRIU00000000_SGD.gff"

# Path of HISAT2 output directory
path_outdir="$path_large/HISAT2"



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

# Check whether GTF file exists
if [ ! -e ${path_genemap} ]; then
	echo "${path_genemap} is missing!"
	exit
fi

# Make output directory
if [ ! -d ${path_outdir} ]; then
	mkdir -p ${path_outdir}
fi

# Make log file
datetime=$(date '+%Y%m%d-%H%M%S')
path_logfile=${path_outdir}/log_star_${datetime}.log
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
		message="[${datetime}] Start HISAT2 ${Sample_ID}..."
		echo $message >> ${path_logfile}
		echo $message

		# # Log command
        echo "hisat2 -q -p 16 --dta -x ${path_index} \
        -1 ${path_indir}/${Sample_ID}_R1_fastp_paired.fastq \
        -2 ${path_indir}/${Sample_ID}_R2_fastp_paired.fastq \
        -S ${path_outdir}/${Sample_ID}.sam"

        echo "hisat2 -q -p 16 --dta -x ${path_index} \
        -1 ${path_indir}/${Sample_ID}_R1_fastp_paired.fastq \
        -2 ${path_indir}/${Sample_ID}_R2_fastp_paired.fastq \
        -S ${path_outdir}/${Sample_ID}.sam" >> ${path_logfile}

		#HISAT2       
        hisat2 -q -p 16 --dta -x ${path_index} \
        -1 ${path_indir}/${Sample_ID}_R1_fastp_paired.fastq \
        -2 ${path_indir}/${Sample_ID}_R2_fastp_paired.fastq \
        -S ${path_outdir}/${Sample_ID}.sam

        # Log end time
		datetime=$(date '+%Y%m%d-%H%M%S')
		message="[${datetime}] Finished STAR ${Sample_ID}!"
		echo $message >> ${path_logfile}
		echo $message
	fi
	COUNT=$(( COUNT + 1 ))
done < ${path_metadata}
