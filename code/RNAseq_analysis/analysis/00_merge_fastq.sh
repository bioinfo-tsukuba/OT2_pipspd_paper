#!/bin/bash

# Path of OT2_yeast
path_ot2_yeast="/mnt/d/200_GitHub_Repository/OT2_yeast"

# Path of large data
path_large="/mnt/k/20240325_pipspd2nd"




# Path of metadata TSV
path_metadata="$path_ot2_yeast/code/RNAseq_analysis/analysis_2nd/metadata.tsv"

# Path of input fastq directory
path_indir="$path_large/LAB_465_FASTQ"

# Path of output directory
path_outdir="$path_large/Merged_FASTQ"

###############################################
# Check whether metadata TSV exists
if [ ! -e ${path_metadata} ]; then
	echo "${path_metadata} is missing!"
	exit
fi

# Check whether metadata TSV exists
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
path_logfile=log_merge_fastq_${datetime}.txt
touch $path_logfile

###############################################
COUNT=0
cd ${path_indir}
while read LINE; do
	if [ $COUNT -gt 0 ]; then
		list=(`echo "$LINE"`)

		# Get Fastq_prefix and Sample_ID
        Fastq_prefix=${list[0]}
		Sample_ID=${list[1]}

		for R in R1 R2; do

			# Define path of merged FASTQ file
			merged_fastq="${path_outdir}/${Sample_ID}_${R}.fastq.gz"

			# Remove existing files
			if [ -e ${merged_fastq} ]; then
				datetime=$(date '+%Y%m%d-%H%M%S')
				message="[${datetime}] Removing existing file: ${merged_fastq}"
				echo $message >> ${path_logfile}
				echo $message
				rm ${merged_fastq}
			fi

			# Log start time
			datetime=$(date '+%Y%m%d-%H%M%S')
			message="[${datetime}] Start generating ${merged_fastq}..."
			echo $message >> ${path_logfile}
			echo $message

			# Log command
			echo "cat ${Fastq_prefix}*_${R}_* >> ${merged_fastq}" >> ${path_logfile}
			echo "cat ${Fastq_prefix}*_${R}_* >> ${merged_fastq}"

			# Concatenate FASTQ files
			cat ${Fastq_prefix}*_${R}_* >> ${merged_fastq}

			# Log end time
			datetime=$(date '+%Y%m%d-%H%M%S')
			message="[${datetime}] Finished generating ${merged_fastq}!"
			echo $message >> ${path_logfile}
			echo $message

		done
	fi
	COUNT=$(( COUNT + 1 ))
done < ${path_metadata}


