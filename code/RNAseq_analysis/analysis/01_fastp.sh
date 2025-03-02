#!/bin/bash

# Path of OT2_yeast
path_ot2_yeast="/mnt/d/200_GitHub_Repository/OT2_yeast"

# Path of large data
path_large="/mnt/k/20240325_pipspd2nd"





# Path of metadata TSV
path_metadata="$path_ot2_yeast/code/RNAseq_analysis/analysis_2nd/metadata.tsv"


# Path of input (separated by *R1.fastq.gz and *R2.fastq.gz)
path_indir="$path_large/Merged_FASTQ"


# Path of trimmomatic output directory
path_outdir="$path_large/fastp"



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
	mkdir -p ${path_outdir}/report_htmls
	mkdir -p ${path_outdir}/report_jsons
fi




# Make log file
datetime=$(date '+%Y%m%d-%H%M%S')
path_logfile=${path_outdir}/log_fastp_${datetime}.log
touch $path_logfile


###############################################

# Activate labopt virtual env
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
		message="[${datetime}] Start fastp ${Sample_ID}..."
		echo $message >> ${path_logfile}
		echo $message



		# # Log command
		echo "fastp -i ${path_indir}/${Sample_ID}.fastq.gz -I ${path_indir}/${Sample_ID}_R2.fastq.gz -3 -o ${path_outdir}/${Sample_ID}_R1_fastp.fastq.gz -O ${path_outdir}/${Sample_ID}_R2_fastp.fastq.gz -h ${path_outdir}/report_htmls/${Sample_ID}_report_fastp.html -j ${path_outdir}/report_jsons/${Sample_ID}_report_fastp.json -q 15 -n 10 -t 1 -T 1 -l 20 -w 16 -A" >> ${path_logfile}
		echo "fastp -i ${path_indir}/${Sample_ID}.fastq.gz -I ${path_indir}/${Sample_ID}_R2.fastq.gz -3 -o ${path_outdir}/${Sample_ID}_R1_fastp.fastq.gz -O ${path_outdir}/${Sample_ID}_R2_fastp.fastq.gz -h ${path_outdir}/report_htmls/${Sample_ID}_report_fastp.html -j ${path_outdir}/report_jsons/${Sample_ID}_report_fastp.json -q 15 -n 10 -t 1 -T 1 -l 20 -w 16 -A"


		# Fastp
        ## parameters are the same as this blog
		## https://kazumaxneo.hatenablog.com/entry/2018/05/21/111947
		fastp -i ${path_indir}/${Sample_ID}_R1.fastq.gz \
        -I ${path_indir}/${Sample_ID}_R2.fastq.gz -3 \
        -o ${path_outdir}/${Sample_ID}_R1_fastp_paired.fastq.gz \
        -O ${path_outdir}/${Sample_ID}_R2_fastp_paired.fastq.gz \
        -h ${path_outdir}/report_htmls/${Sample_ID}_report_fastp.html \
		-j ${path_outdir}/report_jsons/${Sample_ID}_report_fastp.json \
		-q 15 -n 10 -t 1 -T 1 -l 20 -w 16 -A >> ${path_logfile}
		

		# Log end time
		datetime=$(date '+%Y%m%d-%H%M%S')
		message="[${datetime}] Finished fastp ${Sample_ID}!"
		echo $message >> ${path_logfile}
		echo $message
	fi
	COUNT=$(( COUNT + 1 ))
done < ${path_metadata}