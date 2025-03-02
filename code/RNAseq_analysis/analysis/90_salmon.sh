#!/bin/bash

# Path of OT2_yeast
path_ot2_yeast="/mnt/d/200_GitHub_Repository/OT2_yeast"

# Path of large data
path_large="/mnt/k/20240325_pipspd2nd"






# Path of metadata TSV
path_metadata="$path_ot2_yeast/code/RNAseq_analysis/analysis_2nd/metadata.tsv"

# Path of indexed for W303 exon
path_index="$path_large/salmon/index/"

# Path of input for salmon (separated by *R1.fastq.gz and *R2.fastq.gz)
path_indir="$path_large/fastp"

# Path of GTF file
path_genemap="$path_large/reference_genome_W303/W303_JRIU00000000_SGD.gff"

# Path of salmon output directory
path_outdir="$path_large/salmon"



###############################################
# Check whether metadata TSV exists
if [ ! -e ${path_metadata} ]; then
	echo "$path_metadata is missing!"
	exit
fi

# Check whether index exists
if [ ! -d ${path_index} ]; then
	echo "${path_index} is missing!"
	# Make index directory
	mkdir ${path_large}/salmon/index
    salmon index -p 4 \
    -t ${path_large}/reference_genome_W303/W303_JRIU00000000_SGD_orf_genomic.fsa \
    -i ${path_large}/salmon/index
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
path_logfile=${path_large}/log_merge_salmon_after_fastp_${datetime}.txt
touch $path_logfile


###############################################



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
		message="[${datetime}] Start salmon quant ${Sample_ID}..."
		echo $message >> ${path_logfile}
		echo $message

		# # Log command
		echo "salmon quant --index ${path_index} --libType ISR --mates1 ${path_indir}/${Sample_ID}_R1_fastp_paired.fastq --mates2 ${path_indir}/${Sample_ID}_R2_fastp_paired.fastq --writeMappings=${path_outdir}/samout.sam --writeQualities --geneMap ${path_genemap} --seqBias --gcBias --incompatPrior '0.0' --biasSpeedSamp '5' --fldMax '36' --fldMean '36' --fldSD '25' --forgettingFactor '0.65' --maxReadOcc '500'     --numBiasSamples '2000000' --numAuxModelSamples '5000000' --numPreAuxModelSamples '5000' --numGibbsSamples '0'  --numBootstraps '0'  --thinningFactor '16'  --sigDigits '3' --vbPrior '1e-05'   -o ${path_outdir}/${Sample_ID}" >> ${path_logfile}
		echo "salmon quant --index ${path_index} --libType ISR --mates1 ${path_indir}/${Sample_ID}_R1_fastp_paired.fastq --mates2 ${path_indir}/${Sample_ID}_R2_fastp_paired.fastq --writeMappings=${path_outdir}/samout.sam --writeQualities --geneMap ${path_genemap} --seqBias --gcBias --incompatPrior '0.0' --biasSpeedSamp '5' --fldMax '36' --fldMean '36' --fldSD '25' --forgettingFactor '0.65' --maxReadOcc '500'     --numBiasSamples '2000000' --numAuxModelSamples '5000000' --numPreAuxModelSamples '5000' --numGibbsSamples '0'  --numBootstraps '0'  --thinningFactor '16'  --sigDigits '3' --vbPrior '1e-05'   -o ${path_outdir}/${Sample_ID}"


		# Salmon quant
		salmon quant --index ${path_index} --libType ISR --mates1 ${path_indir}/${Sample_ID}_R1_fastp_paired.fastq \
		--mates2 ${path_indir}/${Sample_ID}_R2_fastp_paired.fastq \
		--writeMappings=${path_outdir}/samout.sam --writeQualities \
		--geneMap ${path_genemap} \
		--seqBias --gcBias --incompatPrior '0.0' \
		--biasSpeedSamp '5' --fldMax '36' --fldMean '36' --fldSD '25' --forgettingFactor '0.65'  \
		--maxReadOcc '500'     --numBiasSamples '2000000' --numAuxModelSamples '5000000' --numPreAuxModelSamples '5000'  \
		--numGibbsSamples '0'  --numBootstraps '0'  --thinningFactor '16'  --sigDigits '3' --vbPrior '1e-05'   -o ${path_outdir}/${Sample_ID}


		# Log end time
		datetime=$(date '+%Y%m%d-%H%M%S')
		message="[${datetime}] Finished salmon quant ${Sample_ID}!"
		echo $message >> ${path_logfile}
		echo $message
	fi
	COUNT=$(( COUNT + 1 ))
done < ${path_metadata}

# merging salmon quant.sf by TPM
salmon quantmerge --column TPM \
--quant ${path_outdir}/S* \
--output ${path_large}/salmon/count_by_TPM_salmon_after_fastp.csv

path_index="$path_ot2_yeast/analysis/data/index"

# merging salmon quant.sf by NumReads
salmon quantmerge --column NumReads \
--quant ${path_outdir}/S* \
--output ${path_large}/salmon/count_by_NumRead_salmon_after_fastp.csv