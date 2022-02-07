#!/usr/bin/bash

#SBATCH -t 24:00:00
#SBATCH -J FRiP_score
#SBATCH --mem=250GB
#SBATCH -N 1
#SBATCH --ntasks-per-node=20



peakdir="/projects/fs1_no-backup/ls1-andersson/ANALYSIS/qirui/atac-seq/mouse/callpeak_macs2"
bamdir="/projects/fs1_no-backup/ls1-andersson/ANALYSIS/qirui/atac-seq/mouse/02-BAM/s3_shiftedBAM"
Samples=("m16-8" "m16-11" "m16-12" "m16-14" "m16-15" "m16-16" "m16-17" "m16-19" "m16-23" "m16-25" "m16-26" "m16-27" "m16-30" "m16-31" "m16-33" "m16-35")

#====================================================================================
for sample in ${Samples[@]};do
	time=`date "+%Y-%m-%d %H:%M:%S"`
	echo -e "$time\t${sample}"
	module load GCC/5.4.0-2.26 OpenMPI/1.10.3 BEDTools/2.26.0

	# Reads in peaks
	cut -f1-3 ${peakdir}/${sample}/${sample}_peaks.narrowPeak | bedtools intersect -a ${bamdir}/${sample}.shifted.bam -b stdin -wa -u > ${bamdir}/${sample}_inpeak.bam
	module purge

	module load GCC/7.3.0-2.30 SAMtools/1.9
	reads_in_peaks=$(samtools view -@ 20 -c ${bamdir}/${sample}_inpeak.bam)
	echo -e "Reads in peaks: $reads_in_peaks"
	rm ${bamdir}/${sample}_inpeak.bam

	# Total reads
	total_reads=$(samtools view -@ 20 -c ${bamdir}/${sample}.shifted.bam)
	echo -e "Total reads: $total_reads"
	module purge

	# FRiP score
	FRiP=$(awk "BEGIN {print ${reads_in_peaks}/${total_reads}}")
	echo -e "FRiP score: ${FRiP}\n"

done

