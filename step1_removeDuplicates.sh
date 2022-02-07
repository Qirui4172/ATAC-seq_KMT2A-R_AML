#!/usr/bin/bash

#SBATCH -t 12:00:00
#SBATCH -J removeDuplicate
#SBATCH --mem=250GB
#SBATCH -N 1
#SBATCH --ntasks-per-node=40

WkDir="/projects/fs1_no-backup/ls1-andersson/ANALYSIS/qirui/atac-seq/mouse/02-BAM"
Samples=("m16-8" "m16-11" "m16-12" "m16-14" "m16-15" "m16-16" "m16-17" "m16-19" "m16-23" "m16-25" "m16-26" "m16-27" "m16-30" "m16-31" "m16-33" "m16-35")

#------------------------------------------------------------------------------
# sort BAM files
module load GCC/7.3.0-2.30 SAMtools/1.9

for sample in ${Samples[@]}
do
	samtools sort -@ 16 ${WkDir}/02-1-MERGE_UNSORT/${sample}.bam > ${WkDir}/02-1-MERGE_UNSORT/${sample}.sort.bam
done
module purge

#------------------------------------------------------------------------------
# remove duplicates
module load picard/2.6.0-Java-1.8.0_131

for sample in ${Samples[@]}
do
	java -jar $EBROOTPICARD/picard.jar MarkDuplicates REMOVE_DUPLICATES=true I=${WkDir}/02-1-MERGE_UNSORT/${sample}.sort.bam O=${WkDir}/s1_rmdupBAM/${sample}_rmdup.bam M=${WkDir}/s1_rmdupBAM/${sample}_rmdup.metrics
done
module purge

#------------------------------------------------------------------------------
# index rmdup.bam
module load GCC/7.3.0-2.30 SAMtools/1.9

for sample in ${Samples[@]}
do
#	samtools sort -@ 16 ${WkDir}/s1_rmdupBAM/${sample}_rmdup.bam > ${WkDir}/s1_rmdupBAM/${sample}.sort.bam
	# since BAM was sorted before removing duplicates, it's no need to sort again after removed duplicated, it's in sorted order.

	samtools index -@ 16 ${WkDir}/s1_rmdupBAM/${sample}_rmdup.bam
done
module purge

