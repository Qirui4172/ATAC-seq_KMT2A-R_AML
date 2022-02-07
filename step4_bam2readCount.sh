#!/usr/bin/bash

#SBATCH -t 24:00:00
#SBATCH -J BEDTools
#SBATCH --mem=128GB
#SBATCH -N 1
#SBATCH --ntasks-per-node=20

WkDir="/projects/fs1_no-backup/ls1-andersson/ANALYSIS/qirui/atac-seq/mouse"
Samples=("m16-8" "m16-11" "m16-12" "m16-14" "m16-15" "m16-16" "m16-17" "m16-19" "m16-23" "m16-25" "m16-26" "m16-27" "m16-30" "m16-31" "m16-33" "m16-35")

#------------------------------------------------------------------------------------------------------------------
module load GCC/5.4.0-2.26 OpenMPI/1.10.3 BEDTools/2.26.0

for sample in ${Samples[@]}
do
	time=`date "+%Y-%m-%d %H:%M:%S"`
	echo -e "\t$time\tQuantifying sample $sample ..."

	bedtools multicov -bams ${WkDir}/02-BAM/readyBAM/${sample}.bam -bed ${WkDir}/callpeaks/consensus/consensus_final.bed |sed "1i chr\tstart\tend\tpeakID\t${sample}" > ${WkDir}/readcounts/${sample}_readcounts.txt

	time=`date "+%Y-%m-%d %H:%M:%S"`
	echo -e "\t$time\tDone with sample $sample"
done

time=`date "+%Y-%m-%d %H:%M:%S"`
echo -e "\t$time\tDone with all samples!"

module purge
