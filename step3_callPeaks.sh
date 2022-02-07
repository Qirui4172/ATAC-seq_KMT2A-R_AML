#!/usr/bin/bash

#SBATCH -t 24:00:00
#SBATCH -J MACS2
#SBATCH --mem=250GB
#SBATCH -N 1
#SBATCH --ntasks-per-node=40


WkDir="/projects/fs1_no-backup/ls1-andersson/ANALYSIS/qirui/atac-seq/mouse"
Samples=("m16-8" "m16-11" "m16-12" "m16-14" "m16-15" "m16-16" "m16-17" "m16-19" "m16-23" "m16-25" "m16-26" "m16-27" "m16-30" "m16-31" "m16-33" "m16-35")

#--------------------------------------------------------------------
module load GCC/4.9.3-2.25 OpenMPI/1.10.2 MACS2/2.1.0.20150731-Python-2.7.11

# call peaks for individual samples
for sample in ${Samples[@]}
do
	time=`date "+%Y-%m-%d %H:%M:%S"`
	echo -e "\t$time\tCalling individual peaks of sample ${sample} ..."

	mkdir -p ${WkDir}/callpeaks/${sample}
	macs2 callpeak -t ${WkDir}/02-BAM/s2_readyBAM/${sample}.bam -n ${sample} -g mm -f BAMPE -B --keep-dup all -q 0.05 --outdir ${WkDir}/callpeaks/${sample}/ 2>${WkDir}/callpeaks/${sample}/${sample}_macs2.log

	time=`date "+%Y-%m-%d %H:%M:%S"`
	echo -e "\t$time\tDone with calling peaks of ${sample}\n"
done

#--------------------------------------------------------------------
# call consensus peaks
time=`date "+%Y-%m-%d %H:%M:%S"`
echo -e "\t$time\tCalling consensus peaks ..."

mkdir -p ${WkDir}/callpeaks/consensus
macs2 callpeak -t ${WkDir}/02-BAM/s2_readyBAM/m16-8.bam ${WkDir}/02-BAM/s2_readyBAM/m16-11.bam ${WkDir}/02-BAM/s2_readyBAM/m16-12.bam ${WkDir}/02-BAM/s2_readyBAM/m16-14.bam \
${WkDir}/02-BAM/s2_readyBAM/m16-15.bam ${WkDir}/02-BAM/s2_readyBAM/m16-16.bam ${WkDir}/02-BAM/s2_readyBAM/m16-17.bam ${WkDir}/02-BAM/s2_readyBAM/m16-19.bam \
${WkDir}/02-BAM/s2_readyBAM/m16-23.bam ${WkDir}/02-BAM/s2_readyBAM/m16-25.bam ${WkDir}/02-BAM/s2_readyBAM/m16-26.bam ${WkDir}/02-BAM/s2_readyBAM/m16-27.bam \
${WkDir}/02-BAM/s2_readyBAM/m16-30.bam ${WkDir}/02-BAM/s2_readyBAM/m16-31.bam ${WkDir}/02-BAM/s2_readyBAM/m16-33.bam ${WkDir}/02-BAM/s2_readyBAM/m16-35.bam \
-n consensus -g mm -f BAMPE -B --keep-dup all -q 0.05 --outdir ${WkDir}/callpeaks/consensus/ 2>${WkDir}/callpeaks/consensus/consensusPeaks_macs2.log

time=`date "+%Y-%m-%d %H:%M:%S"`
echo -e "\t$time\tDone with calling consensus peaks\n"

#--------------------------------------------------------------------
# call group-consensus peaks

# Nras-G12D consensus peaks
time=`date "+%Y-%m-%d %H:%M:%S"`
echo -e "\t$time\tCalling Nras-G12D consensus peaks ..."
mkdir -p ${WkDir}/callpeaks/ng
macs2 callpeak -t ${WkDir}/02-BAM/s2_readyBAM/m16-8.bam ${WkDir}/02-BAM/s2_readyBAM/m16-11.bam ${WkDir}/02-BAM/s2_readyBAM/m16-12.bam ${WkDir}/02-BAM/s2_readyBAM/m16-14.bam -n ng -g mm -f BAMPE -B --keep-dup all -q 0.05 --outdir ${WkDir}/callpeaks/ng/ 2>${WkDir}/callpeaks/ng/ng_macs2.log

# FLT3-ITD consensus peaks
time=`date "+%Y-%m-%d %H:%M:%S"`
echo -e "\t$time\tCalling FLT3-ITD consensus peaks ..."
mkdir -p ${WkDir}/callpeaks/fi
macs2 callpeak -t ${WkDir}/02-BAM/s2_readyBAM/m16-15.bam ${WkDir}/02-BAM/s2_readyBAM/m16-16.bam ${WkDir}/02-BAM/s2_readyBAM/m16-17.bam ${WkDir}/02-BAM/s2_readyBAM/m16-19.bam -n fi -g mm -f BAMPE -B --keep-dup all -q 0.05 --outdir ${WkDir}/callpeaks/fi/ 2>${WkDir}/callpeaks/fi/fi_macs2.log

# KMT2A-MLLT3 consensus peaks
time=`date "+%Y-%m-%d %H:%M:%S"`
echo -e "\t$time\tCalling KMT2A-MLLT3 consensus peaks ..."
mkdir -p ${WkDir}/callpeaks/km
macs2 callpeak -t ${WkDir}/02-BAM/s2_readyBAM/m16-23.bam ${WkDir}/02-BAM/s2_readyBAM/m16-25.bam ${WkDir}/02-BAM/s2_readyBAM/m16-26.bam ${WkDir}/02-BAM/s2_readyBAM/m16-27.bam -n km -g mm -f BAMPE -B --keep-dup all -q 0.05 --outdir ${WkDir}/callpeaks/km/ 2>${WkDir}/callpeaks/km/km_macs2.log

# FLT3-N676K consensus peaks
time=`date "+%Y-%m-%d %H:%M:%S"`
echo -e "\t$time\tCalling FLT3-N676K consensus peaks ..."
mkdir -p ${WkDir}/callpeaks/fn
macs2 callpeak -t ${WkDir}/02-BAM/s2_readyBAM/m16-30.bam ${WkDir}/02-BAM/s2_readyBAM/m16-31.bam ${WkDir}/02-BAM/s2_readyBAM/m16-33.bam ${WkDir}/02-BAM/s2_readyBAM/m16-35.bam -n fn -g mm -f BAMPE -B --keep-dup all -q 0.05 --outdir ${WkDir}/callpeaks/fn/ 2>${WkDir}/callpeaks/fn/fn_macs2.log

module purge

