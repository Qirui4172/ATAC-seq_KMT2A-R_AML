#!/usr/bin/bash

#SBATCH -t 12:00:00
#SBATCH -J computePBC_NRF
#SBATCH --mem=250GB
#SBATCH -N 1
#SBATCH --ntasks-per-node=10

WkDir="/projects/fs1_no-backup/ls1-andersson/ANALYSIS/qirui/atac-seq/mouse/02-BAM"
Samples=("m16-8" "m16-11" "m16-12" "m16-14" "m16-15" "m16-16" "m16-17" "m16-19" "m16-23" "m16-25" "m16-26" "m16-27" "m16-30" "m16-31" "m16-33" "m16-35")


#====================================================================================
# Compute PBC1, PBC2, and NRF
module load GCC/5.4.0-2.26 OpenMPI/1.10.3 BEDTools/2.26.0

for sample in ${Samples[@]};do
	echo -e "Computing sample: $sample\n"
	echo -e "TotalReadPairs\tDistinctReadPairs\tOneReadPair\tTwoReadPairs\tNRF\tPBC1\tPBC2" > s6_PBC_NRF/${sample}.pbc
	bedtools bamtobed -i s1_rmdupBAM/${sample}.rmdup.bam |awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$6}' |grep -v 'chrM' |sort |uniq -c |awk 'BEGIN{mt=0;m0=0;m1=0;m2=0} ($1==1){m1=m1+1} ($1==2){m2=m2+1} {m0=m0+1} {mt=mt+$1} END {printf "%d\t%d\t%d\t%d\t%f\t%f\t%f\n",mt,m0,m1,m2,m0/mt,m1/m0,m1/m2}' >> s6_PBC_NRF/${sample}.pbc

done
module purge

# mt: TotalReadPairs
# m0: DistinctReadPairs
# m1: OneReadPair
# m2: TwoReadPairs
# m0/mt: NRF=Distinct/Total
# m1/m0: PBC1=OnePair/Distinct
# m1/m2: PBC2=OnePair/TwoPair


