#!/bin/sh


WkDir="/projects/fs1_no-backup/ls1-andersson/ANALYSIS/qirui/atac-seq/mouse/callpeaks"

#-----------------------------------------------------------------------------------
# unNormalize bedgraph files
module load GCC/5.4.0-2.26 OpenMPI/1.10.3 BEDTools/2.26.0

for sample in "m16-8" "m16-11" "m16-12" "m16-14" "m16-15" "m16-16" "m16-17" "m16-19" "m16-23" "m16-25" "m16-26" "m16-27" "m16-30" "m16-31" "m16-33" "m16-35";do sort -k1,1 -k2,2n ${sample}/${sample}_treat_pileup.bdg > ${sample}/${sample}_treat_sort.bdg; bedtools intersect -v -a ${sample}/${sample}_treat_sort.bdg -b ../mm9_blacklist.bed > ${sample}/${sample}_treat_rmblack.bdg; bedGraphToBigWig.dms ${sample}/${sample}_treat_rmblack.bdg ../mm9_chrom.size ${sample}/${sample}_unNorm.bw;done

module purge

#-----------------------------------------------------------------------------------
# script (from github) normalization
module load GCC/7.3.0-2.30 OpenMPI/3.1.1 wheel/0.31.1-Python-3.6.6

for sample in "m16-8" "m16-11" "m16-12" "m16-14" "m16-15" "m16-16" "m16-17" "m16-19" "m16-23" "m16-25" "m16-26" "m16-27" "m16-30" "m16-31" "m16-33" "m16-35";do python3 ~/.local/bin/normalize_bedgraph.py --to-number-reads 10000000 ${sample}/${sample}_treat_rmblack.bdg > ${sample}/${sample}_scriptNorm.bdg; bedGraphToBigWig.dms ${sample}/${sample}_scriptNorm.bdg ../mm9_chrom.size ${sample}/${sample}_scriptNorm.bw;done

module purge

#-----------------------------------------------------------------------------------
# SPMR normalization
module load GCC/5.4.0-2.26 OpenMPI/1.10.3 BEDTools/2.26.0

for sample in "m16-8" "m16-11" "m16-12" "m16-14" "m16-15" "m16-16" "m16-17" "m16-19" "m16-23" "m16-25" "m16-26" "m16-27" "m16-30" "m16-31" "m16-33" "m16-35";do sort -k1,1 -k2,2n ${sample}_norm/${sample}_treat_pileup.bdg > ${sample}_norm/${sample}_treat_sort.bdg; bedtools intersect -v -a ${sample}_norm/${sample}_treat_sort.bdg -b ../mm9_blacklist.bed > ${sample}_norm/${sample}_treat_rmblack.bdg; bedGraphToBigWig.dms ${sample}_norm/${sample}_treat_rmblack.bdg ../mm9_chrom.size ${sample}_norm/${sample}_spmr.bw;done

module purge


