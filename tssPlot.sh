#!/usr/bin/bash


#SBATCH -t 24:00:00
#SBATCH -J TSS_plot
#SBATCH --mem=250GB
#SBATCH -N 1
#SBATCH --ntasks-per-node=40


wkdir="/projects/fs1_no-backup/ls1-andersson/ANALYSIS/qirui/atac-seq/mouse"
cd ${wkdir}
Samples=("m16-8" "m16-11" "m16-12" "m16-14" "m16-15" "m16-16" "m16-17" "m16-19" "m16-23" "m16-25" "m16-26" "m16-27" "m16-30" "m16-31" "m16-33" "m16-35")

#=================================================================================================
# Normalize bam and compute matrix of TSS regions
module load GCC/7.3.0-2.30 OpenMPI/3.1.1 deepTools/2.5.4-Python-3.6.6
for sample in ${Samples[@]};do
	# Normalize BAM and transfer to bedgraph format
	time=`date "+%Y-%m-%d %H:%M:%S"`
	echo -e "$time\tNormalizing ${sample} using RPKM, transfer to bigwig"
	bamCoverage -b 02-BAM/s3_shiftedBAM/${sample}.shifted.bam -bs 5 -p 40 -bl mm9_blacklist.bed --normalizeUsingRPKM -of bedgraph -o consensus_tssPlot/${sample}_rpkm.bg.tmp

	# fill coordinate gaps
	awk 'BEGIN {FS="\t";OFS="\t"} NR==FNR {array[$1]=$2;next;} NR>FNR && NR==FNR+1 {chr=$1;start=$2;end=$3;score=$4;if(start!=0){print chr,0,start,0};print $0;next;} $2==end {print $0;chr=$1;start=$2;end=$3;score=$4;next;} $2!=end && chr==$1 {print chr,end,$2,0;print $0;chr=$1;start=$2;end=$3;score=$4;next;} $1!=chr {if(end!=array[chr]){print chr,end,array[chr],0}; chr=$1;start=$2;end=$3;score=$4;if(start!=0){print chr,0,start,0};print $0;} END {if(end!=array[chr]){print chr,end,array[chr],0}}' mm9_chr.size consensus_tssPlot/${sample}_rpkm.bg.tmp > consensus_tssPlot/${sample}_rpkm.bg
	rm consensus_tssPlot/${sample}_rpkm.bg.tmp

	# transfer to BIGWIG format
	bedGraphToBigWig.dms consensus_tssPlot/${sample}_rpkm.bg mm9_chr.size consensus_tssPlot/${sample}_rpkm.bw

	#---------------------------------------------------------------------
	# Compute matrix (This part is only for calculating TSS Enrichment Score)
	time=`date "+%Y-%m-%d %H:%M:%S"`
	echo -e "$time\tComputing matrix of ${sample}"
	computeMatrix reference-point --referencePoint center -S consensus_tssPlot/${sample}_rpkm.bw -R consensus_tssPlot/mm9_gencode_tss_unique.bed -b 2100 -a 2100 -bs 10 -p 40 --missingDataAsZero -out consensus_tssPlot/tsse/${sample}_matrix.gz --outFileNameMatrix consensus_tssPlot/tsse/${sample}_matrix.tab --outFileSortedRegions consensus_tssPlot/tsse/${sample}_regions.bed
	#---------------------------------------------------------------------

done


#=================================================================================================
# Plot profile and heatmap

# FLT3-ITD group
time=`date "+%Y-%m-%d %H:%M:%S"`
echo -e "$time\tComputing matrix of FLT3-ITD group ..."
computeMatrix reference-point --referencePoint center -S consensus_tssPlot/m16-15_rpkm.bw consensus_tssPlot/m16-16_rpkm.bw consensus_tssPlot/m16-17_rpkm.bw consensus_tssPlot/m16-19_rpkm.bw -R consensus_tssPlot/mm9_gencode_tss_unique.bed -b 2000 -a 2000 -bs 10 -p 40 --missingDataAsZero -out consensus_tssPlot/fi_matrix.gz --outFileNameMatrix consensus_tssPlot/fi_matrix.tab --outFileSortedRegions consensus_tssPlot/fi_regions.bed

echo -e "Ploting profile and heatmap ..."
plotProfile -m consensus_tssPlot/fi_matrix.gz --perGroup --legendLocation=upper-left --colors '#CC3300' '#FFCC33' '#009900' '#0066cc' -out consensus_tssPlot/fi_profile.pdf
plotHeatmap -m consensus_tssPlot/fi_matrix.gz -out consensus_tssPlot/fi_heatmap.pdf --colorMap=Blues

#---------------------------------------------------------
# FLT3-N676K group
time=`date "+%Y-%m-%d %H:%M:%S"`
echo -e "$time\tComputing matrix of FLT3-N676K group ..."
computeMatrix reference-point --referencePoint center -S consensus_tssPlot/m16-30_rpkm.bw consensus_tssPlot/m16-31_rpkm.bw consensus_tssPlot/m16-33_rpkm.bw consensus_tssPlot/m16-35_rpkm.bw -R consensus_tssPlot/mm9_gencode_tss_unique.bed -b 2000 -a 2000 -bs 10 -p 40 --missingDataAsZero -out consensus_tssPlot/fn_matrix.gz --outFileNameMatrix consensus_tssPlot/fn_matrix.tab --outFileSortedRegions consensus_tssPlot/fn_regions.bed

echo -e "Ploting profile and heatmap ..."
plotProfile -m consensus_tssPlot/fn_matrix.gz --perGroup --legendLocation=upper-left --colors '#CC3300' '#FFCC33' '#009900' '#0066cc' -out consensus_tssPlot/fn_profile.pdf
plotHeatmap -m consensus_tssPlot/fn_matrix.gz -out consensus_tssPlot/fn_heatmap.pdf --colorMap=Blues

#---------------------------------------------------------
# NRAS-G12D
time=`date "+%Y-%m-%d %H:%M:%S"`
echo -e "$time\tComputing matrix of NRAS-G12D group ..."
computeMatrix reference-point --referencePoint center -S consensus_tssPlot/m16-8_rpkm.bw consensus_tssPlot/m16-11_rpkm.bw consensus_tssPlot/m16-12_rpkm.bw consensus_tssPlot/m16-14_rpkm.bw -R consensus_tssPlot/mm9_gencode_tss_unique.bed -b 2000 -a 2000 -bs 10 -p 40 --missingDataAsZero -out consensus_tssPlot/ng_matrix.gz --outFileNameMatrix consensus_tssPlot/ng_matrix.tab --outFileSortedRegions consensus_tssPlot/ng_regions.bed

echo -e "Ploting profile and heatmap ..."
plotProfile -m consensus_tssPlot/ng_matrix.gz --perGroup --legendLocation=upper-left --colors '#CC3300' '#FFCC33' '#009900' '#0066cc' -out consensus_tssPlot/ng_profile.pdf
plotHeatmap -m consensus_tssPlot/ng_matrix.gz -out consensus_tssPlot/ng_heatmap.pdf --colorMap=Blues

#---------------------------------------------------------
# KMT2A-MLLT3
time=`date "+%Y-%m-%d %H:%M:%S"`
echo -e "$time\tComputing matrix of KMT2A-MLLT3 empty GFP group ..."
computeMatrix reference-point --referencePoint center -S consensus_tssPlot/m16-23_rpkm.bw consensus_tssPlot/m16-25_rpkm.bw consensus_tssPlot/m16-26_rpkm.bw consensus_tssPlot/m16-27_rpkm.bw -R consensus_tssPlot/mm9_gencode_tss_unique.bed -b 2000 -a 2000 -bs 10 -p 40 --missingDataAsZero -out consensus_tssPlot/km_matrix.gz --outFileNameMatrix consensus_tssPlot/km_matrix.tab --outFileSortedRegions consensus_tssPlot/km_regions.bed

echo -e "Ploting profile and heatmap ..."
plotProfile -m consensus_tssPlot/km_matrix.gz --perGroup --legendLocation=upper-left --colors '#CC3300' '#FFCC33' '#009900' '#0066cc' -out consensus_tssPlot/km_profile.pdf
plotHeatmap -m consensus_tssPlot/km_matrix.gz -out consensus_tssPlot/km_heatmap.pdf --colorMap=Blues

module purge

