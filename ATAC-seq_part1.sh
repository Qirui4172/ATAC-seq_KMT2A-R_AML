#!/usr/bin/bash

#SBATCH -t 96:00:00
#SBATCH -J ATAC-seq
#SBATCH --mem=250GB
#SBATCH -N 1
#SBATCH --ntasks-per-node=40

WkDir="/projects/fs1_no-backup/ls1-andersson/ANALYSIS/qirui/atac-seq/mouse"
cd ${WkDir}
Samples=("m16-8" "m16-11" "m16-12" "m16-14" "m16-15" "m16-16" "m16-17" "m16-19" "m16-23" "m16-25" "m16-26" "m16-27" "m16-30" "m16-31" "m16-33" "m16-35")

#============================================================================================================
# Step 1. Remove duplicates

time=`date "+%Y-%m-%d %H:%M:%S"`
echo -e "$time Step 1. Removing duplicates\n"

#-------------------------------------------------
# sort BAM files
module load GCC/7.3.0-2.30 SAMtools/1.9
for sample in ${Samples[@]};do
	samtools sort 02-1-MERGE_UNSORT/${sample}.bam > 02-1-MERGE_UNSORT/${sample}.sort.bam
done
module purge

#-------------------------------------------------
# remove duplicates
module load picard/2.6.0-Java-1.8.0_131
for sample in ${Samples[@]};do
	java -jar $EBROOTPICARD/picard.jar MarkDuplicates REMOVE_DUPLICATES=true I=02-1-MERGE_UNSORT/${sample}.sort.bam O=s1_rmdupBAM/${sample}.rmdup.bam M=s1_rmdupBAM/${sample}.rmdup.metrics
done
module purge

#-------------------------------------------------
# create index
module load GCC/7.3.0-2.30 SAMtools/1.9
for sample in ${Samples[@]};do
	samtools index s1_rmdupBAM/${sample}.rmdup.bam
done
module purge


#============================================================================================================
# Step 2. Filter BAM files

time=`date "+%Y-%m-%d %H:%M:%S"`
echo -e "$time Step 2. Filtering BAM files\n"
module load GCC/7.3.0-2.30 SAMtools/1.9

for sample in ${Samples[@]};do
	# extract mapped reads
	samtools view -F 4 -b s1_rmdupBAM/${sample}.rmdup.bam > s2_filteredBAM/${sample}.mapped.bam
	samtools index s2_filteredBAM/${sample}.mapped.bam

	# remove chrM and other abnormal chromosomes (chr_random, *)
	samtools view -h s2_filteredBAM/${sample}.mapped.bam chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chrX chrY |grep -v -E 'chrM|chr.*_random' |samtools view -S -b > s2_filteredBAM/${sample}.cleanChr.bam

	# remove low-quality alignments (<Q30)
	samtools view -b -q 30 s2_filteredBAM/${sample}.cleanChr.bam > s2_filteredBAM/${sample}.highQual.bam

	# remove improper alignments (BAM flag 0x2; then all fragment size are < 1kb)
	samtools view -b -f 2 s2_filteredBAM/${sample}.highQual.bam > s2_filteredBAM/${sample}.properAlign.bam

	# remove remaining single mate of read pairs
	samtools sort -n s2_filteredBAM/${sample}.properAlign.bam |samtools view -h |awk 'BEGIN {FS="\t";OFS="\t";i=0;} $1~/^@/ {print $0;next;} $1!~/^@/ && i==0 {id=$1;previous=$0;i=1;next;} $1==id {print previous;print $0;id=$1;next;} $1!=id {id=$1;previous=$0;next;}' |samtools view -S -b > s2_filteredBAM/${sample}.paired.bam
	samtools sort s2_filteredBAM/${sample}.paired.bam > s2_filteredBAM/${sample}.filtered.bam
	samtools index s2_filteredBAM/${sample}.filtered.bam

	# remove temporary files
	rm s2_filteredBAM/${sample}.cleanChr.bam
	rm s2_filteredBAM/${sample}.highQual.bam
	rm s2_filteredBAM/${sample}.properAlign.bam
	rm s2_filteredBAM/${sample}.paired.bam
done
module purge


#============================================================================================================
# Step 3. Shift reads

time=`date "+%Y-%m-%d %H:%M:%S"`
echo -e "$time Step 3. Shifting reads\n"
module load GCC/9.3.0  OpenMPI/4.0.3 R-bundle-Bioconductor/3.11-R-4.0.0
module load GCC/9.3.0  OpenMPI/4.0.3 R/4.0.0

for sample in ${Samples[@]};do
	Rscript shiftATACseqReads.r s2_filteredBAM/${sample}.filtered.bam ${sample} s3_shiftedBAM/${sample}.shifted.bam
done
module purge


#============================================================================================================
# Step 4. Call peaks

time=`date "+%Y-%m-%d %H:%M:%S"`
echo -e "$time Step 4. Calling peaks\n"
module load GCC/4.9.3-2.25 OpenMPI/1.10.2 MACS2/2.1.0.20150731-Python-2.7.11

#-------------------------------------------------
# call individual sample peaks
for sample in ${Samples[@]};do
	mkdir -p callpeak_macs2/${sample}
	macs2 callpeak -t 02-BAM/s3_shiftedBAM/${sample}.shifted.bam -n ${sample} -g mm -f BAMPE -B --SPMR --keep-dup all -q 0.05 --outdir callpeak_macs2/${sample}/ 2>callpeak_macs2/${sample}/${sample}_macs2.log
done

#-------------------------------------------------
# call genotype consensus peaks
mkdir -p callpeak_macs2/fi callpeak_macs2/fn callpeak_macs2/ng callpeak_macs2/km

# FLT3-ITD
macs2 callpeak -t 02-BAM/s3_shiftedBAM/m16-15.shifted.bam 02-BAM/s3_shiftedBAM/m16-16.shifted.bam 02-BAM/s3_shiftedBAM/m16-17.shifted.bam 02-BAM/s3_shiftedBAM/m16-19.shifted.bam -n fi -g mm -f BAMPE -B --keep-dup all -q 0.05 --outdir callpeak_macs2/fi/ 2>callpeak_macs2/fi/fi_macs2.log

# FLT3-N676K consensus peaks
macs2 callpeak -t 02-BAM/s3_shiftedBAM/m16-30.shifted.bam 02-BAM/s3_shiftedBAM/m16-31.shifted.bam 02-BAM/s3_shiftedBAM/m16-33.shifted.bam 02-BAM/s3_shiftedBAM/m16-35.shifted.bam -n fn -g mm -f BAMPE -B --keep-dup all -q 0.05 --outdir callpeak_macs2/fn/ 2>callpeak_macs2/fn/fn_macs2.log

# NRAS-G12D (sample m16-8 was removed)
macs2 callpeak -t 02-BAM/s3_shiftedBAM/m16-8.shifted.bam 02-BAM/s3_shiftedBAM/m16-12.shifted.bam 02-BAM/s3_shiftedBAM/m16-14.shifted.bam -n ng -g mm -f BAMPE -B --keep-dup all -q 0.05 --outdir callpeak_macs2/ng/ 2>callpeak_macs2/ng/ng_macs2.log

# KMT2A-MLLT3
macs2 callpeak -t 02-BAM/s3_shiftedBAM/m16-23.shifted.bam 02-BAM/s3_shiftedBAM/m16-25.shifted.bam 02-BAM/s3_shiftedBAM/m16-26.shifted.bam 02-BAM/s3_shiftedBAM/m16-27.shifted.bam -n km -g mm -f BAMPE -B --keep-dup all -q 0.05 --outdir callpeak_macs2/km/ 2>callpeak_macs2/km/km_macs2.log

#-------------------------------------------------
# call all-sample consensus peaks (sample m16-8 was removed)
mkdir -p callpeak_macs2/consensus

macs2 callpeak -t 02-BAM/s3_shiftedBAM/m16-8.shifted.bam 02-BAM/s3_shiftedBAM/m16-12.shifted.bam 02-BAM/s3_shiftedBAM/m16-14.shifted.bam \
02-BAM/s3_shiftedBAM/m16-15.shifted.bam 02-BAM/s3_shiftedBAM/m16-16.shifted.bam 02-BAM/s3_shiftedBAM/m16-17.shifted.bam 02-BAM/s3_shiftedBAM/m16-19.shifted.bam \
02-BAM/s3_shiftedBAM/m16-23.shifted.bam 02-BAM/s3_shiftedBAM/m16-25.shifted.bam 02-BAM/s3_shiftedBAM/m16-26.shifted.bam 02-BAM/s3_shiftedBAM/m16-27.shifted.bam \
02-BAM/s3_shiftedBAM/m16-30.shifted.bam 02-BAM/s3_shiftedBAM/m16-31.shifted.bam 02-BAM/s3_shiftedBAM/m16-33.shifted.bam 02-BAM/s3_shiftedBAM/m16-35.shifted.bam \
-n consensus -g mm -f BAMPE -B --keep-dup all -q 0.05 --outdir callpeak_macs2/consensus/ 2>callpeak_macs2/consensus/consensus_macs2.log

module purge


#============================================================================================================
Step 5. Generate consensus peaks

time=`date "+%Y-%m-%d %H:%M:%S"`
echo -e "$time Step 5. Generating consensus peaks\n"
module load GCC/5.4.0-2.26 OpenMPI/1.10.3 BEDTools/2.26.0

#-------------------------------------------------
# remove blacklisted regions

# for individual samples
for sample in ${Samples[@]};do
	cut -f1-4 callpeak_macs2/${sample}/${sample}_peaks.narrowPeak > callpeak_macs2/${sample}/${sample}_raw.bed
	bedtools intersect -v -a callpeak_macs2/${sample}/${sample}_raw.bed -b mm9_blacklist.bed > callpeak_macs2/${sample}/${sample}_rmblacklist.bed
done

# for genotype consensus peaks
for genotype in "fi" "fn" "ng" "km";do
	cut -f1-4 callpeak_macs2/${genotype}/${genotype}_peaks.narrowPeak > callpeak_macs2/${genotype}/${genotype}_raw.bed
	bedtools intersect -v -a callpeak_macs2/${genotype}/${genotype}_raw.bed -b mm9_blacklist.bed > callpeak_macs2/${genotype}/${genotype}_rmblacklist.bed
done

# for all-sample consensus peaks
cut -f1-4 callpeak_macs2/consensus/consensus_peaks.narrowPeak |sed 's/consensus_//g'> callpeak_macs2/consensus/consensus_raw.bed
bedtools intersect -v -a callpeak_macs2/consensus/consensus_raw.bed -b mm9_blacklist.bed > callpeak_macs2/consensus/consensus_rmblacklist.bed

#-------------------------------------------------
# merge genotype peaks
bedtools intersect -a callpeak_macs2/fi/fi_rmblacklist.bed -b callpeak_macs2/m16-15/m16-15_rmblacklist.bed callpeak_macs2/m16-16/m16-16_rmblacklist.bed callpeak_macs2/m16-17/m16-17_rmblacklist.bed callpeak_macs2/m16-19/m16-19_rmblacklist.bed -wa -wb |awk 'BEGIN {FS="\t";OFS="\t";} NR==1 {chr=$1;start=$2;end=$3;peakID=$4;count=1;next;} $4==peakID {count+=1;next;} $4!=peakID {if(count>=3){print chr,start,end,peakID}; chr=$1;start=$2;end=$3;peakID=$4;count=1;} END {if(count>=3){print chr,start,end,peakID}}' > callpeak_macs2/fi/fi_merged.bed	# Keep peaks overlapped in at least 3 out of 4 replicates

bedtools intersect -a callpeak_macs2/fn/fn_rmblacklist.bed -b callpeak_macs2/m16-30/m16-30_rmblacklist.bed callpeak_macs2/m16-31/m16-31_rmblacklist.bed callpeak_macs2/m16-33/m16-33_rmblacklist.bed callpeak_macs2/m16-35/m16-35_rmblacklist.bed -wa -wb |awk 'BEGIN {FS="\t";OFS="\t";} NR==1 {chr=$1;start=$2;end=$3;peakID=$4;count=1;next;} $4==peakID {count+=1;next;} $4!=peakID {if(count>=3){print chr,start,end,peakID}; chr=$1;start=$2;end=$3;peakID=$4;count=1;} END {if(count>=3){print chr,start,end,peakID}}' > callpeak_macs2/fn/fn_merged.bed

bedtools intersect -a callpeak_macs2/ng/ng_rmblacklist.bed -b callpeak_macs2/m16-8/m16-8_rmblacklist.bed callpeak_macs2/m16-12/m16-12_rmblacklist.bed callpeak_macs2/m16-14/m16-14_rmblacklist.bed -wa -wb |awk 'BEGIN {FS="\t";OFS="\t";} NR==1 {chr=$1;start=$2;end=$3;peakID=$4;count=1;next;} $4==peakID {count+=1;next;} $4!=peakID {if(count>=2){print chr,start,end,peakID}; chr=$1;start=$2;end=$3;peakID=$4;count=1;} END {if(count>=3){print chr,start,end,peakID}}' > callpeak_macs2/ng/ng_merged.bed # sample m16-11 was removed, so keep peaks overlapped in at least 2 out of 4 replicates

bedtools intersect -a callpeak_macs2/km/km_rmblacklist.bed -b callpeak_macs2/m16-23/m16-23_rmblacklist.bed callpeak_macs2/m16-25/m16-25_rmblacklist.bed callpeak_macs2/m16-26/m16-26_rmblacklist.bed callpeak_macs2/m16-27/m16-27_rmblacklist.bed -wa -wb |awk 'BEGIN {FS="\t";OFS="\t";} NR==1 {chr=$1;start=$2;end=$3;peakID=$4;count=1;next;} $4==peakID {count+=1;next;} $4!=peakID {if(count>=3){print chr,start,end,peakID}; chr=$1;start=$2;end=$3;peakID=$4;count=1;} END {if(count>=3){print chr,start,end,peakID}}' > callpeak_macs2/km/km_merged.bed

#-------------------------------------------------
# generate genotype consensus peaks
for genotype in "fi" "fn" "ng" "km";do
	bedtools intersect -a callpeak_macs2/consensus/consensus_rmblacklist.bed -b callpeak_macs2/${genotype}/${genotype}_merged.bed -wa |uniq |awk '{print $0"\t""'${genotype}'"}' > callpeak_macs2/${genotype}/${genotype}_consensus.bed
done

#-------------------------------------------------
# generate all-sample consensus peaks
bedtools intersect -a callpeak_macs2/consensus/consensus_rmblacklist.bed -b callpeak_macs2/fi/fi_merged.bed callpeak_macs2/fn/fn_merged.bed callpeak_macs2/ng/ng_merged.bed callpeak_macs2/km/km_merged.bed -wa -wb |cut -f1-4|uniq > callpeak_macs2/consensus/consensus_final.bed

module purge


#============================================================================================================
# Step 6. Generate read counts

time=`date "+%Y-%m-%d %H:%M:%S"`
echo -e "$time Step 6. Generating read counts\n"
module load GCC/5.4.0-2.26 OpenMPI/1.10.3 BEDTools/2.26.0

for sample in ${Samples[@]};do
	bedtools multicov -bams 02-BAM/s3_shiftedBAM/${sample}.shifted.bam -bed callpeak_macs2/consensus/consensus_final.bed |sed "1i chr\tstart\tend\tpeakID\t${sample}" > consensus_readcounts/${sample}_readcounts.txt
	cut -f5 consensus_readcounts/${sample}_readcounts.txt > consensus_readcounts/${sample}.tmp
done
module purge

cd consensus_readcounts/
paste m16-15_readcounts.txt m16-16.tmp m16-17.tmp m16-19.tmp m16-30.tmp m16-31.tmp m16-33.tmp m16-35.tmp m16-8.tmp m16-11.tmp m16-12.tmp m16-14.tmp m16-23.tmp m16-25.tmp m16-26.tmp m16-27.tmp > mmATAC-seq_readcounts.mx
rm m16-*tmp
cd ${WkDir}


#============================================================================================================
# Step 7. Generate BIGWIG

time=`date "+%Y-%m-%d %H:%M:%S"`
echo -e "$time Step 7. Generating normalized BIGWIG files\n"
module load GCC/5.4.0-2.26 OpenMPI/1.10.3 BEDTools/2.26.0

for sample in ${Samples[@]};do
	sort -k1,1 -k2,2n callpeak_macs2/${sample}/${sample}_treat_pileup.bdg > callpeak_macs2/${sample}/${sample}_treat_sort.bdg
	bedtools intersect -v -a callpeak_macs2/${sample}/${sample}_treat_sort.bdg -b mm9_blacklist.bed > callpeak_macs2/${sample}/${sample}_treat_rmblacklist.bdg
	bedGraphToBigWig.dms callpeak_macs2/${sample}/${sample}_treat_rmblacklist.bdg mm9_chr.size norm_bigwig/macs2/${sample}_macs2.bw
done

module purge


#============================================================================================================
# Step 8. Annotate peaks

time=`date "+%Y-%m-%d %H:%M:%S"`
echo -e "$time Step 8. Annotating peaks\n"
module load GCC/7.3.0-2.30 homer/4.10

for group in "fi" "fn" "ng" "km" "consensus";do
	annotatePeaks.pl consensus_annotation/${group}_forHomer.bed mm9 > consensus_annotation/anno_${group}_raw.tsv
	sed -n '2,$p' consensus_annotation/anno_${group}_raw.tsv |cut -f1-4,8,10,16 |sort -k1,1V > consensus_annotation/anno_${group}_short.txt
done

module purge


#============================================================================================================
# Step 9. TSS plot

time=`date "+%Y-%m-%d %H:%M:%S"`
echo -e "$time Step 9. Generating TSS plot\n"
module load GCC/7.3.0-2.30 OpenMPI/3.1.1 deepTools/2.5.4-Python-3.6.6

#-------------------------------------------------
# normalize BAM and compute matrix of TSS regions
for sample in ${Samples[@]};do
	# normalize BAM and transfer to BEDGraph
	bamCoverage -b 02-BAM/s3_shiftedBAM/${sample}.shifted.bam -bs 5 -p 40 -bl mm9_blacklist.bed --normalizeUsingRPKM -of bedgraph -o consensus_tssPlot/${sample}_rpkm.bg.tmp

	# fill coordinate gaps
	awk 'BEGIN {FS="\t";OFS="\t"} NR==FNR {array[$1]=$2;next;} NR>FNR && NR==FNR+1 {chr=$1;start=$2;end=$3;score=$4;if(start!=0){print chr,0,start,0};print $0;next;} $2==end {print $0;chr=$1;start=$2;end=$3;score=$4;next;} $2!=end && chr==$1 {print chr,end,$2,0;print $0;chr=$1;start=$2;end=$3;score=$4;next;} $1!=chr {if(end!=array[chr]){print chr,end,array[chr],0}; chr=$1;start=$2;end=$3;score=$4;if(start!=0){print chr,0,start,0};print $0;} END {if(end!=array[chr]){print chr,end,array[chr],0}}' mm9_chr.size consensus_tssPlot/${sample}_rpkm.bg.tmp > consensus_tssPlot/${sample}_rpkm.bg
	rm consensus_tssPlot/${sample}_rpkm.bg.tmp

	# transfer BEDGraph to BIGWIG format
	bedGraphToBigWig.dms consensus_tssPlot/${sample}_rpkm.bg mm9_chr.size consensus_tssPlot/${sample}_rpkm.bw
done

#-------------------------------------------------
# compute matrix, plot TSS profile and chromatin accessibility heatmap

# FLT3-ITD
computeMatrix reference-point --referencePoint center -S consensus_tssPlot/m16-15_rpkm.bw consensus_tssPlot/m16-16_rpkm.bw consensus_tssPlot/m16-17_rpkm.bw consensus_tssPlot/m16-19_rpkm.bw -R consensus_tssPlot/mm9_gencode_tss_unique.bed -b 2000 -a 2000 -bs 10 -p 40 --missingDataAsZero -out consensus_tssPlot/fi_matrix.gz --outFileNameMatrix consensus_tssPlot/fi_matrix.tab --outFileSortedRegions consensus_tssPlot/fi_regions.bed
plotProfile -m consensus_tssPlot/fi_matrix.gz --perGroup --legendLocation=upper-left --colors '#CC3300' '#FFCC33' '#009900' '#0066cc' -out consensus_tssPlot/fi_profile.pdf
plotHeatmap -m consensus_tssPlot/fi_matrix.gz -out consensus_tssPlot/fi_heatmap.pdf --colorMap=Blues

# FLT3-N676K
computeMatrix reference-point --referencePoint center -S consensus_tssPlot/m16-30_rpkm.bw consensus_tssPlot/m16-31_rpkm.bw consensus_tssPlot/m16-33_rpkm.bw consensus_tssPlot/m16-35_rpkm.bw -R consensus_tssPlot/mm9_gencode_tss_unique.bed -b 2000 -a 2000 -bs 10 -p 40 --missingDataAsZero -out consensus_tssPlot/fn_matrix.gz --outFileNameMatrix consensus_tssPlot/fn_matrix.tab --outFileSortedRegions consensus_tssPlot/fn_regions.bed
plotProfile -m consensus_tssPlot/fn_matrix.gz --perGroup --legendLocation=upper-left --colors '#CC3300' '#FFCC33' '#009900' '#0066cc' -out consensus_tssPlot/fn_profile.pdf
plotHeatmap -m consensus_tssPlot/fn_matrix.gz -out consensus_tssPlot/fn_heatmap.pdf --colorMap=Blues

# NRAS-G12D
computeMatrix reference-point --referencePoint center -S consensus_tssPlot/m16-8_rpkm.bw consensus_tssPlot/m16-11_rpkm.bw consensus_tssPlot/m16-12_rpkm.bw consensus_tssPlot/m16-14_rpkm.bw -R consensus_tssPlot/mm9_gencode_tss_unique.bed -b 2000 -a 2000 -bs 10 -p 40 --missingDataAsZero -out consensus_tssPlot/ng_matrix.gz --outFileNameMatrix consensus_tssPlot/ng_matrix.tab --outFileSortedRegions consensus_tssPlot/ng_regions.bed
plotProfile -m consensus_tssPlot/ng_matrix.gz --perGroup --legendLocation=upper-left --colors '#CC3300' '#FFCC33' '#009900' '#0066cc' -out consensus_tssPlot/ng_profile.pdf
plotHeatmap -m consensus_tssPlot/ng_matrix.gz -out consensus_tssPlot/ng_heatmap.pdf --colorMap=Blues

# KMT2A-MLLT3
computeMatrix reference-point --referencePoint center -S consensus_tssPlot/m16-23_rpkm.bw consensus_tssPlot/m16-25_rpkm.bw consensus_tssPlot/m16-26_rpkm.bw consensus_tssPlot/m16-27_rpkm.bw -R consensus_tssPlot/mm9_gencode_tss_unique.bed -b 2000 -a 2000 -bs 10 -p 40 --missingDataAsZero -out consensus_tssPlot/km_matrix.gz --outFileNameMatrix consensus_tssPlot/km_matrix.tab --outFileSortedRegions consensus_tssPlot/km_regions.bed
plotProfile -m consensus_tssPlot/km_matrix.gz --perGroup --legendLocation=upper-left --colors '#CC3300' '#FFCC33' '#009900' '#0066cc' -out consensus_tssPlot/km_profile.pdf
plotHeatmap -m consensus_tssPlot/km_matrix.gz -out consensus_tssPlot/km_heatmap.pdf --colorMap=Blues

module purge


#====================================================================================
time=`date "+%Y-%m-%d %H:%M:%S"`
echo -e "$time Done with all analysis!\n"

