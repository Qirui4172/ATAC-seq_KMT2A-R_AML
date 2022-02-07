#!/usr/bin/bash

#SBATCH -t 48:00:00
#SBATCH -J TSSEscore
#SBATCH --mem=200GB
#SBATCH -N 1
#SBATCH --ntasks-per-node=20

Wkdir="/projects/fs1_no-backup/ls1-andersson/ANALYSIS/qirui/atac-seq/mouse/02-BAM"
Samples=("m16-8" "m16-11" "m16-12" "m16-14" "m16-15" "m16-16" "m16-17" "m16-19" "m16-23" "m16-25" "m16-26" "m16-27" "m16-30" "m16-31" "m16-33" "m16-35")


#====================================================================================
for sample in ${Samples[@]};do
	echo -e "processing sample ${sample}\n"

	module load GCC/5.4.0-2.26 OpenMPI/1.10.3 BEDTools/2.26.0
	# intersect TSS flanking (-1000/1000) with BAM
	echo -e "intersecting TSS flanking (-1000/1000) with BAM\n"
	bedtools intersect -a ${Wkdir}/s3_shiftedBAM/${sample}.shifted.bam -b ${Wkdir}/s7_tsse/mm9_tss.txt.TSS.unique.2K.bed -ubam -wa -u > ${Wkdir}/s7_tsse/${sample}.sel.bam

	# generate coverage of these regions
	echo -e "generating coverage of these regions\n"
	bedtools genomecov -ibam ${Wkdir}/s7_tsse/${sample}.sel.bam -dz | grep -v 'chrM' > ${Wkdir}/s7_tsse/${sample}.sel.bam.gc

	# transfer coverage file to bed format by +1 of each position
	echo -e "transfering coverage file to bed format by +1 of each position\n"
	awk -v "OFS=\t" '{print $1,$2,$2+1,$3}' ${Wkdir}/s7_tsse/${sample}.sel.bam.gc > ${Wkdir}/s7_tsse/${sample}.sel.bam.gc.tmp

	# intersect again coverage with TSS flanking (maybe some pos out of TSS region after +1 and some located on "-" strand)
	echo -e "intersecting again coverage with TSS flanking\n"
	bedtools intersect -a ${Wkdir}/s7_tsse/mm9_tss.txt.TSS.unique.2K.bed -b ${Wkdir}/s7_tsse/${sample}.sel.bam.gc.tmp -wa -wb > ${Wkdir}/s7_tsse/${sample}.sel.bam.gc.tmp.intersect
	module purge

	# generate TSS coverage matrix
	echo -e "generating TSS coverage matrix\n"
	module load GCC/7.3.0-2.30 OpenMPI/3.1.1 wheel/0.31.1-Python-3.6.6
	python3.6 ${Wkdir}/s7_tsse/getMat.py ${Wkdir}/s7_tsse/${sample}.sel.bam.gc.tmp.intersect ${Wkdir}/s7_tsse/${sample}.tss.mx
	module purge

	# compute TSS enrichment score
	echo -e "${sample} TSSE score:"
	module load GCC/8.2.0-2.31.1  OpenMPI/3.1.3 R/3.6.0
	Rscript ${Wkdir}/s7_tsse/calTSSscore.r ${Wkdir}/s7_tsse/${sample}.tss.mx 50 100 5
	module purge
	echo -e "\n\n"

	# remove temporary files
	rm ${Wkdir}/s7_tsse/${sample}.sel.bam*
	rm ${Wkdir}/s7_tsse/${sample}.tss.mx
done


