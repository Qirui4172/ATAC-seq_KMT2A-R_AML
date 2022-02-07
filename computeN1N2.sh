#!/usr/bin/bash

#SBATCH -t 48:00:00
#SBATCH -J computeN1N2
#SBATCH --mem=250GB
#SBATCH -N 1
#SBATCH --ntasks-per-node=10


#=========================================================================================================
wkdir="/projects/fs1_no-backup/ls1-andersson/ANALYSIS/qirui/atac-seq/mouse/IDR"
cd ${wkdir}
Samples=("m16-15" "m16-16" "m16-17" "m16-19" "m16-30" "m16-31" "m16-33" "m16-35" "m16-8" "m16-11" "m16-12" "m16-14" "m16-23" "m16-25" "m16-26" "m16-27")

# Function: split a BAM to two self pseudo replicates, call peaks, compute IDR, generate number of peaks consistent between self pseudo replicates (N1 or N2).


#=========================================================================================================
for sample in ${Samples[@]};do
	echo -e "processing sample ${sample}\n"
:<<!
	#----------------------------------------------------------------------
	# Split a BAM to two pseudo replicates
	module load GCC/7.3.0-2.30 SAMtools/1.9

	# obtain total and half number of reads
	echo -e "obtaining reads number\n"
	nlines=$(samtools view BAM_true/${sample}.shifted.bam |wc -l)
	nlines=$(( (nlines + 1) / 2 ))

	# shuffle and split into two pseudo replicates with equal size
	echo -e "shuffling and splitting to pseudo replicates\n"
	prefix="BAM_selfpseudo/"${sample}"_"
	samtools view BAM_true/${sample}.shifted.bam | shuf - | split -d -l ${nlines} - ${prefix}
	  # "BAM_selfpseudo/m16-8_" is the prefix, two files (in SAM format) will be generated: BAM_selfpseudo/m16-8_00, and BAM_selfpseudo/m16-8_01

	samtools view -H BAM_true/${sample}.shifted.bam > BAM_selfpseudo/${sample}_header.sam
	cat BAM_selfpseudo/${sample}_header.sam BAM_selfpseudo/${sample}_00 | samtools view -bS - > BAM_selfpseudo/${sample}_00.bam
	cat BAM_selfpseudo/${sample}_header.sam BAM_selfpseudo/${sample}_01 | samtools view -bS - > BAM_selfpseudo/${sample}_01.bam
	module purge
!

	#----------------------------------------------------------------------
	# call peaks
	echo -e "calling peaks of pseudo replicates\n"
	module load GCC/4.9.3-2.25 OpenMPI/1.10.2 MACS2/2.1.0.20150731-Python-2.7.11

	mkdir -p peaks_selfpseudo/${sample}_00
	mkdir -p peaks_selfpseudo/${sample}_01
	macs2 callpeak -t BAM_selfpseudo/${sample}_00.bam -n ${sample}_00 -g mm -f BAMPE -B --SPMR --keep-dup all -q 0.05 --outdir peaks_selfpseudo/${sample}_00/ 2>peaks_selfpseudo/${sample}_00/${sample}_00_macs2.log
	macs2 callpeak -t BAM_selfpseudo/${sample}_01.bam -n ${sample}_01 -g mm -f BAMPE -B --SPMR --keep-dup all -q 0.05 --outdir peaks_selfpseudo/${sample}_01/ 2>peaks_selfpseudo/${sample}_01/${sample}_01_macs2.log

	# sort pvalues in descending order
	echo -e "sorting pvalues of narrowPeak files\n"
	sort -k8,8nr peaks_selfpseudo/${sample}_00/${sample}_00_peaks.narrowPeak > peaks_selfpseudo/${sample}_00/${sample}_00_sorted.narrowPeak
	sort -k8,8nr peaks_selfpseudo/${sample}_01/${sample}_01_peaks.narrowPeak > peaks_selfpseudo/${sample}_01/${sample}_01_sorted.narrowPeak
	module purge

:<<!
	#----------------------------------------------------------------------
	# compute IDR and number of consistent peaks between two self pseudo replicates (N1 or N2)
	echo -e "computing IDR and number of consistent peaks between two self pseudo replicates (N1 or N2)\n"
	module load GCC/8.3.0  OpenMPI/3.1.4 R/3.6.2

	Rscript idr.r peaks_selfpseudo/${sample}_00/${sample}_00_sorted.narrowPeak peaks_selfpseudo/${sample}_01/${sample}_01_sorted.narrowPeak mm9_chr.size idr_selfpseudo/${sample}
	N1=$(awk '$10+0 < 0.01 {$0}' idr_selfpseudo/${sample}_idr.txt |wc -l)
	echo -e "Np: " $N1

	module purge
!

	#----------------------------------------------------------------------
	# Use offical idrCode (recommended!)
	# Compute IDR and number of consistent peaks between two self pseudo replicates (N1 or N2)
	echo -e "computing IDR and number of consistent peaks between two self pseudo replicates (N1 or N2)\n"
	module load GCC/8.3.0  OpenMPI/3.1.4 R/3.6.2

	cd idrCode/  # batch-consistency-analysis.r requires functions-all-clayton-12-13.r and genome_table.txt in the same folder, so it's better to run this script in folder idrCode/

	Rscript batch-consistency-analysis.r ../peaks_selfpseudo/${sample}_00/${sample}_00_peaks.narrowPeak ../peaks_selfpseudo/${sample}_01/${sample}_01_peaks.narrowPeak 200 ../idr_selfpseudo/${sample} 0 F p.value
	cd ${wkdir}

done


