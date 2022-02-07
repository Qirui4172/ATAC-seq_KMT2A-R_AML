#!/usr/bin/bash

#SBATCH -t 48:00:00
#SBATCH -J computeNpNt
#SBATCH --mem=250GB
#SBATCH -N 1
#SBATCH --ntasks-per-node=10


#=========================================================================================================
wkdir="/projects/fs1_no-backup/ls1-andersson/ANALYSIS/qirui/atac-seq/mouse/IDR"
cd ${wkdir}
SamplePairs=("15_16" "15_17" "15_19" "16_17" "16_19" "17_19" "30_31" "30_33" "30_35" "31_33" "31_35" "33_35" "8_11" "8_12" "8_14" "11_12" "11_14" "12_14" "23_25" "23_26" "23_27" "25_26" "25_27" "26_27")

# Function: pool BAMs of two samples, split into two replicates, call peaks, compute IDR, generate number of peaks consistent between pooled pseudoreps (Np), and number of peaks consistent between two true samples (Nt).


#=========================================================================================================
for samplepair in ${SamplePairs[@]};do
	# prepare samples
	samples=(${samplepair/_/ })
	sample1="m16-"${samples[0]}
	sample2="m16-"${samples[1]}
	sample1vs2="m16-"${samples[0]}"vs"${samples[1]}
	echo -e "processing replicates ${sample1} and ${sample2}\n"

:<<!
	#----------------------------------------------------------------------
	# merge BAMs of two samples
	echo -e "merging BAM files of samples ${sample1} and ${sample2}\n"
	module load GCC/7.3.0-2.30 SAMtools/1.9

	samtools merge -u BAM_poolpseudo/${sample1vs2}_merged.bam BAM_true/${sample1}.shifted.bam BAM_true/${sample2}.shifted.bam
	samtools view -H BAM_poolpseudo/${sample1vs2}_merged.bam > BAM_poolpseudo/${sample1vs2}_header.sam

	# obtain total and half number of reads
	echo -e "obtaining reads number\n"
	nlines=$(samtools view BAM_poolpseudo/${sample1vs2}_merged.bam |wc -l)
	nlines=$(( (nlines + 1) / 2 ))

	# shuffle and split into two pseudo samples with equal size
	echo -e "shuffling and splitting pseudo samples\n"
	prefix="BAM_poolpseudo/"${sample1vs2}"_"
	samtools view BAM_poolpseudo/${sample1vs2}_merged.bam | shuf - | split -d -l ${nlines} - ${prefix}
	  # "m16-8vs11_" is the prefix, two files (in SAM format) will be generated: m16-8vs11_00, and m16-8vs11_01

	cat BAM_poolpseudo/${sample1vs2}_header.sam BAM_poolpseudo/${sample1vs2}_00 | samtools view -bS - > BAM_poolpseudo/${sample1vs2}_00.bam
	cat BAM_poolpseudo/${sample1vs2}_header.sam BAM_poolpseudo/${sample1vs2}_01 | samtools view -bS - > BAM_poolpseudo/${sample1vs2}_01.bam
	module purge
!

	#----------------------------------------------------------------------
	# call peaks
	echo -e "calling peaks of pseudo replicates\n"
	module load GCC/4.9.3-2.25 OpenMPI/1.10.2 MACS2/2.1.0.20150731-Python-2.7.11

	mkdir -p peaks_poolpseudo/${sample1vs2}_00
	mkdir -p peaks_poolpseudo/${sample1vs2}_01
	macs2 callpeak -t BAM_poolpseudo/${sample1vs2}_00.bam -n ${sample1vs2}_00 -g mm -f BAMPE -B --SPMR --keep-dup all -q 0.05 --outdir peaks_poolpseudo/${sample1vs2}_00/ 2>peaks_poolpseudo/${sample1vs2}_00/${sample1vs2}_00_macs2.log
	macs2 callpeak -t BAM_poolpseudo/${sample1vs2}_01.bam -n ${sample1vs2}_01 -g mm -f BAMPE -B --SPMR --keep-dup all -q 0.05 --outdir peaks_poolpseudo/${sample1vs2}_01/ 2>peaks_poolpseudo/${sample1vs2}_01/${sample1vs2}_01_macs2.log

	# sort pvalues in descending order
	echo -e "sorting pvalues of narrowPeak files\n"
	sort -k8,8nr peaks_poolpseudo/${sample1vs2}_00/${sample1vs2}_00_peaks.narrowPeak > peaks_poolpseudo/${sample1vs2}_00/${sample1vs2}_00_sorted.narrowPeak
	sort -k8,8nr peaks_poolpseudo/${sample1vs2}_01/${sample1vs2}_01_peaks.narrowPeak > peaks_poolpseudo/${sample1vs2}_01/${sample1vs2}_01_sorted.narrowPeak
	module purge

:<<!
	#----------------------------------------------------------------------
	# compute IDR and number of consistent peaks between two pooled pseudo samples (Np)
	echo -e "computing IDR and number of consistent between two pooled pseudo samples (Np)\n"
	module load GCC/8.3.0  OpenMPI/3.1.4 R/3.6.2

	Rscript idr.r peaks_poolpseudo/${sample1vs2}_00/${sample1vs2}_00_sorted.narrowPeak peaks_poolpseudo/${sample1vs2}_01/${sample1vs2}_01_sorted.narrowPeak mm9_chr.size idr_poolpseudo/${sample1vs2}
	Np=$(awk '$10+0 < 0.01 {print $0}' idr_poolpseudo/${sample1vs2}_idr.txt |wc -l)
	echo -e "Np: " $Np


	#----------------------------------------------------------------------
	# Compute IDR and number of consistent peaks between two true samples (Nt)
	echo -e "computing IDR and number of consistent peaks between two true samples (Nt)\n"

	Rscript idr.r peaks_true/${sample1}_sorted.narrowPeak peaks_true/${sample2}_sorted.narrowPeak mm9_chr.size idr_true/${sample1vs2}
	Nt=$(awk '$10+0 < 0.01 {print $0}' idr_true/${sample1vs2}_idr.txt |wc -l)
	echo -e "Nt: " $Nt

	module purge
!

	#----------------------------------------------------------------------
	# Use offical idrCode (recommended!)
	# compute IDR and number of consistent between two pooled pseudo samples (Np)
	echo -e "computing IDR and number of consistent between two pooled pseudo samples (Np)\n"
	module load GCC/8.3.0  OpenMPI/3.1.4 R/3.6.2

	cd idrCode/  # batch-consistency-analysis.r requires functions-all-clayton-12-13.r and genome_table.txt in the same folder, so it's better to run this script in folder idrCode/
	Rscript batch-consistency-analysis.r ../peaks_poolpseudo/${sample1vs2}_00/${sample1vs2}_00_peaks.narrowPeak ../peaks_poolpseudo/${sample1vs2}_01/${sample1vs2}_01_peaks.narrowPeak 200 ../idr_poolpseudo/${sample1vs2} 0 F p.value
	cd ${wkdir}


	#----------------------------------------------------------------------
	# Compute IDR and number of consistent peaks between two true samples (Nt)
	echo -e "computing IDR and number of consistent peaks between two true samples (Nt)\n"

	cd idrCode/
	Rscript batch-consistency-analysis.r ../peaks_true/${sample1}_peaks.narrowPeak ../peaks_true/${sample2}_peaks.narrowPeak 200 ../idr_true/${sample1vs2} 0 F p.value
	cd ${wkdir}

done


