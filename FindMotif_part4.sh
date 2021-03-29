#!/usr/bin/bash

#SBATCH -t 24:00:00
#SBATCH -J findMotif
#SBATCH --mem=250GB
#SBATCH -N 1
#SBATCH --ntasks-per-node=20


#======================================================================================================
module load GCC/7.3.0-2.30 homer/4.10
for file in `ls *bed`;do
	prefix=$(basename $file .bed)
	outdir=${prefix}_motif
	mkdir -p $outdir
	findMotifsGenome.pl $file /projects/fs1_no-backup/ls1-andersson/ANALYSIS/qirui/ref/mm9/mm9.fa $outdir -p 20
done
module purge


