
#SBATCH -t 12:00:00
#SBATCH -J filterBAM
#SBATCH --mem=250GB
#SBATCH -N 1
#SBATCH --ntasks-per-node=40

WkDir="/projects/fs1_no-backup/ls1-andersson/ANALYSIS/qirui/atac-seq/mouse/02-BAM"
Samples=("m16-8" "m16-11" "m16-12" "m16-14" "m16-15" "m16-16" "m16-17" "m16-19" "m16-23" "m16-25" "m16-26" "m16-27" "m16-30" "m16-31" "m16-33" "m16-35")

#------------------------------------------------------------------------------
module load GCC/7.3.0-2.30 SAMtools/1.9

for sample in ${Samples[@]}
do
	time=`date "+%Y-%m-%d %H:%M:%S"`
c	echo -e "\t$time\t$sample: Removing chrM, inappropriate alignments, and alignments with insert size > 2kb ..."

	# remove chrM and random chromosomes
	samtools view -@ 20 -h ${WkDir}/s1_rmdupBAM/${sample}_rmdup.bam |grep -v -E 'chrM|chr.*_random' |samtools view -@ 20 -S -b > ${WkDir}/s2_readyBAM/${sample}.cleanChr.bam

	# remove inappropriate alignments ($9==0) and alignments with insert size>2kb
	samtools view -@ 20 -h ${WkDir}/s2_readyBAM/${sample}.cleanChr.bam |awk '/^@/ {print $0;next;} $9 > 9 && $9 <= 2009 {print $0;next;} $9 > -2000 && $9 < 0 {print $0;next;}' > ${WkDir}/s2_readyBAM/${sample}.approAligned.sam

	# shift alignments (ONLY AFTER removed inappropriately mapped reads! otherwise perl script gives warning!)
	time=`date "+%Y-%m-%d %H:%M:%S"`
	echo -e "\t$time\t$sample: Shifting alignments ..."

	shiftReads4ATACseq.pl ${WkDir}/s2_readyBAM/${sample}.approAligned.sam ${WkDir}/s2_readyBAM/${sample}.shift.sam

	# remove low-quality alignments (<Q30), sort and index BAMs
	time=`date "+%Y-%m-%d %H:%M:%S"`
	echo -e "\t$time\t$sample: Removing low-quality alignments (<Q30), sort and index BAM (now is real clean) ..."

c	samtools view -@ 16 -S -b -q 30 ${WkDir}/s2_readyBAM/${sample}.shift.sam |samtools sort -@ 16 > ${WkDir}/s2_readyBAM/${sample}.bam
	  # This is real clean BAM/alignments!
	samtools index -@ 16 ${WkDir}/s2_readyBAM/${sample}.bam

	time=`date "+%Y-%m-%d %H:%M:%S"`
	echo -e "\t$time\tDone with sample $sample\n"

	# remove temporary files
	rm ${WkDir}/s2_readyBAM/${sample}.cleanChr.bam
	rm ${WkDir}/s2_readyBAM/${sample}.approAligned.sam
	rm ${WkDir}/s2_readyBAM/${sample}.shift.sam

done
time=`date "+%Y-%m-%d %H:%M:%S"`
echo -e "\t$time\tFinal done!\n"
module purge

