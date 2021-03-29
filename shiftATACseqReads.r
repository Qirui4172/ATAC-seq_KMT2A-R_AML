#!/usr/bin/env Rscript

#=====================================================================================================================================
Usage<-function(){
	cat("Usage: Rscript shiftATACseqReads.r [inBAM] [prefix] [outBAM]\n",

	"Parameters:\n",
	"[inBAM]	Input BAM file, should be sorted by chromosome coordinates and indexed\n",
	"[prefix]	Sample prefix, used for creating plot file, i.e. m16-8, then m16-8_plots.pdf will be generated\n",
	"[outBAM]	Output shifted BAM file, including output dir and file name, i.e. outdir/m16-8_shifted.bam\n\n",

	"Function:\n",
	"Use ATACseqQC package to shift aligned ATAC-seq reads +4 bp for positive and -5 bp for negative strand respectively, to account for the 9-bp duplication created by DNA repair of the nick by Tn5 transposase. Before runing the script, Bioconductor and R4.0 should be loaded:\n",
	"\tmodule load GCC/9.3.0  OpenMPI/4.0.3 R-bundle-Bioconductor/3.11-R-4.0.0\n",
	"\tmodule load GCC/9.3.0  OpenMPI/4.0.3 R/4.0.0\n\n",

	"Example:\n",
	"Rscript shiftATACseqReads.r m16-8.sort.bam shifted/m16-8_shifted.bam\n\n",

	"Author & date:\n",
	"Qirui Zhang (qirui.zhang@med.lu.se)\n",
	"15-11-2020\n\n"
	)
}

args<-commandArgs(TRUE)
if(length(args)!=3){Usage();quit();}

#=====================================================================================================================================
cat("\n=====================================================================================================================================\n")
# Load libraries and data
cat("Start analyze sample ", args[2], "\nLoading libraries and data\n\n")
#local.lib<-"/home/qirui/R/x86_64-pc-linux-gnu-library/4.0"
#.libPaths(unique(c(local.lib, .libPaths())))

library(ATACseqQC)
library(BSgenome.Mmusculus.UCSC.mm9)
library(TxDb.Mmusculus.UCSC.mm9.knownGene)
library(Rsamtools)

bamfile<-args[1]
prefix<-args[2]
outbam<-args[3]

outplot<-paste(dirname(outbam), "/", prefix, "_plots.pdf", sep="")
pdf(outplot)


#=====================================================================================================================================
# Shift reads
time<-format(Sys.time(), format='%H:%M:%S %Y-%m-%d')
cat(time, "\t", prefix, "Shifting reads ...\n\n")
possibleTag<-list("integer"=c("AM", "AS", "CM", "CP", "FI", "H0", "H1", "H2", "HI", "IH", "MQ", "NH", "NM", "OP", "PQ", "SM", "TC", "UQ"), "character"=c("BC", "BQ", "BZ", "CB", "CC", "CO", "CQ", "CR", "CS", "CT", "CY", "E2", "FS", "LB", "MC", "MD", "MI", "OA", "OC", "OQ", "OX", "PG", "PT", "PU", "Q2", "QT", "QX", "R2", "RG", "RX", "SA", "TS", "U2"))
bamTop1million<-scanBam(BamFile(bamfile, yieldSize=1000000), param=ScanBamParam(tag=unlist(possibleTag)))[[1]]$tag
tags<-names(bamTop1million)[lengths(bamTop1million)>0]

chr<-seqinfo(Mmusculus)@seqnames[1:21]
which<-as(seqinfo(Mmusculus)[chr], "GRanges")
gal<-readBamFile(bamFile=bamfile, tag=tags, which=which, asMates=TRUE, bigFile=TRUE)
gal.shift<-shiftGAlignmentsList(gal, outbam=outbam)

time<-format(Sys.time(), format='%H:%M:%S %Y-%m-%d')
cat(time, "\t", "Done with shifting sample ", prefix, "\n\n")

dev.off()
#=====================================================================================================================================
# Fragment size distribution
time<-format(Sys.time(), format='%H:%M:%S %Y-%m-%d')
cat(time, "Ploting fragment size distribution ...\n\n")
fragSize<-fragSizeDist(bamfile, prefix)

