#!/usr/bin/env Rscript


#===================================================================================
# Set parameters
Usage<-function(){
	cat("Usage: Rscript correlateDarGene.r [DARs_file] [atacseqMx] [rnaseqMx] [consen_peaks] [treatment] [control] [out_suffix]\n\n",

	"Parameters:\n",
	"[DARs_file]     DARs_fivskm_up.tsv, output file of \"findDiffpeaks.r\".\n",
	"[atacseqMx]     atac-seq_deseq2norm.tsv, output file of \"findDiffpeaks.r\".\n",
	"[rnaseqMx]      rna-seq_log2fpkm_forCorrelation.mx, it contains 19 columns, col1-17 are gene symbols and read counts matrix, col18-19 are chr & startcodon.\n",
    "[consen_peaks]  consensus_final.bed\n",
	"[treatment]     either of fi, fn, ng\n",
	"[control]       km\n",
	"[out_suffix]    output name suffix, e.g. fivskm_up\n\n",

	"Function:\n",
	"Correlate ATAC-seq peaks and RNA-seq genes within the same chromosome by computing Pearson correlation coefficients, p values, and FDRs (BH method). Output all correlation results including significant and not significant. A further filtration can be performed, i.e. on FDRs, or distance.\n",
	"Note: m16-11 is supposed to be excluded from ATAC-seq matrix (but not RNA-seq matrix) before using this script\n\n",

	"Example:\n",
	"Rscript correlateDarGene.r DARs_fivskm_up.tsv atac-seq_deseq2norm.tsv rna-seq_log2fpkm_forCorrelation.mx consensus_final.bed fi km fivskm_up\n",
	"# the output file \"correlation_fivskm_up.tsv\" will be generted\n\n",

	"Qirui Zhang (qirui.zhang@med.lu.se)\n",
	"Date: 05-10-2020\n\n"
	)
}

args<-commandArgs(TRUE)
if(length(args)!=7){Usage();quit();}


#===================================================================================
# Read parameters, load files, subset mutation samples
cat("loading arguments & reading files...\n\n")
dar.file<-args[1]  # DARs_fivskm_up_fc2.tsv
atac.file<-args[2]  # atac-seq_deseq2norm.tsv
rna.file<-args[3]  # rna-seq_log2fpkm_forCorrelation.mx
peak.file<-args[4]  # consensus_final.bed
treat<-args[5]  # fi
contr<-args[6]  # km
outsuffix<-args[7]  # fivskm_up

dars<-read.table(dar.file, header=T, stringsAsFactors=F)
peaks.sig<-rownames(dars)
consensus.peaks<-read.table(peak.file, header=F, stringsAsFactors=F)
colnames(consensus.peaks)<-c("chr","start","end","peakID")

# create mutation-column table
cat("Creating mutation-column table ...", "\n\n")
mutation.column<-data.frame(start=c(1, 5, 9, 12), end=c(4, 8, 11, 15))
rownames(mutation.column)<-c("fi", "fn", "ng", "km")
start1<-mutation.column[treat,1]
end1<-mutation.column[treat,2]
start2<-mutation.column[contr,1]
end2<-mutation.column[contr,2]

# read ATAC-seq.mx and subset
cat("Reading ATAC-seq matrix and subset mutation samples ...", "\n\n")
atac<-read.table(atac.file, header=T, stringsAsFactors=F)  # m16-11 was already removed from atac-seq matrix
atac$chr<-consensus.peaks$chr
atac$peakcenter<-apply(consensus.peaks[,c("start","end")],1,median)
atac$peakcenter<-as.integer(atac$peakcenter)
atac<-atac[,-16]
colnames(atac)[16:17]<-c("chr","peakcenter")
atac[,c(1:15)]<-log2(atac[,c(1:15)]+0.5)  # atac is DESeq2-normalized, raw reads are only divided by library size factors, and need further transformation before correlation.
atac<-atac[,c(start1:end1, start2:end2, 16:17)]
atac<-atac[peaks.sig,]

# read RNA-seq fpkm.mx and subset
cat("Reading RNA-seq matrix and subset mutation samples ...", "\n\n")
rna<-read.table(rna.file, header=T, stringsAsFactors=F)
rownames(rna)<-rna$Gene
rna<-rna[,-1]
rna<-rna[,-10]	# remove m16-11 (it was not removed from rna-seq matrix beforehand)
colnames(rna)[16:17]<-c("chr","startcodon")
rna<-rna[,c(start1:end1, start2:end2, 16:17)]
if(treat=="ng"){rna<-rna[which(rowSums(rna[,1:7]==-1)<7),]}else{rna<-rna[which(rowSums(rna[,1:8]==-1)<8),]}

#===================================================================================
# Count number of peaks x genes combinations
cat("Counting number of \"peaks x genes\" combinations...","\n\n")
total_row.num=0
for(i in unique(atac$chr)){chr_peak.num=nrow(atac[which(atac$chr==i),]); chr_gene.num=nrow(rna[which(rna$chr==i),]); total_row.num=total_row.num+chr_peak.num*chr_gene.num}
    # "chr_peak.num" is peak numbers on per chr in atac matrix; "chr_gene.num" is gene numbers on the same chr in rna matrix.

#===================================================================================
# Create peak-gene dataframe, coefficient-pvalue-distance matrix, and run Pearson correlation analysis
cat("Creating peak-gene dataframe & coefficient-pvalue-distance matrix...", "\n\n")

k=0
correlation<-data.frame(peakID=character(), geneID=character(), stringsAsFactors=FALSE)  # create peak-gene empty dataframe
correlation.values<-matrix(nrow=total_row.num, ncol=3)  # create coefficient-pvalue-distance empty matrix. The reason why coefficient-pvalue-distance is not created together with peak-gene dataframe is that it takes huge memories and becomes extremely time-consuming when running correlation analysis of so many peak-gene pairs. So it's better to seperately compute Pearson coefficient, pvalue, and distance and combine coefficient-pvalue-distance matrix with peak-gene dataframe in the end after all is done.

cat("Runing Pearson correlation analysis...\n\n")
for(i in 1:nrow(atac)){  # go through each peak
	# Create peak-gene dataframe (only create peak-gene dataframe, no correlation computation)
	peak.id=rownames(atac[i,])  # get peak name
	chr.id=atac[i,"chr"]  # get chr corresponding to each peak
	rna_sub=rna[which(rna$chr==chr.id),]  # extract all genes localized on that chr
	gene.list=rownames(rna_sub)  # these genes' names
	gene.num=nrow(rna_sub)  # how many of these genes
	peak.list=rep(peak.id, gene.num)
	df.tmp=data.frame(peakID=peak.list, geneID=gene.list)  # create a temp DataFrame containing 2 columns: peakID repeated 'num' times, and geneID localized on the same chr
	correlation=rbind(correlation, df.tmp)  # combine this temp DataFrame to the whole correlation DataFrame
	rm(df.tmp)

	# Compute Pearson correlation coefficients and p values
	if(treat=="ng"){atac_for_cor=t(atac[i,c(1:7)])}else{atac_for_cor=t(atac[i,c(1:8)])}  # transfer atac[i,] column to row
	for(j in 1:nrow(rna_sub)){  # go through each gene
		k=k+1
		if(treat=="ng"){rna_for_cor=t(rna_sub[j,c(1:7)])}else{rna_for_cor=t(rna_sub[j,c(1:8)])}  # transfer rna[i,j] column to row
		cor.result=cor.test(atac_for_cor, rna_for_cor)  # correlation computation
		coeff=cor.result$estimate[[1]]  # coefficient
		pvalue=cor.result$p.value  # p value
		dist=atac[i,"peakcenter"]-rna_sub[j,"startcodon"]  # distance between peak and gene, it should NOT use abs() here, because the distance distribution will be plotted in downstream analysis where distance<0 are required.
		correlation.values[k,]=c(coeff, pvalue, dist)  # save result of this peak-gene pair to correlation.values matrix
	}
}

cat("Total row number: ", k, "\n\n")
correlation$peakID=as.character(correlation$peakID)
correlation$geneID=as.character(correlation$geneID)
correlation=cbind(correlation, as.data.frame(correlation.values))  # combine peak-gene correlation dataframe with correlation.values matrix
colnames(correlation)[3:5]=c("coefficient", "pvalue", "distance")


#===================================================================================
# Adjust p values using "BH" method
cat("Adjusting p values using BH method ...", "\n\n")

correlation$BH_padj<-rep(0)
accumulated.num=1
for(i in peaks.sig){
	correlation_sub=correlation[which(correlation$peakID==i),]  # adjust pvalues of one peak at one time, not adjust pvalues of all peaks together.
	num=nrow(correlation_sub)  # how many genes for each peak
	start.num=accumulated.num  # the start row number of the first gene for this peak in the whole correlation dataframe
	end.num=accumulated.num+num-1  # the end row number of the last gene for this peak in the whole correlation dataframe
	padj=p.adjust(correlation_sub$pvalue, method="BH")  # adjust p values
	correlation[start.num:end.num,"BH_padj"]=padj  # add adjusted p value to the correlation dataframe
	accumulated.num=accumulated.num+num  # get accumulated.num ready for the next peak-gene pairs
}


#===================================================================================
# Output results
output<-paste("corr_", outsuffix, ".tsv", sep="")
cat("Outputing results to file", output, "\n\n")
write.table(correlation, output, col.names=T, row.names=F, quote=F, sep="\t")

cat("Done with correlation analysis!", "\n\n")

