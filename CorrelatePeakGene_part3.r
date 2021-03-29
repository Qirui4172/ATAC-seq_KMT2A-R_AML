#!/usr/bin/env Rscript


#===================================================================================
# Set parameters
Usage<-function(){
	cat("Usage: Rscript correlatePeakGene.r [peakList] [atacseqMx] [rnaseqMx] [treatment] [control] [outprefix]\n\n",

	"Parameters:\n",
	"[peakList]    Peak ID list in one column, e.g. peak_1, peak_25\n",
	"[atacseqMx]   ATAC-seq read counts matrix normalized by DESeq2 but not log2-transformed (without m16-11). It contains 18 columns, the first 15 columns are readcounts matrix, the last three columns are mean, chrID, peakCenter\n",
	"[rnaseqMx]    RNA-seq read counts matrixs normalized by log2(fpkm+0.1). It contains 19 columns, the first 17 columns are gene symbols and readcounts matrix, the last two columns are chrID, startCode\n",
	"[treatment]   Treatment short name, either of fi, fn, ng, km\n",
	"[control]     Control short name, either of fi, fn, ng, km\n",
	"[outprefix]   Output name prefix, e.g. fivskm_up\n\n",

	"Function:\n",
	"Correlate ATAC-seq peaks and RNA-seq genes within the same chromosome by computing Pearson correlation coefficients, p values, and FDRs (BH method). Output all correlation results including significant and not significant. A further filtration can be performed, i.e. on FDRs, or distance.\n",
	"Note: m16-11 is supposed to be excluded from ATAC-seq matrix (but not RNA-seq matrix) before using this script\n\n",

	"Example:\n",
	"Rscript correlatePeakGene.r fivskm_down.peaklist atac-seq_deseq2norm_forCorrelation.mx rna-seq_log2fpkm_forCorrelation.mx fi km fivskm_down\n",
	"# the output file \"fivskm_down_correlation.tsv\" will be generted\n\n",

	"Qirui Zhang (qirui.zhang@med.lu.se)\n",
	"Date: 05-10-2020\n\n"
	)
}

args<-commandArgs(TRUE)
if(length(args)!=6){Usage();quit();}


#===================================================================================
# Read parameters, load files, subset mutation samples

# read peak list
cat("Reading peak list ...", "\n\n")
peak.list<-read.table(args[1], header=F, stringsAsFactors=F)
peakID<-peak.list$V1

# create mutation-column table
cat("Creating mutation-column table ...", "\n\n")
mutation.column<-data.frame(start=c(1, 5, 9, 12), end=c(4, 8, 11, 15))
rownames(mutation.column)<-c("fi", "fn", "ng", "km")
treat<-args[4]
start1<-mutation.column[treat,1]
end1<-mutation.column[treat,2]
contr<-args[5]
start2<-mutation.column[contr,1]
end2<-mutation.column[contr,2]

# read ATAC-seq.mx and subset
cat("Reading ATAC-seq matrix and subset mutation samples ...", "\n\n")
atac<-read.table(args[2], header=T, stringsAsFactors=F)
atac<-subset(atac, select=-mean)
atac[,c(1:15)]<-log2(atac[,c(1:15)]+0.5)
atac<-atac[,c(start1:end1, start2:end2, 16, 17)]
atac<-atac[peakID,]

# read RNA-seq fpkm.mx and subset
cat("Reading RNA-seq matrix and subset mutation samples ...", "\n\n")
rna<-read.table(args[3], header=T, stringsAsFactors=F)
rownames(rna)<-rna$Gene
rna<-rna[,-1]
rna<-rna[,-10]	# remove m16-11
rna<-rna[,c(start1:end1, start2:end2, 16, 17)]
if(treat=="ng"){
	rna<-rna[which(rowSums(rna[,1:7]==-1)<7),]
}else{
	rna<-rna[which(rowSums(rna[,1:8]==-1)<8),]
}

#===================================================================================
# Create coefficient-pvalue matrix
cat("Preparing an empty matrix for saving coefficients and pvalues ...","\n")
total_row.num=0
for(i in unique(atac$chrID)){
	chr_peak.num=nrow(atac[which(atac$chrID==i),])	# peak numbers in chr1 (i.e.)
	chr_gene.num=nrow(rna[which(rna$chrID==i),])	# gene numbers in chr1 (i.e.)
	total_row.num=total_row.num+chr_peak.num*chr_gene.num
}
correlation.matrix<-matrix(nrow=total_row.num, ncol=3)

#===================================================================================
# Create peak-gene dataframe and compute Pearson correlation coefficients and p values
cat("Creating peak-gene dataframe and computing Pearson correlation coefficients and p values ...", "\n\n")

k=0
correlation<-data.frame(peakID=character(), geneID=character(), stringsAsFactors=FALSE)
for(i in 1:nrow(atac)){
	# create peak-gene dataframe
	chr=atac[i,"chrID"]
	rna_sub=rna[which(rna$chrID==chr),]
	num=nrow(rna_sub)
	peak=rownames(atac[i,])
	peak.list=rep(peak, num)
	gene.list=rownames(rna_sub)
	df.tmp=data.frame(peakID=peak.list, geneID=gene.list)
	correlation=rbind(correlation, df.tmp)
	rm(df.tmp)

	# compute Pearson correlation coefficients and p values
	if(treat=="ng"){atac_for_cor=t(atac[i,c(1:7)])}else{atac_for_cor=t(atac[i,c(1:8)])}
	for(j in 1:nrow(rna_sub)){
		k=k+1
		if(treat=="ng"){rna_for_cor=t(rna_sub[j,c(1:7)])}else{rna_for_cor=t(rna_sub[j,c(1:8)])}
		cor.result=cor.test(atac_for_cor, rna_for_cor)
		coeff=cor.result$estimate[[1]]
		p.value=cor.result$p.value
		dist=rna_sub[j,"startCode"]-atac[i,"peakCenter"]
		correlation.matrix[k,]=c(coeff, p.value, dist)
	}
}
cat("Total row number: ", k, "\n\n")
correlation$peakID=as.character(correlation$peakID)
correlation$geneID=as.character(correlation$geneID)
correlation=cbind(correlation, as.data.frame(correlation.matrix))
colnames(correlation)[3:5]=c("coefficient", "pvalue", "distance")


#===================================================================================
# Adjust p values using "BH" method
cat("Adjusting p values using BH method ...", "\n\n")

correlation$BH_padj<-rep(0)
accumulated.num=1
for(i in peakID){
	correlation_sub=correlation[which(correlation$peakID==i),]  # adjust pvalues of one peak at one time, not adjust pvalues of all peaks together.
	num=nrow(correlation_sub)
	start.num=accumulated.num
	end.num=accumulated.num+num-1
	padj=p.adjust(correlation_sub$pvalue, method="BH")
	correlation[start.num:end.num,"BH_padj"]=padj
	accumulated.num=accumulated.num+num
}


#===================================================================================
# Output results
output<-paste(args[6], "_correlation.tsv", sep="")
cat("Outputing results to file", output, "\n\n")
write.table(correlation, output, col.names=T, row.names=F, quote=F, sep="\t")

cat("Done with correlation analysis!", "\n\n")

