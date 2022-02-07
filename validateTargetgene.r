#!/usr/bin/env Rscript

args<-commandArgs(TRUE)
if (length(args)!=5){
	cat("Usage: Rscript validateTargetgene.r [corSig] [rnaseqMx] [treatment] [control] [outprefix]","\n\n",

	"Parameters:\n",
	"[corSig]       Significantly correlated peak-gene pairs (pvalue<0.05, distance<=1000000) after filtering result of \"correlatePeakGene.r\", containing 6 columns: peakID, geneID, coefficient, pvalue, distance, BH_padj\n",
	"[rnaseqMx]     rna-seq_log2fpkm_forCorrelation.mx, rna-seq read count metrics normalized by log2(fpkm+0.1), the last columns could be chrID or startSite\n",
	"[treatment]    Treatment genotype name, either of fi, fn, ng, km\n",
	"[control]      Control genotype name, either of fi, fn, ng, km\n",
	"[outprefix]    Output name prefix, e.g. fivskm_up\n\n",

	"Example:\n",
	"Rscript validateTargetgene.r fivsfn_down_correlation-sig.tsv rna-seq_log2fpkm_forCorrelation.mx fi fn fivsfn_down\n",
	"# the file \"fivsfn_down_validation.tsv\" will be generated\n\n",

	"Function:\n",
	"Validate transcript level of target genes between treatment and control using t-test, generate foldChange and p-value for each peak-gene pairs.\n\n",

	"Qirui Zhang (qirui.zhang@med.lu.se)\n",
	"Date: 14-10-2020\n\n"
	)
	quit()
}


#===================================================================================

genotype.column<-data.frame(start=c(1, 5, 9, 13), end=c(4, 8, 12, 16))
rownames(genotype.column)<-c("fi", "fn", "ng", "km")
treat<-args[3]
start1<-genotype.column[treat,1]
end1<-genotype.column[treat,2]
contr<-args[4]
start2<-genotype.column[contr,1]
end2<-genotype.column[contr,2]


correlation<-read.table(args[1], header=T, stringsAsFactors=F)

rna<-read.table(args[2], header=T, stringsAsFactors=F)
rownames(rna)<-rna$Gene
rna<-rna[,2:17]
rna<-2^rna-0.5
rna<-as.data.frame(rna)
correlation$rna_fc<-rep(0)
correlation$rna_pvalue<-rep(0)


#===================================================================================

for(i in 1:nrow(correlation)){
	gene=correlation[i,"geneID"]
	if(!(gene %in% rownames(rna))){correlation[i,c("rna_fc","rna_pvalue")]=c("NA", "NA");next;}
	if(rowSums(rna[gene, c(start1:end1, start2:end2)]==0)==8){next;}
	rna_contr=t(rna[gene,start2:end2])
	rna_treat=t(rna[gene,start1:end1])
	t.result=t.test(rna_contr, rna_treat)
	gene.pvalue=t.result$p.value
	gene.fc=mean(rna_treat)/mean(rna_contr)
	correlation[i,c("rna_fc","rna_pvalue")]=c(gene.fc, gene.pvalue)
}

output<-paste(args[5], "_validation.tsv", sep="")
write.table(correlation, output, col.names=T, row.names=F, quote=F, sep="\t")

