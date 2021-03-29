#!/usr/bin/env Rscript

#===============================================================================================
Usage<-function(){
	cat("Usage: Rscript findDiffPeaks.r [readCounts] [sampleInfo] [removeSampl] [topPeakNum]\n\n",

	"Parameters:\n",
	"[readCounts]		mmATAC_readcounts.mx, including all 16 samples\n",
	"[sampleInfo]		sample.info, containing all 16 samples\n",
	"[removeSampl]		remove which sample, i.e. m16-11, or no for not removing\n",
	"[topPeakNum]		choose topNum or all of peaks for analysis, i.e. 80000, all\n\n",

	"Example (recommended R-3.3.1):\n",
	"Rscript findDiffPeaks.r mmATAC_readcounts.mx sample.info m16-11 all\n",
	"Rscript findDiffPeaks.r mmATAC_readcounts.mx sample.info no 80000\n\n",

	"Function: Use DESeq2 to find differentially accessible regions (DARs) and generate plots.\n",
	"Qirui Zhang (qirui.zhang@med.lu.se)\n",
	"14-05-2020\n\n"
	)
}

args<-commandArgs(TRUE)
if(length(args)!=4){Usage();quit();}


#===============================================================================================
cat("\n=========================================================================================\n")
# Step1. Load libraries
time<-format(Sys.time(), format='%H:%M:%S %Y-%m-%d')
cat(time, "Start differential analysis", "\n\n")
cat("Step1. Load libraries", "\n\n")

library(DESeq2)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
pdf("diffPeaks_plots.pdf")
options(scipen=200)


#===============================================================================================
cat("\n=========================================================================================\n")
# Step2. Read data
cat("Step2. Reading data", "\n\n")

readcounts<-read.table(args[1], header=T, stringsAsFactors=F)
rownames(readcounts)<-readcounts$peakID
readcounts<-readcounts[,5:20]
sampleinfo.origin<-read.table(args[2], header=T, stringsAsFactors=F)
sampleinfo.origin$Mutation<-as.factor(sampleinfo.origin$Mutation)
colnames(readcounts)<-sampleinfo.origin$Sample

rmsample<-args[3]
if(rmsample=="no"){readcounts=readcounts}else{readcounts=readcounts[,!names(readcounts) %in% c(rmsample)]}
sampleinfo<-sampleinfo.origin[which(sampleinfo.origin$Sample!=rmsample),]


#===============================================================================================
cat("\n=========================================================================================\n")
# Step3. Run DESeq2
cat("Step3. Starting to run DESeq2\n\n")

cat("Generating metadata...\n")
dds<-DESeqDataSetFromMatrix(countData=readcounts, colData=sampleinfo, design=~Mutation)
dds<-dds[rowSums(counts(dds))>=1,]
dds<-DESeq(dds)

cat("Normalizing read counts...\n")
normalized.counts<-as.data.frame(counts(dds, normalized=TRUE))
normalized.counts$mean<-apply(normalized.counts, 1, mean)
cat("\nQuantile of mean values of normalized reads count:", "\n")
quantile(normalized.counts$mean)


topNum<-args[4]
if(topNum=="all"){
	cat("\nUsing all peaks...\n")
}else{
	topNum<-as.numeric(topNum)
	topPercent<-round((topNum/nrow(readcounts))*100, 1)
	options(scipen=20)
	cat(paste("\nExtracting top", topNum, " peaks (", topPercent, "% of total peaks) ...", sep=""), "\n")
	keep.rownum<-head(order(normalized.counts$mean, decreasing=T), topNum)
	keep.rownum<-sort(keep.rownum)
	keep.peakid<-rownames(normalized.counts[keep.rownum,])
	normalized.counts<-normalized.counts[keep.peakid,]
	dds<-dds[keep.peakid,]
}
write.table(as.data.frame(normalized.counts), "atac-seq_deseq2norm.tsv", row.names=T, col.names=T, quote=F, sep="\t")

sf<-sizeFactors(dds)
write.table(as.data.frame(sf), "sizefactor.tsv", row.names=T, col.names=T, quote=F, sep="\t")


#===============================================================================================
cat("\n=========================================================================================\n")
# Step4. Transform data with DESeq2::vst
cat("Step4. Transforming data with DESeq2::vst", "\n\n")

vst<-varianceStabilizingTransformation(dds, blind=FALSE)
write.table(as.data.frame(assay(vst)), "atac-seq_vst.tsv", quote=F, row.names=T, col.names=T, sep="\t")


#===============================================================================================
cat("\n=========================================================================================\n")
# Step5. Make plots of sample correlation and variance
cat("Step5. Making plots of sample correlation and variance","\n\n")

# Plot correlation heatmap
cat("Plotting correlation heatmap ...", "\n")
sampleDist<-dist(t(assay(vst)))
sampleDist.mx<-as.matrix(sampleDist)
rownames(sampleDist.mx)<-NULL
colnames(sampleDist.mx)<-sampleinfo$Sample
colors<-colorRampPalette(rev(brewer.pal(9,"Blues")))(255)
pheatmap(as.data.frame(sampleDist.mx), clustering_distance_rows=sampleDist, clustering_distance_cols=sampleDist, col=colors, show_rownames=F)

# Plot PCA using DESeq2::plotPCA
cat("Plotting PCA using DESeq2::plotPCA ...", "\n")
data<-plotPCA(vst, intgroup="Mutation", returnData=TRUE)
percentVar<-round(100*attr(data, "percentVar"))
pca.plot1<-ggplot(data, aes(PC1, PC2, color=Mutation))+geom_point(size=3)+theme_bw()+theme(panel.grid=element_blank(),panel.border=element_rect(size=1),text=element_text(size=18),axis.text=element_text(size=15),legend.text=element_text(size=15))+xlab(paste0("PC1: ",percentVar[1],"% variance"))+ylab(paste0("PC2: ",percentVar[2],"% variance"))+coord_fixed(ratio=2.5)
pca.plot1
pca.plot1+geom_text(aes(label=rownames(data)),size=3)

# Plot PCA using in-house script (recommended!)
cat("Plotting PCA using in-house script ...", "\n")
pca.results<-prcomp(t(as.data.frame(assay(vst))), scale=TRUE)
percentVar<-round(((pca.results$sdev^2)/sum(pca.results$sdev^2))*100, 1)
data<-as.data.frame(pca.results$x)
data$Mutation<-sampleinfo$Mutation
pca.plot2<-ggplot(data, aes(PC1, PC2, colour=Mutation))+geom_point(size=4, shape=17)+xlab(paste("PC1: ", percentVar[1], "% variance", sep=""))+ylab(paste("PC2: ", percentVar[2], "% variance", sep=""))+theme_bw()+theme(panel.grid=element_blank(), panel.border=element_rect(colour="black",size=1), axis.title.x=element_text(size=18), axis.title.y=element_text(size=18), axis.text.x=element_text(size=15), axis.text.y=element_text(size=15), legend.title=element_text(size=15), legend.text=element_text(size=15))+scale_colour_brewer(palette="Set1")+coord_fixed(ratio=1)
pca.plot2
pca.plot2+geom_text(aes(label=rownames(data)),size=3)

#===============================================================================================
cat("\n=========================================================================================\n")
# Step6. Extract results of comparison pairs
cat("Step6. Extract results of comparison pairs", "\n\n")

getDiffPeaks<-function(dds, treatmentString, controlString){
	trt<-treatmentString
	ctl<-controlString
	res<-results(dds, contrast=c("Mutation", trt, ctl), alpha=0.05)
	title<-paste("MA plot of logFC in ", trt, "vs", ctl, sep="")
	plotMA(res, main={title}, ylim=c(-10,10))
	summary(res)
	upfc2.num<-sum(res$padj<0.05 & res$log2FoldChange>=log2(2), na.rm=TRUE)
	downfc2.num<-sum(res$padj<0.05 & -(res$log2FoldChange)>=log2(2), na.rm=TRUE)
	comp<-paste(trt, "vs", ctl, sep="")
	cat("Increased peaks with |foldchange|>2: ", upfc2.num, "\n", sep="")
	cat("Decreased peaks with |foldchange|>2: ", downfc2.num, "\n", sep="")
	dbpeak<-res[!is.na(res$padj) & res$padj<0.05,]
	dbpeak.fc2<-dbpeak[abs(dbpeak$log2FoldChange)>=1,]
	output<-paste("diffPeaks_", comp, ".tsv", sep="")
	write.table(as.data.frame(dbpeak.fc2), output, row.names=T, col.names=T, quote=F, sep="\t")
	out=list(res, dbpeak.fc2)
	return(out)
}

# FLT3-ITD vs KMT2A-MLLT3
cat("--------------------------------------------------\n")
cat("Comparing FLT3-ITD vs KMT2A-MLLT3 ...\n")
diff.results<-getDiffPeaks(dds, "fi", "km")
res.fivskm<-diff.results[[1]]
dbpeak.fivskm.fc2<-diff.results[[2]]

# FLT3-N676K vs KMT2A-MLLT3
cat("--------------------------------------------------\n")
cat("Comparing FLT3-N676K vs KMT2A-MLLT3 ...\n")
diff.results<-getDiffPeaks(dds, "fn", "km")
res.fnvskm<-diff.results[[1]]
dbpeak.fnvskm.fc2<-diff.results[[2]]

# NRAS-G12D vs KMT2A-MLLT3
cat("--------------------------------------------------\n")
cat("Comparing NRAS-G12D vs KMT2A-MLLT3 ...\n")
diff.results<-getDiffPeaks(dds, "ng", "km")
res.ngvskm<-diff.results[[1]]
dbpeak.ngvskm.fc2<-diff.results[[2]]

# FLT3-ITD vs FLT3-N676K
cat("--------------------------------------------------\n")
cat("Comparing FLT3-ITD vs FLT3-N676K ...\n")
diff.results<-getDiffPeaks(dds, "fi", "fn")
res.fivsfn<-diff.results[[1]]
dbpeak.fivsfn.fc2<-diff.results[[2]]

# FLT3-N676K vs FLT3-ITD
cat("--------------------------------------------------\n")
cat("Comparing FLT3-N676K vs FLT3-ITD ...\n")
diff.results<-getDiffPeaks(dds, "fn", "fi")
res.fnvsfi<-diff.results[[1]]
dbpeak.fnvsfi.fc2<-diff.results[[2]]

# FLT3-ITD vs NRAS-G12D
cat("--------------------------------------------------\n")
cat("Comparing FLT3-ITD vs NRAS-G12D ...\n")
diff.results<-getDiffPeaks(dds, "fi", "ng")
res.fivsng<-diff.results[[1]]
dbpeak.fivsng.fc2<-diff.results[[2]]

# NRAS-G12D vs FLT3-ITD
cat("--------------------------------------------------\n")
cat("Comparing NRAS-G12D vs FLT3-ITD ...\n")
diff.results<-getDiffPeaks(dds, "ng", "fi")
res.ngvsfi<-diff.results[[1]]
dbpeak.ngvsfi.fc2<-diff.results[[2]]

# FLT3-N676K vs NRAS-G12D
cat("--------------------------------------------------\n")
cat("Comparing FLT3-N676K vs NRAS-G12D ...\n")
diff.results<-getDiffPeaks(dds, "fn", "ng")
res.fnvsng<-diff.results[[1]]
dbpeak.fnvsng.fc2<-diff.results[[2]]

# NRAS-G12D vs FLT3-N676K
cat("--------------------------------------------------\n")
cat("Comparing NRAS-G12D vs FLT3-N676K ...\n")
diff.results<-getDiffPeaks(dds, "ng", "fn")
res.ngvsfn<-diff.results[[1]]
dbpeak.ngvsfn.fc2<-diff.results[[2]]


#===============================================================================================
cat("\n=========================================================================================\n")
# Step7.Make volcano plots of differentially accessible regions/peaks
cat("Step7. Making volcano plots of differentially accessible regions/peaks", "\n\n")

VolcanoPlot<-function(res, compareString){
	volcano<-as.data.frame(res)
	volcano$significant<-as.factor(ifelse(!is.na(volcano$padj) & volcano$padj < 0.05 & abs(volcano$log2FoldChange)>=1, ifelse(volcano$log2FoldChange>=1, "Up", "Down"), "No"))
	volcano$padj<-ifelse(is.na(volcano$padj),1,volcano$padj)
	volcano$group<-rep(0)
	volcano[which(volcano$padj >= 0.05 | (volcano$padj < 0.05 & abs(volcano$log2FoldChange) < 1)),8]=rep(1)
	volcano[which(volcano$padj < 0.05 & volcano$log2FoldChange < -1),8]=rep(2)
	volcano[which(volcano$padj < 0.05 & volcano$log2FoldChange > 1),8]=rep(3)
	comp=compareString
	title=paste("Volcano plot of ", comp, sep="")
	p<-ggplot(volcano,aes(log2FoldChange,-log10(padj)))+geom_point(data=volcano[which(volcano$group==1),],color="gray",alpha=0.75)+geom_point(data=volcano[which(volcano$group==2),],color="royalblue3",alpha=0.75)+geom_point(data=volcano[which(volcano$group==3),],color=brewer.pal(11,"RdYlBu")[2],alpha=0.75)
	p+labs(title={title},x="log2FoldChange",y="-log10(padj)")+geom_hline(yintercept=-log10(0.05),linetype=2,color="gray80")+geom_vline(xintercept=c(-log2(2),log2(2)),linetype=2,color="gray80")+theme_bw()+theme(text=element_text(size=15),panel.grid=element_blank(),panel.border=element_rect(size=1))
}

# FLT3-ITD vs KMT2A-MLLT3
cat("Ploting FLT3-ITD vs KMT2A-MLLT3 ...", "\n")
VolcanoPlot(res.fivskm, "fivskm")

# FLT3-N676K vs KMT2A-MLLT3
cat("Ploting FLT3-N676K vs KMT2A-MLLT3 ...", "\n")
VolcanoPlot(res.fnvskm, "fnvskm")

# NRAS-G12D vs KMT2A-MLLT3
cat("Ploting NRAS-G12D vs KMT2A-MLLT3 ...", "\n")
VolcanoPlot(res.ngvskm, "ngvskm")

# FLT3-ITD vs FLT3-N676K
cat("Ploting FLT3-ITD vs FLT3-N676K ...", "\n")
VolcanoPlot(res.fivsfn, "fivsfn")

# FLT3-N676K vs FLT3-ITD
cat("Ploting FLT3-N676K vs FLT3-ITD ...", "\n")
VolcanoPlot(res.fnvsfi, "fnvsfi")

# FLT3-ITD vs NRAS-G12D
cat("Ploting FLT3-ITD vs NRAS-G12D ...", "\n")
VolcanoPlot(res.fivsng, "fivsng")

# NRAS-G12D vs FLT3-ITD
cat("Ploting NRAS-G12D vs FLT3-ITD ...", "\n")
VolcanoPlot(res.ngvsfi, "ngvsfi")

# FLT3-N676K vs NRAS-G12D
cat("Ploting FLT3-N676K vs NRAS-G12D ...", "\n")
VolcanoPlot(res.fnvsng, "fnvsng")

# NRAS-G12D vs FLT3-N676K
cat("Ploting NRAS-G12D vs FLT3-N676K ...", "\n")
VolcanoPlot(res.ngvsfn, "ngvsfn")


#===============================================================================================
cat("\n=========================================================================================\n")
# Step8. Make volcano plots in modified version
cat("Step8. Make volcano plots in modified version of differentially accessible regions/peaks", "\n\n")

VolcanoPlotV2<-function(res, compareString){
	volcano<-as.data.frame(res)
	volcano$significant<-as.factor(ifelse(!is.na(volcano$padj) & volcano$padj < 0.05 & abs(volcano$log2FoldChange)>=1, "yes", "no"))
	volcano$padj<-ifelse(is.na(volcano$padj),1,volcano$padj)

	comp<-compareString
	title=paste("Volcano plot in modified version of ", comp, sep="")
	p<-ggplot(volcano, aes(x=log2FoldChange, y=log2(baseMean)))+geom_point(data=subset(volcano, significant=="no"), colour="gray")+geom_point(data=subset(volcano, significant=="yes"), aes(colour=padj), alpha=0.8)+scale_colour_gradientn("Padj", colours=c(rev(brewer.pal(9, "OrRd"))[c(2,5,6,7,8,9)]), limits=c(0, 0.05), breaks=c(0, 0.01, 0.02, 0.03, 0.04, 0.05), labels=c("0", "0.01", "0.02", "0.03", "0.04", "0.05"))
	p+labs(title={title},x="log2FoldChange",y="log2baseMean")+geom_vline(xintercept=c(-1,1),linetype=2,color="black")+theme_bw()+theme(text=element_text(size=15),axis.line=element_line(colour="black"),panel.grid=element_blank(),panel.border=element_blank())+coord_flip(xlim=c(-5, 5))
}

# FLT3-ITD vs KMT2A-MLLT3
cat("Ploting FLT3-ITD vs KMT2A-MLLT3 ...", "\n")
VolcanoPlotV2(res.fivskm, "fivskm")

# FLT3-N676K vs KMT2A-MLLT3
cat("Ploting FLT3-N676K vs KMT2A-MLLT3 ...", "\n")
VolcanoPlotV2(res.fnvskm, "fnvskm")

# NRAS-G12D vs KMT2A-MLLT3
cat("Ploting NRAS-G12D vs KMT2A-MLLT3 ...", "\n")
VolcanoPlotV2(res.ngvskm, "ngvskm")

# FLT3-ITD vs FLT3-N676K
cat("Ploting FLT3-ITD vs FLT3-N676K ...", "\n")
VolcanoPlotV2(res.fivsfn, "fivsfn")

# FLT3-N676K vs FLT3-ITD
cat("Ploting FLT3-N676K vs FLT3-ITD ...", "\n")
VolcanoPlotV2(res.fnvsfi, "fnvsfi")

# FLT3-ITD vs NRAS-G12D
cat("Ploting FLT3-ITD vs NRAS-G12D ...", "\n")
VolcanoPlotV2(res.fivsng, "fivsng")

# NRAS-G12D vs FLT3-ITD
cat("Ploting NRAS-G12D vs FLT3-ITD ...", "\n")
VolcanoPlotV2(res.ngvsfi, "ngvsfi")

# FLT3-N676K vs NRAS-G12D
cat("Ploting FLT3-N676K vs NRAS-G12D ...", "\n")
VolcanoPlotV2(res.fnvsng, "fnvsng")

# NRAS-G12D vs FLT3-N676K
cat("Ploting NRAS-G12D vs FLT3-N676K ...", "\n")
VolcanoPlotV2(res.ngvsfn, "ngvsfn")


#===============================================================================================
cat("\n=========================================================================================\n")
# Step9. Plot log2FoldChange density curves
cat("Step9. Plot log2FoldChange density curves", "\n\n")

DensityPlot<-function(res, compareString){
	density<-as.data.frame(res)
	comp<-compareString
	title=paste("log2FoldChange density plot of ", comp, sep="")

	ggplot(density, aes(x=log2FoldChange))+geom_density(fill="gray")+labs(title={title},x="log2FoldChange",y="Density")+geom_vline(xintercept=c(-1,1),linetype=2,color="black")+theme_bw()+theme(title=element_text(size=12),text=element_text(size=15),axis.line=element_line(colour="black"),panel.grid=element_blank(),panel.border=element_blank(),aspect.ratio=5)+coord_flip(xlim=c(-5, 5))
}

# FLT3-ITD vs KMT2A-MLLT3
cat("Ploting FLT3-ITD vs KMT2A-MLLT3 ...", "\n")
DensityPlot(res.fivskm, "fivskm")

# FLT3-N676K vs KMT2A-MLLT3
cat("Ploting FLT3-N676K vs KMT2A-MLLT3 ...", "\n")
DensityPlot(res.fnvskm, "fnvskm")

# NRAS-G12D vs KMT2A-MLLT3
cat("Ploting NRAS-G12D vs KMT2A-MLLT3 ...", "\n")
DensityPlot(res.ngvskm, "ngvskm")

# FLT3-ITD vs FLT3-N676K
cat("Ploting FLT3-ITD vs FLT3-N676K ...", "\n")
DensityPlot(res.fivsfn, "fivsfn")

# FLT3-N676K vs FLT3-ITD
cat("Ploting FLT3-N676K vs FLT3-ITD ...", "\n")
DensityPlot(res.fnvsfi, "fnvsfi")

# FLT3-ITD vs NRAS-G12D
cat("Ploting FLT3-ITD vs NRAS-G12D ...", "\n")
DensityPlot(res.fivsng, "fivsng")

# NRAS-G12D vs FLT3-ITD
cat("Ploting NRAS-G12D vs FLT3-ITD ...", "\n")
DensityPlot(res.ngvsfi, "ngvsfi")

# FLT3-N676K vs NRAS-G12D
cat("Ploting FLT3-N676K vs NRAS-G12D ...", "\n")
DensityPlot(res.fnvsng, "fnvsng")

# NRAS-G12D vs FLT3-N676K
cat("Ploting NRAS-G12D vs FLT3-N676K ...", "\n")
DensityPlot(res.ngvsfn, "ngvsfn")


#===============================================================================================
cat("\n=========================================================================================\n")
# Step10. Plot heatmap of all differentially accessible regions (DARs)
cat("Step10. Plot heatmap of differentially accessible regions (DARs) generated from all above comparisons with a |FoldChange|>=2", "\n\n")

dbpeak.id<-unique(c(rownames(dbpeak.ngvskm.fc2), rownames(dbpeak.fivskm.fc2), rownames(dbpeak.fnvskm.fc2), rownames(dbpeak.fnvsfi.fc2), rownames(dbpeak.ngvsfi.fc2), rownames(dbpeak.ngvsfn.fc2)))
dbpeak.num<-length(dbpeak.id)
cat("Total number of differentially accessible regions/peaks: ", dbpeak.num, "\n")
vst.dbpeak<-vst[dbpeak.id,]
vst.dbpeak.df<-as.data.frame(assay(vst.dbpeak))
#nb.clust<-NbClust(vst.dbpeak.df, distance="euclidean", min.nc=2, max.nc=15, method="kmeans", index="alllong", alphaBeale=0.1)
  # Give recommended cluster numbers

anno.label<-c(rep("FLT3-ITD",4), rep("FLT3-N676K",4), rep("NRAS-G12D",4), rep("Empty GFP",4))
names(anno.label)<-sampleinfo.origin$Sample
anno.label<-anno.label[which(names(anno.label)!=rmsample)]
anno.label<-as.data.frame(anno.label)
names(anno.label)<-"Mutation"
anno.color<-list(Mutation=c("FLT3-ITD"="black", "FLT3-N676K"="black", "NRAS-G12D"="black", "Empty GFP"="black"))

pheatmap(vst.dbpeak.df, main="Heatmap of differentially accessible regions", scale="row", color=colorRampPalette(c(rev(brewer.pal(9,"Blues")[c(3:9)]), brewer.pal(9, "OrRd")[c(3:9)]))(100), cluster_rows=TRUE, cluster_cols=FALSE, show_rownames=FALSE, annotation_col=anno.label, annotation_colors=anno.color, annotation_names_col=F, border_color=NA)
dev.off()


#===============================================================================================
time<-format(Sys.time(), format='%H:%M:%S %Y-%m-%d')
cat(time, "Done with analysis!", "\n\n")


