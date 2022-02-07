#!/usr/bin/env Rscript

#===============================================================================================
Usage<-function(){
	cat("Usage: Rscript findDiffpeaks.r [reads_count] [sample_info] [adj_Pvalue] [fold_change] [fi_consPeaks] [fn_consPeaks] [ng_consPeaks] [km_consPeaks]\n\n",

	"Parameters:\n",
	"[reads_count]      mmATAC-seq_readscount.mx, including all 16 samples\n",
	"[sample_info]      sample.info, containing all 16 samples\n",
	"[adj_Pvalue]		adjusted Pvalue for defining DARs, i.e. 0.05\n",
    "[fold_change]      fold change for defining DARs, i.e. 2\n",
    "[fi_consPeaks]     fi_consensus.bed, for each comparison only use peaks from fi and km.\n",
    "[fn_consPeaks]     fn_consensus.bed, same as above.\n",
    "[ng_consPeaks]     ng_consensus.bed, same as above.\n",
    "[km_consPeaks]     km_consensus.bed\n\n",

	"Example: (R-3.6.0)\n",
	"Rscript findDiffpeaks.r mmATAC-seq_readscount.mx sample.info 0.05 2 fi_final.bed fn_final.bed ng_final.bed km_final.bed\n\n",

    "The following output files will be generated,\n",
    "  atac-seq_deseq2norm.tsv:    DESeq2-normalized atac-seq reads count matrix\n",
    "  atac-seq_vst.tsv:           vst-transformed atac-seq_deseq2norm.tsv, internally used by DESeq2 for plotting\n",
    "  DARs_fivskm_up.tsv, DARs_fivskm_down.tsv, DARs_fnvskm_up.tsv, DARs_fnvskm_down.tsv, DARs_ngvskm_up.tsv, DARs_ngvskm_down.tsv:\n",
    "                              up-/down-regulated DAR lists of FLT3-ITD/FLT3-N676K/NRAS-G12D vs KM\n",
    "  diffPeaks_plots.pdf:        plots file\n",
    "  sizefactor.tsv:             sizefactor of each library, internally used by DESeq2 for normalizing librarieS\n\n",

	"Function: Use DESeq2 to find differentially accessible regions (DARs) and generate plots.\n",
	"Qirui Zhang (qirui.zhang@med.lu.se)\n",
	"14-05-2020\n\n"
	)
}

args<-commandArgs(TRUE)
if(length(args)!=8){Usage();quit();}


#===============================================================================================
cat("\n=========================================================================================\n")
# load libraries and arguments
time<-format(Sys.time(), format='%H:%M:%S %Y-%m-%d')
cat(time, "Start differential analysis\n\n")
cat("loading libraries and arguments...\n\n")

library(DESeq2)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(gridExtra)
pdf("diffPeaks_plots.pdf", useDingbats=FALSE)
options(scipen=20)

reads.count.file<-args[1]  # mmATAC-seq_readscount.mx
sample.info.file<-args[2]  # sample.info
adjust.pvalue<-as.numeric(args[3])  # 0.05
fold_change<-as.numeric(args[4])  # 2
fi.peaks.file<-args[5]  # fi_final.bed
fn.peaks.file<-args[6]  # fn_final.bed
ng.peaks.file<-args[7]  # ng_final.bed
km.peaks.file<-args[8]  # km_final.bed


# read data
cat("reading files...\n\n")
reads.count<-read.table(reads.count.file, header=T, stringsAsFactors=F)
rownames(reads.count)<-reads.count$peakID
reads.count<-subset(reads.count, select=-c(chr, start, end, peakID, m16.11))

sample.info<-read.table(sample.info.file, header=T, stringsAsFactors=F)
colnames(sample.info)<-c("sample", "mutation", "replicate")
sample.info$mutation<-as.factor(sample.info$mutation)
colnames(reads.count)<-sample.info$sample

fi.peaks<-read.table(fi.peaks.file, header=F, stringsAsFactors=F)
fn.peaks<-read.table(fn.peaks.file, header=F, stringsAsFactors=F)
ng.peaks<-read.table(ng.peaks.file, header=F, stringsAsFactors=F)
km.peaks<-read.table(km.peaks.file, header=F, stringsAsFactors=F)


#===============================================================================================
cat("\n=========================================================================================\n")
# generate metadata and normalize reads count using DESeq2
cat("generating metadata and normalizing reads count using DESeq2...\n\n")
dds.full<-DESeqDataSetFromMatrix(countData=reads.count, colData=sample.info, design= ~ mutation)
dds<-dds.full[rowSums(counts(dds.full))>=10,]
dds<-DESeq(dds)
normalized.counts<-as.data.frame(counts(dds, normalized=TRUE))
normalized.counts$mean<-apply(normalized.counts, 1, mean)
cat("\nnormalized reads count quantile distribution:\n")
quantile(normalized.counts$mean)
write.table(as.data.frame(normalized.counts), "atac-seq_deseq2norm.tsv", row.names=T, col.names=T, quote=F, sep="\t")
sf<-sizeFactors(dds)
write.table(as.data.frame(sf), "sizefactor.tsv", row.names=T, col.names=T, quote=F, sep="\t")


# logarithm transform reads count using DESeq2::vst
cat("\nlogarithm transforming reads count using DESeq2::vst...\n\n")
vst<-varianceStabilizingTransformation(dds.full, blind=FALSE)
write.table(as.data.frame(assay(vst)), "atac-seq_vst.tsv", quote=F, row.names=T, col.names=T, sep="\t")


# plot sample correlation
cat("plotting sample correlation...\n\n")
sampleDist<-dist(t(assay(vst)))
sampleDist.mx<-as.matrix(sampleDist)
rownames(sampleDist.mx)<-NULL
colnames(sampleDist.mx)<-sample.info$Sample
colors<-colorRampPalette(rev(brewer.pal(9,"Blues")))(255)
pheatmap(as.data.frame(sampleDist.mx), clustering_distance_rows=sampleDist, clustering_distance_cols=sampleDist, col=colors, show_rownames=F)


# plot PCA using DESeq2::plotPCA
cat("plotting PCA using DESeq2::plotPCA ...\n\n")
data<-plotPCA(vst, intgroup="mutation", returnData=TRUE)
percentVar<-round(100*attr(data, "percentVar"))

pca.plot1<-ggplot(data, aes(PC1, PC2, color=mutation))+geom_point(size=3)+theme_bw()+theme(panel.grid=element_blank(), panel.border=element_rect(size=1), text=element_text(size=18, family="sans", color="black"), axis.text=element_text(size=15, color="black"), legend.text=element_text(size=15))+xlab(paste0("PC1: ",percentVar[1],"% variance"))+ylab(paste0("PC2: ",percentVar[2],"% variance"))+coord_fixed(ratio=2.5)
pca.plot1
pca.plot1+geom_text(aes(label=rownames(data)),size=3)


# plot PCA using in-house script (recommended!)
cat("plotting PCA using in-house script ...\n")
pca.results<-prcomp(t(as.data.frame(assay(vst))), scale=TRUE)
percentVar<-round(((pca.results$sdev^2)/sum(pca.results$sdev^2))*100, 1)
data<-as.data.frame(pca.results$x)
data$mutation<-sample.info$mutation
pca.plot2<-ggplot(data, aes(PC1, PC2, colour=mutation))+geom_point(size=4, shape=17)+xlab(paste("PC1: ", percentVar[1], "% variance", sep=""))+ylab(paste("PC2: ", percentVar[2], "% variance", sep=""))+theme_bw()+theme(panel.grid=element_blank(), panel.border=element_rect(colour="black",size=1), text=element_text(size=18, family="sans", color="black"), axis.title=element_text(size=18), axis.text=element_text(size=15, color="black"), legend.title=element_text(size=15), legend.text=element_text(size=15))+scale_colour_brewer(palette="Set1")+coord_fixed(ratio=1)
pca.plot2
pca.plot2+geom_text(aes(label=rownames(data)),size=3)


#===============================================================================================
cat("\n=========================================================================================\n")
# run differential comparisons
cat("running differential comparisons...\n\n")

getDiffPeaks<-function(dds.full, mutat, contr, mutat.peaks, contr.peaks){
	mutat=mutat  # "fi"
	contr=contr  # "km"
    mutat.peaks=mutat.peaks  # fi.peaks
    contr.peaks=contr.peaks  # contr.peaks
    all.peaks=unique(rbind(mutat.peaks, contr.peaks))

    # differential comparison
    dds.subset<-dds.full[all.peaks$V4, dds.full$mutation==mutat | dds.full$mutation==contr]
    dds.subset<-dds.subset[rowSums(counts(dds.subset))>=10,]
    dds.subset$mutation<-droplevels(dds.subset$mutation)
    dds.subset$mutation<-relevel(dds.subset$mutation, contr)
    dds.subset<-DESeq(dds.subset)
	res=results(dds.subset, alpha=adjust.pvalue)
    summary(res)

    # extract up- and down-regulated DARs (padj<0.05, |foldchange|>=2)
    dar.up=res[!is.na(res$padj) & res$padj<adjust.pvalue & res$log2FoldChange >= log2(fold_change),]  # up-regulated DARs
    dar.down=res[!is.na(res$padj) & res$padj<adjust.pvalue & res$log2FoldChange <= -log2(fold_change),]  # down-regulated DARs
    dar.all=rbind(dar.up, dar.down)
    dar.up.num=nrow(dar.up)
    dar.down.num=nrow(dar.down)
    dar.all.num=nrow(dar.all)
	cat("up-regulated: ", dar.up.num, "\n")
	cat("down-regulated: ", dar.down.num, "\n")
    cat("total DARs: ", dar.all.num, "\n\n")

    # plotMA plot
    title=paste("MA plot of logFC in ", mutat, "vs", contr, sep="")
    plotMA(res, main={title}, ylim=c(-10,10))

	# output DARs with fc>2
	outname.up=paste("DARs_", mutat, "vs", contr, "_up.tsv", sep="")
	write.table(as.data.frame(dar.up), outname.up, row.names=T, col.names=T, quote=F, sep="\t")
	outname.down=paste("DARs_", mutat, "vs", contr, "_down.tsv", sep="")
	write.table(as.data.frame(dar.down), outname.down, row.names=T, col.names=T, quote=F, sep="\t")

	out=list(res, dar.all)
	return(out)
}

# KM+FLT3-ITD vs KM
cat("--------------------------------------------------\n")
cat("Comparing KM+FLT3-ITD vs KM ...\n")
fivskm.diff.results<-getDiffPeaks(dds.full, "fi", "km", fi.peaks, km.peaks)
res.fivskm<-fivskm.diff.results[[1]]
dar.fivskm<-fivskm.diff.results[[2]]

# KM+FLT3-N676K vs KM
cat("--------------------------------------------------\n")
cat("Comparing KM+FLT3-N676K vs KM ...\n")
fnvskm.diff.results<-getDiffPeaks(dds.full, "fn", "km", fn.peaks, km.peaks)
res.fnvskm<-fnvskm.diff.results[[1]]
dar.fnvskm<-fnvskm.diff.results[[2]]

# KM+NRAS-G12D vs KM
cat("--------------------------------------------------\n")
cat("Comparing KM+NRAS-G12D vs KM ...\n")
ngvskm.diff.results<-getDiffPeaks(dds.full, "ng", "km", ng.peaks, km.peaks)
res.ngvskm<-ngvskm.diff.results[[1]]
dar.ngvskm<-ngvskm.diff.results[[2]]


#===============================================================================================
cat("\n=========================================================================================\n")
# make volcano plots of DARs
cat("making volcano plots of DARs\n\n")

VolcanoPlot<-function(res, compareString){
	volcano<-as.data.frame(res)
	volcano$significant<-as.factor(ifelse(!is.na(volcano$padj) & volcano$padj < adjust.pvalue & abs(volcano$log2FoldChange)>=log2(fold_change), ifelse(volcano$log2FoldChange>=log2(fold_change), "Up", "Down"), "No"))
	volcano$padj<-ifelse(is.na(volcano$padj),1,volcano$padj)
	volcano$group<-rep(0)
	volcano[which(volcano$padj >= adjust.pvalue | (volcano$padj < adjust.pvalue & abs(volcano$log2FoldChange) < log2(fold_change))),8]=rep(1)
	volcano[which(volcano$padj < adjust.pvalue & volcano$log2FoldChange < -log2(fold_change)),8]=rep(2)
	volcano[which(volcano$padj < adjust.pvalue & volcano$log2FoldChange > log2(fold_change)),8]=rep(3)
	comp=compareString
	title=paste("Volcano plot of ", comp, sep="")

	p<-ggplot(volcano,aes(log2FoldChange,-log10(padj)))+geom_point(data=volcano[which(volcano$group==1),], color="gray", alpha=0.75)+geom_point(data=volcano[which(volcano$group==2),], color="royalblue3", alpha=0.75)+geom_point(data=volcano[which(volcano$group==3),], color=brewer.pal(11, "RdYlBu")[2], alpha=0.75)
	p+labs(title={title}, x="log2FoldChange", y="-log10(padj)")+geom_hline(yintercept=-log10(adjust.pvalue), linetype=2, color="gray80")+geom_vline(xintercept=c(-log2(fold_change), log2(fold_change)), linetype=2, color="gray80")+theme_bw()+theme(text=element_text(size=20, family="sans", color="black"), panel.grid=element_blank(), panel.border=element_rect(size=1), axis.ticks=element_line(colour="black"), axis.text=element_text(colour="black"), legend.position = c(.95, .95))
}

# KM+FLT3-ITD vs KM
cat("ploting KM+FLT3-ITD vs KM...", "\n")
VolcanoPlot(res.fivskm, "fivskm")

# KM+FLT3-N676K vs KM
cat("ploting KM+FLT3-N676K vs KM...", "\n")
VolcanoPlot(res.fnvskm, "fnvskm")

# KM+NRAS-G12D vs KM
cat("ploting KM+NRAS-G12D vs KM...", "\n")
VolcanoPlot(res.ngvskm, "ngvskm")


#===============================================================================================
cat("\n=========================================================================================\n")
# make volcano plots in modified version
cat("making volcano plots in modified version\n\n")

VolcanoPlotV2<-function(res, compareString){
	volcano=as.data.frame(res)
	volcano$significant=as.factor(ifelse(!is.na(volcano$padj) & volcano$padj < adjust.pvalue & abs(volcano$log2FoldChange)>=log2(fold_change), "yes", "no"))
	volcano$padj=ifelse(is.na(volcano$padj),1,volcano$padj)

	comp=compareString
	title=paste("Volcano plot of ", comp, sep="")

	p=ggplot(volcano, aes(x=log2FoldChange, y=log2(baseMean)))+geom_point(data=subset(volcano, significant=="no"), colour="gray")+geom_point(data=subset(volcano, significant=="yes"), aes(colour=padj), alpha=0.8)+scale_colour_gradientn("Padj", colours=c(rev(brewer.pal(9, "OrRd"))[c(2,5,6,7,8,9)]), limits=c(0, adjust.pvalue), breaks=c(0.01, 0.03, adjust.pvalue), labels=c("0.01", "0.03", adjust.pvalue))
	p+labs(title={title}, x="Log2foldchange", y="Log2basemean")+geom_vline(xintercept=c(-1,1), linetype=2, color="black")+theme_bw()+theme(text=element_text(family="sans", color="black", size=25), axis.title=element_text(size=30), axis.text=element_text(color="black", size=25), legend.title=element_text(size=25), legend.text=element_text(size=25), axis.line=element_line(color="black", size=1), panel.grid=element_blank(), panel.border=element_blank(), axis.ticks=element_line(colour="black"), legend.position=c(0.9, 0.85))+coord_flip(xlim=c(-5, 5), ylim=c(3,16))
}

# KM+FLT3-ITD vs KM
cat("ploting KM+FLT3-ITD vs KM...", "\n")
VolcanoPlotV2(res.fivskm, "fivskm")

# KM+FLT3-N676K vs KM
cat("ploting KM+FLT3-N676K vs KM...", "\n")
VolcanoPlotV2(res.fnvskm, "fnvskm")

# KM+NRAS-G12D vs KM
cat("ploting KM+NRAS-G12D vs KM...", "\n")
VolcanoPlotV2(res.ngvskm, "ngvskm")


#===============================================================================================
cat("\n=========================================================================================\n")
# plot log2FoldChange density curves
cat("plotting log2FoldChange density curves...\n\n")

DensityPlot<-function(res, compareString){
	density=as.data.frame(res)
	comp=compareString
	title=paste("Density plot of ", comp, sep="")

	ggplot(density, aes(x=log2FoldChange))+geom_density(fill="gray")+labs(title={title}, x="Log2foldchange", y="Density")+geom_vline(xintercept=c(-1,1), linetype=2, color="black")+theme_bw()+theme(title=element_text(size=15), text=element_text(family="sans", color="black", size=25), axis.title=element_text(size=30), axis.text=element_text(color="black", size=25), legend.title=element_text(size=25), legend.text=element_text(size=25), axis.line=element_line(colour="black", size=1), panel.grid=element_blank(), panel.border=element_blank(), aspect.ratio=5)+scale_y_continuous(breaks=c(0, 0.5, 1), labels=c("0", "0.5", "1.0"))+coord_flip(xlim=c(-5, 5), ylim=c(0, 1.4))
}

# KM+FLT3-ITD vs KM
cat("ploting KM+FLT3-ITD vs KM...", "\n")
DensityPlot(res.fivskm, "fivskm")

# KM+FLT3-N676K vs KM
cat("ploting KM+FLT3-N676K vs KM...", "\n")
DensityPlot(res.fnvskm, "fnvskm")

# KM+NRAS-G12D vs KM
cat("ploting KM+NRAS-G12D vs KM...", "\n")
DensityPlot(res.ngvskm, "ngvskm")


#===============================================================================================
cat("\n=========================================================================================\n")
# plot heatmap of all DARs from all above comparisons
cat("plotting heatmap using unique DARs from all three comparisons...\n\n")

dar.id<-unique(c(rownames(dar.fivskm), rownames(dar.fnvskm), rownames(dar.ngvskm)))  # merge all unique DARs from 3 comparisons
dar.num<-length(dar.id)
cat("total number of unique DARs: ", dar.num, "\n")
vst.dar<-vst[dar.id,]
vst.dar.df<-as.data.frame(assay(vst.dar))
#nb.clust<-NbClust(vst.dar.df, distance="euclidean", min.nc=2, max.nc=15, method="kmeans", index="alllong", alphaBeale=0.1)
  # Give recommended cluster numbers

anno.label<-c(rep("KM+FLT3-ITD",4), rep("KM+FLT3-N676K",4), rep("KM+NRAS-G12D",3), rep("KM",4))
names(anno.label)<-sample.info$sample
anno.label<-as.data.frame(anno.label)
names(anno.label)<-"Mutation"
anno.color<-list(Mutation=c("KM+FLT3-ITD"="red", "KM+FLT3-N676K"="orange", "KM+NRAS-G12D"="green", "KM"="blue"))
title=paste("Heatmap of all ", dar.num, " DARs", sep="")

# not cluster columns
pheatmap(vst.dar.df, main={title}, scale="row", color=colorRampPalette(c(rev(brewer.pal(9,"Blues")[c(3:9)]), brewer.pal(9, "OrRd")[c(3:9)]))(100), cluster_rows=TRUE, cluster_cols=FALSE, show_rownames=FALSE, annotation_col=anno.label, annotation_colors=anno.color, annotation_names_col=F, border_color=NA)

# cluster columns
pheatmap(vst.dar.df, main={title}, scale="row", color=colorRampPalette(c(rev(brewer.pal(9,"Blues")[c(3:9)]), brewer.pal(9, "OrRd")[c(3:9)]))(100), cluster_rows=TRUE, cluster_cols=TRUE, show_rownames=FALSE, annotation_col=anno.label, annotation_colors=anno.color, annotation_names_col=F, border_color=NA)

dev.off()


#===============================================================================================
time<-format(Sys.time(), format='%H:%M:%S %Y-%m-%d')
cat(time, "Done with analysis!", "\n\n")


