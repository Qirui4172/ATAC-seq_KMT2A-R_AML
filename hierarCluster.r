#!/usr/bin/env Rscript

#===============================================================================================
Usage<-function(){
    cat("Usage: Rscript hierarCluster.r [atac_mx] [promoter_peak] [distal_peak] [rna_mx]\n\n",

    "Parameters:\n",
    "[atac_mx]          ATAC-seq matrix, mmATAC-seq_quantileNorm.mx\n",
    "[promoter_peak]    promoter peaks in bed format, Promoter-TSS_proximal.bed\n",
    "[distal_peak]      distal peaks in bed format, Distal.bed\n",
    "[rna_mx]           RNA-seq matrix, rna-seq_log2fpkm_forQirui.mx\n\n",

    "Example:\n",
    "Rscript hierarCluster.r mmATAC-seq_quantileNorm.mx Promoter-TSS_proximal.bed Distal.bed rna-seq_log2fpkm_forQirui.mx\n",
      "# file \"clustering_plots.pdf\" will be generated\n\n",

    "Function: hierarchical clustering of samples using quantile-normalized ATAC-seq reads count matrix and log2FPKM-normalized RNA-seq matrix\n",
    "Qirui Zhang (qirui.zhang@med.lu.se)\n",
    "28-09-2021\n\n"
    )
}

args<-commandArgs(TRUE)
if(length(args)!=4){Usage();quit();}


#================================================================================
# load libraries and arguments
cat("loading libraries and arguments ...\n\n")

library(gplots)
library(RColorBrewer)

atac.file<-args[1]  # mmATAC-seq_quantileNorm.mx
promoter.file<-args[2]  # Promoter-TSS_proximal.bed
distal.file<-args[3]  # Distal.bed
rna.file<-args[4]  # rna-seq_log2fpkm_forQirui.mx
output<-"clustering_plots.pdf"
pdf(output)


#==============================================================================
# read files
cat("reading files...\n\n")

# read atac-seq reads count matrix and do log2 transformation
atac.all<-read.table(atac.file, header=T, stringsAsFactors=F)
atac.all<-log2(atac.all+0.5)
atac.cor<-cor(atac.all)

# read promoter-tss proximal peaks
promoter<-read.table(promoter.file, header=F, stringsAsFactors=F)
colnames(promoter)<-c("chr", "start", "end", "peakID")
atac.promoter<-atac.all[promoter$peakID,]
promoter.cor<-cor(atac.promoter)

# read distal peaks
distal<-read.table(distal.file, header=F, stringsAsFactors=F)
colnames(distal)<-c("chr", "start", "end", "peakID")
atac.distal<-atac.all[distal$peakID,]
distal.cor<-cor(atac.distal)


# read rna-seq reads count matrix
rna.mx<-read.table(rna.file, header=T, stringsAsFactors=F)
rownames(rna.mx)<-rna.mx$Gene
rna.mx<-subset(rna.mx, select=-Gene)
rna.cor<-cor(rna.mx)


#==============================================================================
# clustering
clusterPlot<-function(corr.matrix, clust.method, data.type, color.set){
    corr.matrix=corr.matrix
    clust.method=clust.method
    data.type=data.type
    color.set=color.set
    title=paste(data.type, ": ", clust.method, sep="")

    h.dist<-dist(corr.matrix, method="euclidean")  # compute sample correlation distance using "euclidean" method
    h.clust<-hclust(h.dist, method=clust.method)  # compute sample clustering using given method
    plot(h.clust, main={title}, hang=-1, cex=0.8)
    heatmap.2(as.matrix(corr.matrix), scale="none", trace="none", dendrogram="both", Rowv=as.dendrogram(h.clust), Colv=as.dendrogram(h.clust), col=rev(colorRampPalette(brewer.pal(9, color.set))(100)), density.info="none", key.xlab=NA, keysize=1.5, lmat=rbind(c(4,4,5), c(5,5,3), c(5,2,1), c(5,5,5)), lhei=c(1.5, 1.8, 9, 1), lwid=c(2.5, 1.8, 9))
}

# plot
cat("plotting hierarchical clustering...\n\n")

methods<-c("ward.D", "ward.D2", "complete", "average", "single", "mcquitty", "median", "centroid")
  # methods: "ward.D", "ward.D2", "complete", "average", "single", "mcquitty", "median", "centroid".

for(i in methods){clusterPlot(atac.cor, i, "all peaks", "YlGnBu")}
for(i in methods){clusterPlot(promoter.cor, i, "promoter-tss proximal peaks", "YlGnBu")}
for(i in methods){clusterPlot(distal.cor, i, "distal peaks", "YlGnBu")}
for(i in methods){clusterPlot(rna.cor, i, "rna-seq", "RdYlBu")}

dev.off()


#==============================================================================
cat("Done with analysis!\n\n")


