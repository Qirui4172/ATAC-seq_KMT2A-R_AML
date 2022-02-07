#!/usr/bin/env Rscript

args<-commandArgs(TRUE)
if (length(args)!=6){
        cat("Usage: Rscript correlateHeatmap.r [peakgene] [atacseqMx] [rnaseqMx] [treatment] [control] [outprefix]","\n\n",

        "Parameters:\n",
        "[peakgene]    Peak-gene list, having two columns, e.g. peak_2976 Fbxo36\n",
        "[atacseqMx]   ATAC-seq DESeq2-normlaized read counts matrix but not log2-transformed. There may be extra columns after read counts, i.e. mean, chrID, startSite, but they will be excluded in the script\n",
        "[rnaseqMx]    rna-seq read count metrics normalized by log2(fpkm+0.1). There could be some columns after read counts, i.e. chrID, startSite, but they will not be used in the script\n",
        "[treatment]   Treatment short name, either of fi, fn, ng, km\n",
        "[control]     Control short name, either of fi, fn, ng, km\n",
        "[outprefix]   Output name prefix, e.g. fivskm\n\n",

        "Example:\n",
        "Rscript correlateHeatmap.r fivskm_uniquepair.list atac-seq_deseq2norm_forCorrelation.mx rna-seq_log2fpkm_forCorrelation.mx fi km fivskm\n",
        "# the file \"fivskm_cor_heatmap.pdf\" will be generted\n\n",

        "Function:\n",
        "Plot heatmaps for differential peaks (using ATAC-seq data) and significantly correlated target genes (using RNA-seq data), the target genes are kept in the same clustering order as their regulating peaks.\n",

        "Qirui Zhang (qirui.zhang@med.lu.se)\n",
        "Date: 15-10-2020\n\n"
        )
        quit()
}

#===================================================================================
library(pheatmap)
library(RColorBrewer)
out=paste(args[6], "_cor_heatmap.pdf", sep="")
pdf(out)


#===================================================================================
# Read peak-gene list
peakgene.list<-read.table(args[1], header=F, stringsAsFactors=F)

# Create mutation-column table
mutation.column<-data.frame(start=c(1, 5, 9, 13), end=c(4, 8, 12, 16))
rownames(mutation.column)<-c("fi", "fn", "ng", "km")
treat<-args[4]
start1<-mutation.column[treat,1]
end1<-mutation.column[treat,2]
contr<-args[5]
start2<-mutation.column[contr,1]
end2<-mutation.column[contr,2]

# Read ATAC-seq matrix
atac<-read.table(args[2], header=T, stringsAsFactors=F)
atac<-atac[unique(peakgene.list$V1), c(start1:end1, start2:end2)]
atac<-log2(atac+0.5) # log-transform to keep consistent with RNA-seq data

# Read RNA-seq matrix
rna<-read.table(args[3], header=T, stringsAsFactors=F)
rownames(rna)<-rna$Gene
rna<-subset(rna,select=-Gene)
rna<-rna[unique(peakgene.list$V2),c(start1:end1, start2:end2)]


#===================================================================================
# Plot heatmap of differential peaks
p<-pheatmap(as.matrix(atac), scale="row", color=colorRampPalette(brewer.pal(9,"YlGnBu"))(100), cluster_rows=T, cluster_cols=F, show_rownames=F, show_colnames=T, border_color=NA, main="Differential accessible peaks", fontsize=14, cellwidth=20)
p

# Keep peak clustering order and reorder target genes
clust.order<-p$tree_row$order
ordered.peak<-rownames(atac[clust.order,])
rna.orderd<-data.frame()

for(i in ordered.peak){
	targets<-peakgene.list[which(peakgene.list$V1==i),]$V2
	for(j in targets){rna.orderd<-rbind(rna.orderd, rna[j,])}
}

# Plot heatmap of target genes
pheatmap(as.matrix(rna.orderd), scale="row", color=colorRampPalette(brewer.pal(9,"YlOrBr"))(100), cluster_rows=F, cluster_cols=F, show_rownames=F, show_colnames=T, border_color=NA, main="Correlated target genes", fontsize=14, cellwidth=20)

dev.off()


