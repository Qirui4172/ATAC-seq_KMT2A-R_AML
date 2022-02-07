#!/usr/bin/env Rscript



Usage<-function(){
  cat("Usage: Rscript peakgeneCorrPlot.r [peakgene_up] [peakgene_down] [atacseq_mx] [rnaseq_mx] [annotation] [treatment] [control] [out_suffix]","\n\n",

  "Parameters:\n",
  "[peakgene_up]        significantly correlated up-regulated peak-gene list, i.e. \"sigcorr_fivskm_up.tsv\"\n",
  "[peakgene_down]      significantly correlated down-regulated peak-gene list, i.e. \"sigcorr_fivskm_down.tsv\"\n",
  "[atacseq_mx]         atac-seq_deseq2norm_forCorrelation.mx\n",
  "[rnaseq_mx]          rna-seq_log2fpkm_forCorrelation.mx, RNA-seq read count matrix normalized by log2(fpkm+0.1), the extra columns after read counts, i.e. chr, startcondon will not be used.\n",
  "[annotation]         annotation file of consensus peaks in shortmost version, anno_consensus_shortmost.txt\n",
  "[treatment]          either of fi, fn, ng\n",
  "[control]            km\n",
  "[out_suffix]         output suffix, i.e. fivskm, the output file corrPlots_fivskm.pdf will be generated.\n\n",

  "Function:\n",
  "Plot foldchange-distance distribution; gene log2FPKM distribution between km and either of fi, fn, ng; and heatmaps of DARs and significantly correlated genes.\n\n",

  "Example:\n",
  "Rscript peakgeneCorrPlot.r sigcorr_fivskm_up.tsv sigcorr_fivskm_down.tsv atac-seq_deseq2norm_forCorrelation.mx rna-seq_log2fpkm_forCorrelation.mx anno_consensus_shortmost.txt fi km fivskm\n",
  "# the file \"corrPlots_fivskm.pdf\" will be generted\n\n",

  "Qirui Zhang (qirui.zhang@med.lu.se)\n",
  "Date: 15-10-2020\n\n"
  )
}

args<-commandArgs(TRUE)
if (length(args)!=8){Usage();quit();}


#===================================================================================
# Load arguments & libraries
cat("loading arguments & libraries...\n\n")

up.file<-args[1]  # sigcorr_fivskm_up.tsv
down.file<-args[2]  # sigcorr_fivskm_down.tsv
atac.file<-args[3]  # atac-seq_deseq2norm_forCorrelation.mx
rna.file<-args[4]  # rna-seq_log2fpkm_forCorrelation.mx
anno.file<-args[5]  # anno_consensus_shortmost.txt
treat<-args[6]  # fi
contr<-args[7]  # km
out.suffix<-args[8]  # fivskm

library(ggplot2)
library(pheatmap)
library(RColorBrewer)
out<-paste("corrPlots_", out.suffix, ".pdf", sep="")  # corrPlots_fivskm.pdf
pdf(out, useDingbats=FALSE)


#===================================================================================
# Read files
cat("reading files...\n\n")

# Read peak-gene pairs
peakgene.up<-read.table(up.file, header=T, stringsAsFactors=F)
peakgene.up$Regulation<-rep("Up")
peakgene.down<-read.table(down.file, header=T, stringsAsFactors=F)
peakgene.down$Regulation<-rep("Down")
peakgene.all<-rbind(peakgene.up, peakgene.down)
peakgene.all<-peakgene.all[,c(1,2,5,7)]
colnames(peakgene.all)<-c("peakID", "geneID", "distance", "regulation")
peakgene.all$distance<-abs(peakgene.all$distance)/1000

# read annotation file
anno<-read.table(anno.file, header=F, stringsAsFactors=F, sep="\t", quote="")
anno<-anno[,1:2]
colnames(anno)<-c("peakID", "geneID")
rownames(anno)<-anno$peakID

# Read ATAC-seq matrix
atac<-read.table(atac.file, header=T, stringsAsFactors=F)
atac[,1:15]<-log2(atac[,1:15]+0.5) # log-transform to keep consistent with RNA-seq data
colnames(atac)[17:18]<-c("chr", "peakcenter")
atac$peakcenter<-atac$peakcenter/1000
atac$geneID<-anno$geneID

# Read RNA-seq matrix
rna<-read.table(rna.file, header=T, stringsAsFactors=F)
rownames(rna)<-rna$Gene
rna<-subset(rna, select=-Gene)
colnames(rna)[17:18]<-c("chr", "startcodon")
rna$startcodon<-rna$startcodon/1000

# Create mutation-column table
atac.column<-data.frame(start=c(1, 5, 9, 12), end=c(4, 8, 11, 15), fullname=c("FLT3-ITD", "FLT3-N676K", "NRAS-G12D", "KMT2A"))
rownames(atac.column)<-c("fi", "fn", "ng", "km")
atac_treat.start<-atac.column[treat,1]
atac_treat.end<-atac.column[treat,2]
atac_contr.start<-atac.column[contr,1]
atac_contr.end<-atac.column[contr,2]

rna.column<-data.frame(start=c(1, 5, 9, 13), end=c(4, 8, 12, 16), fullname=c("FLT3-ITD", "FLT3-N676K", "NRAS-G12D", "KMT2A"))
rownames(rna.column)<-c("fi", "fn", "ng", "km")
rna_treat.start<-rna.column[treat,1]
rna_treat.end<-rna.column[treat,2]
rna_contr.start<-rna.column[contr,1]
rna_contr.end<-rna.column[contr,2]


#===================================================================================
# Process files
cat("processing files...\n\n")

# Compute log2fc of gene expression
rna$treat_mean<-apply(rna[,rna_treat.start:rna_treat.end], 1, mean)
rna$contr_mean<-apply(rna[,rna_contr.start:rna_contr.end], 1, mean)
rna$lfc<-rna$treat_mean-rna$contr_mean

# extract log2foldchange from rna-seq and add to peakgene.all
peakgene.all$lfc<-rna[peakgene.all$geneID, "lfc"]

# sort geneID & distance and filter
computeMedianDist<-function(peakgene.list){  # compute median distance of the same genes
    peakgene.list=peakgene.list  # the input dataframe should have columns of "geneID", "distance", "regulation", "lfc"
    peakgene.list=peakgene.list[order(peakgene.list$geneID),]  # first sort by gene names

    peakgene.unique=unique(peakgene.list[,c("geneID","regulation","lfc")])
    genes.list=peakgene.unique$geneID
    regulation.list=peakgene.unique$regulation
    lfc.list=peakgene.unique$lfc

    dist.list=peakgene.list$distance
    names(dist.list)=peakgene.list$geneID
    dist.median=sapply(genes.list, function(gene){dist.median=median(dist.list[which(names(dist.list)==gene)])})  # one gene may have multiple distance and take the median distance

    dist_fc=data.frame(geneID=genes.list, regulation=regulation.list, distance=dist.median, lfc=lfc.list)
    return(dist_fc)
}

computeMinimumDist<-function(peakgene.list){  # compute minimum distance of the same genes
    peakgene.list=peakgene.list  # the input dataframe should have columns of "geneID", "distance", "regulation", "lfc"
    peakgene.list=peakgene.list[order(peakgene.list$geneID),]  # first sort by gene names

    peakgene.unique=unique(peakgene.list[,c("geneID","regulation","lfc")])
    genes.list=peakgene.unique$geneID
    regulation.list=peakgene.unique$regulation
    lfc.list=peakgene.unique$lfc

    dist.list=peakgene.list$distance
    names(dist.list)=peakgene.list$geneID
    dist.min=sapply(genes.list, function(gene){dist.min=min(dist.list[which(names(dist.list)==gene)])})  # one gene may have multiple distance and take the median distance

    dist_fc=data.frame(geneID=genes.list, regulation=regulation.list, distance=dist.min, lfc=lfc.list)
    return(dist_fc)
}

computeMinimumDist2<-function(peakgene.list){
    peakgene.list=peakgene.list
    peakgene.list=peakgene.list[order(peakgene.list$geneID, peakgene.list$distance),]
    peakgene.list=peakgene.list[!duplicated(peakgene.list$geneID),]
    dist_fc=peakgene.list[,c("geneID","regulation","distance","lfc")]
    return(dist_fc)
}

dist_fc.sig<-computeMinimumDist2(peakgene.all)

  # one gene might be correlated by multiple peaks, in order to only keep the gene that has the shortest peak-gene distance, it needs to first sort the distance of the same gene in ascending order, then remove other same peak-gene pairs.


#====================================================================================================
# extract non-significantly correlated peak-gene pairs and compute distance
atac.bg<-atac[!rownames(atac) %in% unique(peakgene.all$peakID),]  # for those non-significantly correlated peaks in the background
#genes.bg<-atac.bg$geneID

peakgene.bg<-data.frame(geneID=atac.bg$geneID, peakcenter=atac.bg$peakcenter, startcodon=rna[atac.bg$geneID,"startcodon"], regulation=rep("No"), lfc=rna[atac.bg$geneID,"lfc"], stringsAsFactors=F)
peakgene.bg<-na.omit(peakgene.bg)
peakgene.bg$distance=abs(peakgene.bg$peakcenter-peakgene.bg$startcodon)
peakgene.bg<-peakgene.bg[which(peakgene.bg$distance<=500),]  # remove genes with distance > 500kb
peakgene.bg<-peakgene.bg[which(!peakgene.bg$geneID %in% dist_fc.sig$geneID),]  # remove genes that already exist in dist_fc.sig
dist_fc.bg<-computeMinimumDist2(peakgene.bg)  # compute minimum distance

# combine significantly and non-significantly correlated peak-gene pairs
dist_fc.all=rbind(dist_fc.sig, dist_fc.bg)


#===================================================================================
# plot peak-gene distance & gene expression fold change
cat("ploting distance-foldchange distribution...\n\n")

dist_fc.all$regulation<-factor(dist_fc.all$regulation, levels=c("Up","Down","No"))  # set regulation as factor and set the order of regulations so that the colors can be given in this order 
ggplot(dist_fc.all, aes(distance, lfc))+geom_point(data=dist_fc.all[which(dist_fc.all$regulation=="No"),], color="gray")+geom_point(data=dist_fc.all[which(dist_fc.all$regulation!="No"),], aes(color=regulation))+scale_color_manual("Regulation", breaks=c("Up","Down"), values=c("#F8766D","#00BFC4"))+labs(title="Distance-log2foldchange distribution", x="Distance (kb)", y="Log2foldchange of gene expression")+ylim(-4.5, 7.5)+theme_bw()+theme(panel.grid=element_blank(), panel.border=element_blank(), axis.text=element_text(color="black", size=20), axis.title=element_text(size=25), axis.line=element_line(color="black", size=0.5), legend.title=element_text(size=20), legend.text=element_text(size=20))


#===================================================================================
# plot gene log2fpkm distribution between control and treatment
cat("ploting gene log2fpkm distribution between control and treatment...\n\n")

# add "Up/Down" information to rna matrix
gene_regul.list<-unique(peakgene.all[,c("geneID","regulation")])
rownames(gene_regul.list)<-gene_regul.list$geneID
rna$regulation<-gene_regul.list[rownames(rna),"regulation"]
rna[which(is.na(rna$regulation)),"regulation"]<-"No"

y.title<-paste("KMT3 + ", atac.column[treat,"fullname"], " log2FPKM", sep="")
rna$regulation<-factor(rna$regulation, levels=c("Up","Down","No"))
ggplot(rna, aes(contr_mean, treat_mean))+geom_point(data=rna[which(rna$regulation=="No"),], color="gray")+geom_point(data=rna[which(rna$regulation!="No"),], aes(color=regulation))+scale_color_manual("Regulation", breaks=c("Up","Down"), values=c("#F8766D","#00BFC4"))+labs(title="Gene log2FPKM distribution between control and treatment", x="KMT3 log2FPKM", y={y.title})+theme_bw()+theme(panel.grid=element_blank(), panel.border=element_blank(), axis.text=element_text(color="black", size=20), axis.title=element_text(size=25), axis.line=element_line(color="black", size=0.5), legend.title=element_text(size=20), legend.text=element_text(size=20))  # "#CCCCCC" is gray, "#F8766D" is red, "#00BFC4" is green.


#===================================================================================
# Plot heatmap of differential peaks
cat("ploting peak-gene heatmaps...\n\n")

atac.sig<-atac[unique(peakgene.all$peakID), c(atac_treat.start:atac_treat.end, atac_contr.start:atac_contr.end)]
p<-pheatmap(as.matrix(atac.sig), scale="row", color=colorRampPalette(brewer.pal(9,"YlGnBu"))(100), cluster_rows=T, cluster_cols=F, show_rownames=F, show_colnames=T, border_color=NA, main="DARs", fontsize=14, cellwidth=20)
p

# Keep peak clustering order and reorder target genes
clust.order<-p$tree_row$order
ordered.peaks<-rownames(atac.sig[clust.order,])

# keep nearest genes and plot
peak_gene.uniq<-peakgene.all[order(peakgene.all$peakID, peakgene.all$distance),]
peak_gene.uniq<-peak_gene.uniq[!duplicated(peak_gene.uniq$peakID),]  # one peak may correlate with multiple genes and only keep the nearest one
order.numbers<-sapply(peak_gene.uniq$peakID, function(peak){order.num=which(ordered.peaks==peak)})
peak_gene.uniq$order<-order.numbers
peak_gene.uniq<-peak_gene.uniq[order(peak_gene.uniq$order),]  # sort peaks to the same order as ordered.peaks
rna.sig<-rna[peak_gene.uniq$geneID, c(rna_treat.start:rna_treat.end, rna_contr.start:rna_contr.end)]
pheatmap(as.matrix(rna.sig), scale="row", color=colorRampPalette(brewer.pal(9,"YlOrBr"))(100), cluster_rows=F, cluster_cols=F, show_rownames=F, show_colnames=T, border_color=NA, main="Target genes (nearest one if multiple)", fontsize=14, cellwidth=20)


dev.off()

#===================================================================================
cat("Done with analysis!\n\n")


