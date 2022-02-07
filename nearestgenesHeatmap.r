#!/usr/bin/env Rscript



Usage<-function(){
  cat("Usage: Rscript nearestgenesHeatmap.r [contr_sort_mx] [mutat_sort_mx] [DARs_list] [rnaseq_mx] [contr] [mutat]","\n\n",

  "Parameters:\n",
  "[contr_sort_mx]      control sorted peaks matrix file, \"km_matrix_sort\"\n",
  "[mutat_sort_mx]      mutant sorted peaks matrix file, i.e. \"fi_matrix_sort\"\n",
  "[DARs_list]          DARs list in bed format, i.e. \"DARs_fivskm.bed\"\n",
  "[rnaseq_mx]          RNA-seq read count matrix prepared in peak-gene correlation analysis, \"rna-seq_log2fpkm_forCorrelation.mx\"\n",
  "[contr]              control, should be km\n",
  "[mutat]              mutat, either of fi, fn, ng\n\n",

  "Function:\n",
  "Extract genes nearest to DARs, and plot heatmap by increasing ratio of mean(mutant) to mean (km).\n\n",

  "Example:\n",
  "Rscript nearestgenesHeatmap.r km_matrix_sort fi_matrix_sort DARs_fivskm.bed rna-seq_log2fpkm_forCorrelation.mx km fi\n",
  "# the file \"nearestgenesHeatmap_fivskm.pdf\" will be generted\n\n",

  "Qirui Zhang (qirui.zhang@med.lu.se)\n",
  "Date: 19-09-2021\n\n"
  )
}

args<-commandArgs(TRUE)
if (length(args)!=6){Usage();quit();}


#=========================================================================================================
# Load libraries & arguments
cat("loading libraries & arguments...\n\n")

library(RColorBrewer)
library(pheatmap)

contr_peak.sort.file<-args[1]  # "km_matrix_sort"
mutat_peak.sort.file<-args[2]  # "fi_matrix_sort"
dars.list.file<-args[3]  # "DARs_fivskm.bed"
rna.mx.file<-args[4]  # "rna-seq_log2fpkm_forCorrelation.mx"
contr<-args[5]  # "km"
mutat<-args[6]  # either of "fi", "fn", "ng"
out.file<-paste("nearestgenesHeatmap_", mutat, "vs", contr, ".pdf", sep="")  # nearestgenesHeatmap_fivskm.pdf
pdf(out.file)


#=========================================================================================================
# Read and process files
cat("reading and processing files...\n\n")

# read sorted peaks matrix files with binned values in each group (i.e. "km_matrix_sort", generated from computeMatrix of deepTools)
contr_peak.sort<-read.table(contr_peak.sort.file, skip=1, header=F, sep="\t", stringsAsFactors=F)
mutat_peak.sort<-read.table(mutat_peak.sort.file, skip=1, header=F, sep="\t", stringsAsFactors=F)

# read DARs full list file ("DARs_fivskm.bed")
dars.list<-read.table(dars.list.file, header=F, stringsAsFactors=F)
colnames(dars.list)<-c("chr", "start", "end", "peakID")
dars.list$peakcenter<-apply(dars.list[,c("start","end")], 1, median)
rownames(dars.list)<-paste(dars.list$chr, ":", dars.list$start, "-", dars.list$end, sep="")
dars.list<-dars.list[contr_peak.sort$V4,]  # sort DARs to the same order as km_matrix_sort or fi_matrix_sort

# read RNA-seq matrix file ("rna-seq_log2fpkm_forCorrelation.mx")
rna.mx<-read.table(rna.mx.file, header=T, stringsAsFactors=F)
rownames(rna.mx)<-rna.mx$Gene
rna.mx<-subset(rna.mx, select=-Gene)

# prepare RNA-seq start-end column dataframe
rna.column<-data.frame(start=c(1, 5, 9, 13), end=c(4, 8, 12, 16), fullname=c("FLT3-ITD", "FLT3-N676K", "NRAS-G12D", "KMT2A"))
rownames(rna.column)<-c("fi", "fn", "ng", "km")


#=========================================================================================================
# Plot sorted peaks heatmap (abandoned! use deeptools::plotHeatmap instead!)
#cat("plotting sorted peaks heatmap...\n\n")

#color.stock<-c(brewer.pal(9, "Reds")[1:3], brewer.pal(9, "Reds")[7:9])
#colors=c(colorRampPalette(color.stock[1:3])(200), colorRampPalette(color.stock[3:4])(50), colorRampPalette(color.stock[4:5])(50), colorRampPalette(color.stock[5:6])(150))
    # there're many values between brewer.pal(9, "Reds")[3] and brewer.pal(9, "Reds")[7] which are light red and dark red respectively, to avoid showing too many mid-red colors in the background only 50 pieces were set between color.stock[3] and [4] (which are in fact brewer.pal(9, "Reds")[3] and [7]).

#bks<-c(seq(0,600,length=200), seq(600.00001,800,length=50), seq(800.00001,17000,length=201))
    # the signal values should correspond to colors, and the pieces of values should also same as the pieces of colors, that is 200 pieces (0-600), 50 pieces (600-800), 201 pieces (800-17000).
    # to make the values corresponding to colors perfectly and gives ideal visible results, it first should be clear with the values distribution by some methods, i.e. summary(km), 1300 is the max(Q3) of km/fi/fn/ng, 17000 is the lowest max values of km/fi/fn/ng, so here set values range 0-17000 corresponding to colors, values above 17000 will have same color as 17000.

#peakHeatmapSortbyRatio<-function(data, group){
#    data=data
#    title=group
#    data=data[,7:46]
#    pheatmap(data, main={title}, scale="none", color=colors, cluster_rows=FALSE, cluster_cols=FALSE, show_rownames=FALSE, show_colnames=FALSE, border_color=NA, cellwidth=2, fontsize=20, legend_breaks=seq(0,15000,5000), legend_labels=c("0","5000","10000","15000"), breaks=bks)
#}

#peakHeatmapSortbyRatio(km_peak.sort, "km")
#peakHeatmapSortbyRatio(fi_peak.sort, "fi")
#peakHeatmapSortbyRatio(fn_peak.sort, "fn")
#peakHeatmapSortbyRatio(ng_peak.sort, "ng")


#=========================================================================================================
# extract genes that nearest to DARs
cat("extracting nearest genes...\n\n")

nearestgenes<-data.frame()
for(i in 1:nrow(dars.list)){chr_id=dars.list[i,"chr"]; peak_center=dars.list[i,"peakcenter"]; rna.tmp=rna.mx[which(rna.mx$chr==chr_id),]; rna.tmp$dist=abs(rna.tmp$startcodon-peak_center); rna.tmp=rna.tmp[order(rna.tmp$dist),]; nearestgenes=rbind(nearestgenes, rna.tmp[1,]);}

#FindNearestGenes<-function(peak){
#    chr.id=dars.list[which(dars.list$peakID==peak),"chr"]
#    peakcenter=dars.list[which(dars.list$peakID==peak),"peakcenter"]
#    rna.tmp=rna.mx[which(rna.mx$chr==chr_id),]
#    rna.tmp$dist=abs(rna.tmp$startcodon-peakcenter)
#    rna.tmp$gene=rownames(rna.tmp)
#    rna.tmp=rna.tmp[order(rna.tmp$dist),]
#    return(rna.tmp[1,])
#}
#peak.list<-dars.list$peakID
#nearestgenes<-sapply(peak.list, FindNearestGenes)  # extract nearest genes
#nearestgenes<-as.data.frame(t(nearestgenes))  # switch columns/rows and change matrix to dataframe
#nearestgenes<-as.data.frame(sapply(nearestgenes, as.numeric))  # modify data to numeric format


# compute mean values of gene expression levels in each group
for(group in c(contr, mutat)){start=rna.column[group,"start"]; end=rna.column[group,"end"]; nearestgenes[,group]=apply(nearestgenes[,start:end], 1, mean)}
newcol.name<-paste(mutat, "2", contr, sep="")
nearestgenes[,newcol.name]<-nearestgenes[,mutat]-nearestgenes[,contr]
    # because there're only two columns, plotting two columns with scale="row" will gives only two colors "red" or "green" but not "black", so it's not a good idea to plot two columns. the good way is to compute the mutat/contr ratio then plot ratio.


#=========================================================================================================
# plot
cat("ploting...\n\n")

#bks<-c(seq(-1, 0, length=100), seq(0.00001, 1, length=100))
bk1<-quantile(nearestgenes[,newcol.name], probs=0.2)
bk2<-quantile(nearestgenes[,newcol.name], probs=0.5)
bk3<-quantile(nearestgenes[,newcol.name], probs=0.5)+0.0001
bk4<-quantile(nearestgenes[,newcol.name], probs=0.8)
bks<-c(seq(bk1, bk2, length=100), seq(bk3, bk4, length=100))  # set min->Q1 as green; Q1->median as green->black; median->Q3 as black->red; Q3->max as red.

color.stock<-c("green", "black", "red")
colors<-c(colorRampPalette(color.stock[1:2])(100), colorRampPalette(color.stock[2:3])(100))
title<-paste("Nearest Genes ratio of ", mutat, "/", contr, sep="")
pheatmap(as.matrix(nearestgenes[,newcol.name]), main={title}, scale="none", breaks=bks, color=colors, cluster_rows=F, cluster_cols=F, show_rownames=F, show_colnames=T, angle_col=45, border_color=NA, fontsize=12, cellwidth=100)

dev.off()


#=========================================================================================================
cat("Done with analysis!\n\n")

