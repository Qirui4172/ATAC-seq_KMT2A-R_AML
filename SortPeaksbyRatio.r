#!/usr/bin/env Rscript


Usage<-function(){
  cat("Usage: Rscript SortPeaksbyRatio.r [km_mx] [mut_mx] [mut_name]\n\n",

  "Parameters:\n",
  "[km_mx]       km matrix file, \"km_matrix\"\n",
  "[mut_mx]      activating mutation matrix file, i.e. \"fi_matrix\"\n",
  "[mut_name]    activating mutation name, either of \"fi\", \"fn\", or \"ng\".\n\n",

  "Function:\n",
  "Sort matrix by mean(km_matrix)/mean(mut_matrix) ratio, km_matrix is \"km_matrix\", mut_matrix is either of \"fi_matrix\", \"fn_matrix\", \"ng_matrix\".\n\n",

  "Example:\n",
  "Rscript SortPeaksbyRatio.r km_matrix fi_matrix fi\n",
  "\t# the file \"sortDARs_fi2km_ratio.bed\" will be generted\n\n",

  "Qirui Zhang (qirui.zhang@med.lu.se)\n",
  "Date: 03-10-2021\n\n"
  )
}

args<-commandArgs(TRUE)
if (length(args)!=3){Usage();quit();}


#=========================================================================================================
km_matrix.file<-args[1]  # "km_matrix"
mut_matrix.file<-args[2]  # "fi_matrix"
mut<-as.character(args[3])  # "fi"

# read files and compute mean values
km_matrix<-read.table(km_matrix.file, skip=1, header=F, stringsAsFactors=F, sep="\t")
km_matrix$mean<-apply(km_matrix[,7:86], 1, mean)
km_matrix<-km_matrix[order(km_matrix$V1, km_matrix$V2),]

mut_matrix<-read.table(mut_matrix.file, skip=1, header=F, stringsAsFactors=F, sep="\t")
mut_matrix$mean<-apply(mut_matrix[,7:86], 1, mean)
mut_matrix<-mut_matrix[order(mut_matrix$V1, mut_matrix$V2),]

# compute ratio
mut2km<-data.frame(chr=km_matrix$V1, start=km_matrix$V2, end=km_matrix$V3, km_mean=km_matrix$mean, mut_mean=mut_matrix$mean)
mut2km$mut2km_ratio<-mut2km$mut_mean/mut2km$km_mean

# sort by increasing ratio
mut2km<-mut2km[order(mut2km$mut2km_ratio),]
colnames(mut2km)[5:6]<-c(paste0(mut, "_mean"), paste0(mut, "2km", "_ratio"))

# output
out.name<-paste0("sortDARs_", mut, "2km_ratio.bed")
write.table(mut2km, out.name, quote=F, sep="\t", col.names=T, row.names=F)


