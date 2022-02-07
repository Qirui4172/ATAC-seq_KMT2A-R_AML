#!/usr/bin/env Rscript


#===================================================================================
# Set parameters
Usage<-function(){
	cat("Usage: Rscript replicatesCorrelation.r [atacseqMx]\n\n",

	"Parameters:\n",
	"[atacseqMx]     atac-seq_deseq2norm.tsv, is the output file of \"findDiffpeaks.r\".\n\n",

	"Function:\n",
	"Compute correlatin of pairwise replicates within groups and plot.\n\n",

	"Example:\n",
	"Rscript replicatesCorrelation.r atac-seq_deseq2norm.tsv\n",
	"  # the file \"replicatesCorr_plots.pdf\" will be generted\n\n",

	"Qirui Zhang (qirui.zhang@med.lu.se)\n",
	"Date: 16-09-2021\n\n"
	)
}

args<-commandArgs(TRUE)
if(length(args)!=1){Usage();quit();}


#===================================================================================
# load libraries, arguments and read files
cat("loading libraries, arguments and reading files...\n\n")

library(ggplot2)
atac.file<-args[1]  # atac-seq_deseq2norm.tsv
atac.mx<-read.table(atac.file, header=T, stringsAsFactors=F)
atac.mx<-subset(atac.mx, select=-mean)
atac.mx<-log2(atac.mx+0.5)
pdf("replicatesCorr_plots.pdf")


# create position dataframe
atac.column<-data.frame(start=c(1, 5, 9, 12), end=c(4, 8, 11, 15), fullname=c("FLT3-ITD", "FLT3-N676K", "NRAS-G12D", "KMT2A"))
rownames(atac.column)<-c("fi", "fn", "ng", "km")


# compute replicates correlation and plot
cat("computing replicates correlation and ploting...\n\n")
for(i in rownames(atac.column)){
    start=atac.column[i,"start"]
    end=atac.column[i,"end"]
    for(j in start:(end-1)){
        for(k in (j+1):end){
            cor.result=cor.test(atac.mx[,j], atac.mx[,k])
            coeff=cor.result$estimate[[1]]
            pvalue=cor.result$p.value

            sample_name.x=colnames(atac.mx)[j]
            sample_name.y=colnames(atac.mx)[k]
            title=paste(sample_name.x, "-", sample_name.y, " cor=", round(coeff,2), sep="")
            cor.plot=ggplot(atac.mx, aes(atac.mx[,j], atac.mx[,k]))+geom_point(alpha=0.3)+geom_smooth(method="lm", color="red", se=FALSE)+theme_bw()+theme(panel.grid=element_blank(), panel.border=element_blank(), axis.text=element_text(size=14, color="black"), axis.title=element_text(size=16), axis.line=element_line(color="black", size=0.5), legend.title=element_text(size=14))+labs(x={sample_name.x}, y={sample_name.y}, title={title})
            print(cor.plot)  # if not use print() it will not generate plots!
        }
    }
}

dev.off()

cat("Done with analysis!\n\n")


