#!/usr/bin/env Rscript


#===================================================================================
# Set parameters
Usage<-function(){
	cat("Usage: Rscript peakgeneLinkPlot.r [peakgene_up] [peakgene_down] [mm9_genes] [consen_peaks] [out_suffix]\n\n",

	"Parameters:\n",
	"[peakgene_up]      significantly correlated up-regulated peak-gene list, i.e. \"fivskm_up_corSig.tsv\"\n",
	"[peakgene_down]    significantly correlated down-regulated peak-gene list, i.e. \"fivskm_down_corSig.tsv\"\n",
	"[mm9_genes]        mm9 genes in bed format, \"mm9_ensembl-symbol.bed\"\n",
	"[consen_peaks]     consensus peaks file, \"consensus_final.bed\"\n",
	"[out_suffix]       output suffix, e.g. fivskm, the file linkPlots_fivskm.pdf will be generated\n\n",

	"Function:\n",
	"Analyse peak-gene links, including distribution of peak-gene distance, genes per peak, peaks per gene, and gene numbers covered by peak-gene links, and plot them in bars.\n\n",

	"Example:\n",
	"Rscript peakgeneLinkPlot.r fivskm_up_corSig.tsv fivskm_down_corSig.tsv mm9_ensembl-symbol.bed consensus_final.bed fivskm\n",
	"# the output file \"linkPlots_fivskm.pdf\" will be generted\n\n",

	"Qirui Zhang (qirui.zhang@med.lu.se)\n",
	"Date: 15-09-2021\n\n"
	)
}

args<-commandArgs(TRUE)
if(length(args)!=5){Usage();quit();}


#=======================================================================================
# Load libraries & arguments
cat("loading libraries & arguments...\n\n")

library(ggplot2)
library(RColorBrewer)

up.file<-args[1]  # fivskm_up_corSig.tsv
down.file<-args[2]  # fivskm_down_corSig.tsv
mm9.file<-args[3]  # mm9_ensembl-symbol.bed
peak.file<-args[4]  # consensus_final.bed
out.suffix<-args[5]  # fivskm
out<-paste("linkPlots_", out.suffix, ".pdf", sep="")  # linkPlots_fivskm.pdf
pdf(out)


#=======================================================================================
# Read and process files
cat("reading and processing files...\n\n")

peakgene.up<-read.table(up.file, header=T, stringsAsFactors=F)
peakgene.up$distance<-peakgene.up$distance/1000
peakgene.down<-read.table(down.file, header=T, stringsAsFactors=F)
peakgene.down$distance<-peakgene.down$distance/1000
peakgene.all<-rbind(peakgene.up, peakgene.down)

mm9.genes<-read.table(mm9.file, header=F, stringsAsFactors=F)
colnames(mm9.genes)<-c("chr","start","end","gene")
consensus.peaks<-read.table(peak.file, header=F, stringsAsFactors=F)
colnames(consensus.peaks)<-c("chr", "start", "end", "peakID")


#=======================================================================================
# Analysis Part I: Plot peak-gene link distance distribution
cat("ploting peak-gene link distance distribution...\n\n")

peakgeneDistPlot<-function(peakgene, regulation){
    regul=regulation

#    title.count=paste("Peak-gene distance distribution in count", regul, sep=" ")
#    ggplot(peakgene, aes(x=distance))+geom_histogram(bins=30, fill=brewer.pal(9,"Greens")[6])+labs(title={title.count}, x="Distance to TSS (kb)", y="Count")+theme_bw()+theme(panel.grid=element_blank(), axis.text=element_text(size=14, color="black"), axis.title=element_text(size=16))+scale_x_continuous(limits=c(-500,500), breaks=c(-500,-250,0,250,500), labels=c(-500,-250,0,250,500))  # the real bins num is 29

    #---------------------------------------------
    total.count=nrow(peakgene)
    peakgeneDist.summary<-data.frame(Count=integer(), Count_sum=integer(), Frequency=numeric())
    bins.range<-seq(-500, 500, length.out=30)  # divide 1000kb (from -500kb to 500kb) into 29 bins
    bin.width=1000/29
    for(i in 1:29){bin.start=bins.range[i]; bin.end=bins.range[i+1]; bin.mid=mean(c(bin.start, bin.end)); if(i<29){bin.count=nrow(peakgene[which(peakgene$distance>=bin.start & peakgene$distance<bin.end),]);}else{bin.count=nrow(peakgene[which(peakgene$distance>=bin.start & peakgene$distance<=bin.end),])}; bin.freq=(bin.count/total.count)*100; peakgeneDist.summary[i,]=as.numeric(c(bin.mid, bin.count, bin.freq));}

    title.freq=paste("Peak-gene distance distribution in frequency", regul, sep=" ")
    ggplot(peakgeneDist.summary, aes(Count, Frequency))+geom_bar(stat="identity", width={bin.width}, fill=brewer.pal(9,"Greens")[6])+labs(title={title.freq}, x="Distance to TSS (kb)", y="Frequency")+theme_bw()+theme(panel.grid=element_blank(), axis.text=element_text(size=14, color="black"), axis.title=element_text(size=16))+scale_x_continuous(limits=c(-501,501), breaks=c(-500,-250,0,250,500), labels=c("-500","-250","0","250","500"))  # show frequency (in self-defined function it only can generate one plot); x-axis limits should include from -500 to 500, so the range is better a little larger than that.
}

peakgeneDistPlot(peakgene.up, "up")
peakgeneDistPlot(peakgene.down, "down")
peakgeneDistPlot(peakgene.all, "all")


#=======================================================================================
# Analysis Part II: Plot genes linked by per peak distribution
cat("plotting genes linked by per peak distribution...\n\n")

genesperPeakPlot<-function(peakgene, regulation){
    # sort peaks
    genesperPeak=peakgene
    genesperPeak=genesperPeak[order(genesperPeak$peakID),]

    # counts of genes for each peak
	for(i in 1:nrow(genesperPeak)){if(i==1){peak=genesperPeak[i,"peakID"]; genesperPeak[i,7]=1;next}; if(genesperPeak[i,"peakID"]==peak){genesperPeak[i,7]=genesperPeak[i-1,7]+1; genesperPeak[i-1,8]="remove";}else{genesperPeak[i-1,8]="keep"; genesperPeak[i,7]=1; peak=genesperPeak[i,"peakID"]}; if(i==nrow(genesperPeak)){genesperPeak[i,8]="keep";}}

    # unique peaks
	genesperPeak<-genesperPeak[which(genesperPeak$V8=="keep"),c("peakID","V7")]
	colnames(genesperPeak)[2]<-"Count"

    # compute frequency of each count in total counts
    total.count=nrow(genesperPeak)
    genesperPeak.summary<-data.frame(Count=integer(), Count_sum=integer(), Frequency=numeric())
    for(i in 1:10){  # only show 10 bars, so set 1:10. all rest counts >= 10 will be shown in the last bar together.
        if(i<10){count.sum=nrow(genesperPeak[which(genesperPeak$Count==i),])}else{count.sum=nrow(genesperPeak[which(genesperPeak$Count>=i),])}
        freq=(count.sum/total.count)*100
        genesperPeak.summary[i,]=as.numeric(c(i, count.sum, freq))  # save count, count.sum, frequency to the new dataframe genesperPeak.summary.
        if(i==10){rownames(genesperPeak.summary)[i]=">=10"}  # modify the last bar's name
    }

    # plot
    regul=regulation
#    title.count=paste("Genes per peak distribution in counts", regul, sep=" ")
#    ggplot(genesperPeak, aes(x=Count))+geom_histogram(bins=12, fill=brewer.pal(9,"Blues")[6])+labs(title={title.count}, x="Number of genes linked by per peak", y="Count")+theme_bw()+theme(panel.grid=element_blank(), axis.text=element_text(size=14, color="black"), axis.title=element_text(size=16))+scale_x_continuous(limits=c(0,11), breaks=c(1,5,10), labels=c(1,5,">=10"))  # show counts

    title.freq=paste("Genes per peak distribution in frequency", regul, sep=" ")
    ggplot(genesperPeak.summary, aes(Count, Frequency))+geom_bar(stat="identity", width=1.0, fill=brewer.pal(9,"Blues")[6])+labs(title={title.freq}, x="Number of genes linked by per peak", y="Frequency")+theme_bw()+theme(panel.grid=element_blank(), axis.text=element_text(size=14, color="black"), axis.title=element_text(size=16))+scale_x_continuous(limits=c(0,11), breaks=c(1,5,10), labels=c(1,5,">=10"))  # show frequency (in self-defined function it only can generate one plot)
}

genesperPeakPlot(peakgene.up, "up")
genesperPeakPlot(peakgene.down, "down")
genesperPeakPlot(peakgene.all, "all")


#=======================================================================================
# Analysis Part III:Plot peaks linked by per gene distribution
cat("plotting peaks linked by per gene distribution...\n\n")

peaksperGenePlot<-function(peakgene, regulation){
    # sort genes
    peaksperGene=peakgene
    peaksperGene=peaksperGene[order(peaksperGene$geneID),]

    # counts of peaks for each gene
    for(i in 1:nrow(peaksperGene)){if(i==1){gene=peaksperGene[i,"geneID"]; peaksperGene[i,7]=1;next}; if(peaksperGene[i,"geneID"]==gene){peaksperGene[i,7]=peaksperGene[i-1,7]+1; peaksperGene[i-1,8]="remove";}else{peaksperGene[i-1,8]="keep"; peaksperGene[i,7]=1; gene=peaksperGene[i,"geneID"]}; if(i==nrow(peaksperGene)){peaksperGene[i,8]="keep";}}

    # unique peaks
    peaksperGene<-peaksperGene[which(peaksperGene$V8=="keep"),c("geneID","V7")]
    colnames(peaksperGene)[2]<-"Count"

    # compute frequency of each count in total counts
    total.count=nrow(peaksperGene)
    peaksperGene.summary<-data.frame(Count=integer(), Count_sum=integer(), Frequency=numeric())
    for(i in 1:10){  # only show 10 bars, so set 1:10. all rest counts >= 10 will be shown in the last bar together.
        if(i<10){count.sum=nrow(peaksperGene[which(peaksperGene$Count==i),])}else{count.sum=nrow(peaksperGene[which(peaksperGene$Count>=i),])}
        freq=(count.sum/total.count)*100
        peaksperGene.summary[i,]=as.numeric(c(i, count.sum, freq))  # save count, count.sum, frequency to the new dataframe peaksperGene.summary.
        if(i==10){rownames(peaksperGene.summary)[i]=">=10"}  # modify the last bar's name
    }

    # plot
    regul=regulation
#    title.count=paste("Peaks per gene distribution in counts", regul, sep=" ")
#    ggplot(peaksperGene, aes(x=Count))+geom_histogram(bins=12, fill=brewer.pal(9,"Reds")[6])+labs(title={title.count}, x="Number of peaks linked by per gene", y="Count")+theme_bw()+theme(panel.grid=element_blank(), axis.text=element_text(size=14, color="black"), axis.title=element_text(size=16))+scale_x_continuous(limits=c(0,11), breaks=c(1,5,10), labels=c(1,5,">=10")) # show counts

    title.freq=paste("Peaks per gene distribution in frequency", regul, sep=" ")
    ggplot(peaksperGene.summary, aes(Count, Frequency))+geom_bar(stat="identity", width=1.0, fill=brewer.pal(9,"Reds")[6])+labs(title={title.freq}, x="Number of peaks linked by per gene", y="Frequency")+theme_bw()+theme(panel.grid=element_blank(), axis.text=element_text(size=14, color="black"), axis.title=element_text(size=16))+scale_x_continuous(limits=c(0,11), breaks=c(1,5,10), labels=c(1,5,">=10"))  # show frequency (in self-defined function it only can generate one plot)
}

peaksperGenePlot(peakgene.up, "up")
peaksperGenePlot(peakgene.down, "down")
peaksperGenePlot(peakgene.all, "all")


#=======================================================================================
# Analysis Part IV: Count number of genes covered by peak-gene links and plot
cat("ploting num of genes covered by peak-gene links distribution...\n\n")

genesinPeakGenePlot<-function(peakgene, regulation){
    mm9.genes<-mm9.genes[order(mm9.genes$chr, mm9.genes$gene, mm9.genes$start),]  # sort genes
    peakgene.interval<-peakgene[,c("peakID","geneID")]  # create peak-gene interval dataframe

    #--------------------------------------------------------------------------------
    # extract peak coordinates, gene coordinates, and peak-gene interval coordinates
    for(i in 1:nrow(peakgene)){
        peak_id=peakgene[i,"peakID"]
        gene_id=peakgene[i,"geneID"]

        peak.coord=consensus.peaks[which(consensus.peaks$peakID==peak_id),]
        chr_id=peak.coord[1,"chr"]
        peak.start=peak.coord[1,"start"]
        peak.end=peak.coord[1,"end"]

        gene.coord=mm9.genes[which(mm9.genes$chr==chr_id & mm9.genes$gene==gene_id),]
        if(nrow(gene.coord)>1){  # some genes have more than one ensembl_ID and coordinates
            gene.coord$center=apply(gene.coord[,2:3],1,mean)  # gene coordinate center
            peak.coord$center=apply(peak.coord[,2:3],1,mean)  # peak coordinate center
            gene.coord$dist=abs(gene.coord$center-peak.coord$center)  # distance between peak & gene center
            gene.coord=gene.coord[which(gene.coord$dist==min(gene.coord$dist)),]  # keep the coordinate that is most closest to the peak coordinate
        }
        gene.start=gene.coord[1,"start"]
        gene.end=gene.coord[1,"end"]

        if(peak.end<gene.start){  # peak in the left side of gene
            interval.start=peak.end
            interval.end=gene.start
        }else if(peak.start>gene.end){  # peak in the right side of gene
            interval.start=gene.end
            interval.end=peak.start
        }else{interval.start=0; interval.end=0}  # peak cross with gene

        peakgene.interval[i,3:9]=c(chr_id, peak.start, peak.end, gene.start, gene.end, interval.start, interval.end)
    }
    colnames(peakgene.interval)[3:9]<-c("chr", "peak_start", "peak_end", "gene_start", "gene_end", "interval_start", "interval_end")
    for(i in 4:ncol(peakgene.interval)){peakgene.interval[,i]=as.integer(peakgene.interval[,i])}  # modify coordinate columns to integer format

    #--------------------------------------------------------------------------------
    # count how many genes are covered by (or in between) peak-gene links (including target gene itself)
    for(i in 1:nrow(peakgene.interval)){
        chr_id=peakgene.interval[i,"chr"]
        interval.start=peakgene.interval[i,"interval_start"]
        interval.end=peakgene.interval[i,"interval_end"]

        intersect=mm9.genes[which((mm9.genes$chr==chr_id & mm9.genes$start>interval.start & mm9.genes$start<interval.end) | (mm9.genes$chr==chr_id & mm9.genes$end>interval.start & mm9.genes$end<interval.end)),]  # extract all genes that intersect with this peak-gene link

        gene.num=length(unique(intersect[,"gene"]))+1  # "+1" to include target gene itself
        peakgene.interval[i,10]=gene.num
    }
    colnames(peakgene.interval)[10]<-"gene_count"

    # generate gene_count (num of genes covered by peak-gene links) summary dataframe
    total.count=nrow(peakgene.interval)
    genesinPeakGene.summary=data.frame(Count=integer(), Count_sum=integer(), Frequency=numeric())
    for(i in 1:25){if(i<25){count.sum=nrow(peakgene.interval[which(peakgene.interval$gene_count==i),])}else{count.sum=nrow(peakgene.interval[which(peakgene.interval$gene_count>=i),])}; freq=(count.sum/total.count)*100; genesinPeakGene.summary[i,]=as.numeric(c(i, count.sum, freq)); if(i==25){rownames(genesinPeakGene.summary)[i]=">=25"}}

    #--------------------------------------------------------------------------------
    # plot
    regul=regulation
    #title.count=paste("Genes covered by peak-gene links distribution in count", regul, sep=" ")
    #ggplot(peakgene.interval, aes(x=gene_count))+geom_histogram(bins=27, fill=brewer.pal(9,"Purples")[6])+labs(title={title.count}, x="Number of genes covered by peak-gene links", y="Count")+theme_bw()+theme(panel.grid=element_blank(), axis.text=element_text(size=14, color="black"), axis.title=element_text(size=16))+scale_x_continuous(limits=c(0,26), breaks=c(1,5,10,15,20,25), labels=c(1,5,10,15,20,">=25")) # show counts

    title.freq=paste("Genes covered by peak-gene links distribution in frequency", regul, sep=" ")
    ggplot(genesinPeakGene.summary, aes(Count, Frequency))+geom_bar(stat="identity", width=1.0, fill=brewer.pal(9,"Purples")[6])+labs(title={title.freq}, x="Number of genes covered by peak-gene links", y="Frequency")+theme_bw()+theme(panel.grid=element_blank(), axis.text=element_text(size=14, color="black"), axis.title=element_text(size=16))+scale_x_continuous(limits=c(0,26), breaks=c(1,5,10,15,20,25), labels=c(1,5,10,15,20,">=25"))  # show frequency (in self-defined function it only can generate one plot)
}

genesinPeakGenePlot(peakgene.up, "up")
genesinPeakGenePlot(peakgene.down, "down")
genesinPeakGenePlot(peakgene.all, "all")


cat("done with analysis!\n\n")


