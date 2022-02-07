#!/usr/bin/env Rscript


#===============================================================================================
Usage<-function(){cat("Usage: Rscript TSSEscore.r [matrix]\n\n",

	"Parameters:\n",
	"[matrix]	Modified matrix that was originally generated from deepTools computeMatrix, removed headers (first 2 rows)\n\n",

	"Function:\n",
	"Calculate TSS Enrichment Score of provided matrix following the definition by ENCODE. TSS region includes 2kb in each side (upstream and downstream), and another extended 100bp. The only information generated is a TSSE score, not saving any files. \n\n",

	"Example:\n",
	"Rscript TSSEscore.r m16-8_modified.mx\n\n",

	"Author and date:\n",
	"Qirui Zhang (qirui.zhang@med.lu.se)\n",
	"01/12/2020\n\n"
	)
}

args<-commandArgs(TRUE)
if(length(args)!=1){Usage();quit();}


#================================================================================
cat("\n==========================================================================\n")
cat("Calculating sample", args[1], "\n")

# Read input.matrix
mx<-read.table(args[1], header=T)
bins.mean<-colMeans(mx)

# Calculate background noise mean value
bgnoise.avg<-mean(c(head(bins.mean,10), tail(bins.mean,10)))

# Divide each bins by background noise
bins.enrich<-bins.mean/bgnoise.avg

# Use maximum value as TSSE Score (defined by ENCODE)
tsse.score<-max(bins.enrich)

cat("TSSE score:", tsse.score, "\n")
cat("bin location:\n")
which.max(bins.enrich)



