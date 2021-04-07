### Plot het. in sliding windows across the genome, plus histograms of window-based het.
# Input data: output of SlidingWindowHet.py with 1 Mb win size and 1 Mb step (non-overlapping windows)

library(plyr)

setwd("~/condor/analysis/winHet/condorRef/1Mbwins")
hetdatafiles=list.files(pattern="step.txt")

winsize=1e6
dkcol="#005a32"
ltcol="#259349"

pdf("plot_het_wins_hist.pdf", width=6.85, height=5.5, pointsize=8)

split.screen(rbind(
c(0.0, 0.75, 0.75, 1.0),
c(0.75, 1.0, 0.75, 1.0),
c(0.0, 0.75, 0.50, 0.75),
c(0.75, 1.0, 0.50, 0.75),
c(0.0, 0.75, 0.25, 0.5),
c(0.75, 1.0, 0.25, 0.5),
c(0.0, 0.75, 0.0, 0.25),
c(0.75, 1.0, 0.0, 0.25)))


plotwinhet=function(samplename, labelx, panellabel){
	# Calculate mean het. and generate plot title
	meanhet=sum(allhet[,5], na.rm=T)/sum(allhet[,4], na.rm=T)
	plotname=paste(samplename, "\nmean het. = ", sprintf("%.3f", 1000*meanhet), " per kb", sep="")
	
	# Get positions of chromosome boundaries, midpoints
	pos=as.numeric(rownames(unique(data.frame(temp$chrom)[1])))
	pos=append(pos,length(temp$chrom))
	numpos=NULL
	for (i in 1:length(pos)-1){numpos[i]=(pos[i]+pos[i+1])/2}

	# Define colors for each chr
	colset=c(dkcol,ltcol)
	mycols=NULL
	mycols[1]=colset[1]
	for (i in 2:length(temp$chrom)){
		if (temp$chrom[i]==temp$chrom[i-1]){mycols[i]=mycols[i-1]}
		else {mycols[i]=colset[which(colset!=mycols[i-1])]}
	}
	
	# Plot
	par(mar=c(3,3.75,2,0))
	par(xpd=T)
	b=barplot(1000*temp[,5]/temp[,4], ylim=c(0,6), border=mycols, col=mycols, ylab="", main="", axes=F)
	axis(side=2, line=-0.75, labels=F)
	axis(side=2, line=-1, tick=F, labels=T)
	title(main=plotname, font.main=1, cex.main=1)
	title(ylab="Heterozygosity per kb", line=1.25)
	axis(side=1, at=b[pos], labels=F)
	mylabels=as.character(unique(temp$chrom))
	mylabels[c(11,13,15,17,19,21,22,23,25,26,27)]=""
	axis(side=1, at=b[numpos], tick=F, labels=mylabels, las=3, line=-.25, cex.axis=1)
}


plotwinhethist=function(samplename, labelx, panellabel){
	par(mar=c(3,3.75,2,0.5))
	par(xpd=T)
	mycol=dkcol
	
	# Get bin counts
	h=hist(1000*temp[,5]/temp[,4], breaks=seq(0,10, by=.1), plot=F)
	
	# Plot
	barplot(c(h$counts[1:59],sum(h$counts[60:length(h$counts)])), ylim=c(0,250), space=0, border=mycol, col=mycol, xlab="", ylab="")
	title(ylab="# of windows", line=2)
	axis(1, at=seq(0,60,10), labels=F)
	axis(1, at=seq(0,60,10), tick=F, labels=c(0:5,"6+"), line=-0.4)
	if (labelx){title(xlab="Het. per kb", line=1.75)}
}


searchstring="CYW1141"
hfiles=hetdatafiles[grep(pattern=searchstring, hetdatafiles)]
allhet=ldply(hfiles, read.table, header=TRUE, sep="\t")
temp=allhet[which(allhet[,4]>=(0.5*winsize)),]
screen(1)
plotwinhet("CA condor 1")
screen(2)
plotwinhethist("CA condor 1", labelx=F)

searchstring="CRW1112"
hfiles=hetdatafiles[grep(pattern=searchstring, hetdatafiles)]
allhet=ldply(hfiles, read.table, header=TRUE, sep="\t")
temp=allhet[which(allhet[,4]>=(0.5*winsize)),]
screen(3)
plotwinhet("CA condor 2")
screen(4)
plotwinhethist("CA condor 2", labelx=F)

searchstring="VulGry1"
hfiles=hetdatafiles[grep(pattern=searchstring, hetdatafiles)]
allhet=ldply(hfiles, read.table, header=TRUE, sep="\t")
temp=allhet[which(allhet[,4]>=(0.5*winsize)),]
screen(5)
plotwinhet("Andean condor")
screen(6)
plotwinhethist("Andean condor", labelx=F)

searchstring="BGI_N323"
hfiles=hetdatafiles[grep(pattern=searchstring, hetdatafiles)]
allhet=ldply(hfiles, read.table, header=TRUE, sep="\t")
temp=allhet[which(allhet[,4]>=(0.5*winsize)),]
screen(7)
plotwinhet("Turkey vulture")
title(xlab="Chromosome", line=1.75)
screen(8)
plotwinhethist("Turkey vulture", labelx=T)


dev.off()


