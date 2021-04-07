### Plots and stats for het. vs. chr. length and distance from end of chromosome

library(plyr)

# Significance symbols to add to plots
starfn=function(pval){
	if (pval<0.001) {return('***')}
	else if (pval<0.01) {return('**')}
	else if (pval<0.05) {return('*')}
	else {return('NS')}
}

setwd("~/condor/analysis/winHet/condorRef")
chrkey=read.table("condor_chr_scaff_length.txt", header=T, sep="\t", stringsAsFactors=F)

pdf("chromosome_end_het.pdf", width=6.5, height=6, pointsize=10)
par(mfrow=c(3,4))

setwd("~/condor/analysis/winHet/condorRef/1Mbwins")
hetdatafiles=list.files(pattern="step.txt")

winsize=1e6

# Get data files, change chromosome names, exclude chrZ and small scaffolds not in key
get_data=function(){
	hfiles=hetdatafiles[grep(pattern=searchstring, hetdatafiles)]
	allhet=ldply(hfiles, read.table, header=TRUE, sep="\t")
	allhet=allhet[which(as.character(allhet$chrom) %in% chrkey$SCAFF),]
	chrom=NULL
	for (i in 1:length(allhet$chrom)){ chrom[i]=chrkey[which(chrkey$SCAFF==allhet[i,]$chrom),]$CHR }
	allhet=cbind(chrom,allhet[,2:length(names(allhet))])
	allhet=allhet[which(allhet$chrom %in% 1:100),]
	allhet=allhet[with(allhet, order(as.numeric(as.character(allhet$chrom)), as.numeric(allhet$window_start))),]
	temp=allhet[which(allhet[,4]>=(0.5*winsize)),]
}


################################################################################
# Plot mean chromosome het as a function of chromosome length

chet_scatter=function(samplename, labely){
	# Get chromosome lengths and mean het. per chromosome
	cc=NULL
	cl=NULL
	ch=NULL
	for (i in 1:length(unique(temp$chrom))){
		cc[i]=as.character(unique(temp$chrom)[i])
		cl[i]=chrkey[which(chrkey$CHR==cc[i]),]$LENGTH
		tt=temp[which(temp$chrom==cc[i]),]
		ch[i]=sum(tt[,5], na.rm=T)/sum(tt[,4], na.rm=T)
	}
	
	# Plot
	par(mar=c(5,3.5,2,1))
	plot(log(cl, base=10), 1000*ch, xlim=c(6.25, 8.5), ylim=c(0,3.5), xlab="", ylab="")
	title(main=samplename, font.main=1, cex.main=1, line=0.5)
	title(xlab=expression("Chromosome length (log"[10]*")"), line=2.5)
	if (labely){title(ylab="Heterozygosity per kb", line=2.5)}
	
	# Add regression line and Pearson correlation
	reg=lm(1000*ch ~ log(cl, base=10))
	abline(reg)
	pcorr=cor.test(log(cl, base=10), 1000*ch, method="pearson")
	text(x=7.95, y=3.45, paste("R: ", sprintf("%.3f", pcorr$estimate), starfn(pcorr$p.value), sep=""))
}

# Plots
searchstring="CYW1141"
samplename="CA condor 1"
get_data()
chet_scatter(samplename, labely=T)

searchstring="CRW1112"
samplename="CA condor 2"
get_data()
chet_scatter(samplename, labely=F)

searchstring="VulGry1"
samplename="Andean condor"
get_data()
chet_scatter(samplename, labely=F)

searchstring="BGI_N323"
samplename="Turkey vulture"
get_data()
chet_scatter(samplename, labely=F)


################################################################################
# Boxplots for het in distal versus proximal windows

end_het_box=function(samplename, labely){
	# Divide distal and proximal windows
	dist=10e6
	ends=NULL
	mids=NULL
	for (i in 1:length(unique(temp$chrom))){
		chrom=as.character(unique(temp$chrom)[i])
		ll=chrkey[which(chrkey$CHR==chrom),]$LENGTH
		tt=temp[which(temp$chrom==chrom),]
		if (dim(tt)[1]>0){
			ends=rbind(ends, tt[which((tt$window_start+winsize-1)<=dist),])
			ends=rbind(ends, tt[which((tt$window_start-1)>=(ll-dist)),])
			mids=rbind(mids, tt[which((tt$window_start+winsize-1)>dist & (tt$window_start-1)<(ll-dist)),])
		}
	}
	
	# Calculate het and statistical significance
	aa=1000*ends[,5]/ends[,4]
	bb=1000*mids[,5]/mids[,4]
	wt=wilcox.test(aa, bb, alternative="greater")
	print(wt)
	
	# Plot
	ymin=0
	ymax=6.5
	yrange=ymax-ymin
	
	par(mar=c(5,3.5,2,1))
	par(xpd=T)
	boxplot(aa, bb, ylim=c(ymin, ymax), boxwex=0.25, names=c("Distal", "Proximal"), xlab="")
	title(main=samplename, font.main=1, cex.main=1, line=0.5)
	if (labely){title(ylab="Heterozygosity per kb", line=2.5)}
	arrows(1, ymax-(.075*yrange), 2, ymax-(.075*yrange), angle=90, code=3, length=.05)
	text(1.5, ymax-(.025*yrange), starfn(wt$p.value))
	}

# Plots
searchstring="CYW1141"
samplename="CA condor 1"
get_data()
end_het_box(samplename, labely=T)

searchstring="CRW1112"
samplename="CA condor 2"
get_data()
end_het_box(samplename, labely=F)

searchstring="VulGry1"
samplename="Andean condor"
get_data()
end_het_box(samplename, labely=F)

searchstring="BGI_N323"
samplename="Turkey vulture"
get_data()
end_het_box(samplename, labely=F)


################################################################################
# Boxplots for calls in distal versus proximal windows

end_call_box=function(samplename, labely){
	# Divide distal and proximal windows
	dist=10e6
	ends=NULL
	mids=NULL
	for (i in 1:length(unique(temp$chrom))){
		chrom=as.character(unique(temp$chrom)[i])
		ll=chrkey[which(chrkey$CHR==chrom),]$LENGTH
		tt=temp[which(temp$chrom==chrom),]
		if (dim(tt)[1]>0){
			ends=rbind(ends, tt[which((tt$window_start+winsize-1)<=dist),])
			ends=rbind(ends, tt[which((tt$window_start-1)>=(ll-dist)),])
			mids=rbind(mids, tt[which((tt$window_start+winsize-1)>dist & (tt$window_start-1)<(ll-dist)),])
		}
	}
	
	# Calculate call rate and statistical significance
	aa=1e-6*ends[,4]
	bb=1e-6*mids[,4]
	wt=wilcox.test(aa, bb, alternative="less")
	print(wt)
	
	# Plot
	ymin=0.5
	ymax=1
	yrange=ymax-ymin
	
	par(mar=c(5,3.5,2,1))
	par(xpd=T)
	boxplot(aa, bb, ylim=c(ymin, ymax), boxwex=0.25, names=c("Distal", "Proximal"), xlab="")
	title(main=samplename, font.main=1, cex.main=1, line=0.5)
	if (labely){title(ylab="% of window called", line=2.5)}
	arrows(1, ymax-(.075*yrange), 2, ymax-(.075*yrange), angle=90, code=3, length=.05)
	text(1.5, ymax-(.025*yrange), starfn(wt$p.value))
	}

# Plots
searchstring="CYW1141"
samplename="CA condor 1"
get_data()
end_call_box(samplename, labely=T)

searchstring="CRW1112"
samplename="CA condor 2"
get_data()
end_call_box(samplename, labely=F)

searchstring="VulGry1"
samplename="Andean condor"
get_data()
end_call_box(samplename, labely=F)

searchstring="BGI_N323"
samplename="Turkey vulture"
get_data()
end_call_box(samplename, labely=F)


dev.off()

