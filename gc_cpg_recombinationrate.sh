# GC content + CpG island identification and recombination rate analysis

### Get GC content in non-overlapping 1 Mb windows
# Note: same coordinates as iSMC 1 Mb windows

TWOBIT=gc_PacBio_HiC.2bit
DATA=CYW1141_ismc.rho.1Mb.bedgraph_scaffnames_sorted.bed

echo -e "chromosome\tstart\tend\tbase_count\tgc_proportion\trho" > ${TWOBIT%.2bit}_CYW1141_ismc.rho.1Mb_GC.txt
while read -r CHR START END RHO; do
twoBitToFa ${TWOBIT}:${CHR}:${START}-${END} temp.fa
countgc.sh in=temp.fa out=temp.fa_gc.txt format=4
C=$(cut -f2 temp.fa_gc.txt)
GC=$(cut -f3 temp.fa_gc.txt)
echo -e "${CHR}\t${START}\t${END}\t${C}\t${GC}\t${RHO}" >> ${TWOBIT%.2bit}_CYW1141_ismc.rho.1Mb_GC.txt
done < ${DATA}


### Get CpG islands (not within repeats)
# cpg_lh part of UCSC tools (http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/

FASTA=gc_PacBio_HiC.fasta
MASK=gc_PacBio_HiC_repeats_TRF_WMdust.bed

maskOutFa ${FASTA} ${MASK} ${FASTA%.fa*}_hardMask.fa
cpg_lh ${FASTA%.fa*}_hardMask.fa > ${FASTA%.fa*}_cpgIslands_norepeats.tmp
awk '{printf "%s\t%s\t%s\n", $1, $2-1, $3}' ${FASTA%.fa*}_cpgIslands_norepeats.tmp \
| sort -k1,1 -k2,2n > ${FASTA%.fa*}_cpgIslands_norepeats.bed


### Get distance to nearest CpG island for each 1 kb window of iSMC output
# Note: both inputs must already be sorted

CGI=gc_PacBio_HiC_cpgIslands_norepeats.bed
IN=CYW1141_ismc.rho.1kb.bedgraph_sorted.bed
OUT=CYW1141_ismc.rho.1kb.bedgraph_sorted_distToNearestCgi.bed

bedtools closest -d -t first -a ${IN} -b ${CGI} > ${OUT}


################################################################################
### Plotting in R

### Plot GC content across the genome

# Read in data and rename chromosomes
df=read.table("gc_PacBio_HiC_CYW1141_ismc.rho.1Mb_GC.txt", header=T)
chrkey=read.table("condor_chr_scaff_length.txt", header=T, sep="\t", stringsAsFactors=F)
df=df[which(df$chromosome %in% chrkey$SCAFF),]
chromosome=NULL
for (i in 1:length(df$chromosome)){
	chromosome[i]=chrkey[which(chrkey$SCAFF==df[i,]$chromosome),]$CHR
}
df[,1]=chromosome
temp=df[which(df$chromosome %in% 1:100),]
temp=temp[with(temp, order(as.numeric(as.character(temp$chromosome)), as.numeric(temp$start))),]

# Get chromosome boundaries, midpoints
pos=as.numeric(rownames(unique(data.frame(temp$chromosome)[1])))
pos=append(pos, length(temp$chromosome))
numpos=NULL
for (i in 1:length(pos)-1){numpos[i]=(pos[i]+pos[i+1])/2}

# Define gray bg color and chromosome colors
bgcol="#ebebeb"
c1=c(0, 0, 0)/255
c2=c(171, 0, 0)/255
colset=c(rgb(c1[1], c1[2], c1[3], alpha=1), rgb(c2[1], c2[2], c2[3], alpha=1))

# Define chromosome label vector
mylabels=as.character(unique(temp$chromosome))
mylabels[c(15,17,19,21,22,23,25,26,27,28)]=""

# Plot (using no internal padding with xaxs="i" but then adding some padding (10) manually on either side of data)
ymin=0.35
ymax=0.625

par(mar=c(3,4,1.5,0.5))
plot(0, 0, type="n", xlim=c(0, nrow(temp)+20), ylim=c(ymin, ymax), xlab="", ylab="", frame.plot=F, axes=F, xaxs="i")
title(ylab="GC content (%)", line=2.5)
axis(side=2)
axis(side=1, at=pos+10, labels=F)
axis(side=1, at=numpos+10, tick=F, labels=mylabels, las=3, line=-.25, cex.axis=1)
title(xlab="Chromosome", line=1.75)
# Add gray stripes to background
for (i in seq(ymin, ymax, by=.1)){rect(-1000, i, nrow(temp)+1000, i+.05, col=bgcol, border=NA)}

# Add chromosome lines
chrs=as.numeric(unique(temp$chromosome))
# Add padding
xx=10
mycol=colset[1]
for (i in 1:length(chrs)){
  tt=temp[which(temp$chromosome==chrs[i]),]
  lines(1:nrow(tt)+xx, tt$gc_proportion, col=mycol)
  xx=xx+nrow(tt)
  mycol=colset[which(colset!=mycol)]
}

box()


################################################################################
### Boxplot of GC content in distal versus proximal regions
# Note: Uses same input data as above (GC content across the genome)

# Significance symbols to add to plot
starfn=function(pval){
	if (pval<0.001) {return('***')}
	else if (pval<0.01) {return('**')}
	else if (pval<0.05) {return('*')}
	else {return('NS')}
}

# Including or excluding windows with fewer sites does not affect the results
temp=temp[which(temp$base_count>=1e6),]

# Divide windows into ends (distal, within 10 Mb of ends) and mids (proximal, >10 Mb from ends)
dist=10e6
ends=NULL
mids=NULL

for (i in 1:length(unique(temp$chromosome))){
	chrom=unique(temp$chromosome)[i]
	ll=chrkey[which(chrkey$CHR==chrom),]$LENGTH
	tt=temp[which(temp$chromosome==chrom),]
	if (dim(tt)[1]>0){
		ends=rbind(ends, tt[which(tt$end<=dist),])
		ends=rbind(ends, tt[which(tt$start>=(ll-dist)),])
		mids=rbind(mids, tt[which(tt$end>dist & tt$start<(ll-dist)),])
	}
}

# Test if GC content is higher in distal regions
wt=wilcox.test(ends$gc_proportion, mids$gc_proportion, alternative="greater")
wt
# 	Wilcoxon rank sum test with continuity correction
# 
# data:  ends$gc_proportion and mids$gc_proportion
# W = 258390, p-value < 2.2e-16
# alternative hypothesis: true location shift is greater than 0

# Plot
ymin=0.35
ymax=0.625
yrange=ymax-ymin

par(mar=c(3,2.5,1.5,0.5))
b=boxplot(ends$gc_proportion, mids$gc_proportion, ylim=c(ymin, ymax), boxwex=0.25, names=mynames, xlab="")

# Add significance of wt
arrows(1, ymax-(.075*yrange), 2, ymax-(.075*yrange), angle=90, code=3, length=.05)
text(1.5, ymax-(.025*yrange), starfn(wt$p.value))


################################################################################
### Recombination rate by distance from CpG islands plus bootstrap

# Read in rho data annotated with distance to nearest CGI
# Columns: chromosome, window start, window end, mean rho/bp, distance to nearest CGI
df=read.table("CYW1141_ismc.rho.1kb.bedgraph_sorted_distToNearestCgi.bed", header=F, sep='\t')
df=df[,c(1:4,8)]
colnames(df)=c("chrom","start_pos","end_pos","rho","cgi_dist")

# Keep only windows within a particular distance from CGI (in case input was not already filtered for this)
max_dist_to_CGI=100e3
df=df[which(df$cgi_dist <= max_dist_to_CGI),]

# Discard truncated windows (found at the ends of chromosomes)
window_size=1e3
df=df[which(df$end_pos-df$start_pos == window_size),]

### Get mean rho across windows binned by distance to nearest CGI

# Generate bins 
bin_size=1e3
bin_starts=seq(0, max_dist_to_CGI-bin_size, by=1e3)
bin_ends=bin_starts+bin_size

# Calculate mean rho per distance bin
mean_rho_per_dist_bin=NULL

for (i in 1:length(bin_starts)){
	# Get windows within distance bin
	tt=df[which(df$cgi_dist >= bin_starts[i] & df$cgi_dist < bin_ends[i]),]
	# Get mean rho (weighted)
	if (dim(tt)[1] > 0){
		n_sites=tt$end_pos-tt$start_pos
		total_rho=tt$rho*n_sites
		mean_rho_per_dist_bin[i]=sum(total_rho)/sum(n_sites)
	} 
	else {
		mean_rho_per_dist_bin[i]=NA
	}
}

# Boootstrap (resample rows from df with replacement)
n_boot=100
boot_mat=matrix(nrow=n_boot, ncol=length(bin_starts))

for (b in 1:n_boot){
	print(b)
	newrows=sample(1:dim(df)[1], size=dim(df)[1], replace=T)
	df_boot=df[newrows,]
	mean_rho_per_dist_bin_boot=NULL
	for (i in 1:length(bin_starts)){
		# Get windows within distance bin
		tt=df_boot[which(df_boot$cgi_dist >= bin_starts[i] & df_boot$cgi_dist < bin_ends[i]),]
		# Get mean rho (weighted)
		if (dim(tt)[1] > 0){
			n_sites=tt$end_pos-tt$start_pos
			total_rho=tt$rho*n_sites
			mean_rho_per_dist_bin_boot[i]=sum(total_rho)/sum(n_sites)
		} 
		else {
			mean_rho_per_dist_bin_boot[i]=NA
		}
	}
	boot_mat[b,]=mean_rho_per_dist_bin_boot
}	
	
scale=1e3
myylab=expression(paste("Mean ", rho, "/bp x 10"^-3, sep=""))

# Set plot y-axis limits by getting span of mean rho values and adding 5% padding
yrange=max(rbind(mean_rho_per_dist_bin, boot_mat)) - min(rbind(mean_rho_per_dist_bin, boot_mat))
ymax=max(rbind(mean_rho_per_dist_bin, boot_mat)) + .05*yrange
ymin=min(rbind(mean_rho_per_dist_bin, boot_mat)) - .05*yrange

# Initialize plot
par(mar=c(5,4,2,0.5))
plot(1:length(mean_rho_per_dist_bin), scale*mean_rho_per_dist_bin, type="n", ylim=scale*c(ymin, ymax), xlab="", ylab="")
title(xlab="Distance to nearest\nCpG island (kb)", line=3.5)
title(ylab=myylab, line=2.5)

# Add ribbon showing range of bootstrap min and max rho values per bin
ribbon_col="#c8c8c8"
ribbon_xs=c(1:dim(boot_mat)[2], dim(boot_mat)[2]:1)
ribbon_ys=c(apply(boot_mat, 2, function(x) max(x, na.rm = TRUE)), rev(apply(boot_mat, 2, function(x) min(x, na.rm = TRUE))))
polygon(ribbon_xs, scale*ribbon_ys, col=ribbon_col, border=NA)

# Add the line for the actual mean rho per distance bin
lines(1:length(mean_rho_per_dist_bin), scale*mean_rho_per_dist_bin, lwd=2)
