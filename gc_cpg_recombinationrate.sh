# GC content + CpG island identification and recombination rate analysis (Fig. 3C-F)

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


### Get het. data for same window coordinates as in iSMC for correlogram

# Get coordinates of iSMC 1 Mb windows
BED=CYW1141_ismc.rho.1Mb.bedgraph_scaffnames_sorted.bed
awk '{printf "%s\t%s\t%s\n", $1, $2, $3-1}' ${BED} > condor_ismc_coords_1Mb.list

# Get numbers of hets, calls in coordinates from list
SCRIPT=WindowHetCoords.py
for NAME in BGI_N323 CRW1112 CYW1141 VulGry1 ; do
for i in {1..32} ; do
VCF=${NAME}*_$(printf %03d ${i})_*_Filter.vcf.gz
awk -v var=HiC_scaffold_${i} '$1==var' condor_ismc_coords_1Mb.list > tempcoords_${NAME}_${i}
OUT=${VCF%.vcf.gz}_het_ismc_coords_1Mb.txt
python ${SCRIPT} ${VCF} tempcoords_${NAME}_${i} ${OUT}
rm tempcoords_${NAME}_${i}
done & done

# Concatenate to make final file (note: scaffolds renamed because iSMC assigns chromosome names sequentially)
for NAME in BGI_N323 CRW1112 CYW1141 VulGry1 ; do
echo -e "chrom\twindow_start\tsites_total\thet_${NAME}" > ${NAME}_het_ismc_coords_1Mb.txt
for i in {001..032} ; do
tail -n +2 ${NAME}*_${i}_*_het_ismc_coords_1Mb.txt \
| awk '{if ($4 > 0) {printf "%s\t%s\t%s\t%s\n", $1, $2, $3, $5/$4} else {printf "%s\t%s\t%s\t%s\n", $1, $2, $3, "NA"}}' >> ${NAME}_het_ismc_coords_1Mb.txt
done ; done


################################################################################
### Plotting in R

### Plot GC content across the genome (Fig. 3C)

# Read in data files and rename chromosomes (from scaffold to proper chromosome names, also exclude chrZ)
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
### Boxplot of GC content in distal versus proximal regions (Fig. 3D)
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
### Recombination rate by distance from CpG islands plus bootstrap (Fig. 3E)

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


################################################################################
### Correlogram (Fig. 3F)
# Adapted from https://stackoverflow.com/questions/19012529/correlation-corrplot-configuration

library(corrgram)

# Pearson correlations
panel.shadeNtext <- function (x, y, corr = NULL, col.regions, ...) 
{
  corr <- cor(x, y, use = "pair")
  results <- cor.test(x, y, alternative = "two.sided")
  est <- results$p.value
  stars <- ifelse(est < 0.001, "***", 
                  ifelse(est < 0.01, "**", 
                         ifelse(est < 0.05, "*", "")))
  ncol <- 14
  pal <- col.regions(ncol)
  col.ind <- as.numeric(cut(corr, breaks = seq(from = -1, to = 1, 
                                               length = ncol + 1), include.lowest = TRUE))
  usr <- par("usr")
  rect(usr[1], usr[3], usr[2], usr[4], col = pal[col.ind], 
       border = NA)
  #box(col = "lightgray")
  box(col = NULL, lwd=.5)
  on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- formatC(corr, digits = 2, format = "f")
  cex.cor <- .6/strwidth("-X.xx")
  fonts <- ifelse(stars != "", 2,1)
  # option 1: stars:
  text(0.5, 0.4, paste0(r,"\n", stars), cex = cex.cor)
  # option 2: bolding:
  #text(0.5, 0.5, r, cex = cex.cor, font=fonts)
}

# Read in data files and rename chromosomes (from scaffold to proper chromosome names, also exclude chrZ)
aa=read.table("gc_PacBio_HiC_CYW1141_ismc_20200720_a.rho.1Mb_GC.txt", header=T, sep='\t')
aa=aa[with(aa, order(aa$chromosome, as.numeric(aa$start))),]

hetfiles=list.files(pattern="_het_ismc_coords_1Mb.txt")
chrkey=read.table("condor_chr_scaff_length.txt", header=T, sep="\t", stringsAsFactors=F)
df=aa
for (i in 1:length(hetfiles)){
	bb=read.table(hetfiles[i], header=T, sep="\t")
	bb=bb[with(bb, order(bb$chrom, as.numeric(bb$window_start))),]
	newcolnames=c(names(df), names(bb)[4])
	df=cbind(df, bb[,4])
	colnames(df)=newcolnames
}
df=df[which(df$chromosome %in% chrkey$SCAFF),]

# Get distance from chromosome end for each window
pos_pct=NULL
pos_dist=NULL
for (i in 1:length(df$chrom)){
	c=df[i,]$chromosome
	s=df[i,]$start
	l=chrkey[which(chrkey$SCAFF==c),]$LENGTH
	pos_pct[i]=s/l
	pos_dist[i]=ifelse(pos_pct[i] <= 0.5, pos_pct[i], 1-pos_pct[i])
}
df=cbind(df, pos_pct, pos_dist)

# Make final data table
tt=cbind(df$pos_dist, df$gc_proportion, df$rho, df$het_CYW1141, df$het_CRW1112, df$het_VulGry1, df$het_BGI_N323)

# Plot
mylabels=c(
"Distance\nfrom\nchr. end",
"GC (%)",
"rho",
"Het.\nCA condor\n1",
"Het.\nCA condor\n2",
"Het.\nAndean\ncondor",
"Het.\nTurkey\nvulture")
mycols=colorRampPalette(rev(c("red", "salmon", "white", "royalblue", "navy")))

pdf("correlogram_rho_gc_het.pdf", width=4.2, height=2.8, pointsize=8)

par(xpd=T)
corrgram(tt, type="data", lower.panel=NULL, upper.panel=panel.shadeNtext, col.regions=mycols, labels=mylabels, cex.labels=1.225)

dev.off()


### Function to plot color bar (added in post with image editor)
# From https://stackoverflow.com/questions/9314658/colorbar-from-custom-colorramppalette

color.bar <- function(lut, min, max=-min, nticks=11, ticks=seq(min, max, len=nticks), title='') {
    scale = (length(lut)-1)/(max-min)
    #dev.new(width=5, height=2)
    plot(c(min,max), c(0,10), type='n', bty='n', xaxt='n', xlab='', yaxt='n', ylab='', main=title)
    axis(1, ticks, las=1)
    for (i in 1:(length(lut)-1)) {
     x = (i-1)/scale + min
     rect(x,0,x+1/scale,10, col=lut[i], border=NA)
    }
    lines(c(-1,1), c(10,10))
    lines(c(-1,1), c(0,0))
	lines(c(-1,-1), c(0,10))
	lines(c(1,1), c(0,10))
}

mycols=colorRampPalette(c("red", "salmon", "white", "royalblue", "navy"))

pdf("scale.pdf", width=2, height=.5, pointsize=8)
par(mar=c(2,1,1,1))
color.bar(mycols(100), 1)
dev.off()



