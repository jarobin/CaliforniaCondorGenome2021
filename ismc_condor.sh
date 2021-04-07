### iSMC recombination rate analysis 

### Use filtered VCF file (sites failing filters excluded, missing genotypes excluded, only scaffolds 1-4, 6-32 included)
VCF=CYW1141_simple_PASS_nomissing_scaffold1-32_autos.vcf.gz

### Generate .tab file
# Format is described here, from https://github.com/gvbarroso/iSMC/blob/master/doc/ismc_opt.bpp:

# // The 1st column is the block ID (eg chr1)
# // The 2nd is the start coordinate of the block mapped to the reference genome
# // The 3rd is the end coordinate of the block mapped to the reference genome
# // NOTE: start and end coordinates are relative to the block, meaning it restarts
# // e.g. at every chromosome. 
# // The 4th is 0 for all blocks and the 5th is the difference betwee 3rd & 2nd columns
# // The 6th and 7th columns are bottom and top cut-offs that can be convenient to 'synchronise'
# // coordinates when using ismc_mapper. This can happen when you want to consider a range of sites that
# // matches that of another file. For example, if your sequence data for chromosome 1 starts at position
# // 5000 and goes until position 250000000, and you want to 'synchronise' it with an experimental
# // genetic map that extends from 10000 to 230000000, the 6th and 7th columns of the first line of
# // your tab_file should be  9999 and 229999999.
# // Otherwise, if you want to include all sites, they should be the same as the 2nd and 3rd columns.

### Tab format uses start and end positions of chromosomes, based on sites that are in the VCF file

zcat ${VCF} | grep -v "^#" | \
awk '{ if (NR==1) {c=$1; s=$2; e=$2} else {if ($1!=c) {printf "%s\t%s\t%s\t0\t%s\t%s\t%s\n", c,s,e,e-s,s,e; c=$1; s=$2; e=$2} else {e=$2} } }END{printf "%s\t%s\t%s\t0\t%s\t%s\t%s\n", c,s,e,e-s,s,e}' \
> ${VCF}_ismc.tab


### Generate iSMC options file
# Most options left at default values (except for filenames, etc.)

cat > CYW1141_ismc.bpp <<EOF
dataset_label = CYW1141_ismc
input_file_type = cVCF 
sequence_file_path = CYW1141_simple_PASS_nomissing_scaffold1-32_autos.vcf.gz
seq_compression_type = gzip 
tab_file_path = CYW1141_simple_PASS_nomissing_scaffold1-32_autos.vcf.gz_ismc.tab
number_threads = 31
optimize = true
decode = true
decode_breakpoints_parallel = false
number_rho_categories = 5 
number_intervals = 40
function_tolerance = 1e-6
fragment_size = 2000000
EOF


### Run iSMC (~1 week)

ISMC_DIR=~/project/programs/ismc/src
PAR=CYW1141_ismc.bpp
LOG=CYW1141_ismc.log
/usr/bin/time -v  ${ISMC_DIR}/ismc params=${PAR} |& tee ${LOG}


### Generate iSMC mapping options file

cat > CYW1141_ismc.map.bpp <<EOF
dataset_label = CYW1141_ismc
tab_file_path = CYW1141_simple_PASS_nomissing_scaffold1-32_autos.vcf.gz_ismc.tab
bin_sizes = 1000,1000000 
bin_rate = rho 
tmrca = true 
EOF


### Run iSMC mapper (~5 min)

ISMC_DIR=~/project/programs/ismc/src
PAR=CYW1141_ismc.map.bpp
LOG=CYW1141_ismc.map.log
/usr/bin/time -v  ${ISMC_DIR}/ismc_mapper params=${PAR} |& tee ${LOG}


### Rename chromosomes in iSMC output (iSMC gives sequential numeric chromosome names)

for DATA in CYW1141_ismc.rho.*.bedgraph ; do 
grep -v "^chrom" ${DATA} | sed 's/^chr//g' \
| awk '{if ($1>4) {printf "HiC_scaffold_%s\t%s\t%s\t%s\n", $1+1, $2, $3, $4} else {printf "HiC_scaffold_%s\t%s\t%s\t%s\n", $1, $2, $3, $4}}' > ${DATA}_scaffnames.bed
sort -k1,1 -k2,2n ${DATA}_scaffnames.bed > ${DATA}_scaffnames_sorted.bed
done


################################################################################
### Plot in R

### Recombination rate across the genome

# Read in data and rename chromosomes
df=read.table("CYW1141_ismc.rho.1Mb.bedgraph_scaffnames_sorted.bed", header=F)
colnames(df)=c("chromosome","start","end","rho")
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
ymin=0
ymax=2.25

par(mar=c(3,4,1.5,0.5))
plot(0, 0, type="n", xlim=c(0, nrow(temp)+20), ylim=c(ymin, ymax), xlab="", ylab="", frame.plot=F, axes=F, xaxs="i")
title(ylab=expression(paste(rho, "/bp x 10"^-3, sep="")), line=2.5)
axis(side=2)
axis(side=1, at=pos+10, labels=F)
axis(side=1, at=numpos+10, tick=F, labels=mylabels, las=3, line=-.25, cex.axis=1)
title(xlab="Chromosome", line=1.75)
# Add gray stripes to background
for (i in seq(ymin, ymax, by=1)){rect(-1000, i, nrow(temp)+1000, i+.5, col=bgcol, border=NA)}

# Add chromosome lines
chrs=as.numeric(unique(temp$chromosome))
# Add padding
xx=10
mycol=colset[1]
for (i in 1:length(chrs)){
  tt=temp[which(temp$chromosome==chrs[i]),]
  lines(1:nrow(tt)+xx, 1e3*tt$rho, col=mycol)
  xx=xx+nrow(tt)
  mycol=colset[which(colset!=mycol)]
}

box()


################################################################################
### Boxplot of recombination rate in distal versus proximal regions
# Note: Uses same input data as above (recombination rate across the genome)

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

# Test if rho is higher in distal regions
wt=wilcox.test(ends$rho, mids$rho, alternative="greater")
wt
# 	Wilcoxon rank sum test with continuity correction
# 
# data:  ends$rho and mids$rho
# W = 272980, p-value < 2.2e-16
# alternative hypothesis: true location shift is greater than 0

# Plot
ymin=0
ymax=2.25
yrange=ymax-ymin

par(mar=c(3,2.5,1.5,0.5))
b=boxplot(1e3*ends$rho, 1e3*mids$rho, ylim=c(ymin, ymax), boxwex=0.25, names=c("Distal","Proximal"), xlab="")

# Add significance of wt
arrows(1, ymax-(.075*yrange), 2, ymax-(.075*yrange), angle=90, code=3, length=.05)
text(1.5, ymax-(.025*yrange), starfn(wt$p.value))
