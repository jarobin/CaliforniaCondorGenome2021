### PSMC

### Generate complement masks for making consensus sequence for PSMC
# Note: Generates the complement of masks created while making input files for MSMC

GENOME=~/project/condor/reference/gc_PacBio_HiC/gc_PacBio_HiC.sizes
for NAME in BGI_N323 CRW1112 CYW1141 VulGry1 ; do
zcat ${NAME}_simple_PASS_nomissing_HiC_scaffold_{1..542}_sitescovered.bed.gz > ${NAME}_simple_PASS_nomissing_sitescovered.bed ; done
for NAME in BGI_N323 CRW1112 CYW1141 VulGry1 ; do
bedtools complement -i ${NAME}_simple_PASS_nomissing_sitescovered.bed -g ${GENOME} > ${NAME}_simple_PASS_nomissing_sitescovered_comp.bed
done


### Generate consensus sequences for PSMC

REF=~/project/condor/reference/gc_PacBio_HiC/gc_PacBio_HiC.fasta
for NAME in BGI_N323 CRW1112 CYW1141 VulGry1 ; do
MASK=${NAME}_simple_PASS_nomissing_sitescovered_comp.bed
VCF=${NAME}_simple_PASS_nomissing.vcf.gz
OUT=${NAME}_simple_PASS_nomissing.fa.gz
qsub -N bcftools_consensus -V -cwd -l h_rt=24:00:00,mem_free=4G <<EOF
bcftools consensus -f ${REF} -I -M N -m ${MASK} -s ${NAME} ${VCF} | gzip > ${OUT}
EOF
done


### Generate PSMC input files

PSMCDIR=~/project/programs/psmc/
for i in *.fa.gz ; do ${PSMCDIR}/utils/fq2psmcfa ${i} > ${i%.fa.gz}.psmcfa ; done

# Index PSMC input files
for i in *nomissing.psmcfa ; do samtools faidx ${i} ; done

# Retrieve scaffolds for analysis
for i in *nomissing.psmcfa; do for c in {1..4} {6..32} ; do 
samtools faidx ${i} HiC_scaffold_${c} > ${i%.psmcfa}_HiC_scaffold_${c}.psmcfa
done ; done

# Remove files with no data
wc -l *HiC_scaffold*psmcfa | sed -e 's/^\ *//g' | awk '$1==1 {print $2}' | xargs rm

# Concatenate PSMC input files
for i in BGI_N323 CRW1112 CYW1141 VulGry1 ; do
cat ${i}_simple_PASS_nomissing_HiC_scaffold_{1..32}.psmcfa > ${i}_simple_PASS_nomissing_subset.psmcfa
done
mkdir -p temp
mv *_HiC_scaffold_*.psmcfa temp


### Run PSMC on scaffolds 1-4, 6-32 ("subset")
PSMCDIR=~/project/programs/psmc/
for i in *subset.psmcfa ; do 
qsub -N psmc -V -cwd -l h_rt=24:00:00,mem_free=4G <<EOF
${PSMCDIR}/psmc -N25 -t15 -r5 -p "4+25*2+4+6" -o ${i%.psmcfa}.out.psmc ${i}
EOF
done

# Check likelihood scores to see plateau
for i in *psmc ; do 
grep "^LK" ${i} | tail -n +2 | cut -f2 | \
gnuplot -p -e "set terminal dumb size 120, 30; set autoscale; plot '-' using 1 with lines notitle"
done

# Check for intervals with <10 recombination events
for i in *psmc ; do 
echo ${i}
tail -n +$(grep -P -n "^RD\t21" ${i} | cut -d':' -f1) ${i} | \
grep "^RS" | awk '$5<10'
done

# Note: Andean condor had intervals with too few recombinations, ended up using -p "4+20*2+6"


### Bootstrap

PSMCDIR=~/project/programs/psmc/
for i in {B,C}*subset.psmcfa ; do
${PSMCDIR}/utils/splitfa ${i} > ${i%.psmcfa}_split.psmcfa
qsub -N psmcboot -t 1-100 -V -cwd -l h_rt=24:00:00,mem_free=4G <<EOF
${PSMCDIR}/psmc -b -N25 -t15 -r5 -p "4+25*2+4+6" -o ${i%.psmcfa}_split_boot\${SGE_TASK_ID}.out.psmc ${i%.psmcfa}_split.psmcfa
EOF
done
for i in V*subset.psmcfa ; do
${PSMCDIR}/utils/splitfa ${i} > ${i%.psmcfa}_split.psmcfa
qsub -N psmcboot -t 1-100 -V -cwd -l h_rt=24:00:00,mem_free=4G <<EOF
${PSMCDIR}/psmc -b -N25 -t15 -r5 -p "4+20*2+6" -o ${i%.psmcfa}_split_boot\${SGE_TASK_ID}.out.psmc ${i%.psmcfa}_split.psmcfa
EOF
done

cat BGI_N323_simple_PASS_nomissing_subset.out.psmc BGI_N323_*_boot* > BGI_N323_simple_PASS_nomissing_subset_split_boot.combined.out.psmc
cat CRW1112_simple_PASS_nomissing_subset.out.psmc CRW1112_*_boot* > CRW1112_simple_PASS_nomissing_subset_split_boot.combined.out.psmc
cat CYW1141_simple_PASS_nomissing_subset.out.psmc CYW1141_*_boot* > CYW1141_simple_PASS_nomissing_subset_split_boot.combined.out.psmc
cat VulGry1_simple_PASS_nomissing_subset.out8.psmc VulGry1_*_boot* > VulGry1_simple_PASS_nomissing_subset_split_boot.combined.out.psmc


### Generate plot files
# Use -R to keep the files generated during plotting - the .plot.*.txt files contain the rescaled coordinates
# M is the per-year mutation rate, MU is the per-generation mutation rate
# Main results use the mutation rate of 1.4e-8 per generation and generation time of 10 years

PSMCDIR=~/project/programs/psmc/

for M in 1.4e-9 0.65e-9 ; do
for G in 5 10 15 ; do
MU=$(echo ${M} | awk -v var=${G} '{print $0*var}')
mkdir plot_mu${MU}_g${G}
for INPUT in *subset.out*.psmc ; do
PLOT=plot_mu${MU}_g${G}/${INPUT}_mu${MU}_g${G}.plot
perl ${PSMCDIR}/utils/psmc_plot.pl -R -p -u ${MU} -g ${G} ${PLOT} ${INPUT}
done ; done ; done


### PSMC with high, low diversity regions of condor (distal and proximal regions)

SIZES=gc_PacBio_HiC.sizes
LIM=10000000

# Get end coordinates (distal)
ENDS=${SIZES}_ends_10Mb
head -n 21 ${SIZES} | grep -v "HiC_scaffold_5" | awk -v var=${LIM} '{printf "%s:0-%s\n%s:%s-%s\n", $1,var,$1,$2-var,$2}' > ${ENDS}
head -n 32 ${SIZES} | tail -n 11 | awk '{printf "%s:0-%s\n", $1,$2}' >> ${ENDS}

# Get mid coordinates (proximal)
MIDS=${SIZES}_mids_10Mb
head -n 21 ${SIZES} | grep -v "HiC_scaffold_5" | awk -v var=${LIM} '{printf "%s:%s-%s\n", $1,var+1,$2-var-1}' > ${MIDS}

# Index fastas
for i in *_simple_PASS_nomissing.fa.gz ; do samtools faidx ${i} & done

# Extract ends and mids
for i in *_simple_PASS_nomissing.fa.gz ; do
samtools faidx -r ${ENDS} ${i} | bgzip > ${i%.fa.gz}_ends.fa.gz
samtools faidx -r ${MIDS} ${i} | bgzip > ${i%.fa.gz}_mids.fa.gz
done

# Generate PSMC input files
PSMCDIR=~/project/programs/psmc/
for i in *ends.fa.gz ; do ${PSMCDIR}/utils/fq2psmcfa ${i} > ${i%.fa.gz}.psmcfa ; done
for i in *mids.fa.gz ; do ${PSMCDIR}/utils/fq2psmcfa ${i} > ${i%.fa.gz}.psmcfa ; done

# Run PSMC (use same parameters as for full genome)
PSMCDIR=~/project/programs/psmc/
for i in {B,C}*.psmcfa ; do 
qsub -N psmc -V -cwd -l h_rt=24:00:00,mem_free=4G <<EOF
${PSMCDIR}/psmc -N25 -t15 -r5 -p "4+25*2+4+6" -o ${i%.psmcfa}.out.psmc ${i}
EOF
done
for i in V*.psmcfa ; do 
qsub -N psmc -V -cwd -l h_rt=24:00:00,mem_free=4G <<EOF
${PSMCDIR}/psmc -N25 -t15 -r5 -p "4+20*2+6" -o ${i%.psmcfa}.out.psmc ${i}
EOF
done

# Check likelihood scores to see plateau
for i in *psmc ; do 
grep "^LK" ${i} | tail -n +2 | cut -f2 | \
gnuplot -p -e "set terminal dumb size 120, 30; set autoscale; plot '-' using 1 with lines notitle"
done

# Check for intervals with <10 recombination events
for i in *psmc ; do 
echo ${i}
tail -n +$(grep -P -n "^RD\t21" ${i} | cut -d':' -f1) ${i} | \
grep "^RS" | awk '$5<10'
done


### Generate plot files with mu=1.4e-8, g=10
# Use -R to keep the files generated during the plotting - the .plot.*.txt files contain the rescaled coordinates

MU=1.4e-8
G=10

PSMCDIR=~/project/programs/psmc/
for INPUT in *.out.psmc ; do
PLOT=${INPUT}_mu${MU}_g${G}.plot
perl ${PSMCDIR}/utils/psmc_plot.pl -R -p -u ${MU} -g ${G} ${PLOT} ${INPUT}
done


### Make final plots in R

setwd("~/condor/analysis/psmc")

pdf("psmc_main_5panels_mu1.4e-08_g10.pdf", width=4.488, height=5.25, pointsize=8)

split.screen(rbind(
c(0.0, 1.0, 0.5, 1.0),
c(0.0, 0.5, 0.25, 0.5),
c(0.5, 1.0, 0.25, 0.5),
c(0.0, 0.5, 0.0, 0.25),
c(0.5, 1.0, 0.0, 0.25)))


# A: Full genome results with bootstrap replicates
setwd("main_with_boot/plot_mu1.4e-08_g10")
allfiles=list.files(pattern="txt")

xmin=1e3
xmax=4e6
ymin=0
ymax=5.5

# Turkey vulture color
c1=c(0, 192, 0)/255
# CA condor 1 color
c2=c(171, 0, 0)/255
# CA condor 2 color
c3=c(240, 105, 20)/255
# Andean condor color
c4=c(0, 128, 255)/255

# Main colors (no transparency)
mycols1=c(
rgb(c1[1], c1[2], c1[3], alpha=1),
rgb(c2[1], c2[2], c2[3], alpha=1),
rgb(c3[1], c3[2], c3[3], alpha=1),
rgb(c4[1], c4[2], c4[3], alpha=1)
)

# Bootstrap replicate colors (with transparency)
transp=.075
mycols2=c(
rgb(c1[1], c1[2], c1[3], alpha=transp),
rgb(c2[1], c2[2], c2[3], alpha=transp),
rgb(c3[1], c3[2], c3[3], alpha=transp),
rgb(c4[1], c4[2], c4[3], alpha=transp)
)

# Plot
screen(1)
par(mar=c(3.5,3.75,0.5,0.5))
plot(1, 1, type="n", log="x", axes=F, xlim=c(xmin, xmax), ylim=c(ymin, ymax), xlab="", ylab="")

title(xlab="Years before present", line=2)
title(ylab=expression("Effective population size (x10"^4*")"), line=2.25)
	
axis(side=2, line=0, labels=F)
axis(side=2, line=-.25, labels=T, tick=F)
	
at.x=outer(1:9, 10^(3:8))
lab.x=NULL
for (i in 1:length(at.x)){
	p=log10(at.x[i])
	if (p %% 1 == 0) {lab.x[i]=as.expression(bquote(10^ .(p)))}
	else {lab.x[i]=""}
}
axis(1, at=at.x, labels=lab.x, las=1)
	
legend("topleft", lwd=2, col=mycols1[c(2, 3, 4, 1)], legend=c("CA condor 1", "CA condor 2", "Andean condor", "Turkey vulture"), bty="n")
	
box()

# Add lines for each sample
psmc_plot_fill=function(){
	dfiles=allfiles[grep(pattern=samplename, allfiles)]
	aa=read.table(dfiles[1])
	for (i in 2:length(dfiles)){
		aa=read.table(dfiles[i])
		lines(aa$V1, aa$V2, type="s", col=mycols2[nn], lwd=1)
	}
	lines(aa$V1, aa$V2, type="s", col=mycols1[nn], lwd=2)
}

samplename="BGI_N323"
nn=1
psmc_plot_fill()

samplename="CYW1141"
nn=2
psmc_plot_fill()

samplename="CRW1112"
nn=3
psmc_plot_fill()

samplename="VulGry1"
nn=4
psmc_plot_fill()


# B: Ends & mids (note: x- and y-axis labels added in post with image editor)

setwd("~/condor/analysis/psmc")
endfiles=list.files(path="ends_mids/plot_mu1.4e-8_g10", pattern="plot.0.txt", full.names=T)[grep(pattern="ends", list.files(path="ends_mids/plot_mu1.4e-8_g10", pattern="plot.0.txt"))]
midfiles=list.files(path="ends_mids/plot_mu1.4e-8_g10", pattern="plot.0.txt", full.names=T)[grep(pattern="mids", list.files(path="ends_mids/plot_mu1.4e-8_g10", pattern="plot.0.txt"))]
allfiles=list.files(path="main_with_boot/plot_mu1.4e-08_g10", pattern="plot.0.txt", full.names=T)

xmin=1e3
xmax=4e6
ymin=0
ymax=5.5

# Turkey vulture color
c1=c(0, 192, 0)/255
# CA condor 1 color
c2=c(171, 0, 0)/255
# CA condor 2 color
c3=c(240, 105, 20)/255
# Andean condor color
c4=c(0, 128, 255)/255

mycols1=c(
rgb(c1[1], c1[2], c1[3], alpha=1),
rgb(c2[1], c2[2], c2[3], alpha=1),
rgb(c3[1], c3[2], c3[3], alpha=1),
rgb(c4[1], c4[2], c4[3], alpha=1)
)

psmc_base_plot=function(){
	par(mar=c(2,3.75,0.5,0.5))
	plot(1, 1, type="n", log="x", axes=F, xlim=c(xmin, xmax), ylim=c(ymin, ymax), xlab="", ylab="")	
	axis(side=2, line=0, labels=F)
	axis(side=2, line=-.25, labels=T, tick=F)
	at.x=outer(1:9, 10^(3:8))
	lab.x=NULL
	for (i in 1:length(at.x)){
		p=log10(at.x[i])
		if (p %% 1 == 0) {lab.x[i]=as.expression(bquote(10^ .(p)))}
		else {lab.x[i]=""}
	}
	axis(1, at=at.x, labels=lab.x, las=1)		
	box()
}

psmc_plot_fill=function(){
	aa=read.table(endfiles[grep(pattern=samplename, endfiles)])
	bb=read.table(midfiles[grep(pattern=samplename, midfiles)])
	cc=read.table(allfiles[grep(pattern=samplename, allfiles)])
	lines(aa$V1, aa$V2, type="s", col=mycols1[nn], lty=1)
	lines(bb$V1, bb$V2, type="s", col=mycols1[nn], lty=2)
	lines(cc$V1, cc$V2, type="s", col="black", lty=1)
	text(x=6e3, y=5, plottitle)
}

screen(2)
samplename="CYW1141"
plottitle="CA condor 1"
nn=2
psmc_base_plot()
psmc_plot_fill()

screen(3)
samplename="CRW1112"
plottitle="CA condor 2"
nn=3
psmc_base_plot()
psmc_plot_fill()

screen(4)
samplename="VulGry1"
plottitle="Andean condor"
nn=4
psmc_base_plot()
psmc_plot_fill()

screen(5)
samplename="BGI_N323"
plottitle="Turkey vulture"
nn=1
psmc_base_plot()
psmc_plot_fill()

close.screen(all = TRUE)
dev.off()


### Plots with varying mu, g (in R)

# Turkey vulture color
c1=c(0, 192, 0)/255
# CA condor 1 color
c2=c(171, 0, 0)/255
# CA condor 2 color
c3=c(240, 105, 20)/255
# Andean condor color
c4=c(0, 128, 255)/255

mycols=c(
rgb(c1[1], c1[2], c1[3], alpha=1),
rgb(c2[1], c2[2], c2[3], alpha=1),
rgb(c3[1], c3[2], c3[3], alpha=1),
rgb(c4[1], c4[2], c4[3], alpha=1)
)

psmc_base_plot=function(){
	par(mar=c(5,4,1.5,1))
	plot(1, 1, type="n", log="x", axes=F, xlim=c(xmin, xmax), ylim=c(ymin, ymax), xlab="", ylab="")	
	axis(side=2)	
	at.x=outer(1:9, 10^(3:8))
	lab.x=NULL
	for (i in 1:length(at.x)){
		p=log10(at.x[i])
		if (p %% 1 == 0) {lab.x[i]=as.expression(bquote(10^ .(p)))}
		else {lab.x[i]=""}
	}
	axis(1, at=at.x, labels=lab.x, las=1)		
	if(labx){title(xlab="Years before present", line=3.75)}
	title(xlab=myxlab_sub, line=2.5)
	if(laby){title(ylab=expression("Effective population size (x10"^4*")"), line=2.5)}
	box()
}

psmc_plot_fill=function(){
	allfiles=list.files(pattern="txt")
	samplename="BGI_N323"
	nn=1
	dfiles=allfiles[grep(pattern=samplename, allfiles)]
	aa=read.table(dfiles[1])
	lines(aa$V1, aa$V2, type="s", col=mycols[nn])
	samplename="CYW1141"
	nn=2
	dfiles=allfiles[grep(pattern=samplename, allfiles)]
	aa=read.table(dfiles[1])
	lines(aa$V1, aa$V2, type="s", col=mycols[nn])
	samplename="CRW1112"
	nn=3
	dfiles=allfiles[grep(pattern=samplename, allfiles)]
	aa=read.table(dfiles[1])
	lines(aa$V1, aa$V2, type="s", col=mycols[nn])
	samplename="VulGry1"
	nn=4
	dfiles=allfiles[grep(pattern=samplename, allfiles)]
	aa=read.table(dfiles[1])
	lines(aa$V1, aa$V2, type="s", col=mycols[nn])
}

setwd("~/condor/analysis/psmc/variable_g/")
pdf("psmc_6panel_variable_mu_g.pdf", width=6, height=5, pointsize=8)
par(mfrow=c(3,2))

xmin=1e3
xmax=1e7
ymin=0

setwd("~/condor/analysis/psmc/variable_g/")
setwd("plot_mu3.25e-09_g5")
myxlab_sub=expression(paste(mu, "=3.25x10"^-9, "/generation"))
ymax=25
labx=F
laby=T
psmc_base_plot()
psmc_plot_fill()
title(main=expression(paste("M=6.5x10"^-10, "/year")), font=1)
legend("topleft", col=mycols[c(2, 3, 4, 1)], legend=c("CA condor 1", "CA condor 2", "Andean condor", "Turkey vulture"), bty="n")

setwd("~/condor/analysis/psmc/variable_g/")
setwd("plot_mu7e-09_g5")
myxlab_sub=expression(paste(mu, "=7.0x10"^-9, "/generation"))
ymax=25
labx=F
laby=F
psmc_base_plot()
psmc_plot_fill()
title(main=expression(paste("M=1.4x10"^-9, "/year")), font=1)

setwd("~/condor/analysis/psmc/variable_g/")
setwd("plot_mu6.5e-09_g10")
myxlab_sub=expression(paste(mu, "=6.5x10"^-9, "/generation"))
ymax=12
labx=F
laby=T
psmc_base_plot()
psmc_plot_fill()

setwd("~/condor/analysis/psmc/variable_g/")
setwd("plot_mu1.4e-08_g10")
myxlab_sub=expression(paste(mu, "=1.4x10"^-8, "/generation"))
ymax=12
labx=F
laby=F
psmc_base_plot()
psmc_plot_fill()

setwd("~/condor/analysis/psmc/variable_g/")
setwd("plot_mu9.75e-09_g15")
myxlab_sub=expression(paste(mu, "=9.75x10"^-9, "/generation"))
ymax=8
labx=T
laby=T
psmc_base_plot()
psmc_plot_fill()

setwd("~/condor/analysis/psmc/variable_g/")
setwd("plot_mu2.1e-08_g15")
myxlab_sub=expression(paste(mu, "=2.1x10"^-8, "/generation"))
ymax=8
labx=T
laby=F
psmc_base_plot()
psmc_plot_fill()

dev.off()


