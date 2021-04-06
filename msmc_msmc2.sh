### MSMC/MSMC2

### Generate VCF files and masks for use with MSMC with vcfAllSiteParser.py

MSMCUTILS=~/project/programs/msmc/msmc-tools-master
for NAME in BGI_N323 CRW1112 CYW1141 VulGry1 ; do
qsub -t 1-542 -N vcfAllSiteParser -V -cwd -l h_rt=24:00:00,mem_free=4G <<EOF
VCF=${NAME}_simple_PASS_nomissing.vcf.gz
CHR=HiC_scaffold_\${SGE_TASK_ID}
MASK_OUT=${NAME}_simple_PASS_nomissing_\${CHR}_sitescovered.bed.gz
VCF_OUT=${NAME}_simple_PASS_nomissing_\${CHR}_msmc.vcf.gz
tabix -h \${VCF} \${CHR} | ${MSMCUTILS}/vcfAllSiteParser.py \${CHR} \${MASK_OUT} | bgzip > \${VCF_OUT}
tabix -p vcf \${VCF_OUT}
EOF
done


### Generate MSMC input files with generate_multihetsep.py

MSMCUTILS=~/project/programs/msmc/msmc-tools-master
for NAME in BGI_N323 CRW1112 CYW1141 VulGry1 ; do
qsub -t 1-542 -N generate_multihetsep -V -cwd -l h_rt=24:00:00,mem_free=4G <<EOF
CHR=HiC_scaffold_\${SGE_TASK_ID}
MASK=${NAME}_simple_PASS_nomissing_\${CHR}_sitescovered.bed.gz
VCF=${NAME}_simple_PASS_nomissing_\${CHR}_msmc.vcf.gz
${MSMCUTILS}/generate_multihetsep.py --mask=\${MASK} \${VCF} > \${VCF%.vcf.gz}.txt
EOF
done


### Run MSMC
# Note: Same -p patterns as for PSMC
# Only use scaffolds 1-4, 6-32

MSMC=~/project/programs/msmc/msmc-1.1.0/build/msmc
NTHREADS=8
for NAME in BGI_N323 CRW1112 CYW1141 ; do
qsub -N msmc -V -cwd -l h_rt=24:00:00,mem_free=4G -pe smp ${NTHREADS} <<EOF
${MSMC} -t ${NTHREADS} -i 25 -p "1*4+25*2+1*4+1*6" -o ./results/${NAME}_msmc_subset.out ./inputs/${NAME}_simple_PASS_nomissing_HiC_scaffold_{1..4}_msmc.txt ./inputs/${NAME}_simple_PASS_nomissing_HiC_scaffold_{6..32}_msmc.txt
EOF
done
for NAME in VulGry1 ; do
qsub -N msmc -V -cwd -l h_rt=24:00:00,mem_free=4G -pe smp ${NTHREADS} <<EOF
${MSMC} -t ${NTHREADS} -i 25 -p "1*4+20*2+1*6" -o ./results/${NAME}_msmc_subset.out ./inputs/${NAME}_simple_PASS_nomissing_HiC_scaffold_{1..4}_msmc.txt ./inputs/${NAME}_simple_PASS_nomissing_HiC_scaffold_{6..32}_msmc.txt
EOF
done


### Run MSMC2
# Note: Same -p patterns as for PSMC
# Only use scaffolds 1-4, 6-32

MSMC2=/wynton/home/walllab/robinsonj/project/programs/msmc2/msmc2_linux64bit
NTHREADS=8
for NAME in BGI_N323 CRW1112 CYW1141 ; do
qsub -N msmc2 -V -cwd -l h_rt=24:00:00,mem_free=4G -pe smp ${NTHREADS} <<EOF
${MSMC2} -t ${NTHREADS} -i 25 -p "1*4+25*2+1*4+1*6" -o ./results/${NAME}_msmc_subset.out ./inputs/${NAME}_simple_PASS_nomissing_HiC_scaffold_{1..4}_msmc.txt ./inputs/${NAME}_simple_PASS_nomissing_HiC_scaffold_{6..32}_msmc.txt
EOF
done
for NAME in VulGry1 ; do
qsub -N msmc2 -V -cwd -l h_rt=24:00:00,mem_free=4G -pe smp ${NTHREADS} <<EOF
${MSMC2} -t ${NTHREADS} -i 25 -p "1*4+20*2+1*6" -o ./results/${NAME}_msmc_subset.out ./inputs/${NAME}_simple_PASS_nomissing_HiC_scaffold_{1..4}_msmc.txt ./inputs/${NAME}_simple_PASS_nomissing_HiC_scaffold_{6..32}_msmc.txt
EOF
done


### Plot MSMC/MSMC2/PSMC results in R

setwd("~/condor/analysis/msmc")
DIR1="results_msmc"
DIR2="results_msmc2"
DIR3="~/condor/analysis/psmc/main_with_boot/plot_mu1.4e-08_g10"

results1=list.files(path=DIR1, pattern="final.txt")
results2=list.files(path=DIR2, pattern="final.txt")
results3=list.files(path=DIR3, pattern="plot.0.txt")

mu=1.4e-8
xmin=1e3
xmax=1e7
ymin=0
ymax=12

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

msmc_plot=function(sampletitle, mycol){
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

	if(labx){title(xlab="Years before present", line=2.5)}
	if(laby){title(ylab=expression("Effective population size (x 10"^4*")"), line=2.5)}

	xx=gen*((dat1$left_time_boundary+dat1$right_time_boundary)/2)/mu
	yy=(1/dat1$lambda)/(2*mu)
	lines(xx, 1e-4*yy, type="s", col=mycol, lty=1)
	xx=gen*((dat2$left_time_boundary+dat2$right_time_boundary)/2)/mu
	yy=(1/dat2$lambda)/(2*mu)
	lines(xx, 1e-4*yy, type="s", col=mycol, lty=2)
	lines(dat3$V1, dat3$V2, type="s", col="black")

	legend("topleft", legend=sampletitle, bty="n")
	
	box()
}

pdf("msmc_msmc2_plot.pdf", width=6, height=4, pointsize=8)

par(mfrow=c(2,2))
	par(mar=c(5,4,1.5,1))

# Set generation time
gen=10

# Plot
samplename="CYW1141"
sampletitle="CA condor 1"
dat1=read.table(paste(DIR1, results1[grep(samplename, results1)], sep="/"), header=T)
dat2=read.table(paste(DIR2, results2[grep(samplename, results2)], sep="/"), header=T)
dat3=read.table(paste(DIR3, results3[grep(samplename, results3)], sep="/"), header=F)
labx=F
laby=T
msmc_plot(sampletitle, mycols[2])

samplename="CRW1112"
sampletitle="CA condor 2"
dat1=read.table(paste(DIR1, results1[grep(samplename, results1)], sep="/"), header=T)
dat2=read.table(paste(DIR2, results2[grep(samplename, results2)], sep="/"), header=T)
dat3=read.table(paste(DIR3, results3[grep(samplename, results3)], sep="/"), header=F)
labx=F
laby=F
msmc_plot(sampletitle, mycols[3])

samplename="VulGry1"
sampletitle="Andean condor"
dat1=read.table(paste(DIR1, results1[grep(samplename, results1)], sep="/"), header=T)
dat2=read.table(paste(DIR2, results2[grep(samplename, results2)], sep="/"), header=T)
dat3=read.table(paste(DIR3, results3[grep(samplename, results3)], sep="/"), header=F)
labx=T
laby=T
msmc_plot(sampletitle, mycols[4])

samplename="BGI_N323"
sampletitle="Turkey vulture"
dat1=read.table(paste(DIR1, results1[grep(samplename, results1)], sep="/"), header=T)
dat2=read.table(paste(DIR2, results2[grep(samplename, results2)], sep="/"), header=T)
dat3=read.table(paste(DIR3, results3[grep(samplename, results3)], sep="/"), header=F)
labx=T
laby=F
msmc_plot(sampletitle, mycols[1])

dev.off()

