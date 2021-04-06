### MSMC/MSMC2

### Generate VCF files and masks for use with MSMC with vcfAllSiteParser.py

MSMCUTILS=~/project/programs/msmc/msmc-tools-master
REPORTDIR=~/project/reports/condor/msmc/vcfAllSiteParser
for NAME in BGI_N323 CRW1112 CYW1141 VulGry1 ; do
qsub -t 1-542 -N vcfAllSiteParser -V -cwd -l h_rt=24:00:00,mem_free=4G -o ${REPORTDIR} -e ${REPORTDIR} <<EOF
VCF=${NAME}_simple_PASS_nomissing.vcf.gz
CHR=HiC_scaffold_\${SGE_TASK_ID}
MASK_OUT=${NAME}_simple_PASS_nomissing_\${CHR}_sitescovered.bed.gz
VCF_OUT=${NAME}_simple_PASS_nomissing_\${CHR}_msmc.vcf.gz
cd \${TMPDIR}
mkdir mydir_\${SGE_TASK_ID} ; cd mydir_\${SGE_TASK_ID}
cp /wynton/group/wall/jacqueline/condor/condorRef/vcfs/Mask/${NAME}/\${VCF} .
cp /wynton/group/wall/jacqueline/condor/condorRef/vcfs/Mask/${NAME}/\${VCF}.tbi .
tabix -h \${VCF} \${CHR} | ${MSMCUTILS}/vcfAllSiteParser.py \${CHR} \${MASK_OUT} | bgzip > \${VCF_OUT}
tabix -p vcf \${VCF_OUT}
cp \${MASK_OUT} /wynton/scratch/robinsonj/condor/msmc
cp \${VCF_OUT} /wynton/scratch/robinsonj/condor/msmc
cp \${VCF_OUT}.tbi /wynton/scratch/robinsonj/condor/msmc
cd ..
rm -rf mydir_\${SGE_TASK_ID}
EOF
done


### Generate MSMC input files with generate_multihetsep.py

MSMCUTILS=~/project/programs/msmc/msmc-tools-master
REPORTDIR=~/project/reports/condor/msmc/generate_multihetsep
cd ~/myscratch/condor/msmc
for NAME in BGI_N323 CRW1112 CYW1141 VulGry1 ; do
qsub -t 1-542 -N generate_multihetsep -V -cwd -l h_rt=24:00:00,mem_free=4G -o ${REPORTDIR} -e ${REPORTDIR} <<EOF
CHR=HiC_scaffold_\${SGE_TASK_ID}
MASK=${NAME}_simple_PASS_nomissing_\${CHR}_sitescovered.bed.gz
VCF=${NAME}_simple_PASS_nomissing_\${CHR}_msmc.vcf.gz
${MSMCUTILS}/generate_multihetsep.py --mask=\${MASK} \${VCF} > \${VCF%.vcf.gz}.txt
EOF
done


### Run MSMC
# Note: Using the same patterns as for PSMC

MSMC=/wynton/home/walllab/robinsonj/project/programs/msmc/msmc-1.1.0/build/msmc
NTHREADS=8
REPORTDIR=~/project/reports/condor/msmc
cd ~/mygroup/condor/condorRef/msmc
OUTDIR=$(date +%Y%m%d)
mkdir -p results/${OUTDIR}
for NAME in BGI_N323 CRW1112 CYW1141 ; do
qsub -N msmc -V -cwd -l h_rt=24:00:00,mem_free=4G -pe smp ${NTHREADS} -o ${REPORTDIR} -e ${REPORTDIR} <<EOF
${MSMC} -t ${NTHREADS} -i 25 -p "1*4+25*2+1*4+1*6" -o ./results/${OUTDIR}/${NAME}_msmc_subset.out ./inputs/${NAME}_simple_PASS_nomissing_HiC_scaffold_{1..4}_msmc.txt ./inputs/${NAME}_simple_PASS_nomissing_HiC_scaffold_{6..32}_msmc.txt
EOF
done
for NAME in VulGry1 ; do
qsub -N msmc -V -cwd -l h_rt=24:00:00,mem_free=4G -pe smp ${NTHREADS} -o ${REPORTDIR} -e ${REPORTDIR} <<EOF
${MSMC} -t ${NTHREADS} -i 25 -p "1*4+20*2+1*6" -o ./results/${OUTDIR}/${NAME}_msmc_subset.out ./inputs/${NAME}_simple_PASS_nomissing_HiC_scaffold_{1..4}_msmc.txt ./inputs/${NAME}_simple_PASS_nomissing_HiC_scaffold_{6..32}_msmc.txt
EOF
done


### Run MSMC2
# Note: Using the same patterns as for PSMC

MSMC2=/wynton/home/walllab/robinsonj/project/programs/msmc2/msmc2_linux64bit
NTHREADS=8
REPORTDIR=~/project/reports/condor/msmc
cd ~/mygroup/condor/condorRef/msmc
OUTDIR=$(date +%Y%m%d)
mkdir -p results/${OUTDIR}
for NAME in BGI_N323 CRW1112 CYW1141 ; do
qsub -N msmc2 -V -cwd -l h_rt=24:00:00,mem_free=4G -pe smp ${NTHREADS} -o ${REPORTDIR} -e ${REPORTDIR} <<EOF
${MSMC2} -t ${NTHREADS} -i 25 -p "1*4+25*2+1*4+1*6" -o ./results/${OUTDIR}/${NAME}_msmc_subset.out ./inputs/${NAME}_simple_PASS_nomissing_HiC_scaffold_{1..4}_msmc.txt ./inputs/${NAME}_simple_PASS_nomissing_HiC_scaffold_{6..32}_msmc.txt
EOF
done
for NAME in VulGry1 ; do
qsub -N msmc2 -V -cwd -l h_rt=24:00:00,mem_free=4G -pe smp ${NTHREADS} -o ${REPORTDIR} -e ${REPORTDIR} <<EOF
${MSMC2} -t ${NTHREADS} -i 25 -p "1*4+20*2+1*6" -o ./results/${OUTDIR}/${NAME}_msmc_subset.out ./inputs/${NAME}_simple_PASS_nomissing_HiC_scaffold_{1..4}_msmc.txt ./inputs/${NAME}_simple_PASS_nomissing_HiC_scaffold_{6..32}_msmc.txt
EOF
done

