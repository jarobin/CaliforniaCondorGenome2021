### PSMC

### Generate complement masks for making consensus sequence for PSMC
# Note: Generates the complement of masks generated while making input files for MSMC

GENOME=~/project/condor/reference/gc_PacBio_HiC/gc_PacBio_HiC.sizes
cd ~/myscratch/condor/msmc/masks
for NAME in BGI_N323 CRW1112 CYW1141 VulGry1 ; do
zcat ${NAME}_simple_PASS_nomissing_HiC_scaffold_{1..542}_sitescovered.bed.gz > ${NAME}_simple_PASS_nomissing_sitescovered.bed ; done
mv *bed ~/myscratch/condor/psmc
cd ~/myscratch/condor/psmc
for NAME in BGI_N323 CRW1112 CYW1141 VulGry1 ; do
bedtools complement -i ${NAME}_simple_PASS_nomissing_sitescovered.bed -g ${GENOME} > ${NAME}_simple_PASS_nomissing_sitescovered_comp.bed
done


### Generate consensus sequences for PSMC

cd ~/myscratch/condor/psmc
REF=/wynton/home/walllab/robinsonj/project/condor/reference/gc_PacBio_HiC/gc_PacBio_HiC.fasta
REPORTDIR=~/project/reports/condor/psmc
for NAME in BGI_N323 CRW1112 CYW1141 VulGry1 ; do
MASK=${NAME}_simple_PASS_nomissing_sitescovered_comp.bed
VCF=${NAME}_simple_PASS_nomissing.vcf.gz
OUT=${NAME}_simple_PASS_nomissing.fa.gz
qsub -N bcftools_consensus -V -cwd -l h_rt=24:00:00,mem_free=4G -o ${REPORTDIR} -e ${REPORTDIR} <<EOF
bcftools consensus -f ${REF} -I -M N -m ${MASK} -s ${NAME} ${VCF} | gzip > ${OUT}
EOF
done


### Generate PSMC input files (quick)

cd ~/myscratch/condor/psmc
PSMCDIR=~/project/programs/psmc/
for i in *.fa.gz ; do ${PSMCDIR}/utils/fq2psmcfa ${i} > ${i%.fa.gz}.psmcfa ; done

# Index PSMC input files
for i in *nomissing.psmcfa ; do samtools faidx ${i} ; done

# Retrieve scaffolds
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


### Run PSMC on scaffolds 1-4, 6-32
PSMCDIR=~/project/programs/psmc/
REPORTDIR=~/project/reports/condor/psmc
cd ~/myscratch/condor/psmc
for i in *subset.psmcfa ; do 
qsub -N psmc -V -cwd -l h_rt=24:00:00,mem_free=4G -o ${REPORTDIR} -e ${REPORTDIR} <<EOF
${PSMCDIR}/psmc -N25 -t15 -r5 -p "4+25*2+4+6" -o ${i%.psmcfa}.out.psmc ${i}
EOF
done

# Check likelihood scores
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
REPORTDIR=~/project/reports/condor/psmc
for i in {B,C}*.psmcfa ; do 
qsub -N psmc -V -cwd -l h_rt=24:00:00,mem_free=4G -o ${REPORTDIR} -e ${REPORTDIR} <<EOF
${PSMCDIR}/psmc -N25 -t15 -r5 -p "4+25*2+4+6" -o ${i%.psmcfa}.out.psmc ${i}
EOF
done

for i in V*.psmcfa ; do 
qsub -N psmc -V -cwd -l h_rt=24:00:00,mem_free=4G -o ${REPORTDIR} -e ${REPORTDIR} <<EOF
${PSMCDIR}/psmc -N25 -t15 -r5 -p "4+20*2+6" -o ${i%.psmcfa}.out.psmc ${i}
EOF
done

# Check likelihood scores
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

