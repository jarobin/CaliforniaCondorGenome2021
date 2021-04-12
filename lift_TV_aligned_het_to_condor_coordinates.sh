# Lift over turkey vulture-aligned VCF files to get het. in sliding windows on condor coordinate system (for Fig. S2)

# Workflow 
# - Turn coordinates of calls and hets from VCF files into bed files
# - LiftOver bed files
# - Use bedtools intersect to get counts per window
# - Concatenate output


### Generate bed files for hets and calls

#! /bin/bash
# USAGE:
# SCRIPT=~/project/scripts/SlidingWindowHet/condor_vcf2LiftedBed_20200221.sh
# REPORTDIR=~/project/reports/condor/vcf2LiftedBed
# mkdir -p ${REPORTDIR}
# NUMJOBS=105
# for i in BGI_N323 CRW1112 CYW1141 VulGry1; do
# export NAME=${i}
# qsub -N vcf2LiftedBed -t 1-${NUMJOBS} -V -cwd -l h_rt=96:00:00 -o ${REPORTDIR} -e ${REPORTDIR} ${SCRIPT}
# done

source /wynton/home/walllab/robinsonj/anaconda3/etc/profile.d/conda.sh
conda activate gentools

CHAIN=~/project/condor/reference/gc_PacBio_HiC/liftover/ASM69994v1.gc_PacBio_HiC.rbest.chain.gz

cd ~/mygroup/condor/turkeyVultureRef/vcfs/${NAME}
IDX=$(printf %03d ${SGE_TASK_ID})
VCF=$(ls *_${IDX}_*Mask.vcf.gz)

# Print coordinates from VCF file to make hets and calls bed files
zcat ${VCF} | grep -v "^#" | grep -v "FAIL\|REP" | grep -v "\./\." | \
awk -v hets=${VCF%.vcf.gz}_hets.bed -v calls=${VCF%.vcf.gz}_calls.bed '{if ($10~"0/1") {printf ("%s\t%s\t%s\n", $1, $2-1, $2) >> hets ; printf ("%s\t%s\t%s\n", $1, $2-1, $2) >> calls } else {printf ("%s\t%s\t%s\n", $1, $2-1, $2) >> calls } }'

mv ${VCF%.vcf.gz}_hets.bed ~/mygroup/condor/turkeyVultureRef/winhet/${NAME}
mv ${VCF%.vcf.gz}_calls.bed ~/mygroup/condor/turkeyVultureRef/winhet/${NAME}

# LiftOver
cd ~/mygroup/condor/turkeyVultureRef/winhet/${NAME}
liftOver ${VCF%.vcf.gz}_hets.bed ${CHAIN} ${VCF%.vcf.gz}_hets_lifted.bed ${VCF%.vcf.gz}_hets_unMapped.bed
liftOver ${VCF%.vcf.gz}_calls.bed ${CHAIN} ${VCF%.vcf.gz}_calls_lifted.bed ${VCF%.vcf.gz}_calls_unMapped.bed


##########

### Concatenate and sort lifted bed files

# Sort hets (small files)

for NAME in BGI_N323 CRW1112 CYW1141 VulGry1; do
cd ~/mygroup/condor/turkeyVultureRef/winhet/${NAME}
cat *_hets_lifted.bed | sort -k1,1 -k2,2n > ${NAME}_hets_lifted_sorted.bed
done

# Sort calls (big files)

REPORTDIR=~/project/reports/condor/sortLiftedBed
mkdir -p ${REPORTDIR}
NUMJOBS=105
for i in BGI_N323 CRW1112 CYW1141 VulGry1; do
export NAME=${i}
cd ~/mygroup/condor/turkeyVultureRef/winhet/${i}
qsub -N sortLiftedBed -t 1-${NUMJOBS} -V -cwd -l h_rt=48:00:00,mem_free=16G -o ${REPORTDIR} -e ${REPORTDIR} <<EOF
cat *_calls_lifted.bed | awk -v var="HiC_scaffold_\${SGE_TASK_ID}" '\$1==var' > ${NAME}_calls_lifted_HiC_scaffold_\${SGE_TASK_ID}.bed
sort -k1,1 -k2,2n ${NAME}_calls_lifted_HiC_scaffold_\${SGE_TASK_ID}.bed > ${NAME}_calls_lifted_sorted_HiC_scaffold_\${SGE_TASK_ID}.bed
EOF
done


### Get overlap with 1-Mb windows to tally numbers of hets, calls per window
# NOTE: lifted coordinates will not always be single coverage on the new reference
#  - Made one set of files with overlapping coordinates merged, and one set with no merging
#  - Ultimately, results were virtually the same between the two sets
#  - Unmerged data used in final figure

# Coordinates of 1 Mb windows made with: 
# bedtools makewindows -g gc_PacBio_HiC.sizes -w 1000000 > gc_PacBio_HiC.sizes_1Mbwins.bed

WINS=~/project/condor/reference/gc_PacBio_HiC/gc_PacBio_HiC.sizes_1Mbwins.bed
cd ~/mygroup/condor/turkeyVultureRef/winhet

# Only use scaffolds 1-32 (scaffolds that are 1 Mb or larger)
for i in {1..32} ; do awk -v var=HiC_scaffold_${i} '$1==var' ${WINS} ; done > wins.temp
NUMJOBS=$(wc -l wins.temp | cut -d' ' -f1)

REPORTDIR=~/project/reports/condor/winOverlap
mkdir -p ${REPORTDIR}
for i in BGI_N323 CRW1112 CYW1141 VulGry1; do
export NAME=${i}
cd ~/mygroup/condor/turkeyVultureRef/winhet/${i}
qsub -N winOverlap -t 1-${NUMJOBS} -V -cwd -l h_rt=96:00:00,mem_free=4G -e ${REPORTDIR} -o ${REPORTDIR} <<EOF
MYTEMP=\${TMPDIR}/myworkdir_\${SGE_TASK_ID}
mkdir -p \${MYTEMP}
tWIN=win_\${SGE_TASK_ID}.bed
head -n \${SGE_TASK_ID} ~/mygroup/condor/turkeyVultureRef/winhet/wins.temp | tail -n 1 > \${MYTEMP}/\${tWIN}
cp hets_lifted_sorted/${NAME}_hets_lifted_sorted.bed \${MYTEMP}
CHR=\$(cut -f1 \${MYTEMP}/\${tWIN})
cBED=${NAME}_calls_lifted_sorted_\${CHR}.bed
cp calls_lifted_sorted/\${cBED} \${MYTEMP}
cd \${MYTEMP}
bedtools intersect -sorted -a \${tWIN} -b ${NAME}_hets_lifted_sorted.bed | awk '{sum+=\$3-\$2}END{print sum}' > hets_\${SGE_TASK_ID}.temp
bedtools intersect -sorted -a \${tWIN} -b ${NAME}_hets_lifted_sorted.bed | bedtools merge -i stdin | awk '{sum+=\$3-\$2}END{print sum}' > hets_merged_\${SGE_TASK_ID}.temp
bedtools intersect -sorted -a \${tWIN} -b \${cBED} | awk '{sum+=\$3-\$2}END{print sum}' > calls_\${SGE_TASK_ID}.temp
bedtools intersect -sorted -a \${tWIN} -b \${cBED} | bedtools merge -i stdin | awk '{sum+=\$3-\$2}END{print sum}' > calls_merged_\${SGE_TASK_ID}.temp
mv hets_\${SGE_TASK_ID}.temp ~/mygroup/condor/turkeyVultureRef/winhet/${NAME}
mv hets_merged_\${SGE_TASK_ID}.temp ~/mygroup/condor/turkeyVultureRef/winhet/${NAME}
mv calls_\${SGE_TASK_ID}.temp ~/mygroup/condor/turkeyVultureRef/winhet/${NAME}
mv calls_merged_\${SGE_TASK_ID}.temp ~/mygroup/condor/turkeyVultureRef/winhet/${NAME}
cd \${TMPDIR}
rm -rf myworkdir_\${SGE_TASK_ID}
EOF
done


### Organize files

cd ~/mygroup/condor/turkeyVultureRef/winhet
for i in BGI_N323 CRW1112 CYW1141 VulGry1 ; do 
mkdir -p ${i}/calls_lifted_sorted ; mv ${i}/*calls_lifted_sorted*bed ${i}/calls_lifted_sorted
mkdir -p ${i}/calls ; mv ${i}/*_calls.bed ${i}/calls
mkdir -p ${i}/calls_lifted_unsorted ; mv ${i}/*calls_lifted*bed ${i}/calls_lifted_unsorted
mkdir -p ${i}/calls_lifted_unmapped ; mv ${i}/*_calls_unMapped.bed ${i}/calls_lifted_unmapped
mkdir -p ${i}/hets ; mv ${i}/*_hets.bed ${i}/hets
mkdir -p ${i}/hets_lifted_unsorted ; mv ${i}/*hets_lifted*bed ${i}/hets_lifted_unsorted
mkdir -p ${i}/hets_lifted_unmapped ; mv ${i}/*_hets_unMapped.bed ${i}/hets_lifted_unmapped
mkdir -p ${i}/counts_hets_merged ; mv ${i}/hets_merged*temp ${i}/counts_hets_merged
mkdir -p ${i}/counts_calls_merged ; mv ${i}/calls_merged*temp ${i}/counts_calls_merged
mkdir -p ${i}/counts_hets_unmerged ; mv ${i}/hets_*temp ${i}/counts_hets_unmerged
mkdir -p ${i}/counts_calls_unmerged ; mv ${i}/calls_*temp ${i}/counts_calls_unmerged
done


### Concatenate counts to produce final files

for NAME in BGI_N323 CRW1112 CYW1141 VulGry1 ; do
cd ~/mygroup/condor/turkeyVultureRef/winhet/${NAME}
cd counts_hets_unmerged; cat $(ls -v) > ../counts_hets_unmerged.temp; cd ..
cd counts_hets_merged; cat $(ls -v) > ../counts_hets_merged.temp; cd ..
cd counts_calls_unmerged; cat $(ls -v) > ../counts_calls_unmerged.temp; cd ..
cd counts_calls_merged; cat $(ls -v) > ../counts_calls_merged.temp; cd ..
echo -e "chrom\twindow_start\tsites_total\tcalls_${NAME}\thets_${NAME}" > ${NAME}_het_1000000win_1000000step.txt
echo -e "chrom\twindow_start\tsites_total\tcalls_${NAME}\thets_${NAME}" > ${NAME}_het_1000000win_1000000step_merged.txt
paste ~/mygroup/condor/turkeyVultureRef/winhet/wins.temp counts_calls_unmerged.temp counts_hets_unmerged.temp | awk '{printf "%s\t%s\tNA\t%s\t%s\n", $1, $2+1, $4, $5}' >> ${NAME}_het_1000000win_1000000step.txt
paste ~/mygroup/condor/turkeyVultureRef/winhet/wins.temp counts_calls_merged.temp counts_hets_merged.temp | awk '{printf "%s\t%s\tNA\t%s\t%s\n", $1, $2+1, $4, $5}' >> ${NAME}_het_1000000win_1000000step_merged.txt
done


### Compress large bed files to save space

REPORTDIR=~/project/reports/condor/compress
mkdir -p ${REPORTDIR}
cd ~/mygroup/condor/turkeyVultureRef/winhet
for i in BGI_N323 CRW1112 CYW1141 VulGry1 ; do 
cd ${i}
for j in calls calls_lifted_unmapped calls_lifted_unsorted calls_lifted_sorted ; do 
cd ${j}
qsub -N compress -V -cwd -l h_rt=24:00:00,mem_free=4G -e ${REPORTDIR} -o ${REPORTDIR} <<EOF
for b in *.bed ; do bgzip \${b} ; done
EOF
cd ..
done
cd ..
done
