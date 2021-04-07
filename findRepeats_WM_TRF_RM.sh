#!/bin/bash
# Repeat masking with WindowMasker, RepeatMasker, and Tandem Repeat Finder
#
# USAGE EXAMPLE:
# Make working directory (MAINDIR) and make this script executable from within MAINDIR 
# Define variables:
# export FASTA=gc_PacBio_HiC.fasta
# export SPECIES="Gymnogyps"
#
# export NAME=${FASTA%.fa*}
# export MAINDIR=~/myscratch/${NAME}/repeats_${NAME}
# export SCRIPT=findRepeats_WM_TRF_RM.sh
# export CHUNKSIZE=1000000
# export QSUB=/opt/sge/bin/lx-amd64/qsub
# export DATE=$(date +%Y%m%d)
# export REPORTDIR=${MAINDIR}/reports_${DATE}
# mkdir -p ${REPORTDIR}
# Execute (~10 mins):
# cd ${MAINDIR} ; chmod u+x ${SCRIPT}; ./${SCRIPT}

source ~/anaconda3/etc/profile.d/conda.sh
conda activate repeatmask

# Initialize MAINDIR
mkdir -p WM ; mkdir -p interval_lists ; mkdir -p parts
touch fail.log; touch success.log

# Generate 2bit and sizes files, removing any soft masking (if present)
faToTwoBit -noMask ${FASTA} ${NAME}.2bit
twoBitInfo ${NAME}.2bit ${NAME}.sizes

# Prep for Tandem Repeats Finder and RepeatMasker
# Get lists of large and small scaffolds 
awk -v var=${CHUNKSIZE} '$2>=var' ${NAME}.sizes > ${NAME}.sizes_big
awk -v var=${CHUNKSIZE} '$2<var' ${NAME}.sizes > ${NAME}.sizes_small

# Split large scaffolds into chunks
cat ${NAME}.sizes_big | while read -r CHR SIZE ; do
    mkdir -p ${MAINDIR}/parts/${CHR} ; cd ${MAINDIR}/parts/${CHR}
    mkdir -p TRF; mkdir -p RepeatMasker
    touch fail.log
    twoBitToFa ${MAINDIR}/${NAME}.2bit:${CHR} ${CHR}.fa
    if [ ${?} -ne 0 ]; then echo "FAIL twoBitToFa ${CHR}" >> ${MAINDIR}/fail.log ; exit ; fi    
    faSplit size ${CHR}.fa ${CHUNKSIZE} ${CHR}_ -lift=${CHR}.lft
    if [ ${?} -ne 0 ]; then echo "FAIL faSplit ${CHR}" >> ${MAINDIR}/fail.log ; exit ; fi    
    rm ${CHR}.fa
done
cd ${MAINDIR}

# Group small scaffolds into lists with total length <=CHUNKSIZE
cat ${NAME}.sizes_small | awk -v c=${CHUNKSIZE} -v NAME=${NAME} 'BEGIN{i=1 ; s=0}{ s+=$2 ; if (s<=c) {file=NAME".sizes_small_"i ; print $0 >> file} else {i+=1 ; s=$2 ; file=NAME".sizes_small_"i ; print $0 >> file} }' 
mv ${NAME}.sizes_small_* interval_lists
n=$(ls interval_lists/${NAME}.sizes_small_* | wc -l)
i=1
while [ ${i} -le ${n} ] ; do 
    mkdir -p ${MAINDIR}/parts/small_${i} ; cd ${MAINDIR}/parts/small_${i}
    mkdir -p TRF; mkdir -p RepeatMasker
    touch fail.log
    LIST=${MAINDIR}/interval_lists/${NAME}.sizes_small_${i}
    cut -f1 ${LIST} > small_${i}_seq.list
    twoBitToFa -seqList=small_${i}_seq.list ${MAINDIR}/${NAME}.2bit small_${i}.fa
    if [ ${?} -ne 0 ]; then echo "FAIL twoBitToFa small_${i}" >> ${MAINDIR}/fail.log ; exit ; fi
    ((i=i+1))
done
cd ${MAINDIR}


# Submit WindowMasker job
cat > repeatmask_1_${NAME}_${DATE}_WM.sh <<EOF
source ~/anaconda3/etc/profile.d/conda.sh
conda activate repeatmask

# Get counts and intervals
windowmasker -mk_counts -in ${FASTA} -out ${NAME}_repeats_WM.counts
if [ \${?} -ne 0 ]; then echo "FAIL WM -mk_counts" >> fail.log ; exit ; fi
windowmasker -ustat ${NAME}_repeats_WM.counts -in ${FASTA} -out ${NAME}_repeats_WM.intervals
if [ \${?} -ne 0 ]; then echo "FAIL WM -ustat" >> fail.log ; exit ; fi
windowmasker -ustat ${NAME}_repeats_WM.counts -dust true -in ${FASTA} -out ${NAME}_repeats_WMdust.intervals
if [ \${?} -ne 0 ]; then echo "FAIL WM -ustat with DUST" >> fail.log ; exit ; fi

# Reformat intervals to bed format and sort
perl -wpe 'if (s/^>(.*)//) { \$chr = \$1; } if (/^(\d+) - (\d+)/) { \$s=\$1; \$e=\$2+1; s/(\d+) - (\d+)/\$chr\t\$s\t\$e/; }' ${NAME}_repeats_WM.intervals | sed '/^\$/d' | sort -k1,1 -k2,2n > ${NAME}_repeats_WM.bed_temp1
if [ \${?} -ne 0 ]; then echo "FAIL reformat WM intervals (temp1)" >> fail.log ; exit ; fi
cut -d' ' -f1 ${NAME}_repeats_WM.bed_temp1 > ${NAME}_repeats_WM.bed_temp2a
if [ \${?} -ne 0 ]; then echo "FAIL reformat WM intervals (temp2a)" >> fail.log ; exit ; fi
cut -f2-3 ${NAME}_repeats_WM.bed_temp1 > ${NAME}_repeats_WM.bed_temp2b
if [ \${?} -ne 0 ]; then echo "FAIL reformat WM intervals (temp2b)" >> fail.log ; exit ; fi
paste ${NAME}_repeats_WM.bed_temp2a ${NAME}_repeats_WM.bed_temp2b > ${NAME}_repeats_WM.bed_temp3 
if [ \${?} -ne 0 ]; then echo "FAIL reformat WM intervals (temp3)" >> fail.log ; exit ; fi

perl -wpe 'if (s/^>(.*)//) { \$chr = \$1; } if (/^(\d+) - (\d+)/) { \$s=\$1; \$e=\$2+1; s/(\d+) - (\d+)/\$chr\t\$s\t\$e/; }' ${NAME}_repeats_WMdust.intervals | sed '/^\$/d' | sort -k1,1 -k2,2n > ${NAME}_repeats_WMdust.bed_temp1
if [ \${?} -ne 0 ]; then echo "FAIL reformat WM intervals with DUST" >> fail.log ; exit ; fi
cut -d' ' -f1 ${NAME}_repeats_WMdust.bed_temp1 > ${NAME}_repeats_WMdust.bed_temp2a
if [ \${?} -ne 0 ]; then echo "FAIL reformat WMdust intervals (temp2a)" >> fail.log ; exit ; fi
cut -f2-3 ${NAME}_repeats_WMdust.bed_temp1 > ${NAME}_repeats_WMdust.bed_temp2b
if [ \${?} -ne 0 ]; then echo "FAIL reformat WMdust intervals (temp2b)" >> fail.log ; exit ; fi
paste ${NAME}_repeats_WMdust.bed_temp2a ${NAME}_repeats_WMdust.bed_temp2b > ${NAME}_repeats_WMdust.bed_temp3 
if [ \${?} -ne 0 ]; then echo "FAIL reformat WMdust intervals (temp3)" >> fail.log ; exit ; fi

# Merge overlapping intervals and sort by genome order
sort -k1,1 -k2,2n ${NAME}_repeats_WM.bed_temp3 | bedtools merge -i stdin | bedtools sort -g ${NAME}.sizes -i stdin > ${NAME}_repeats_WM.bed
if [ \${?} -ne 0 ]; then echo "FAIL bedtools merge & sort WM intervals" >> fail.log ; exit ; fi
sort -k1,1 -k2,2n ${NAME}_repeats_WMdust.bed_temp3 | bedtools merge -i stdin | bedtools sort -g ${NAME}.sizes -i stdin > ${NAME}_repeats_WMdust.bed
if [ \${?} -ne 0 ]; then echo "FAIL bedtools merge & sort WM intervals with DUST" >> fail.log ; exit ; fi

# Clean up
mv ${NAME}_repeats_WM.counts WM
mv ${NAME}_repeats_WM*.intervals WM
mv ${NAME}_repeats_WM*.bed_temp WM

echo "SUCCESS WM" >> success.log

echo "DONE"
EOF

${QSUB} -N repeatmask_1_${NAME}_WM -V -wd ${MAINDIR} -l h_rt=72:00:00,mem_free=4G -e ${REPORTDIR} -o ${REPORTDIR} repeatmask_1_${NAME}_${DATE}_WM.sh


# Submit an array job for TRF/RM for each large scaffold
cat > repeatmask_2_${NAME}_${DATE}_TRF_RM_big.sh <<EOF
source ~/anaconda3/etc/profile.d/conda.sh
conda activate repeatmask

cd ${MAINDIR}/parts/\${CHR}
f=\$(ls *.fa | head -n \${SGE_TASK_ID} | tail -n 1)

# Make a temp file to be deleted upon job completion
echo "\${f}" > \${f}_fail.temp

# Run TRF
trf \${f} 2 7 7 80 10 50 2000 -m -d
if [ \${?} -ne 0 ]; then echo "FAIL trf \${f}" >> fail.log ; exit ; fi
TRFdat_to_bed_with_period.py --dat \${f}.*.dat --bed \${f}.trf.bed
if [ \${?} -ne 0 ]; then echo "FAIL TRFdat_to_bed_with_period.py \${f}" >> fail.log ; exit ; fi
mv \${f}*.html TRF; mv \${f}*.mask TRF; mv \${f}*.dat TRF

# Run RepeatMasker
RepeatMasker -engine crossmatch -s -align -species "${SPECIES}" \${f}
if [ \${?} -ne 0 ]; then echo "FAIL RepeatMasker \${f}" >> fail.log ; exit ; fi
mv \${f}*.align RepeatMasker; mv \${f}*.cat RepeatMasker; mv \${f}*.log RepeatMasker; mv \${f}*.masked RepeatMasker; mv \${f}*.tbl RepeatMasker

# Delete temp file and create temp success file
rm \${f}_fail.temp
echo "\${f}" > \${f}_success.temp

echo "DONE"
EOF

cat ${NAME}.sizes_big | while read -r CHR SIZE ; do
export CHR
NUMJOBS=$(ls parts/${CHR}/*.fa | wc -l)
${QSUB} -N repeatmask_2_${NAME}_TRF_RM_big_${CHR} -t 1-${NUMJOBS} -V -wd ${MAINDIR} -l h_rt=72:00:00,mem_free=4G -e ${REPORTDIR} -o ${REPORTDIR} repeatmask_2_${NAME}_${DATE}_TRF_RM_big.sh
done


# Submit an array job for TRF/RM for smaller chunks
cat > repeatmask_2_${NAME}_${DATE}_TRF_RM_small.sh <<EOF
source ~/anaconda3/etc/profile.d/conda.sh
conda activate repeatmask

cd ${MAINDIR}/parts/small_\${SGE_TASK_ID}
f=small_\${SGE_TASK_ID}.fa

# Make a temp file to be deleted upon job completion
echo "\${f}" > \${f}_fail.temp

# Run TRF
trf \${f} 2 7 7 80 10 50 2000 -m -d
if [ \${?} -ne 0 ]; then echo "FAIL trf \${f}" >> fail.log ; exit ; fi
TRFdat_to_bed_with_period.py --dat \${f}.*.dat --bed \${f}.trf.bed
if [ \${?} -ne 0 ]; then echo "FAIL TRFdat_to_bed_with_period.py \${f}" >> fail.log ; exit ; fi
mv \${f}*.html TRF; mv \${f}*.mask TRF; mv \${f}*.dat TRF

# Run RepeatMasker
RepeatMasker -engine crossmatch -s -align -species "${SPECIES}" \${f}
if [ \${?} -ne 0 ]; then echo "FAIL RepeatMasker \${f}" >> fail.log ; exit ; fi
mv \${f}*.align RepeatMasker; mv \${f}*.cat RepeatMasker; mv \${f}*.log RepeatMasker; mv \${f}*.masked RepeatMasker; mv \${f}*.tbl RepeatMasker

# Delete temp file and create temp success file
rm \${f}_fail.temp
echo "\${f}" > \${f}_success.temp

echo "DONE"
EOF

NUMJOBS=$(ls ${MAINDIR}/interval_lists/${NAME}.sizes_small_* | wc -l)
${QSUB} -N repeatmask_2_${NAME}_TRF_RM_small -t 1-${NUMJOBS} -V -wd ${MAINDIR} -l h_rt=72:00:00,mem_free=4G -e ${REPORTDIR} -o ${REPORTDIR} repeatmask_2_${NAME}_${DATE}_TRF_RM_small.sh


# Submit a job to concatenate/flatten/sort TRF/RM output, archive
cat > repeatmask_3_${NAME}_${DATE}.sh <<EOF
source ~/anaconda3/etc/profile.d/conda.sh
conda activate repeatmask

# Make sure all jobs completed successfully
aa=\$(ls parts/*/*.fa | grep -v "^parts/small_" | wc -l)
bb=\$(ls interval_lists/${NAME}.sizes_small_* | wc -l)
NUMJ=\$((\$aa + \$bb))

NUMFa=\$(cat parts/*/fail.log | wc -l)
NUMFb=\$(ls parts/*/fail.temp | wc -l)
if [ \${NUMFa} -gt 0 ]; then echo "FAIL one or more TRF/RM jobs failed" >> fail.log ; exit ; fi
if [ \${NUMFb} -gt 0 ]; then echo "FAIL one or more TRF/RM jobs failed" >> fail.log ; exit ; fi

NUMS=\$(cat parts/*/*success.temp | wc -l)
if [ \${NUMS} -ne \${NUMJ} ]; then echo "FAIL not all TRF/RM jobs completed successfully" >> fail.log ; exit ; fi
rm parts/*/*success.temp

# Lift, sort and concatenate big scaffold repeats
cat ${NAME}.sizes_big | while read -r CHR SIZE ; do
    cd ${MAINDIR}/parts/\${CHR}
    liftUp -type=.bed stdout \${CHR}.lft error *.fa.trf.bed | sort -k1,1 -k2,2n | awk '\$4<=12 {printf "%s\t%s\t%s\n", \$1, \$2, \$3}' >> ${MAINDIR}/${NAME}_repeats_TRF.bed_temp
    if [ \${?} -ne 0 ]; then echo "FAIL TRF lift & sort \${f}" >> ${MAINDIR}/fail.log ; exit ; fi
    liftUp -type=.out stdout \${CHR}.lft error *.fa.out | tail -n+4 | sed -e 's/  */\t/g' | sort -k5,5 -k6,6n | awk '{printf "%s\t%s\t%s\n", \$5, \$6-1, \$7}' >> ${MAINDIR}/${NAME}_repeats_RM.bed_temp
    if [ \${?} -ne 0 ]; then echo "FAIL RM lift & sort \${f}" >> ${MAINDIR}/fail.log ; exit ; fi
done
cd ${MAINDIR}

# Sort and concatenate small scaffold repeats
n=\$(ls ${MAINDIR}/interval_lists/${NAME}.sizes_small_* | wc -l)
i=1
while [ \${i} -le \${n} ] ; do 
    sort -k1,1 -k2,2n  parts/small_\${i}/small_\${i}.fa.trf.bed | awk '\$4<=12 {printf "%s\t%s\t%s\n", \$1, \$2, \$3}' >> ${MAINDIR}/${NAME}_repeats_TRF.bed_temp
    if [ \${?} -ne 0 ]; then echo "FAIL TRF sort & concatenate small_\${i}" >> fail.log ; exit ; fi
    grep -v "There were no repetitive sequences detected" parts/small_\${i}/small_\${i}.fa.out | tail -n+4 | sed -e 's/  */\t/g' | sort -k5,5 -k6,6n | awk '{printf "%s\t%s\t%s\n", \$5, \$6-1, \$7}' >> ${MAINDIR}/${NAME}_repeats_RM.bed_temp
    if [ \${?} -ne 0 ]; then echo "FAIL RM sort & concatenate small_\${i}" >> fail.log ; exit ; fi
    ((i=i+1))
done

# Merge overlapping intervals and sort by genome order
bedtools merge -i ${NAME}_repeats_TRF.bed_temp | bedtools sort -g ${NAME}.sizes -i stdin > ${NAME}_repeats_TRF.bed
if [ \${?} -ne 0 ]; then echo "FAIL bedtools merge & sort TRF" >> fail.log ; exit ; fi
echo "SUCCESS TRF" >> success.log
bedtools merge -i ${NAME}_repeats_RM.bed_temp | bedtools sort -g ${NAME}.sizes -i stdin > ${NAME}_repeats_RM.bed
if [ \${?} -ne 0 ]; then echo "FAIL bedtools merge & sort RM" >> fail.log ; exit ; fi
echo "SUCCESS RM" >> success.log

# Check if all jobs completed successfully and archive
NUMF=\$(wc -l fail.log | cut -d' ' -f1)
NUMS=\$(wc -l success.log | cut -d' ' -f1)
if [ \${NUMF} -gt 0 ]; then echo "FAIL WM/TRF/RM" >> fail.log ; exit ; fi
if [ \${NUMS} -ne 3 ]; then echo "FAIL one of WM/TRF/RM did not complete successfully" >> fail.log ; exit ; fi
cd ${MAINDIR}/..
tar -czf repeats_${NAME}.tar.gz ${MAINDIR}

echo "DONE"
EOF

${QSUB} -N repeatmask_3_${NAME} -hold_jid "repeatmask_1_${NAME}_*" -hold_jid "repeatmask_2_${NAME}_*" -V -wd ${MAINDIR} -l h_rt=72:00:00,mem_free=4G -e ${REPORTDIR} -o ${REPORTDIR} repeatmask_3_${NAME}_${DATE}.sh

