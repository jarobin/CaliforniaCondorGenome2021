#!/bin/bash

# Repeat identification in turkey vulture genome, ASM69994v1
# Runs RepeatMasker, TRF, and WindowMasker
# Note: Ultimately used WindowMasker and TRF repeats for masking

source ~/anaconda3/etc/profile.d/conda.sh
conda activate repeatmask

mkdir -p ~/myscratch/ASM69994v1_repeats

# Get fasta for repeat masking and remove previous masking
cd ~/project/condor/reference/turkeyVulture/ASM69994v1/
twoBitToFa -noMask ASM69994v1.2bit ~/myscratch/ASM69994v1_repeats/ASM69994v1.fa

cd ~/myscratch/ASM69994v1_repeats
faToTwoBit ASM69994v1.fa ASM69994v1.2bit
twoBitInfo ASM69994v1.2bit ASM69994v1.sizes
cut -f1 ASM69994v1.sizes > ASM69994v1.seq_list

mkdir parts
split --suffix-length=2 --numeric-suffixes=1 --number=l/99 ASM69994v1.seq_list parts/ASM69994v1.seq_list_

export FASTA=ASM69994v1.fa
export SPECIES="Cathartes"
export MAINDIR=~/myscratch/ASM69994v1_repeats

cd ${MAINDIR}

export QSUB=/opt/sge/bin/lx-amd64/qsub
export DATE=$(date +%Y%m%d)
export SEQDIR=${MAINDIR}/parts
export NUMJOBS=$(ls ${SEQDIR}/${FASTA%.fa*}.seq_list_* | wc -l)
export REPORTDIR=~/project/reports/repeatmask/${FASTA%.fa*}_${DATE}
mkdir -p ${REPORTDIR}


cat > repeatmask_1_${FASTA%.fa*}_${DATE}.sh <<EOF
source ~/anaconda3/etc/profile.d/conda.sh
conda activate repeatmask

NUM=\$(printf %02d \${SGE_TASK_ID})
cd ${SEQDIR}
SEQLIST=${FASTA%.fa*}.seq_list_\${NUM}
mkdir -p seq\${NUM}
mv \${SEQLIST} seq\${NUM}
cd seq\${NUM}
mkdir -p TRF; mkdir -p RepeatMasker

sed '/^[ \t]*$/d' \${SEQLIST} | while read n ; do
twoBitToFa ${MAINDIR}/${FASTA%.fa*}.2bit:\${n} ${FASTA%.fa*}_\${n}.fa
if [ \${?} -ne 0 ]; then echo "FAIL twoBitToFa \${n}" >> fail.log ; continue ; fi
f=${FASTA%.fa*}_\${n}.fa
trf \${f} 2 7 7 80 10 50 2000 -m -d
if [ \${?} -ne 0 ]; then echo "FAIL trf \${n}" >> fail.log ; continue ; fi
TRFdat_to_bed_with_period.py --dat \${f}.*.dat --bed \${f}.trf.bed
if [ \${?} -ne 0 ]; then echo "FAIL TRFdat_to_bed_with_period.py \${n}" >> fail.log ; continue ; fi
mv *.html TRF; mv *.mask TRF; mv *.dat TRF
RepeatMasker -engine crossmatch -s -align -species "${SPECIES}" \${f}
if [ \${?} -ne 0 ]; then echo "FAIL RepeatMasker \${n}" >> fail.log ; continue ; fi
mv *.align RepeatMasker; mv *.cat RepeatMasker; mv *.log RepeatMasker; mv *.masked RepeatMasker; mv *.tbl RepeatMasker
rm \${f}
done

if [ -f fail.log ]; then echo "FAIL" ; exit; fi

i=${FASTA%.fa*}_\${NUM}
sed '/^[ \t]*$/d' \${SEQLIST} | while read n ; do
f=${FASTA%.fa*}_\${n}.fa
sort -k1,1 -k2,2n \${f}.trf.bed | awk '{printf "%s\t%s\t%s\t%s\n", \$1, \$2, \$3, \$4}' >> \${i}_trf.merged.bed
if [ \${?} -ne 0 ]; then echo "FAIL trf merge \${n}" >> fail.log ; continue ; fi
grep -v "There were no repetitive sequences detected" \${f}.out | tail -n+4 | sed -e 's/  */\t/g' | sort -k5,5 -k6,6n | awk '{printf "%s\t%s\t%s\n", \$5, \$6-1, \$7}' >> \${i}_RM.merged.bed
if [ \${?} -ne 0 ]; then echo "FAIL RM merge \${n}" >> fail.log ; continue ; fi
done

if [ -f fail.log ]; then echo "FAIL" ; exit; fi

awk '\$4<=12' \${i}_trf.merged.bed > \${i}_trf.merged.maxperiod12.bed
if [ \${?} -ne 0 ]; then echo "FAIL trf filter \${n}" >> fail.log ; exit ; fi

echo "SUCCESS" > success.log

echo "DONE"

EOF

${QSUB} -N rmask1 -t 1-${NUMJOBS} -V -wd ${MAINDIR} -l h_rt=72:00:00,mem_free=4G -e ${REPORTDIR} -o ${REPORTDIR} repeatmask_1_${FASTA%.fa*}_${DATE}.sh


cat > repeatmask_2_${FASTA%.fa*}_${DATE}.sh <<EOF
source ~/anaconda3/etc/profile.d/conda.sh
conda activate repeatmask

cd ${SEQDIR}
NUMS=\$(cat */success.log | grep "SUCCESS" | wc -l)
if [ \${NUMS} -ne ${NUMJOBS} ]; then echo "FAIL" ; exit ; fi

for i in * ; do
cat \${i}/*trf.merged.maxperiod12.bed >> ${MAINDIR}/${FASTA%.fa*}_repeats_TRF.bed
cat \${i}/*RM.merged.bed >> ${MAINDIR}/${FASTA%.fa*}_repeats_RM.bed
done

cd ${MAINDIR}
windowmasker -mk_counts -in ${FASTA} -out ${FASTA%.fa*}_repeats_WM.counts
windowmasker -ustat ${FASTA%.fa*}_repeats_WM.counts -in ${FASTA} -out ${FASTA%.fa*}_repeats_WM.intervals
windowmasker -ustat ${FASTA%.fa*}_repeats_WM.counts -dust true -in ${FASTA} -out ${FASTA%.fa*}_repeats_WMdust.intervals

perl -wpe 'if (s/^>(.*)//) { \$chr = \$1; } if (/^(\d+) - (\d+)/) { \$s=\$1; \$e=\$2+1; s/(\d+) - (\d+)/\$chr\t\$s\t\$e/; }' ${FASTA%.fa*}_repeats_WM.intervals | sed '/^\$/d' > ${FASTA%.fa*}_repeats_WM.bed
perl -wpe 'if (s/^>(.*)//) { \$chr = \$1; } if (/^(\d+) - (\d+)/) { \$s=\$1; \$e=\$2+1; s/(\d+) - (\d+)/\$chr\t\$s\t\$e/; }' ${FASTA%.fa*}_repeats_WMdust.intervals | sed '/^\$/d' > ${FASTA%.fa*}_repeats_WMdust.bed

mkdir -p WM
mv ${FASTA%.fa*}_repeats_WM.counts WM
mv ${FASTA%.fa*}_repeats_WM*.intervals WM

echo "DONE"

cd ${MAINDIR}
cd ..
tar -czf ${FASTA%.fa*}_repeats.tar.gz ${MAINDIR}

EOF

${QSUB} -N rmask2 -hold_jid rmask1 -V -wd ${MAINDIR} -l h_rt=24:00:00,mem_free=4G -e ${REPORTDIR} -o ${REPORTDIR} repeatmask_2_${FASTA%.fa*}_${DATE}.sh



