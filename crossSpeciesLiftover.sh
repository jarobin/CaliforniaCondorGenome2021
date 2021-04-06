#!/bin/bash
#
# Script for generating cross-species reciprocal best chain files for liftOver
# Adapted from UCSC scripts by Jacqueline Robinson, 18 Nov. 2019
#
# To be run on a HPC cluster with SGE/UGE queuing system
#
# Starting from repeat-masked target and query genomes in .2bit format, produces 
# all.chain, over.chain, and rbest.chain for target to query, plus reciprocal rbest.chain 
# for query to target
# 
#
# USAGE:
# - Create a working directory ($BASEDIR)
# - Put .2bit and .sizes files for target and query genomes in $BASEDIR
# - Set variables below
# - In an interactive session: ./crossSpeciesLiftover.sh
#
#
# REQUIREMENTS:
# - UCSC utilities (make accessible via PATH)
#   - http://hgdownload.cse.ucsc.edu/admin/exe/
#
# - UCSC perl scripts (make accessible via PATH)
#   - partitionSequence.pl
#     - https://github.com/ucscGenomeBrowser/kent/blob/master/src/hg/utils/automation/partitionSequence.pl
#   - blastz-normalizeLav
#     - https://github.com/ucscGenomeBrowser/kent/blob/master/src/hg/utils/automation/blastz-normalizeLav
#   - blastz-run-ucsc
#   - blastz-run-ucsc (comment out lines 559-568, then comment out lines 85-91 and instead use: my $defaultPath = $ENV{'PATH'}; )
#     - https://github.com/ucscGenomeBrowser/kent/blob/master/src/hg/utils/automation/blastz-run-ucsc
#
# - .2bit and .sizes files for repeat-masked target and query genomes
#   - Repeats should be soft-masked (in lower-case letters)
#   - .2bit files can be generated from fasta files with faToTwoBit
#   - .sizes files can be generated from fasta files with faSize
#
# - Values for DEF parameter file
#   - Look to http://genomewiki.ucsc.edu/index.php/DoBlastzChainNet.pl and UCSC logs for guidance on setting parameter values
#
#
# NOTES:
# - Target genome is the genome you want to convert coordinates FROM, query is the genome you want to convert coordinates TO
#   - "Target" = "reference" = "from" = "old"
#   - "Query" = "to" = "new"
#
# - Outline of the procedure:
#   - Split the target genome into overlapping chunks (tLAP>0), split the query genome into non-overlapping chunks (qLAP=0)
#   - Align each query sequence against each target sequence (all * all) with lastz
#   - Convert alignments (.psl) to chains
#   - Convert chains (all.chain) to nets
#   - Extract chains for liftOver (over.chain)
#   - Generate reciprocal best chains in each direction (.rbest.chain)
#
# - Adjust the partitionining parameters (SEQ?_CHUNK, SEQ?_LAP, SEQ?_LIMIT) to yield a reasonable number of lastz jobs (<=100,000)
#   - There will be N_target_chunks lastz jobs, each of which is an array of N_query_chunks jobs
#   - There will be 1 lastzCheck job, which submits a single array of N_chains pslToChain jobs (N_chains < N_target_chunks)
#   - There will be 1 reciprocalBest job

 
################################################################################

# PIPELINE
# - Note: Ensure .2bit and .sizes files for both target and query are present in $BASEDIR folder!

# Set variables:
QSUB=/opt/sge/bin/lx-amd64/qsub

BASEDIR=/wynton/scratch/robinsonj/gc_PacBio_HiC_to_galGal6
TMPDIR=${BASEDIR}/temp
REPORTDIR=/wynton/home/walllab/robinsonj/project/reports/liftover/gc_PacBio_HiC_to_galGal6_$(date +"%Y%m%d")

tNAME=gc_PacBio_HiC
tTWOBIT=${BASEDIR}/${tNAME}.2bit
tSIZES=${BASEDIR}/${tNAME}.sizes
tCHUNK=20000000
tLAP=10000
tLIMIT=40

qNAME=galGal6
qTWOBIT=${BASEDIR}/${qNAME}.2bit
qSIZES=${BASEDIR}/${qNAME}.sizes
qCHUNK=1000000
qLAP=0
qLIMIT=2000

BLASTZ=/wynton/home/walllab/robinsonj/bin/lastz-1.04.00
BLASTZ_Q=/wynton/home/walllab/robinsonj/project/liftover/HoxD55.q
BLASTZ_O=400
BLASTZ_E=30
BLASTZ_K=2200
BLASTZ_L=4000
BLASTZ_H=2000
BLASTZ_M=254
BLASTZ_Y=3400
MIN_SCORE=1000
LINEAR_GAP=loose

# Note: ensure there is a blank line at the end of BLASTZ_Q
# HoxD55.q:
#    A    C    G    T
#   91  -90  -25 -100
#  -90  100 -100  -25
#  -25 -100  100  -90
# -100  -25  -90   91
#

################################################################################

# Create DEF

mkdir -p ${TMPDIR}
mkdir -p ${REPORTDIR}

cd ${BASEDIR}
cat <<EOF > DEF
BASEDIR=${BASEDIR}
TMPDIR=${TMPDIR}
# TARGET
tTWOBIT=${tTWOBIT}
tSIZES=${tSIZES}
tCHUNK=${tCHUNK}
tLIMIT=${tLIMIT}
tLAP=${tLAP}
# QUERY
qTWOBIT=${qTWOBIT}
qSIZES=${qSIZES}
qCHUNK=${qCHUNK}
qLIMIT=${qLIMIT}
qLAP=${qLAP}
# Lastz/axtChain parameters
BLASTZ=${BLASTZ}
BLASTZ_Q=${BLASTZ_Q}
BLASTZ_O=${BLASTZ_O}
BLASTZ_E=${BLASTZ_E}
BLASTZ_K=${BLASTZ_K}
BLASTZ_L=${BLASTZ_L}
BLASTZ_H=${BLASTZ_H}
BLASTZ_M=${BLASTZ_M}
BLASTZ_Y=${BLASTZ_Y}
MIN_SCORE=${MIN_SCORE}
LINEAR_GAP=${LINEAR_GAP}
EOF


# Partition
partitionSequence.pl ${tCHUNK} ${tLAP} ${tTWOBIT} ${tSIZES} ${tLIMIT} -xdir xdir.sh -rawDir ./psl -lstDir tParts > target.lst
partitionSequence.pl ${qCHUNK} ${qLAP} ${qTWOBIT} ${qSIZES} ${qLIMIT} -lstDir qParts > query.lst

if [ -d tParts ]; then
   echo 'constructing tParts/*.2bit files'
   ls tParts/*.lst | sed -e 's#tParts/##; s#.lst##;' | while read tPart
   do
     sed -e 's#.*.2bit:##;' tParts/${tPart}.lst \
       | twoBitToFa -seqList=stdin ${tTWOBIT} stdout \
         | faToTwoBit stdin tParts/${tPart}.2bit
   done
fi

if [ -d qParts ]; then
   echo 'constructing qParts/*.2bit files'
   ls qParts/*.lst | sed -e 's#qParts/##; s#.lst##;' | while read qPart
   do
     sed -e 's#.*.2bit:##;' qParts/${qPart}.lst \
       | twoBitToFa -seqList=stdin ${qTWOBIT} stdout \
         | faToTwoBit stdin qParts/${qPart}.2bit
   done
fi

chmod a+x xdir.sh
./xdir.sh


# Write and submit lastz.sh per target
NUMT=$(wc -l target.lst | cut -d' ' -f1)
NUMQ=$(wc -l query.lst | cut -d' ' -f1)

for T_path in $(cat target.lst); do
T=$(basename ${T_path})
Tsub=$(echo ${T} | sed 's/:/_/g')
mkdir -p ${REPORTDIR}/${Tsub}
cat <<EOF > lastz.sh 
Q_path=\$(head -n \${SGE_TASK_ID} query.lst | tail -n 1)
Q=\$(basename \${Q_path})
#
blastz-run-ucsc -outFormat psl ${T_path} \${Q_path} DEF ./psl/${T}/${T}_\${Q}.psl
if [ \${?} -eq 0 ]; then echo "SUCCESS" ; fi
EOF
${QSUB} -N lastz -t 1-${NUMQ} -V -wd ${BASEDIR} -l h_rt=24:00:00,mem_free=4G -o ${REPORTDIR}/${Tsub} -e ${REPORTDIR}/${Tsub} lastz.sh
done


# Write lastzCheck.sh
cat <<EOF > lastzCheck.sh 
set -o pipefail
#
# Check lastz jobs finished successfully
touch lastz.fail
touch lastz.success
#
for T_path in \$(cat target.lst); do
T=\$(basename \${T_path})
Tsub=\$(echo \${T} | sed 's/:/_/g')
TSCOUNT=0
TSCOUNT=\$(cat ${REPORTDIR}/\${Tsub}/*.o* | grep "SUCCESS" | wc -l) 
if [ \${TSCOUNT} -ne ${NUMQ} ] ; then echo "FAIL \${T}" >> lastz.fail ; else echo "SUCCESS \${T}" >> lastz.success ; fi
done
#
FCOUNT=\$(wc -l lastz.fail | cut -d' ' -f1)
SCOUNT=\$(wc -l lastz.success | cut -d' ' -f1)
if [ \${FCOUNT} -ne 0 ] ; then echo "Lastz job failure! Check lastz.fail" >> error.log ; exit 1 ; fi
if [ \${SCOUNT} -ne ${NUMT} ] ; then echo "Lastz job failure! Check target.lst and lastz.success" >> error.log ; exit 1 ; fi
#
# Concatenate .psl files per chunk
(cd ./psl; find . -maxdepth 1 -type d | grep '^./') | sed -e 's#/\$##; s#^./##' > tParts.lst
mkdir -p pslParts
for i in \$(cat tParts.lst) ; do find ./psl/\${i} -name "*.psl" | xargs cat | gzip -c > pslParts/\${i}.psl.gz ; done
if [ \${?} -ne 0 ] ; then echo "Concatenation of .psl files failed!" >> error.log ; exit 1 ; fi
#
# Submit pslToChain
echo "O=${BLASTZ_O} E=${BLASTZ_E}" | cat ${BLASTZ_Q} - | sed '/^$/d' > scores.matrix
mkdir -p chain
grep -v "^part" tParts.lst | cut -d':' -f1-2 | sort | uniq | sed 's/$/:/g' > chainParts.lst
grep "^part" tParts.lst >> chainParts.lst
NUMCHAINS=\$(wc -l chainParts.lst | cut -d' ' -f1)
mkdir -p ${REPORTDIR}/pslToChain
${QSUB} -N pslToChain -hold_jid lastz,lastzCheck -t 1-\${NUMCHAINS} -V -wd ${BASEDIR} -l h_rt=24:00:00,mem_free=4G -o ${REPORTDIR}/pslToChain -e ${REPORTDIR}/pslToChain pslToChain.sh
#
# Submit reciprocalBest
${QSUB} -N reciprocalBest -hold_jid lastz,lastzCheck,pslToChain -V -wd ${BASEDIR} -l h_rt=48:00:00,mem_free=4G -o ${REPORTDIR} -e ${REPORTDIR} reciprocalBest.sh 
#
echo "SUCCESS"
EOF


# Write pslToChain.sh
cat <<EOF >> pslToChain.sh
p=\$(head -n \${SGE_TASK_ID} chainParts.lst | tail -n 1)
zcat ./pslParts/\${p}*.psl.gz \\
| axtChain -verbose=0 -scoreScheme=scores.matrix -minScore=${MIN_SCORE} -linearGap=${LINEAR_GAP} -psl stdin ${tTWOBIT} ${qTWOBIT} stdout \\
| chainAntiRepeat ${tTWOBIT} ${qTWOBIT} stdin ./chain/\${p}.chain
if [ \${?} -eq 0 ]; then echo "SUCCESS" ; fi
EOF


# Write reciprocalBest.sh
cat <<EOF > reciprocalBest.sh 
set -o pipefail
#
# Check pslToChain jobs finished successfully
NUMCHAINS=\$(wc -l chainParts.lst | cut -d' ' -f1)
SCOUNT=0
SCOUNT=\$(cat ${REPORTDIR}/pslToChain/*.o* | grep "SUCCESS" | wc -l)
if [ \${NUMCHAINS} -ne \${SCOUNT} ] ; then echo "pslToChain failure! Check chainParts.lst, ${REPORTDIR}/pslToChain, ${BASEDIR}/chain" >> error.log ; exit 1 ; fi
#
# Make .all.chain (target to query)
find ./chain -name "*.chain" | chainMergeSort -inputList=stdin | gzip -c > ${tNAME}.${qNAME}.all.chain.gz
if [ \${?} -ne 0 ] ; then echo "Make .all.chain failure! Check ${REPORTDIR}/reciprocalBest*" >> error.log ; exit 1 ; fi
#
# Net chains
chainPreNet ${tNAME}.${qNAME}.all.chain.gz ${tSIZES} ${qSIZES} stdout \\
| chainNet stdin -minSpace=1 ${tSIZES} ${qSIZES} stdout /dev/null \\
| netSyntenic stdin noClass.net
if [ \${?} -ne 0 ] ; then echo "Netting .all.chain failure! Check ${REPORTDIR}/reciprocalBest*" >> error.log ; exit 1 ; fi
#
# Make .over.chain (target to query)
netChainSubset -verbose=0 noClass.net ${tNAME}.${qNAME}.all.chain.gz stdout \\
| chainStitchId stdin stdout | gzip -c > ${tNAME}.${qNAME}.over.chain.gz
if [ \${?} -ne 0 ] ; then echo "Make .over.chain failure! Check ${REPORTDIR}/reciprocalBest*" >> error.log ; exit 1 ; fi
gzip noClass.net
#
# Swap tNAME and qNAME
chainStitchId ${tNAME}.${qNAME}.over.chain.gz stdout \\
| chainSwap stdin stdout \\
| chainSort stdin ${qNAME}.${tNAME}.tBest.chain
if [ \${?} -ne 0 ] ; then echo "Swap tNAME and qNAME failure! Check ${REPORTDIR}/reciprocalBest*" >> error.log ; exit 1 ; fi
#
# Generate reciprocal best net (query to target)
chainPreNet ${qNAME}.${tNAME}.tBest.chain ${qSIZES} ${tSIZES} stdout \\
| chainNet -minSpace=1 -minScore=0 stdin ${qSIZES} ${tSIZES} stdout /dev/null \\
| netSyntenic stdin stdout \\
| gzip -c > ${qNAME}.${tNAME}.rbest.net.gz
if [ \${?} -ne 0 ] ; then echo "Make ${qNAME}.${tNAME}.rbest.net.gz failure! Check ${REPORTDIR}/reciprocalBest*" >> error.log ; exit 1 ; fi
#
# Extract reciprocal best chains (query to target)
netChainSubset ${qNAME}.${tNAME}.rbest.net.gz ${qNAME}.${tNAME}.tBest.chain stdout \\
| chainStitchId stdin stdout \\
| gzip -c > ${qNAME}.${tNAME}.rbest.chain.gz
if [ \${?} -ne 0 ] ; then echo "Make ${qNAME}.${tNAME}.rbest.chain.gz failure! Check ${REPORTDIR}/reciprocalBest*" >> error.log ; exit 1 ; fi
gzip ${qNAME}.${tNAME}.tBest.chain
#
# Swap reciprocal best chains (target to query)
chainSwap ${qNAME}.${tNAME}.rbest.chain.gz stdout \\
| chainSort stdin stdout \\
| gzip -c > ${tNAME}.${qNAME}.rbest.chain.gz
if [ \${?} -ne 0 ] ; then echo "Make ${tNAME}.${qNAME}.rbest.chain.gz failure! Check ${REPORTDIR}/reciprocalBest*" >> error.log ; exit 1 ; fi
#
# Generate swapped reciprocal best net (target to query)
chainPreNet ${tNAME}.${qNAME}.rbest.chain.gz ${tSIZES} ${qSIZES} stdout \\
| chainNet stdin -minSpace=1 -minScore=0 ${tSIZES} ${qSIZES} stdout /dev/null \\
| netSyntenic stdin stdout \\
| gzip -c > ${tNAME}.${qNAME}.rbest.net.gz
if [ \${?} -ne 0 ] ; then echo "Make ${tNAME}.${qNAME}.rbest.net.gz failure! Check ${REPORTDIR}/reciprocalBest*" >> error.log ; exit 1 ; fi
#
# Generate psl files from reciprocal best chains to test coverage
chainToPsl ${tNAME}.${qNAME}.rbest.chain.gz ${tSIZES} ${qSIZES} ${tTWOBIT} ${qTWOBIT} ${tNAME}.${qNAME}.rbest.chain.psl
if [ \${?} -ne 0 ] ; then echo "Make ${tNAME}.${qNAME}.rbest.chain.psl failure! Check ${REPORTDIR}/reciprocalBest*" >> error.log ; exit 1 ; fi
chainToPsl ${qNAME}.${tNAME}.rbest.chain.gz ${qSIZES} ${tSIZES} ${qTWOBIT} ${tTWOBIT} ${qNAME}.${tNAME}.rbest.chain.psl
if [ \${?} -ne 0 ] ; then echo "Make ${qNAME}.${tNAME}.rbest.chain.psl failure! Check ${REPORTDIR}/reciprocalBest*" >> error.log ; exit 1 ; fi
#
# Generate bed files from reciprocal best nets to test coverage
netToBed -maxGap=1 ${tNAME}.${qNAME}.rbest.net.gz ${tNAME}.${qNAME}.rbest.net.bed
if [ \${?} -ne 0 ] ; then echo "Make ${tNAME}.${qNAME}.rbest.net.bed failure! Check ${REPORTDIR}/reciprocalBest*" >> error.log ; exit 1 ; fi
netToBed -maxGap=1 ${qNAME}.${tNAME}.rbest.net.gz ${qNAME}.${tNAME}.rbest.net.bed
if [ \${?} -ne 0 ] ; then echo "Make ${qNAME}.${tNAME}.rbest.net.bed failure! Check ${REPORTDIR}/reciprocalBest*" >> error.log ; exit 1 ; fi
#
# Get coverage of reciprocal best chains and nets
tChCov=\$(awk '{print \$19}' ${tNAME}.${qNAME}.rbest.chain.psl | sed -e 's/,/\n/g' | awk '{N += \$1} END {printf "%d\n", N}')
qChCov=\$(awk '{print \$19}' ${qNAME}.${tNAME}.rbest.chain.psl | sed -e 's/,/\n/g' | awk '{N += \$1} END {printf "%d\n", N}')
tNetCov=\$(awk '{N += (\$3 - \$2)} END {printf "%d\n", N}' ${tNAME}.${qNAME}.rbest.net.bed)
qNetCov=\$(awk '{N += (\$3 - \$2)} END {printf "%d\n", N}' ${qNAME}.${tNAME}.rbest.net.bed)
echo -e "tChCov = \${tChCov}\nqChCov = \${qChCov}\ntNetCov = \${tNetCov}\nqNetCov = \${qNetCov}" > reciprocal_best_coverage.txt
gzip ${tNAME}.${qNAME}.rbest.chain.psl
gzip ${qNAME}.${tNAME}.rbest.chain.psl
gzip ${tNAME}.${qNAME}.rbest.net.bed
gzip ${qNAME}.${tNAME}.rbest.net.bed
#
# Generate reciprocal best axt files (optional)
netToAxt ${tNAME}.${qNAME}.rbest.net.gz ${tNAME}.${qNAME}.rbest.chain.gz ${tTWOBIT} ${qTWOBIT} stdout \\
| axtSort stdin stdout \\
| gzip -c > ${tNAME}.${qNAME}.rbest.axt.gz
if [ \${?} -ne 0 ] ; then echo "Make ${tNAME}.${qNAME}.rbest.axt.gz failure! Check ${REPORTDIR}/reciprocalBest*" >> error.log ; exit 1 ; fi
netToAxt ${qNAME}.${tNAME}.rbest.net.gz ${qNAME}.${tNAME}.rbest.chain.gz ${qTWOBIT} ${tTWOBIT} stdout \\
| axtSort stdin stdout \\
| gzip -c > ${qNAME}.${tNAME}.rbest.axt.gz
if [ \${?} -ne 0 ] ; then echo "Make ${qNAME}.${tNAME}.rbest.axt.gz failure! Check ${REPORTDIR}/reciprocalBest*" >> error.log ; exit 1 ; fi
#
# Generate reciprocal best maf files (optional)
axtToMaf -tPrefix=${tNAME}. -qPrefix=${qNAME}. ${tNAME}.${qNAME}.rbest.axt.gz ${tSIZES} ${qSIZES} stdout \\
| gzip -c > ${tNAME}.${qNAME}.rbest.maf.gz
if [ \${?} -ne 0 ] ; then echo "Make ${tNAME}.${qNAME}.rbest.maf.gz failure! Check ${REPORTDIR}/reciprocalBest*" >> error.log ; exit 1 ; fi
axtToMaf -tPrefix=${qNAME}. -qPrefix=${tNAME}. ${qNAME}.${tNAME}.rbest.axt.gz ${qSIZES} ${tSIZES} stdout \\
| gzip -c > ${qNAME}.${tNAME}.rbest.maf.gz
if [ \${?} -ne 0 ] ; then echo "Make ${qNAME}.${tNAME}.rbest.maf.gz failure! Check ${REPORTDIR}/reciprocalBest*" >> error.log ; exit 1 ; fi
#
echo "SUCCESS"
EOF


# Submit lastzCheck.sh
${QSUB} -N lastzCheck -hold_jid lastz -V -wd ${BASEDIR} -l h_rt=24:00:00,mem_free=4G -o ${REPORTDIR} -e ${REPORTDIR} lastzCheck.sh


echo "Done"

