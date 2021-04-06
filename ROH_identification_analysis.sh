### Condor ROH identification and analysis
# Input VCF files filtered to remove all sites without PASS in filter column or without
# a called genotype. Only kept scaffolds 1-32 (scaffolds >=1 Mb)

### Convert to plink format and reformat scaffold names

REPORTDIR=~/project/reports/condor/plink
for i in $(ls *_simple_PASS_nomissing_scaffold1-32.vcf.gz); do
export VCF=${i}
qsub -N plink_convert -V -cwd -l h_rt=24:00:00,mem_free=8G -o ${REPORTDIR} -e ${REPORTDIR} << EOF
plink --vcf ${VCF} --keep-allele-order --allow-extra-chr --chr-set 32 --out ${VCF%.vcf.gz}
sed 's/HiC_scaffold_//g' ${VCF%.vcf.gz}.bim > ${VCF%.vcf.gz}.bim_rename
mv ${VCF%.vcf.gz}.bim ${VCF%.vcf.gz}.bim_orig
mv ${VCF%.vcf.gz}.bim_rename ${VCF%.vcf.gz}.bim
EOF
done


### ROH identification
# Exclude scaffold 5 (chrZ)

# Default parameters
REPORTDIR=~/project/reports/condor/plink
for i in $(ls *fam); do
export FILE=${i%.fam}
qsub -N plink_ROH -V -cwd -l h_rt=24:00:00,mem_free=48G -o ${REPORTDIR} -e ${REPORTDIR} << EOF
plink --bfile ${FILE} --out ${FILE}_ROH --homozyg --allow-extra-chr --chr-set 32 --chr 1,2,3,4,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32
EOF
done


# FID  IID  PHE     NSEG       KB    KBAVG
# BGI_N323	BGI_N323   -9       33  46958.6  1422.99 
# CRW1112   CRW1112   -9      117   244252  2087.62 
# CYW1141   CYW1141   -9      122   288291  2363.04 
# VulGry1   VulGry1   -9       45  63405.7  1409.01 


### ROH length sums
# BGI_N323	46958632
# CRW1112	244251515
# CYW1141	288291227
# VulGry1	63405667

# Percentages (divide summed length by 1106731338, the sum of scaffolds 1-4, 6-32)
# BGI_N323	0.04243002
# CRW1112	0.22069630
# CYW1141	0.26048890
# VulGry1	0.05729093


### Shared ROH

BED_A=CRW1112_simple_PASS_nomissing_scaffold1-32_ROH.hom.bed
BED_B=CYW1141_simple_PASS_nomissing_scaffold1-32_ROH.hom.bed
OUT=CRW1112_CYW1141_simple_PASS_nomissing_scaffold1-32_ROH.hom_shared.bed
bedtools intersect -a ${BED_A} -b ${BED_B} > ${OUT}

# CRW1112/CYW1141 shared ROH sum: 103406110 (103.4 Mb)


# Jaccard intersection (ROH sharing statistic)

cp ~/project/condor/reference/condor/gc_PacBio_HiC/*sizes .
head -n 32 g*sizes | grep -v "HiC_scaffold_5" | sed 's/HiC_scaffold_//g' > autos.sizes
bedtools jaccard -g autos.sizes -a CRW1112_simple_PASS_nomissing_scaffold1-32_ROH.hom.bed -b CYW1141_simple_PASS_nomissing_scaffold1-32_ROH.hom.bed
# intersection	union	jaccard	n_intersections
# 103406110	429136632	0.240963	73

# 24% ROH sharing


### ROH sharing and identity

# Reformat shared ROH bed to get proper scaffold names to match with VCF files
sed 's/^/HiC_scaffold_/g' CRW1112_CYW1141_simple_PASS_nomissing_scaffold1-32_ROH.hom_shared.bed > CRW1112_CYW1141_simple_PASS_nomissing_scaffold1-32_ROH.hom_shared_scaffnames.bed

# Extract regions in shared ROH
VCF_A=/mnt/Walllab/jrobinson/condor/condorRef/vcfs/CRW1112_simple.vcf.gz
VCF_B=/mnt/Walllab/jrobinson/condor/condorRef/vcfs/CYW1141_simple.vcf.gz
BED=CRW1112_CYW1141_simple_PASS_nomissing_scaffold1-32_ROH.hom_shared_scaffnames.bed
OUT_A=CRW1112_simple_sharedROH.vcf.gz
OUT_B=CYW1141_simple_sharedROH.vcf.gz

bcftools view -R ${BED} -Oz -o ${OUT_A} ${VCF_A} &
bcftools view -R ${BED} -Oz -o ${OUT_B} ${VCF_B} &
tabix -p vcf ${OUT_A}
tabix -p vcf ${OUT_B}

# Bcftools complained in subsequent merge step about "AD" in FORMAT, so remove all FORMAT fields besides genotypes
bcftools annotate -x FORMAT --threads 4 -Oz -o ${OUT_A%.vcf.gz}_recode.vcf.gz ${OUT_A} &
bcftools annotate -x FORMAT --threads 4 -Oz -o ${OUT_B%.vcf.gz}_recode.vcf.gz ${OUT_B} &
tabix -p vcf ${OUT_A%.vcf.gz}_recode.vcf.gz
tabix -p vcf ${OUT_B%.vcf.gz}_recode.vcf.gz

# Make joint VCF file
bcftools merge -m all --threads 4 -Oz -o CRW1112_CYW1141_simple_sharedROH.vcf.gz ${OUT_A%.vcf.gz}_recode.vcf.gz ${OUT_B%.vcf.gz}_recode.vcf.gz
tabix -p vcf CRW1112_CYW1141_simple_sharedROH.vcf.gz 

# Count number of fixed differences per shared ROH region
# Note: Checked that all hets at PASS sites encoded as 0/1
VCF=CRW1112_CYW1141_simple_sharedROH.vcf.gz
A=CRW1112_CYW1141_simple_PASS_nomissing_scaffold1-32_ROH.hom_shared_scaffnames.bed
B=CRW1112_CYW1141_simple_PASS_nomissing_scaffold1-32_ROH.hom_shared_scaffnames_fixeddiff.txt

while read -r CHR START END ; do
tabix ${VCF} ${CHR}:${START}-${END} | awk '$7=="PASS" && $10!=$11 && $10!="./." && $11!="./."' | grep -v "0/1" | wc -l
done < ${A} > ${B}

# Total shared ROH length (should equal what is given above)
paste ${A} ${B} | awk '{printf "%s\t%s\t%s\t%s\t%s\n", $1, $2, $3, $4, $4/($3-$2)}' | awk '{sum+=$3-$2}END{print sum}'
# 103406110

# Shared ROH and same haplotype (frequency of fixed differences <5e-5 over ROH length)
paste ${A} ${B} | awk '{printf "%s\t%s\t%s\t%s\t%s\n", $1, $2, $3, $4, $4/($3-$2)}' | awk '$5<5e-5{sum+=$3-$2}END{print sum}'
# 85234082

# Shared ROH and different haplotype (frequency of fixed differences >=5e-5 over ROH length)
paste ${A} ${B} | awk '{printf "%s\t%s\t%s\t%s\t%s\n", $1, $2, $3, $4, $4/($3-$2)}' | awk '$5>=5e-5{sum+=$3-$2}END{print sum}'
# 18172028

# Percent identical: 100*85234082/103406110 = 82.4%
# Percent not identical: 100*18172028/103406110 = 17.6%

