# Reads to variants pipeline, example commands

### Prep reference genome by building indexes

export FASTA=gc_PacBio_HiC.fasta

bwa index -a bwtsw ${FASTA} 
samtools faidx ${FASTA}
java -jar picard.jar CreateSequenceDictionary \
REFERENCE=${FASTA} OUTPUT=${FASTA%.fa*}.dict 


### Trim adapters from reads with BBDUK
# Note: Andean condor reads previously trimmed to remove the 10 leftmost bases from R1 and R2

export BBDUK_REF=~/anaconda3/envs/gentools/opt/bbmap-38.58-0/resources/adapters.fa
export BBDUK_RAM=32G
export BBDUK_THREADS=4

export FQ1=CYW1141.1_R1.fastq.gz
export FQ2=CYW1141.1_R2.fastq.gz
export RGID=CYW1141.1_A

bbduk.sh -Xmx${BBDUK_RAM} threads=${BBDUK_THREADS} ref=${BBDUK_REF} \
ktrim=r k=23 mink=11 hdist=1 tpe tbo mlf=0.5 \
in1=${FQ1} \
in2=${FQ2} \
out1=${RGID}_R1_trim.fastq.gz \
out2=${RGID}_R2_trim.fastq.gz \
stats=${RGID}_bbduk_report.txt 


### Align trimmed reads to reference and sort by query name

export BWAMEM_THREADS=4

export REFDIR=~/project/condor/reference/condor/gc_PacBio_HiC
export REFERENCE=gc_PacBio_HiC.fasta
export RGID=CYW1141.1_A
export READ1=${RGID}_R1_trim.fastq.gz
export READ2=${RGID}_R2_trim.fastq.gz

bwa mem -M -t ${BWAMEM_THREADS} ${REFDIR}/${REFERENCE} ${READ1} ${READ2} | \
samtools sort -n -@ ${BWAMEM_THREADS} -O BAM -T ${TMPDIR} -o ${RGID}_querysort.bam - 


### Fix mate information

export FIXMATE_THREADS=4

export RGID=CYW1141.1_A

samtools fixmate -@ ${FIXMATE_THREADS} -O BAM ${RGID}_querysort.bam ${RGID}_fixmate.bam


### Sort by coordinate and add read group info

export ADDRG_RAM=32G

export NAME=CYW1141
export RGID=CYW1141.1_A
export RGLB=Lib1
export RGPU=HYJLYCCXY.8.CYW1141
export RGPL=illumina

picard -Xmx${ADDRG_RAM} -Djava.io.tmpdir=${TMPDIR} AddOrReplaceReadGroups \
INPUT=${RGID}_fixmate.bam \
OUTPUT=${RGID}_addRG.bam \
RGID=${RGID} RGSM=${NAME} RGLB=${RGLB} RGPU=${RGPU} RGPL=${RGPL} \
SORT_ORDER=coordinate CREATE_INDEX=true VALIDATION_STRINGENCY=LENIENT \
TMP_DIR=${TMPDIR} 


### Remove duplicates

export MARKDUP_RAM=32G

export NAME=CYW1141

picard -Xmx${MARKDUP_RAM} -Djava.io.tmpdir=${TMPDIR} MarkDuplicates \
$(for i in *addRG.bam; do echo "INPUT=${i} "; done) \
OUTPUT=${NAME}_rmDup.bam \
METRICS_FILE=${NAME}_rmDup_metrics.txt \
MAX_RECORDS_IN_RAM=250000 MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000 \
REMOVE_DUPLICATES=true \
CREATE_INDEX=true TMP_DIR=${TMPDIR} 


### Indel realignment (can safely skip if using HaplotypeCaller to genotype)

# A: RealignerTargetCreator (multithreaded)

export INDEL_RTC_THREADS=4
export INDEL_RTC_RAM=32G

export REFDIR=~/project/condor/reference/condor/gc_PacBio_HiC
export REFERENCE=gc_PacBio_HiC.fasta
export NAME=CYW1141

gatk3 -Xmx${INDEL_RTC_RAM} -Djava.io.tmpdir=${TMPDIR} -T RealignerTargetCreator \
-nt ${INDEL_RTC_THREADS} \
-R ${REFDIR}/${REFERENCE} \
-I ${NAME}_rmDup.bam \
-o ${NAME}_realign.bam.intervals 


# B: IndelRealigner (not multithreaded)

export INDEL_IR_RAM=32G

export REFDIR=~/project/condor/reference/condor/gc_PacBio_HiC
export REFERENCE=gc_PacBio_HiC.fasta
export NAME=CYW1141

gatk3 -Xmx${INDEL_IR_RAM} -Djava.io.tmpdir=${TMPDIR} -T IndelRealigner \
-R ${REFDIR}/${REFERENCE} \
-I ${NAME}_rmDup.bam \
-o ${NAME}_realign.bam \
-targetIntervals ${NAME}_realign.bam.intervals


### Qualimap (for very large bams, increase memory to prevent job failure)

export QMAP_THREADS=4
export QMAP_RAM=32G

export NAME=CYW1141

qualimap bamqc -bam ${NAME}_realign.bam -c -nt ${QMAP_THREADS} --java-mem-size=${QMAP_RAM} 


### HaplotypeCaller
# SGE_TASK_ID is the job number within the array (SGE job scheduling)
# Lists of chromosome(s) provided in ${REFDIR}/interval_lists/*.list

export REFDIR=~/project/condor/reference/condor/gc_PacBio_HiC
export REFERENCE=gc_PacBio_HiC.fasta
export NAME=CYW1141

IDX=$(printf %03d ${SGE_TASK_ID})
SCAFFLIST=$(ls ${REFDIR}/interval_lists/*.list | head -n ${SGE_TASK_ID} | tail -n 1 )

gatk3 -Xmx32g -Djava.io.tmpdir=${TMPDIR} \
-T HaplotypeCaller \
-nct 4 \
-R ${REFDIR}/${REFERENCE} \
-ERC BP_RESOLUTION \
-out_mode EMIT_ALL_SITES \
--dontUseSoftClippedBases \
-mbq 20 \
-L ${SCAFFLIST} \
-I ${NAME}_realign.bam \
-o ${NAME}_${IDX}.g.vcf.gz 


# gVCF to VCF conversion and other post-processing
# As in HaplotypeCaller, jobs submitted as array

export REFDIR=~/project/condor/reference/condor/gc_PacBio_HiC
export REFERENCE=gc_PacBio_HiC.fasta
export NAME=CYW1141

SCAFF=$(head -n ${SGE_TASK_ID} ${REFDIR}/${SCAFFLIST} | tail -n 1 )
IDX=$(printf %03d ${SGE_TASK_ID})

gatk3 -Xmx32g -Djava.io.tmpdir=${TMPDIR} \
-T GenotypeGVCFs \
-R ${REFDIR}/${REFERENCE} \
-allSites \
-stand_call_conf 0 \
-L ${SCAFF} \
-V ${NAME}_${IDX}.g.vcf.gz \
-o ${NAME}_${IDX}.vcf.gz

gatk3 -Xmx32g -Djava.io.tmpdir=${TMPDIR} \
-T SelectVariants \
-R ${REFDIR}/${REFERENCE} \
-trimAlternates \
-V ${NAME}_${IDX}.vcf.gz \
-o ${NAME}_${IDX}_TrimAlt.vcf.gz

gatk3 -Xmx32g -Djava.io.tmpdir=${TMPDIR} \
-T VariantAnnotator \
-R ${REFDIR}/${REFERENCE} \
-G StandardAnnotation \
-A VariantType \
-A AlleleBalance \
-V ${NAME}_${IDX}_TrimAlt.vcf.gz \
-o ${NAME}_${IDX}_TrimAlt_Annot.vcf.gz


### Filtering

# A: Mask repeats and apply site-level filters

export REFDIR=~/project/condor/reference/condor/gc_PacBio_HiC
export REFERENCE=gc_PacBio_HiC.fasta
export MASK=gc_PacBio_HiC_repeats_TRF_WMdust.bed
export NAME=CYW1141

IDX=$(printf %03d ${SGE_TASK_ID})

gatk3 -Xmx16g -Djava.io.tmpdir=${TMPDIR} \
-T VariantFiltration \
-R ${REFDIR}/${REFERENCE} \
-mask ${REFDIR}/${MASK} -maskName "REP_TRF_WMdust" \
-filter "QD < 2.0" -filterName "FAIL_QD" \
-filter "FS > 60.0" -filterName "FAIL_FS" \
-filter "MQ < 40.0" -filterName "FAIL_MQ" \
-filter "MQRankSum < -12.5" -filterName "FAIL_MQRankSum" \
-filter "ReadPosRankSum < -8.0" -filterName "FAIL_ReadPosRankSum" \
-filter "ReadPosRankSum > 8.0" -filterName "FAIL_ReadPosRankSum" \
-filter "SOR > 4.0" -filterName "FAIL_SOR" \
-l ERROR \
-V ${NAME}_${IDX}_TrimAlt_Annot.vcf.gz \
-o ${NAME}_${IDX}_TrimAlt_Annot_Mask.vcf.gz 


# B: Custom filtering

export SCRIPT=~/project/condor/scripts/customFilterVCF.py

export NAME=CYW1141

IDX=$(printf %03d ${SGE_TASK_ID})

IDX=$(printf %03d ${SGE_TASK_ID})
VCF_IN=${NAME}_${IDX}_TrimAlt_Annot_Mask.vcf.gz
VCF_OUT=${NAME}_${IDX}_TrimAlt_Annot_Mask_Filter.vcf.gz

python2.7 ${SCRIPT} ${VCF_IN} | bgzip > ${VCF_OUT}
tabix -p vcf ${VCF_OUT}


### Simplify VCF files to reduce size

export NAME=CYW1141

bcftools concat -O u $(ls *Filter.vcf.gz) | bcftools annotate -x INFO,FORMAT -O z -o ${NAME}_simple.vcf.gz -
tabix -p vcf ${NAME}_simple.vcf.gz


