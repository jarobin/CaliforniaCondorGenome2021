'''
Input = raw VCF
Output = filtered VCF prints to screen
- Filtered sites are marked as FAIL_? in the 7th column
- Sites that pass go on to genotype filtering
- Filtered out genotypes are changed to './.', all others reported

NOTE: Min and max genotype depth hard-coded here (see below)

Possible usage:

SCRIPT=customFilterVCF.py
python ${SCRIPT} myfile.vcf.gz | bgzip > myfile_filtered.vcf.gz
tabix -p vcf myfile_filtered.vcf.gz

'''

import sys
import gzip
import re
import scipy.stats as ss


vcf_file = sys.argv[1]
VCF = gzip.open(vcf_file, 'r')


# Min depth (1/3x mean) and max depth (2x mean)
### Mean depth CYW1141: 68.5129
### Mean depth VulGry1: 38.2535
### Mean depth BGI_N323: 25.2678
### Mean depth CRW1112: 29.2387
minD={'CYW1141':22,'VulGry1':12,'BGI_N323':8,'CRW1112':9}
maxD={'CYW1141':137,'VulGry1':76,'BGI_N323':50,'CRW1112':58}


# Get list of samples in VCF file
samples=[]
for line in VCF:
    if line.startswith('##'):
        pass
    else:
        for i in line.split()[9:]: samples.append(i)
        break


# Go back to beginning of file
VCF.seek(0)


# Filter to be applied to individual genotypes
### sample is the sample name
### GT_entry is the entire genotype entry for that individual
### ADpos is the position of the AD in FORMAT (typically GT:AD:DP:GQ)
### DPpos is the position of the DP in FORMAT
### GQpos is the position of the GQ in FORMAT
def GTfilter(sample, GT_entry, ADpos, DPpos, GQpos):
    if GT_entry[:1]=='.' : return GT_entry
    else:
        gt=GT_entry.split(':')
        if gt[0] in ('0/0','0/1','1/1') and gt[GQpos]!='.' and gt[DPpos]!='.':
            DP=int(gt[DPpos])
            GQ=float(gt[GQpos])
            if GQ>=0.0 and minD[sample]<=DP<=maxD[sample]:
                REF=float(gt[ADpos].split(',')[0])
                AB=float(REF/DP)
                if gt[0]=='0/0':
                    if AB>=0.9: return GT_entry
                    else: return './.:' + ':'.join(gt[1:])
                elif gt[0]=='0/1':
                    if 0.2<=AB<=0.8: return GT_entry
                    else: return './.:' + ':'.join(gt[1:])
                elif gt[0]=='1/1':
                    if AB<=0.1: return GT_entry
                    else: return './.:' + ':'.join(gt[1:])
                else: './.:' + ':'.join(gt[1:])
            else: return './.:' + ':'.join(gt[1:])
        else: return './.:' + ':'.join(gt[1:])


# Write header lines
### Add new header lines for filters being added - for GATK compatibility
for line0 in VCF:
    if line0.startswith('#'):
        if line0.startswith('##FORMAT'):
            sys.stdout.write('##FILTER=<ID=FAIL_refN,Description="Low quality">\n')
            sys.stdout.write('##FILTER=<ID=FAIL_multiAlt,Description="Low quality">\n')
            sys.stdout.write('##FILTER=<ID=FAIL_badAlt,Description="Low quality">\n')
            sys.stdout.write('##FILTER=<ID=FAIL_noQUAL,Description="Low quality">\n')
            sys.stdout.write('##FILTER=<ID=FAIL_noINFO,Description="Low quality">\n')
            sys.stdout.write('##FILTER=<ID=FAIL_mutType,Description="Low quality">\n')
            sys.stdout.write('##FILTER=<ID=FAIL_noGQi,Description="Low quality">\n')
            sys.stdout.write('##FILTER=<ID=FAIL_noDPi,Description="Low quality">\n')
            sys.stdout.write('##FILTER=<ID=FAIL_noADi,Description="Low quality">\n')
            sys.stdout.write(line0)
            break
        else: sys.stdout.write(line0)


# Go through VCF file line by line to apply filters
for line0 in VCF:
    if line0.startswith('#'):
        sys.stdout.write(line0); continue

### For all other lines:
    line=line0.strip().split('\t')

### Site filtering:
### Keep any filters that have already been applied
    filter=[]
    if line[6] not in ('.', 'PASS'):
        filter.append(line[6])

### Reference must not be N
    if line[3]=='N':
        filter.append('FAIL_refN')
        sys.stdout.write('%s\t%s\t%s\n' % ('\t'.join(line[0:6]), ';'.join(filter), '\t'.join(line[7:])) ) ; continue

### Alternate allele must not be multiallelic or <NON_REF>
    if ',' in line[4]:
        filter.append('FAIL_multiAlt')

    if '<NON_REF>' in line[4]:
        filter.append('FAIL_badAlt')

### Must have a valid QUAL
    if line[5]=='.':
        filter.append('FAIL_noQUAL')

### Access INFO field annotations
    if ';' in line[7]:
        INFO=line[7].split(';')
        d=dict(x.split('=') for x in INFO)
    else:
        INFO=line[7]
        if '=' in INFO:
            d={INFO.split('=')[0]:INFO.split('=')[1]}
        else: filter.append('FAIL_noINFO')

### Only accept sites that are monomorphic or simple SNPs
    if 'VariantType' not in d or d['VariantType'] not in ('NO_VARIATION', 'SNP'):
        filter.append('FAIL_mutType')

### Get the position of AD, DP, GQ value in genotype fields
    if 'AD' in line[8]:
        ADpos=line[8].split(':').index('AD')
    else: filter.append('FAIL_noADi')

    if 'DP' in line[8]:
        DPpos=line[8].split(':').index('DP')
    else: filter.append('FAIL_noDPi')

    if 'GQ' in line[8]:
        ff=line[8].split(':')
        GQpos=[ff.index(x) for x in ff if 'GQ' in x][0]
    else: filter.append('FAIL_noGQi')

### If any filters failed, write out line and continue
    if filter!=[]:
        sys.stdout.write('%s\t%s\t%s\n' % ('\t'.join(line[0:6]), ';'.join(filter), '\t'.join(line[7:])) ) ; continue

### Genotype filtering:
    GT_list=[]
    for i in range(0,len(samples)):
        GT=GTfilter(samples[i],line[i+9],ADpos,DPpos,GQpos)
        GT_list.append(GT)
    if filter==[]:
        filter.append('PASS')

### Write out new line
    sys.stdout.write('%s\t%s\t%s\t%s\t%s\n' % ('\t'.join(line[0:6]), ';'.join(filter), ';'.join('{0}={1}'.format(key, val) for key, val in sorted(d.items())), line[8], '\t'.join(GT_list)) )


# Close files and exit
VCF.close()
exit()


