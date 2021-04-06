### Calculate absolute divergence
# See calc_dxy.py below

python ./calc_dxy.py BGI_N323_simple.vcf.gz CRW1112_simple.vcf.gz
# BGI_N323_simple.vcf.gz CRW1112_simple.vcf.gz Dxy = 30427420 / (2 * 829591817) = 0.01833879

python ./calc_dxy.py BGI_N323_simple.vcf.gz CYW1141_simple.vcf.gz
# BGI_N323_simple.vcf.gz CYW1141_simple.vcf.gz Dxy = 30476647 / (2 * 828107741) = 0.01840138

python ./calc_dxy.py VulGry1_simple.vcf.gz CRW1112_simple.vcf.gz
# VulGry1_simple.vcf.gz CRW1112_simple.vcf.gz Dxy = 23964940 / (2 * 938061047) = 0.01277366

python ./calc_dxy.py VulGry1_simple.vcf.gz CYW1141_simple.vcf.gz
# VulGry1_simple.vcf.gz CYW1141_simple.vcf.gz Dxy = 23982822 / (2 * 930737550) = 0.01288377

python ./calc_dxy.py VulGry1_simple.vcf.gz BGI_N323_simple.vcf.gz
# VulGry1_simple.vcf.gz BGI_N323_simple.vcf.gz Dxy = 29880274 / (2 * 820073349) = 0.01821805

python ./calc_dxy.py CYW1141_simple.vcf.gz CRW1112_simple.vcf.gz 
# CYW1141_simple.vcf.gz CRW1112_simple.vcf.gz Dxy = 2338458 / (2 * 948508593) = 0.00123270


################################################################################
# calc_dxy.py

import sys
import gzip

VCF1 = gzip.open(sys.argv[1], 'rt')
VCF2 = gzip.open(sys.argv[2], 'rt')

sitecount=0
dcount=0

for line1 in VCF1:
    line2=VCF2.readline()
    if line1.startswith('#'): continue
    line1=line1.strip().split('\t')
    line2=line2.strip().split('\t')
    if line1[6]!="PASS" or line2[6]!="PASS" : continue
    gt1=line1[9][:3]
    gt2=line2[9][:3]
    if gt1=='./.' or gt2=='./.' : continue
    sitecount+=1
    if gt1=='0/1' or gt2=='0/1' : dcount+=1 ; continue
    if gt1!=gt2: dcount+=2 ; continue

dxy=float(dcount)/(2*sitecount)

print("%s %s Dxy = %s / (2 * %s) = %s\n" % (sys.argv[1], sys.argv[2], dcount, sitecount, format(dxy, '.8f')))

VCF1.close()
VCF2.close()
exit()
