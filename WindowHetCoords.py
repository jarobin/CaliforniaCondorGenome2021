# Script to count number of called genotypes and number of heterozygotes per sample in 
# regions given by a list of coordinates (1-based, inclusive of the start and end positions).
#
# Example:
# HiC_scaffold_1	4652	1004651
# HiC_scaffold_1	1004652	2004651
# HiC_scaffold_1	2004652	3004651
# HiC_scaffold_1	3004652	4004651
# HiC_scaffold_1	4004652	5004651
#
# Used to get heterozygosity in windows with the same coordinates as those from iSMC for 
# computing correlations shown in Fig. 3F.
#
# Executed with Python2.7. Not sure if it works with Python 3+.
#
# Input file is a single- or multi-sample VCF file that has been filtered (passing sites 
# have "PASS" in the FILTER column) and compressed with gzip/bgzip. Sites without PASS in
# filter column are skipped.
#
# Usage: 
# python ./WindowHetCoords.py [vcf] [coordinate_list] [outfile]
#
# Example: 
# python ./WindowHetCoords.py input.vcf.gz my_coords.list output_het_coords.txt

import sys
import pysam
import os
import gzip


# Open input file and make sure the VCF file is indexed (if not, create index)
filename = sys.argv[1]
VCF = gzip.open(filename, 'r')

if not os.path.exists("%s.tbi" % filename):
    pysam.tabix_index(filename, preset="vcf")
parsevcf = pysam.Tabixfile(filename)


# Set variables
coords=open(sys.argv[2], 'r')


# Get list of samples from VCF file header
samples=[]
for line in VCF:
    if line.startswith('##'):
        pass
    else:
        for i in line.split()[9:]: samples.append(i)
        break


# Create output file
output = open(sys.argv[3], 'w')
output.write('chrom\twindow_start\tsites_total\tcalls_%s\thets_%s\n' % ('\tcalls_'.join(samples), '\thets_'.join(samples)) )


# Fetch a region, ignore sites that fail filters, tally genotype calls and heterozygotes        
def snp_cal(chrom,window_start,window_end):
    print("%s:%s-%s" % (chrom, window_start, window_end))
    rows = tuple(parsevcf.fetch(region="%s:%s-%s" % (chrom, window_start, window_end), parser=pysam.asTuple()))    
    sites_total=0
    if rows==():
        calls=['NA']*len(samples)
        hets=['NA']*len(samples)
        output.write('%s\t%s\t%s\t%s\t%s\n' % (chrom, window_start, sites_total, '\t'.join(map(str,calls)), '\t'.join(map(str,hets))) )
    else:
        calls=[0]*len(samples)
        hets=[0]*len(samples)
        for line in rows:
            if line[6]!="PASS": continue
            sites_total+=1
            for i in range(0,len(samples)):
                if line[i+9][:1]=='.': continue
                calls[i]+=1
                GT=line[i+9].split(':')[0]
                if '/' in GT: sp='/'
                if '|' in GT: sp='|'
                if GT.split(sp)[0]!=GT.split(sp)[1]: hets[i]+=1
        output.write('%s\t%s\t%s\t%s\t%s\n' % (chrom, window_start, sites_total, '\t'.join(map(str,calls)), '\t'.join(map(str,hets))) )


# Get tallies for each region in coordinates list
for line in coords:
    line=line.strip().split('\t')
    chrom=line[0]
    window_start=line[1]
    window_end=line[2]
    snp_cal(chrom,window_start,window_end)


# Close files and exit
VCF.close()
coords.close()
output.close()

exit()

