### iSMC recombination rate analysis 

### Use filtered VCF file (sites failing filters excluded, missing genotypes excluded, only scaffolds 1-4, 6-32 included)
VCF=CYW1141_simple_PASS_nomissing_scaffold1-32_autos.vcf.gz

### Generate .tab file
# Format is described here, from https://github.com/gvbarroso/iSMC/blob/master/doc/ismc_opt.bpp:

# // The 1st column is the block ID (eg chr1)
# // The 2nd is the start coordinate of the block mapped to the reference genome
# // The 3rd is the end coordinate of the block mapped to the reference genome
# // NOTE: start and end coordinates are relative to the block, meaning it restarts
# // e.g. at every chromosome. 
# // The 4th is 0 for all blocks and the 5th is the difference betwee 3rd & 2nd columns
# // The 6th and 7th columns are bottom and top cut-offs that can be convenient to 'synchronise'
# // coordinates when using ismc_mapper. This can happen when you want to consider a range of sites that
# // matches that of another file. For example, if your sequence data for chromosome 1 starts at position
# // 5000 and goes until position 250000000, and you want to 'synchronise' it with an experimental
# // genetic map that extends from 10000 to 230000000, the 6th and 7th columns of the first line of
# // your tab_file should be  9999 and 229999999.
# // Otherwise, if you want to include all sites, they should be the same as the 2nd and 3rd columns.

### Tab format uses start and end positions of chromosomes, based on sites that are in the VCF file

zcat ${VCF} | grep -v "^#" | \
awk '{ if (NR==1) {c=$1; s=$2; e=$2} else {if ($1!=c) {printf "%s\t%s\t%s\t0\t%s\t%s\t%s\n", c,s,e,e-s,s,e; c=$1; s=$2; e=$2} else {e=$2} } }END{printf "%s\t%s\t%s\t0\t%s\t%s\t%s\n", c,s,e,e-s,s,e}' \
> ${VCF}_ismc.tab


### Generate iSMC options file
# Most options left at default values (except for filenames, etc.)

cat > CYW1141_ismc.bpp <<EOF
dataset_label = CYW1141_ismc
input_file_type = cVCF 
sequence_file_path = CYW1141_simple_PASS_nomissing_scaffold1-32_autos.vcf.gz
seq_compression_type = gzip 
tab_file_path = CYW1141_simple_PASS_nomissing_scaffold1-32_autos.vcf.gz_ismc.tab
number_threads = 31
optimize = true
decode = true
decode_breakpoints_parallel = false
number_rho_categories = 5 
number_intervals = 40
function_tolerance = 1e-6
fragment_size = 2000000
EOF


### Run iSMC (~1 week)

ISMC_DIR=~/project/programs/ismc/src
PAR=CYW1141_ismc.bpp
LOG=CYW1141_ismc.log
/usr/bin/time -v  ${ISMC_DIR}/ismc params=${PAR} |& tee ${LOG}


### Generate iSMC mapping options file

cat > CYW1141_ismc.map.bpp <<EOF
dataset_label = CYW1141_ismc
tab_file_path = CYW1141_simple_PASS_nomissing_scaffold1-32_autos.vcf.gz_ismc.tab
bin_sizes = 1000,1000000 
bin_rate = rho 
tmrca = true 
EOF


### Run iSMC mapper (~5 min)

ISMC_DIR=~/project/programs/ismc/src
PAR=CYW1141_ismc.map.bpp
LOG=CYW1141_ismc.map.log
/usr/bin/time -v  ${ISMC_DIR}/ismc_mapper params=${PAR} |& tee ${LOG}


### Rename chromosomes in iSMC output

for DATA in CYW1141_ismc.rho.*.bedgraph ; do 
grep -v "^chrom" ${DATA} | sed 's/^chr//g' \
| awk '{if ($1>4) {printf "HiC_scaffold_%s\t%s\t%s\t%s\n", $1+1, $2, $3, $4} else {printf "HiC_scaffold_%s\t%s\t%s\t%s\n", $1, $2, $3, $4}}' > ${DATA}_scaffnames.bed
sort -k1,1 -k2,2n ${DATA}_scaffnames.bed > ${DATA}_scaffnames_sorted.bed
done

