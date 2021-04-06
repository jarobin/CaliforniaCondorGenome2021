# Condor-chicken Circos plot

cd ~/project/condor/circos_condor

mkdir 20201220
cd 20201220

### Get .sizes files (tab-delimited: chr length) and reciprocal best net from liftover
cp ~/project/condor/reference/c*/*/*sizes .
cp ~/project/condor/reference/condor/gc_PacBio_HiC/liftover/gc_PacBio_HiC_to_galGal6/gc_PacBio_HiC.galGal6.rbest.net.gz .

### Make karyotype files
# Chicken karyotype, keep chromosomes only
while read NAME SIZE ; do COL=$(echo col_gg_${NAME} | tr -s '[:upper:]' '[:lower:]') ; echo "chr - gg_${NAME} ${NAME#chr} 0 ${SIZE} ${COL}" ; done < galGal6.sizes > galGal6.karyotype_raw
grep -v "chrUn_" galGal6.karyotype_raw | grep -v "_random" | grep -v M | grep -v Z | grep -v W | sort -V > galGal6.karyotype_temp
grep -v "chrUn_" galGal6.karyotype_raw | grep -v "_random" | grep Z >> galGal6.karyotype_temp
grep -v "chrUn_" galGal6.karyotype_raw | grep -v "_random" | grep W >> galGal6.karyotype_temp
sed 's/col_g_/col_gg_/g' galGal6.karyotype_temp > galGal6.karyotype

# Condor karyotype
# Need to rename chromosomes first (scaffold names in net and sizes files)
cat > rename_chr_1.py <<EOF
import sys

infile = open(sys.argv[1], 'r')

old2new={
'HiC_scaffold_1':'chr1',
'HiC_scaffold_2':'chr2',
'HiC_scaffold_3':'chr3',
'HiC_scaffold_4':'chr4',
'HiC_scaffold_6':'chr5',
'HiC_scaffold_8':'chr6',
'HiC_scaffold_7':'chr7',
'HiC_scaffold_9':'chr8',
'HiC_scaffold_11':'chr9',
'HiC_scaffold_10':'chr10',
'HiC_scaffold_13':'chr11',
'HiC_scaffold_14':'chr12',
'HiC_scaffold_12':'chr13',
'HiC_scaffold_15':'chr14',
'HiC_scaffold_16':'chr15',
'HiC_scaffold_18':'chr16',
'HiC_scaffold_17':'chr17',
'HiC_scaffold_21':'chr18',
'HiC_scaffold_19':'chr19',
'HiC_scaffold_20':'chr20',
'HiC_scaffold_22':'chr21',
'HiC_scaffold_23':'chr22',
'HiC_scaffold_28':'chr23',
'HiC_scaffold_27':'chr24',
'HiC_scaffold_24':'chr25',
'HiC_scaffold_29':'chr26',
'HiC_scaffold_26':'chr27',
'HiC_scaffold_25':'chr28',
'HiC_scaffold_30':'chr29',
'HiC_scaffold_5':'chrZ'
}

for line in infile:
    line_list=line.strip().split('\t')
    old=line_list[0]
    if old in old2new.keys():
        new=old2new[old]
        newline=line.replace(old, new)
        sys.stdout.write(newline)

infile.close()
exit()

EOF

chmod u+x rename_chr_1.py
python ./rename_chr_1.py gc_PacBio_HiC.sizes > gc_PacBio_HiC_rename.sizes

for i in {1..29} Z ; do 
awk -v var=chr${i} '$1==var' gc_PacBio_HiC_rename.sizes >> gc_PacBio_HiC_rename_sort.sizes
done

while read NAME SIZE ; do 
COL=$(echo col_gc_${NAME} | tr -s '[:upper:]' '[:lower:]') ; echo "chr - gc_${NAME} ${NAME#chr} 0 ${SIZE} ${COL}"
done > gc_PacBio_HiC.karyotype < gc_PacBio_HiC_rename_sort.sizes

### Make .conf file for chromosome order
cut -d' ' -f3 galGal6.karyotype | tr '\n' ',' > galGal6.karyotype_chrlist_temp
cut -d' ' -f3 gc_PacBio_HiC.karyotype | tac | tr '\n' ',' | rev | cut -d',' -f2- | rev > gc_PacBio_HiC.karyotype_chrlist_temp
echo "chromosomes_order = "$(cat galGal6.karyotype_chrlist_temp gc_PacBio_HiC.karyotype_chrlist_temp) > segment.order.conf

### Make data files for chromosome/scaffold labels
KARYO=galGal6.karyotype
cut -d' ' -f3,6 ${KARYO} | awk -v var="gg_chr" '{gsub(var,"",$1) ; printf "%s%s\t%.0f\t%.0f\t%s\n", var, $1, $2/2-1, $2/2, $1}' | sed -n 'p;n' > ${KARYO}.labels1
cut -d' ' -f3,6 ${KARYO} | awk -v var="gg_chr" '{gsub(var,"",$1) ; printf "%s%s\t%.0f\t%.0f\t%s\n", var, $1, $2/2-1, $2/2, $1}' | sed -n 'n;p' > ${KARYO}.labels2

KARYO=gc_PacBio_HiC.karyotype
cut -d' ' -f3,6 ${KARYO} | awk -v var="gc_chr" '{gsub(var,"",$1) ; printf "%s%s\t%.0f\t%.0f\t%s\n", var, $1, $2/2-1, $2/2, $1}' | sed -n 'p;n' > ${KARYO}.labels1
cut -d' ' -f3,6 ${KARYO} | awk -v var="gc_chr" '{gsub(var,"",$1) ; printf "%s%s\t%.0f\t%.0f\t%s\n", var, $1, $2/2-1, $2/2, $1}' | sed -n 'n;p' > ${KARYO}.labels2

### Make links from net file(s)
# Extract top-level chains
# If orientation is reversed (column 5 is - instead of +), print end coordinate first and then start
zcat gc_PacBio_HiC.galGal6.rbest.net.gz | awk '{if (/net/) {tChrom=$2} ; if (/top/) { if ($5=="+") {printf "%s\t%s\t%s\t%s\t%s\t%s\n", tChrom, $2, $2+$3, $4, $6, $6+$7} else {printf "%s\t%s\t%s\t%s\t%s\t%s\n", tChrom, $2, $2+$3, $4, $6+$7, $6} }}' > gc_PacBio_HiC.galGal6.rbest.net.links_temp
python ./rename_chr_1.py gc_PacBio_HiC.galGal6.rbest.net.links_temp | awk '{printf "gc_%s\t%s\t%s\tgg_%s\t%s\t%s\n", $1, $2, $3, $4, $5, $6}' > gc_PacBio_HiC.galGal6.rbest.net.links_temp_rename

# Must be >=1 Mb in both condor and chicken
awk '$3-$2>=1000000 { if ($6-$5>0 && $6-$5>=1000000) {printf "%s\tcolor=col_%s_link\n", $0, tolower($1)} else { if ($6-$5<=-1000000) {printf "%s\tcolor=col_%s_link\n", $0, tolower($1)} } }' gc_PacBio_HiC.galGal6.rbest.net.links_temp_rename > gc_PacBio_HiC.galGal6.rbest.net.links_1Mb_tcol

# Define my colors
cat>cols.temp<<EOF
193,33,31
221,58,37
227,95,44
229,123,50
229,146,54
223,164,58
209,180,64
185,188,74
154,189,93
125,184,119
105,177,143
94,169,160
86,160,174
78,152,185
77,147,189
77,146,191
77,144,192
77,141,195
76,139,197
76,136,197
77,133,197
77,129,197
77,127,197
78,123,197
77,120,196
80,116,193
81,114,191
82,110,189
83,107,186
84,103,184
87,101,180
89,97,178
92,95,175
94,93,172
EOF

### Make .conf file for chicken chr colors (transparency=0.6)
KARYO=galGal6.karyotype
cut -d' ' -f7 ${KARYO} | paste -d'=' - cols.temp | grep "col_" | sed 's/$/,0.6/g' > col_${KARYO}.conf

### Make .conf file for condor scaff colors
KARYO=gc_PacBio_HiC.karyotype
cut -d' ' -f7 ${KARYO} | paste -d'=' - cols.temp | grep "col_" | sed 's/$/,0.6/g' | grep -v "chrz" > col_${KARYO}.conf
grep "chrz" col_galGal6.karyotype.conf | sed 's/gg/gc/g' >> col_${KARYO}.conf

### Make .conf file for link colors (with transparency, as above)
# Here, KARYO determines whether link colors match condor or chicken chromosome color
# KARYO=galGal6.karyotype
KARYO=gc_PacBio_HiC.karyotype
sed 's/=/_link=/g' col_${KARYO}.conf > col_links.conf

### Execute Circos (see circos_condor.conf below)
circos -conf circos_condor.conf


################################################################################
# circos_condor.conf

# Condor Circos plot

# Chromosome name, size and color definition
karyotype = galGal6.karyotype,gc_PacBio_HiC.karyotype

# Define segment order and reverse the order of condor chromosomes
<<include segment.order.conf>>
chromosomes_reverse = /gc_chr*/

# The <ideogram> block defines the position, size, labels and other
# properties of the segments on which data are drawn. 

<ideogram>

# <<include ideogram.label.conf>>
show_label       = no
label_radius     = dims(ideogram,radius) + 0.05r
label_size       = 35
label_parallel   = yes

label_font       = light

# Radial position within the image of the ideograms.
radius           = 0.80r

# Thickness of ideograms, which can be absolute (e.g. pixels, "p"
# suffix) or relative ("r" suffix). 
thickness        = 50p

# Ideograms can be drawn as filled, outlined, or both. When filled,
# the color will be taken from the last field in the karyotype file,
# or set by chromosomes_colors. When stroke_thickness=0p or if the 
# parameter is missing, the ideogram has no outline and the value of 
# stroke_color is not used.
fill             = yes  
stroke_color     = dgrey
stroke_thickness = 2p   

<spacing>
# Spacing between ideograms. 
default = 0.004r

# Leave extra space in between condor and chicken ideograms at top and bottom
<pairwise gc_chr1 gg_chr1>
spacing = 30r
</pairwise>

<pairwise gg_chrW gc_chrZ>
spacing = 30r
</pairwise>

</spacing>

</ideogram>

# Links
<links>

<link>

ribbon = yes
file   = gc_PacBio_HiC.galGal6.rbest.net.links_1Mb_tcol
radius        = 0.99r
bezier_radius = 0r
stroke_color     = dgrey
stroke_thickness = 2

</link>

</links>


<image>

# Included from Circos distribution.
# <<include image.conf>>  
# image.conf just includes image.generic.conf and background.white.conf, contents copied here:
dir   = . 
file  = condor_circos.png
png   = yes
svg   = no

# radius of inscribed circle in image
radius         = 1500p

auto_alpha_colors = yes
auto_alpha_steps  = 5

background = white

# Override angle_offset defined in etc/image.conf
# Between -70 and -80
angle_offset* = -75
              
</image>

<plots>

<plot>
type = text
file = galGal6.karyotype.labels1
r0 = 1.02r
r1 = 1.5r
label_size = 55p
color = black
label_font* = sans_serif
label_rotate = no
label_font = condensed
padding = 0p
rpadding = 0p
label_snuggle = no
</plot>

<plot>
type = text
file = galGal6.karyotype.labels2
r0 = 1.1r
r1 = 1.5r
label_size = 55p
color = black
label_font* = sans_serif
label_rotate = no
label_font = condensed
padding = 0p
rpadding = 0p
label_snuggle = no
</plot>

<plot>
type = text
file = gc_PacBio_HiC.karyotype.labels1
r0 = 1.02r
r1 = 1.5r
label_size = 55p
color = black
label_font* = sans_serif
label_rotate = yes
label_font = condensed
padding = 0p
rpadding = 0p
label_snuggle = no
</plot>

<plot>
type = text
label_rotate = yes
file = gc_PacBio_HiC.karyotype.labels2
r0 = 1.1r
r1 = 1.5r
label_size = 55p
color = black
label_font* = sans_serif
label_font = condensed
padding = 0p
rpadding = 0p
label_snuggle = no
</plot>

</plots>


# RGB/HSV/LCH color definitions, color lists, location of fonts, fill
# patterns. Included from Circos distribution.

<<include etc/colors_fonts_patterns.conf>> 

<colors>
<<include etc/colors.conf>>

# My custom colors
<<include col_galGal6.karyotype.conf>>
<<include col_gc_PacBio_HiC.karyotype.conf>>
<<include col_links.conf>>

</colors>


<fonts>
<<include etc/fonts.conf>>
</fonts>

<patterns>
<<include etc/patterns.conf>>
</patterns>


# Debugging, I/O and other system parameters
# Included from Circos distribution.
<<include etc/housekeeping.conf>> 
track_defaults* = undef


