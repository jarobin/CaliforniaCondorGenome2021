# Generate condor-chicken dotplot with MashMap (https://github.com/marbl/MashMap)
# Fig. S1B

Q=gc_PacBio_HiC.fasta
R=galGal6_ordered.fa
OUT=mashmap_${R%.fa*}_${Q%.fa*}.out
mashmap -t 16 --perc_identity 95 --filter_mode one-to-one -r ${R} -q ${Q} -o ${OUT}
generateDotPlot png medium ${OUT}


# Note: manually edited out.gp for style (see below), then ran: gnuplot out.gp ; epstopdf out.eps

# out.gp
set term eps size 5,5 font "Arial,12"
set output "out.eps"
set bmargin 5
set xtics rotate ( \
 "chr1" 1.0, \
 "chr2" 197608386.0, \
 "chr3" 347290434.0, \
 "chr4" 458128851.0, \
 "chrZ" 549444095.0, \
 "chr5" 631974015.0, \
 "chr7" 691783112.0, \
 "chr6" 728525419.0, \
 "chr8" 764900119.0, \
 "chr9" 795119564.0, \
 "chr10" 819272649.0, \
 "chr12" 840392488.0, \
 "chr11" 860779765.0, \
 "chr13" 880979806.0, \
 "chr14" 900146519.0, \
 "chr20" 916365826.0, \
 "chr15" 930263112.0, \
 "chr18" 943325295.0, \
 "chr17" 954698434.0, \
 "chr19" 965460945.0, \
 "chr27" 975784156.0, \
 "chr21" 983864587.0, \
 "chr24" 990709565.0, \
 "chr23" 997200786.0, \
 "chr26" 1003350365.0, \
 "chr22" 1009406074.0, \
 "chr28" 1014865535.0, \
 "chr31" 1019982416.0, \
 "chrW" 1026135449.0, \
 "" 1032948590.0 \
)
set lmargin 15
set ytics ( \
 "*scaffold1" 1.0, \
 "scaffold32" 219006665.0, \
 "*scaffold2" 220977927.0, \
 "scaffold3" 389192215.0, \
 "*scaffold11" 517350168.0, \
 "*scaffold4" 542453957.0, \
 "*scaffold35" 627988238.0, \
 "scaffold5" 629100663.0, \
 "scaffold6" 714447480.0, \
 "scaffold7" 786722089.0, \
 "*scaffold8" 831753112.0, \
 "scaffold9" 873721020.0, \
 "scaffold10" 910805487.0, \
 "scaffold13" 940982696.0, \
 "*scaffold12" 965662022.0, \
 "*scaffold14" 990612078.0, \
 "*scaffold15" 1014420804.0, \
 "*scaffold16" 1036540231.0, \
 "scaffold17" 1056177006.0, \
 "scaffold18" 1073934376.0, \
 "*scaffold19" 1091080625.0, \
 "scaffold21" 1105373541.0, \
 "scaffold20" 1118429736.0, \
 "scaffold25" 1131569303.0, \
 "scaffold22" 1139408872.0, \
 "*scaffold24" 1149047815.0, \
 "scaffold23" 1157021302.0, \
 "*scaffold26" 1165319823.0, \
 "*scaffold28" 1172959752.0, \
 "scaffold27" 1178828917.0, \
 "" 1184958256.0 \
)
set size 1,1
set grid
unset key
set border 0
set tics scale 0 
set xlabel "Chicken (galGal6)" offset 0,0,0
set ylabel "California condor" offset 1,0,0
set format "%.0f"
set xrange [1.0:1032948590.0]
set yrange [1.0:1184958256.0]
set style line 1  lt 1 lw 1 pt 7 ps .5
set style line 2  lt 3 lw 1 pt 7 ps .5
set style line 3  lt 2 lw 1 pt 7 ps .5
plot \
 "out.fplot" title "FWD" w lp ls 1, \
 "out.rplot" title "REV" w lp ls 2


