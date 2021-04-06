# Genome statistics (given in Table S1)
# - statswrapper.sh is part of the BBTools suite (http://sourceforge.net/projects/bbmap)
# - Output: filename, scaf_bp, gap_pct, gc_avg, n_scaffolds, scaf_max, scaf_L50, scaf_L90, scaf_N50, scaf_N90, n_contigs, ctg_max, ctg_L50, ctg_L90, ctg_N50, ctg_N90

statswrapper.sh n=100 in=$(ls *.fa* | tr '\n' ',' | sed 's/,$//g') > genome_stats.txt_temp
cat genome_stats.txt_temp | tr "_N_L" "_L_N" | sed 's/\/home\/jacqueline\/project\/condor\/assembly_stats\///g' | awk '{printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" , $20,$3,$5,$18,$1,$14,$6,$10,$7,$11,$2,$15,$8,$12,$9,$13}' > genome_stats.txt
rm genome_stats.txt_temp
