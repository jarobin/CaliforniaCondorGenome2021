### Genome statistics (Table S1)
# - statswrapper.sh is part of the BBTools suite (http://sourceforge.net/projects/bbmap)
# - Output: filename, scaf_bp, gap_pct, gc_avg, n_scaffolds, scaf_max, scaf_L50, scaf_L90, scaf_N50, scaf_N90, n_contigs, ctg_max, ctg_L50, ctg_L90, ctg_N50, ctg_N90

statswrapper.sh n=100 in=$(ls *.fa* | tr '\n' ',' | sed 's/,$//g') > genome_stats.txt_temp
cat genome_stats.txt_temp | tr "_N_L" "_L_N" | sed 's/\/home\/jacqueline\/project\/condor\/assembly_stats\///g' | awk '{printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" , $20,$3,$5,$18,$1,$14,$6,$10,$7,$11,$2,$15,$8,$12,$9,$13}' > genome_stats.txt
rm genome_stats.txt_temp


### Plot cumulative genome proportion by number of scaffolds in R (Fig. S1A)
# Commands for getting input files:
# for i in *.fa* ; do faToTwoBit ${i} ${i%.fa*}.2bit ; done
# for i in *.2bit ; do twoBitInfo ${i} ${i%.2bit}.sizes ; done

library(RColorBrewer)

fs=list.files()
fss=fs[c(
grep(pattern="gc_PacBio_HiC", fs),
grep(pattern="ASM69994v1", fs),
grep(pattern="galGal6", fs),
grep(pattern="FicAlb1.5", fs),
grep(pattern="Taeniopygia_guttata", fs),
grep(pattern="Haliaeetus_leucocephalus", fs),
grep(pattern="Aquila_chrysaetos", fs),
grep(pattern="ASM69197v1", fs)
)]

pdf("assembly_comp.pdf", width=5, height=2.25, pointsize=8)

par(mar=c(4,4,1,3))

plot(0, 0, type="n", xlim=c(0,5.5), ylim=c(0,1), xlab="", ylab="", axes=F, xaxs="i")
title(xlab="Scaffold number", line=2.5)
title(ylab="Proportion of assembly (cumulative)", line=2.75)
axis(1, at=0:5, labels=c(1, 10, 100, "1,000", "10,000", "100,000"), las=1)
axis(2)

mycols=brewer.pal(n=8, name="Dark2")
mycols=mycols[c(1:5,7,6,8)]

for (i in 1:length(fss)){
	f=fss[i]
	df=read.table(f, header=F, sep="\t")
	colnames(df)=c("chr","len")
	df=df[order(df$len, decreasing=T),]
	yy=cumsum(df$len)/sum(df$len)
	xx=1:length(yy)
	lines(log(xx, base=10), yy, col=mycols[i], lwd=2)
}

legend("bottomright", c("CA condor", "turkey vulture", "chicken", "collared flycatcher", "zebra finch", "bald eagle", "golden eagle", "American crow"), col=mycols, bty="n", lwd=2)

dev.off()

