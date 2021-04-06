# BUSCO
# Used v4.0.5 and v4.0.6 through conda
# Dataset used for analysis: https://busco-data.ezlab.org/v4/data/lineages/aves_odb10.2019-11-20.tar.gz

export FASTA=galGal6.fa
export LINEAGE=aves_odb10
busco --cpu 8 -m genome --augustus_species=chicken -i ${FASTA} -o busco_${FASTA%.fa*} -l ${LINEAGE}
