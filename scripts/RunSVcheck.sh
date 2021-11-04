RDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"


#call is FINAL.vcf.gz <your human ref>.gff3
perl $RDIR/VCFtoSVbed.pl <(zcat $1) > $1.SV.bed

bedtools intersect -a $1.SV.bed -b $2 -wb > $1.intersect.out

bash $RDIR/processGFFintersect.sh $1.intersect.out 
