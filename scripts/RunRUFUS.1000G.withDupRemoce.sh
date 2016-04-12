date
Parent1=$1
Parent2=$2
Parent3=$3
MutantGenerator=$4
MutantBam=$5
K=$6
Threads=$7
Out=$8

echo "You gave
Parent1=$1
Parent2=$2
Parent3=$2
MutantGenerator=$4
MutantBam=$5
K=$6
Threads=$7
Out=$8
"

if [ -z "$Out" ]
then 
	echo "out file not specified"
	exit
fi

RUFUSmodel=~/d1/home/farrelac/bin/RUFUSAwesome/ModelDist7.7.withDip
RUFUSbuild=~/d1/home/farrelac/bin/RUFUSAwesome/RUFUS.Build.8.cutoff 
RUFUSfilter=~/d1/home/farrelac/bin/RUFUS_git/RUFUS/bin/RUFUSv5.Filter
RUFUSOverlap=~/d1/home/farrelac/bin/RUFUSAwesome/OverlapBashMultiThread.sh
samtools=~/bin/samtools-1.2/samtools
RUFUSfilter2=/uufs/chpc.utah.edu/common/home/u0991464/d1/home/farrelac/bin/RUFUSAwesome/RUFUSv6.FilterHash.memMap.multiTrhead

perl -ni -e 's/ /\t/;print' $MutantGenerator.histo
perl -ni -e 's/ /\t/;print' $Parent1.histo
perl -ni -e 's/ /\t/;print' $Parent2.histo
perl -ni -e 's/ /\t/;print' $Parent3.histo

ParentMaxE=0 
MutantMinCov=$(head -2 $MutantGenerator.histo.7.7.model | tail -1 )
echo "$ParentMaxE \n $MutantMinCov \n"


date

echo "$samtools view $Out".k$K"_m"$ParentMaxE"_c"$MutantMinCov".overlap.hashcount.fastq.bam | ~/d1/home/farrelac/bin/RUFUSAwesome/RUFUS.interpret -r ~/d1/data/project_fasta/human_g1k_v37_decoy_phix.fasta -hf $Out".k$K"_m"$ParentMaxE"_c"$MutantMinCov".HashList -o $Out".k$K"_m"$ParentMaxE"_c"$MutantMinCov".Interpret2 -c $MutantGenerator.sorted.min2.tab -c $Parent1.sorted.min2.tab -c $Parent2.sorted.min2.tab -c $Parent3.sorted.min2.tab #> $Out".k$K"_m"$ParentMaxE"_c"$MutantMinCov".Interpret2.out "

$samtools view $Out".k$K"_m"$ParentMaxE"_c"$MutantMinCov".overlap.hashcount.fastq.bam | ~/d1/home/farrelac/bin/RUFUSAwesome/RUFUS.interpret -r ~/d1/data/project_fasta/human_g1k_v37_decoy_phix.fasta -hf $Out".k$K"_m"$ParentMaxE"_c"$MutantMinCov".HashList -o $Out".k$K"_m"$ParentMaxE"_c"$MutantMinCov".Interpret2 -c $MutantGenerator.sorted.min2.tab -c $Parent1.sorted.min2.tab -c $Parent2.sorted.min2.tab -c $Parent3.sorted.min2.tab #> $Out".k$K"_m"$ParentMaxE"_c"$MutantMinCov".Interpret2.out  
 
