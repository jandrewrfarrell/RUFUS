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

RDIR=/uufs/chpc.utah.edu/common/home/u0991464/d1/home/farrelac/bin/RUFUS_git/RUFUS/

RUFUSmodel=$RDIR/bin/ModelDist
RUFUSbuild=$RDIR/bin/RUFUS.Build 
RUFUSfilter=$RDIR/bin/RUFUS.Filter
RUFUSOverlap=$RDIR/scripts/OverlapBashMultiThread.sh
DeDupDump=$RDIR/scripts/HumanDedup.grenrator.tenplate


perl -ni -e 's/ /\t/;print' $MutantGenerator.Jhash.histo
perl -ni -e 's/ /\t/;print' $Parent1.Jhash.histo
perl -ni -e 's/ /\t/;print' $Parent2.Jhash.histo
perl -ni -e 's/ /\t/;print' $Parent3.Jhash.histo

echo "staring model"
if [ -e "$MutantGenerator.Jhash.histo.7.7.model" ]
then 
	echo "skipping model"
else
	/usr/bin/time -v $RUFUSmodel $MutantGenerator.Jhash.histo $K 150 $Threads > $Out.Run.out 
	/usr/bin/time -v $RUFUSmodel $Parent1.Jhash.histo $K 150 $Threads > $Out.Run.out
	/usr/bin/time -v $RUFUSmodel $Parent2.Jhash.histo $K 150 $Threads > $Out.Run.out
	/usr/bin/time -v $RUFUSmodel $Parent3.Jhash.histo $K 150 $Threads > $Out.Run.out
	echo "done with model "
fi

ParentMaxE=0 
MutantMinCov=$(head -2 $MutantGenerator.Jhash.histo.7.7.model | tail -1 )
echo "$ParentMaxE \n $MutantMinCov \n"


date
echo "starting RUFUS build "
let "Max= $MutantMinCov*100"
if [ -e "$Out".k$K"_m"$ParentMaxE"_c"$MutantMinCov".HashList.prefilter" ]
then 
	echo "Skipping build"
else
	/usr/bin/time -v $RUFUSbuild -c $Parent1.Jhash.sorted.min2.tab -c $Parent2.Jhash.sorted.min2.tab -c $Parent3.Jhash.sorted.min2.tab  -s $MutantGenerator.Jhash.sorted.min2.tab -o $Out".k$K"_m"$ParentMaxE"_c"$MutantMinCov".HashList.prefilter -hs $K -mS $MutantMinCov -mC $ParentMaxE  -max $Max -t 1  >> $Out.Run.out
fi 


if [ -e "$Out".k$K"_m"$ParentMaxE"_c"$MutantMinCov".HashList" ]
then
        echo "Skipping filter"
else

awk '{print $4 "\t" $3}' $Out".k$K"_m"$ParentMaxE"_c"$MutantMinCov".HashList.prefilter > $Out".k$K"_m"$ParentMaxE"_c"$MutantMinCov".HashList.prefilter.rearrange

cp $Out".k$K"_m"$ParentMaxE"_c"$MutantMinCov".HashList.prefilter $Out".k$K"_m"$ParentMaxE"_c"$MutantMinCov".HashList

#$RUFUSfilter -d ' ' -c /uufs/chpc.utah.edu/common/home/u0991464/lustr/RUFUS.1000g.reference/1000G.RUFUSreference.sorted.min45.tab  -s $Out".k$K"_m"$ParentMaxE"_c"$MutantMinCov".HashList.prefilter.rearrange  -o $Out".k$K"_m"$ParentMaxE"_c"$MutantMinCov".HashList -hs $k -mS $MutantMinCov -mC 0 -max $Max -t 20
fi



echo "done with RUFUS build "
echo "startin RUFUS filter"
rm  $MutantGenerator.temp
mkfifo $MutantGenerator.temp 
/usr/bin/time -v  bash $DeDupDump $MutantBam >  $MutantGenerator.temp &
/usr/bin/time -v   $RUFUSfilter  $Out".k$K"_m"$ParentMaxE"_c"$MutantMinCov".HashList $MutantGenerator.temp $Out".k$K"_m"$ParentMaxE"_c"$MutantMinCov".filtered.fq $K 0 5 10 $Threads >> $Out.Run.out   & 
wait
exit
echo "startin RUFUS overlap"
/usr/bin/time -v bash $RUFUSOverlap $Out".k$K"_m"$ParentMaxE"_c"$MutantMinCov".filtered.fq.Mutations.fastq 5 $Out".k$K"_m"$ParentMaxE"_c"$MutantMinCov" $Out".k$K"_m"$ParentMaxE"_c"$MutantMinCov".HashList $Threads 
echo "done with everything "
