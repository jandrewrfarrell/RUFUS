date
ParentGenerator=$1
MutantGenerator=$2
K=$3
Threads=$4
Out=$5

echo "You gave
ParentGenerator=$1
MutantGenerator=$2
K=$3
Threads=$4
Out=$5
"

if [ -z "$Out" ]
then
        echo "out file not specified"
        exit
fi

RDIR=/home/ubuntu/work/RUFUS
RUFUSmodel=$RDIR/bin/ModelDist
RUFUSbuild=$RDIR/bin/RUFUS.Build
RUFUSfilter=$RDIR/bin/RUFUS.Filter
RUFUSOverlap=$RDIR/scripts/OverlapBashMultiThread.sh
DeDupDump=$RDIR/scripts/HumanDedup.grenrator.tenplate


perl -ni -e 's/ /\t/;print' $MutantGenerator.Jhash.histo
perl -ni -e 's/ /\t/;print' $ParentGenerator.Jhash.histo

echo "staring model"
if [ -e "$MutantGenerator.Jhash.histo.7.7.model" ]
then
        echo "skipping model"
else
        /usr/bin/time -v $RUFUSmodel $MutantGenerator.Jhash.histo $K 150 $Threads > $Out.Run.out
#        /usr/bin/time -v $RUFUSmodel $ParentGenerator.Jhash.histo $K 150 $Threads > $Out.Run.out
 #       /usr/bin/time -v $RUFUSmodel $Parent2.Jhash.histo $K 150 $Threads > $Out.Run.out
  #      /usr/bin/time -v $RUFUSmodel $Parent3.Jhash.histo $K 150 $Threads > $Out.Run.out
        echo "done with model "
fi

ParentMaxE=0
MutantMinCov=$(head -2 $MutantGenerator.Jhash.histo.7.7.model | tail -1 )
echo "$ParentMaxE \n $MutantMinCov \n"


date
echo "starting RUFUS build "
let "Max= $MutantMinCov*100"
if [ -e "$Out".k$K"_m"$ParentMaxE"_c"$MutantMinCov".HashList" ]
then
        echo "Skipping build"
else
        /usr/bin/time -v $RUFUSbuild  -c $ParentGenerator.Jhash.sorted.min2.tab -s $MutantGenerator.Jhash.sorted.min2.tab -o $Out".k$K"_m"$ParentMaxE"_c"$MutantMinCov".HashList -hs $K -mS $MutantMinCov -mC $ParentMaxE  -max $Max -t 1  >> $Out.Run.out
        #/usr/bin/time -v $RUFUSbuild -c <(s3cmd  get --no-progress s3://rufus.marth.lab/1000G.RUFUSreference.sorted.min45.tab.gz - | zcat) -c $ParentGenerator.Jhash.sorted.min2.tab -c $Parent2.Jhash.sorted.min2.tab -c $Parent3.Jhash.sorted.min2.tab  -s $MutantGenerator.Jhash.sorted.min2.tab -o $Out".k$K"_m"$ParentMaxE"_c"$MutantMinCov".HashList -hs $K -mS $MutantMinCov -mC $ParentMaxE  -max $Max -t 1  >> $Out.Run.out
fi




echo "done with RUFUS build "
echo "startin RUFUS filter"
rm  $MutantGenerator.temp
mkfifo $MutantGenerator.temp
/usr/bin/time -v  bash $MutantGenerator >  $MutantGenerator.temp &
/usr/bin/time -v   $RUFUSfilter  $Out".k$K"_m"$ParentMaxE"_c"$MutantMinCov".HashList $MutantGenerator.temp $Out".k$K"_m"$ParentMaxE"_c"$MutantMinCov".filtered.fq $K 0 5 10 $Threads >> $Out.Run.out   &
wait

echo "startin RUFUS overlap"
/usr/bin/time -v bash $RUFUSOverlap $Out".k$K"_m"$ParentMaxE"_c"$MutantMinCov".filtered.fq.Mutations.fastq 5 $Out".k$K"_m"$ParentMaxE"_c"$MutantMinCov" $Out".k$K"_m"$ParentMaxE"_c"$MutantMinCov".HashList $Threads
echo "done with everything "

