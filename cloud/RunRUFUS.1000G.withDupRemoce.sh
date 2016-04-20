date
Parent1Generator=$1
Parent2Generator=$2
SiblingGenerator=$3
ProbandGenerator=$4
K=$5
Threads=$6
Out=$7

echo "You gave
Parent1Generator=$1
Parent2Generator=$2
SiblingGenerator=$2
ProbandGenerator=$4
K=$5
Threads=$6
Out=$7
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
PullSampleHashes=$RDIR/cloud/C

perl -ni -e 's/ /\t/;print' $ProbandGenerator.Jhash.histo
perl -ni -e 's/ /\t/;print' $Parent1Generator.Jhash.histo
perl -ni -e 's/ /\t/;print' $Parent2Generator.Jhash.histo
perl -ni -e 's/ /\t/;print' $SiblingGenerator.Jhash.histo

echo "staring model"
if [ -e "$ProbandGenerator.Jhash.histo.7.7.model" ]
then
        echo "skipping model"
else
        /usr/bin/time -v $RUFUSmodel $ProbandGenerator.Jhash.histo $K 150 $Threads > $Out.Run.out
#        /usr/bin/time -v $RUFUSmodel $Parent1Generator.Jhash.histo $K 150 $Threads > $Out.Run.out
 #       /usr/bin/time -v $RUFUSmodel $Parent2Generator.Jhash.histo $K 150 $Threads > $Out.Run.out
  #      /usr/bin/time -v $RUFUSmodel $SiblingGenerator.Jhash.histo $K 150 $Threads > $Out.Run.out
        echo "done with model "
fi

ParentMaxE=0
MutantMinCov=$(head -2 $ProbandGenerator.Jhash.histo.7.7.model | tail -1 )
echo "$ParentMaxE \n $MutantMinCov \n"


date
echo "starting RUFUS build "
let "Max= $MutantMinCov*100"
if [ -e "Family.Unique.HashList.prefilter" ]
then
        echo "Skipping build"
else
	/usr/bin/time -v ../RUFUS/bin/jellyfish-MODIFIED/bin/jellyfish merge $Parent1Generator.Jhash  $Parent2Generator.Jhash $SiblingGenerator.Jhash $ProbandGenerator.Jhash >  Family.Unique.HashList.prefilter
fi

if [ -d "Family.Unique.HashList" ]
then 
 	echo "skipping 1kg filter"
else
	##add 1kg filter  
fi

if [ -e $ProbandGenerator.HashList ]
then 
	echo "skipping $ProbandGenerator.HashList pull "
else

	$PullSampleHashes  $ProbandGenerator.Jhash Family.Unique.HashList > $ProbandGenerator.HashList
fi 

if [ -e $SiblingGenerator.HashList ]
then 
	echo "skipping $SiblingGenerator.HashList pull"
else
	$PullSampleHashes $SiblingGenerator.Jhash Family.Unique.HashList > $SiblingGenerator.HashList 
fi 




echo "done with RUFUS build "
echo "startin RUFUS filter"
if [ -e $ProbandGenerator".k$K"_m"$ParentMaxE"_c"$MutantMinCov".Mutations.fq ]
rm  $ProbandGenerator.temp
mkfifo $ProbandGenerator.temp
/usr/bin/time -v  bash $ProbandGenerator >  $ProbandGenerator.temp &
/usr/bin/time -v   $RUFUSfilter  $Out".k$K"_m"$ParentMaxE"_c"$MutantMinCov".HashList $ProbandGenerator.temp $ProbandGenerator".k$K"_m"$ParentMaxE"_c"$MutantMinCov" $K 0 5 10 $Threads >> $Out.Run.out   &
wait

echo "startin RUFUS overlap"
/usr/bin/time -v bash $RUFUSOverlap $Out".k$K"_m"$ParentMaxE"_c"$MutantMinCov".filtered.fq.Mutations.fastq 5 $Out".k$K"_m"$ParentMaxE"_c"$MutantMinCov" $Out".k$K"_m"$ParentMaxE"_c"$MutantMinCov".HashList $Threads
echo "done with everything "

