date
Parent1Generator=$1
Parent2Generator=$2
ProbandGenerator=$3
K=$4
Threads=$5
Out=$6

echo "You gave
Parent1Generator=$1
Parent2Generator=$2
ProbandGenerator=$3
K=$4
Threads=$5
Out=$6
"

if [ -z "$Out" ]
then
        echo "out file not specified"
        exit
fi

RDIR=/uufs/chpc.utah.edu/common/home/u0991464/d1/home/farrelac/RUFUS
RUFUSmodel=$RDIR/bin/ModelDist
RUFUSbuild=$RDIR/bin/RUFUS.Build
RUFUSfilter=$RDIR/bin/RUFUS.Filter
RUFUSOverlap=$RDIR/scripts/OverlapBashMultiThread.trio.sh
DeDupDump=$RDIR/scripts/HumanDedup.grenrator.tenplate
PullSampleHashes=$RDIR/cloud/CheckJellyHashList.sh
RUFUS1kgFilter=$RDIR/bin/RUFUS.1kg.filter
RunJelly=$RDIR/cloud/RunJellyForRUFUS.UGP



/usr/bin/time -v bash $RunJelly $Parent1Generator $K $(echo $Threads -2 | bc) 2
/usr/bin/time -v bash $RunJelly $Parent2Generator $K $(echo $Threads -2 | bc) 2
/usr/bin/time -v bash $RunJelly $ProbandGenerator $K $(echo $Threads -2 | bc) 2

perl -ni -e 's/ /\t/;print' $ProbandGenerator.Jhash.histo
perl -ni -e 's/ /\t/;print' $Parent1Generator.Jhash.histo
perl -ni -e 's/ /\t/;print' $Parent2Generator.Jhash.histo

if [ -e "$ProbandGenerator.Jhash.histo.7.7.model" ]
then
        echo "skipping model"
else
	echo "staring model"
        /usr/bin/time -v $RUFUSmodel $ProbandGenerator.Jhash.histo $K 150 $Threads 
        echo "done with model "
fi

ParentMaxE=0
MutantMinCov=$(head -2 $ProbandGenerator.Jhash.histo.7.7.model | tail -1 )

date
echo "starting RUFUS build "
let "Max= $MutantMinCov*100"
if [ -s "Family.Unique.HashList" ]
then
        echo "Skipping build"
else
	/usr/bin/time -v $RDIR/cloud/jellyfish-MODIFIED-merge/bin/jellyfish merge $Parent1Generator.Jhash  $Parent2Generator.Jhash $ProbandGenerator.Jhash >  Family.Unique.HashList
fi

echo "Mut cov = $MutantMinCov "
if [ -s $ProbandGenerator.k$K_c$MutantMinCov.HashList.prefilter ]
then 
	echo "skipping $ProbandGenerator.HashList pull "
else
	/usr/bin/time -v bash $PullSampleHashes $ProbandGenerator.Jhash Family.Unique.HashList $MutantMinCov > $ProbandGenerator.k$K_c$MutantMinCov.HashList.prefilter
fi 

if [ -s $ProbandGenerator.k$K_c$MutantMinCov.HashList ]
then 
	echo "skipping 1kg filter"
else
	 /usr/bin/time -v $RDIR/cloud/RUFUS.search.1kg -hf <(awk '{print $1 "\t" $2}' $ProbandGenerator.k$K_c$MutantMinCov.HashList.prefilter ) -o $ProbandGenerator.k$K_c$MutantMinCov.HashList  -c /scratch/ucgd/lustre/work/u0991464/RUFUS.1000g.reference/1000G.RUFUSreference.sorted.min45.tab -hs 25
fi 

echo "done with RUFUS build "

echo "startin RUFUS filter"
if [ -s $ProbandGenerator.Mutations.fastq ]
then 
	echo "skipping filter"
else 
echo "crap"
	rm  $ProbandGenerator.temp
	mkfifo $ProbandGenerator.temp
	/usr/bin/time -v  bash $ProbandGenerator | $RDIR/cloud/PassThroughSamCheck $ProbandGenerator.filter.chr >  $ProbandGenerator.temp &
	/usr/bin/time -v   $RUFUSfilter  $ProbandGenerator.k$K_c$MutantMinCov.HashList $ProbandGenerator.temp $ProbandGenerator $K 5 5 10 $(echo $Threads -2 | bc) &
	wait

fi 

if [ -e $ProbandGenerator.V2.overlap.hashcount.fastq.bam.vcf ]
then 
	echo "skipping overlap"
else
	echo "startin RUFUS overlap"
	/usr/bin/time -v bash $RUFUSOverlap $ProbandGenerator.Mutations.fastq 5 $ProbandGenerator $ProbandGenerator.k$MutantMinCov.HashList $Threads $ProbandGenerator.Jhash $Parent1Generator.Jhash $Parent2Generator.Jhash 
fi


echo "done with everything "



