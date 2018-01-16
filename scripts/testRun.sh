#!/bin/bash
args=("$@")

echo "attempting to print out command line args\n"

echo "Values of \"\$@\":"
for arg in "$@"
do
  echo "Arg #$cnt= $arg"
  let "cnt+=1"
done


numArgs=$#

Parents=("${args[@]:0:$((numArgs-4))}")
Parents="${@:0:numArgs-3}"


echo "Parents: "
echo $Parents


ProbandGenerator=${args[$((numArgs-4))]}

echo "ProbandGenerator: "
echo $ProbandGenerator

K=${args[$((numArgs-3))]}
Threads=${args[$((numArgs-2))]}
Out=${args[$((numArgs-1))]}

echo "K: "
echo $K

echo "Threads: "
echo $Threads

echo "Out:"

echo $Out

K=$args[$(($numArgs-3))]
Threads=$args[$(($numArgs-2))]
Out=$args[$(($numArgs-1))]

parentsString=""
space=" "
jhash=".Jhash"

for parent in "${parents[@]}"
do
 parentsString=$parentsString$space$parent$jhash
 echo "parents string equals " $parentsString
done


echo "you gave"

for arg in "${args[@]}"
do
  echo $arg
done
if [ -z "$Out" ]
then
        echo "out file not specified"
        exit
fi
RDIR=/uufs/chpc.utah.edu/common/home/u0991464/d1/home/farrelac/RUFUS
RDIR=$PWD
RUFUSmodel=$RDIR/bin/ModelDist
RUFUSbuild=$RDIR/bin/RUFUS.Build
RUFUSfilter=$RDIR/bin/RUFUS.Filter
RUFUSOverlap=$RDIR/scripts/OverlapBashMultiThread.trio.sh
DeDupDump=$RDIR/scripts/HumanDedup.grenrator.tenplate
PullSampleHashes=$RDIR/cloud/CheckJellyHashList.sh
RUFUS1kgFilter=$RDIR/bin/RUFUS.1kg.filter
RunJelly=$RDIR/cloud/RunJellyForRUFUS
for parent in "${parents[@]}"
do
    /usr/bin/time -v bash $RunJelly $parent $K $(echo $Threads -2 | bc) 2 &
done
/usr/bin/time -v bash $RunJelly $ProbandGenerator $K $(echo $Threads -2 | bc) 2 &
wait
perl -ni -e 's/ /\t/;print' $ProbandGenerator.Jhash.histo
for parent in "${parents[@]}"
do
  perl -ni -e 's/ /\t/;print' $parent.Jhash.histo
done
if [ -e "$ProbandGenerator.Jhash.histo.7.7.model" ]
then
        echo "skipping model"
else
    echo "staring model"
        /usr/bin/time -v $RUFUSmodel $ProbandGenerator.Jhash.histo $K 150 $Threads
        echo "done with model"
fi
ParentMaxE=1
MutantMinCov=$(head -2 $ProbandGenerator.Jhash.histo.7.7.model | tail -1 )

echo "starting RUFUS build"
let "Max= $MutantMinCov*100"
if [ -e "$Out.Family.Unique.HashList" ]
then
        echo "Skipping build"
else
    /usr/bin/time -v $RDIR/cloud/jellyfish-MODIFIED-merge/bin/jellyfish merge $parentsString $ProbandGenerator.Jhash >  $Out.Family.Unique.HashList
fi
echo "Mut cov = $MutantMinCov"
if [ -e $ProbandGenerator.k$K_c$MutantMinCov.HashList ]
then
    echo "skipping $ProbandGenerator.HashList pull "
else
    /usr/bin/time -v bash $PullSampleHashes $ProbandGenerator.Jhash $Out.Family.Unique.HashList $MutantMinCov > $ProbandGenerator.k$K_c$MutantMinCov.HashList
fi
echo "done with RUFUS build "
echo "starting RUFUS filter"
if [ -e $ProbandGenerator.Mutations.fastq ]
then
    echo "skipping filter"
else
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
    echo "starting RUFUS overlap"
    /usr/bin/time -v bash $RUFUSOverlap $ProbandGenerator.Mutations.fastq 5 $ProbandGenerator $ProbandGenerator.k$MutantMinCov.HashList $K $Threads $ProbandGenerator.Jhash $parentsString
fi
echo "done with everything"
