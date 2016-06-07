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
PullSampleHashes=$RDIR/cloud/CheckJellyHashList.sh
RUFUS1kgFilter=$RDIR/bin/RUFUS.1kg.filter
RunJelly=$RDIR/cloud/RunJellyForRUFUS


aws s3 --region us-east-1 cp s3://marthlab.rufus/ASC.scripts/$Parent1Generator ./
aws s3 --region us-east-1 cp s3://marthlab.rufus/ASC.scripts/$Parent2Generator ./
aws s3 --region us-east-1 cp s3://marthlab.rufus/ASC.scripts/$SiblingGenerator ./
aws s3 --region us-east-1 cp s3://marthlab.rufus/ASC.scripts/$ProbandGenerator ./

/usr/bin/time -v bash $RunJelly $Parent1Generator $K $(echo $Threads -2 | bc)
if [ "$( tail -n 2 $Parent1Generator.Jelly.chr | head -1)" == "*" ]
then 
	echo "Jelly on $Parent1Generator successfull"
else
	echo "ReRunning jelly on $Parent1Generator"
	rm $Parent1Generator.Jhash
	rm $Parent1Generator.Jelly.chr
	/usr/bin/time -v bash $RunJelly $Parent1Generator $K $(echo $Threads -2 | bc)
fi 

/usr/bin/time -v bash $RunJelly $Parent2Generator $K $(echo $Threads -2 | bc)
if [ "$( tail -n 2 $Parent2Generator.Jelly.chr | head -1)" == "*" ]
then
        echo "Jelly on $Parent2Generator successfull"
else
        echo "ReRunning jelly on $Parent2Generator"
        rm $Parent2Generator.Jhash
	rm $Parent2Generator.Jelly.chr
       /usr/bin/time -v bash $RunJelly $Parent2Generator $K $(echo $Threads -2 | bc)
fi




/usr/bin/time -v bash $RunJelly $SiblingGenerator $K $(echo $Threads -2 | bc)
if [ "$( tail -n 2 $SiblingGenerator.Jelly.chr | head -1)" == "*" ]
then
        echo "Jelly on $SiblingGenerator successfull"
else
        echo "ReRunning jelly on $SiblingGenerator"
        rm $SiblingGenerator.Jhash
	rm  $SiblingGenerator.Jelly.chr
       /usr/bin/time -v bash $RunJelly $SiblingGenerator $K $(echo $Threads -2 | bc)
fi



/usr/bin/time -v bash $RunJelly $ProbandGenerator $K $(echo $Threads -2 | bc)
if [ "$( tail -n 2 $ProbandGenerator.Jelly.chr | head -1)" == "*" ]
then
        echo "Jelly on $ProbandGenerator successfull"
else
        echo "ReRunning jelly on $ProbandGenerator"
        rm $ProbandGenerator.Jhash
	rm $ProbandGenerator.Jelly.chr
        /usr/bin/time -v bash $RunJelly $ProbandGenerator $K $(echo $Threads -2 | bc)
fi


perl -ni -e 's/ /\t/;print' $ProbandGenerator.Jhash.histo
perl -ni -e 's/ /\t/;print' $Parent1Generator.Jhash.histo
perl -ni -e 's/ /\t/;print' $Parent2Generator.Jhash.histo
perl -ni -e 's/ /\t/;print' $SiblingGenerator.Jhash.histo

if [ -e "$ProbandGenerator.Jhash.histo.7.7.model" ]
then
        echo "skipping model"
else
	echo "staring model"
        /usr/bin/time -v $RUFUSmodel $ProbandGenerator.Jhash.histo $K 150 $Threads 
        echo "done with model "
fi

if [ -e "$SiblingGenerator.Jhash.histo.7.7.model" ]
then 
	 echo "skipping model"
else
	 echo "staring model"
	/usr/bin/time -v $RUFUSmodel $SiblingGenerator.Jhash.histo $K 150 $Threads
	echo "done with model "
fi 

ParentMaxE=0
MutantMinCov=$(head -2 $ProbandGenerator.Jhash.histo.7.7.model | tail -1 )
SiblingMinCov=$(head -2 $SiblingGenerator.Jhash.histo.7.7.model | tail -1 )
echo "$ParentMaxE \n $MutantMinCov \n"

date
echo "starting RUFUS build "
let "Max= $MutantMinCov*100"
if [ -e "Family.Unique.HashList" ]
then
        echo "Skipping build"
else
	/usr/bin/time -v ../RUFUS/cloud/jellyfish-MODIFIED-merge/bin/jellyfish merge $Parent1Generator.Jhash  $Parent2Generator.Jhash $SiblingGenerator.Jhash $ProbandGenerator.Jhash >  Family.Unique.HashList
fi

echo "Mut cov = $MutantMinCov and SiblingMinCov = $SiblingMinCov"
if [ -e $ProbandGenerator.k$K_c$MutantMinCov.HashList.prefilter ]
then 
	echo "skipping $ProbandGenerator.HashList pull "
else
echo "crap"
	/usr/bin/time -v bash $PullSampleHashes $ProbandGenerator.Jhash Family.Unique.HashList $MutantMinCov > $ProbandGenerator.k$K_c$MutantMinCov.HashList.prefilter
fi 

if [ -e $SiblingGenerator.k$K_c$SiblingMinCov.HashList.prefilter ]
then 
	echo "skipping $SiblingGenerator.HashList pull"
else
	/usr/bin/time -v bash $PullSampleHashes $SiblingGenerator.Jhash Family.Unique.HashList $SiblingMinCov > $SiblingGenerator.k$K_c$SiblingMinCov.HashList.prefilter
fi 


if [ -e $ProbandGenerator.k$K_c$MutantMinCov.HashList ]
then 
	echo "skipping 1kg filter"
else
echo "crap"
	 /usr/bin/time -v ../RUFUS/cloud/RUFUS.search.1kg -hf <(awk '{print $1 "\t" $2}' $ProbandGenerator.k$K_c$MutantMinCov.HashList.prefilter ) -o $ProbandGenerator.k$K_c$MutantMinCov.HashList  -c $RDIR/cloud/1000G.RUFUSreference.sorted.min45.tab -hs 25
fi 

if [ -e $SiblingGenerator.k$K_c$SiblingMinCov.HashList ]
then
        echo "skipping 1kg filter"
else
         /usr/bin/time -v  ../RUFUS/cloud/RUFUS.search.1kg -hf <(awk '{print $1 "\t" $2}' $SiblingGenerator.k$K_c$SiblingMinCov.HashList.prefilter ) -o $SiblingGenerator.k$K_c$SiblingMinCov.HashList  -c $RDIR/cloud/1000G.RUFUSreference.sorted.min45.tab -hs 25
fi

echo "done with RUFUS build "

echo "startin RUFUS filter"
if [ -e $ProbandGenerator.Mutations.fastq ]
then 
	echo "skipping filter"
else 
echo "crap"
	rm  $ProbandGenerator.temp
	mkfifo $ProbandGenerator.temp
	/usr/bin/time -v  bash $ProbandGenerator | /home/ubuntu/work/RUFUS/cloud/PassThroughSamCheck $ProbandGenerator.filter.chr >  $ProbandGenerator.temp &
	/usr/bin/time -v   $RUFUSfilter  $ProbandGenerator.k$K_c$MutantMinCov.HashList $ProbandGenerator.temp $ProbandGenerator $K 5 5 10 $(echo $Threads -2 | bc) &
	wait

	if [ "$( tail -n 2 $ProbandGenerator.filter.chr | head -1)" == "*" ]
        then
                echo "Filter on $ProbandGenerator successfull"
        else
                echo "ReRunning filter on $ProbandGenerator"
		rm  $ProbandGenerator.temp
	        mkfifo $ProbandGenerator.temp
	        /usr/bin/time -v  bash $ProbandGenerator | /home/ubuntu/work/RUFUS/cloud/PassThroughSamCheck $ProbandGenerator.filter.chr >  $ProbandGenerator.temp &
	        /usr/bin/time -v   $RUFUSfilter  $ProbandGenerator.k$K_c$MutantMinCov.HashList $ProbandGenerator.temp $ProbandGenerator $K 5 5 10 $(echo $Threads -2 | bc) &
		wait
	fi
fi 

if [ -e $ProbandGenerator.V2.overlap.hashcount.fastq.bam.vcf ]
then 
	echo "skipping overlap"
else
	echo "startin RUFUS overlap"
	/usr/bin/time -v bash $RUFUSOverlap $ProbandGenerator.Mutations.fastq 5 $ProbandGenerator $ProbandGenerator.k$MutantMinCov.HashList $Threads $ProbandGenerator.Jhash $SiblingGenerator.Jhash $Parent1Generator.Jhash $Parent2Generator.Jhash 
fi

echo "starting RUFUS filter"
if [ -e $SiblingGenerator.Mutations.fastq ]
then 
	echo "skipping filter" 
else
	rm  $SiblingGenerator.temp
	mkfifo $SiblingGenerator.temp
	/usr/bin/time -v  bash $SiblingGenerator | /home/ubuntu/work/RUFUS/cloud/PassThroughSamCheck $SiblingGenerator.filter.chr >  $SiblingGenerator.temp &
	/usr/bin/time -v   $RUFUSfilter  $SiblingGenerator.k$K_c$SiblingMinCov.HashList $SiblingGenerator.temp $SiblingGenerator $K 5 5 10 $(echo $Threads -2 | bc) &
	wait

	if [ "$( tail -n 2 $SiblingGenerator.filter.chr | head -1)" == "*" ]
	then
	        echo "Filter on $SiblingGenerator successfull"
	else
	        echo "ReRunning filter on $SiblingGenerator"
		rm  $SiblingGenerator.temp
        	mkfifo $SiblingGenerator.temp
        	/usr/bin/time -v  bash $SiblingGenerator | /home/ubuntu/work/RUFUS/cloud/PassThroughSamCheck $SiblingGenerator.filter.chr >  $SiblingGenerator.temp &
        	/usr/bin/time -v   $RUFUSfilter  $SiblingGenerator.k$K_c$SiblingMinCov.HashList $SiblingGenerator.temp $SiblingGenerator $K 5 5 10 $(echo $Threads -2 | bc) &
        	wait
	fi

fi

if [ -e $SiblingGenerator.V2.overlap.hashcount.fastq.bam.vcf ]
then
	echo "skipping overlap"
else
	echo "startin RUFUS overlap"
	/usr/bin/time -v bash $RUFUSOverlap $SiblingGenerator.Mutations.fastq 5 $SiblingGenerator $SiblingGenerator.k$SiblingMinCov.HashList $Threads $SiblingGenerator.Jhash $ProbandGenerator.Jhash $Parent1Generator.Jhash $Parent2Generator.Jhash

fi

for i in *vcf*; do aws s3  --region us-east-1 cp $i s3://marthlab.rufus/ASC.out/$Out/ ; done 
for i in *Mutations.fastq ; do aws s3  --region us-east-1 cp $i s3://marthlab.rufus/ASC.out/$Out/ ; done
for i in *HashList;  do aws s3  --region us-east-1 cp $i s3://marthlab.rufus/ASC.out/$Out/ ; done
for i in *bam; do aws s3  --region us-east-1 cp $i s3://marthlab.rufus/ASC.out/$Out/ ; done
for i in *chr; do aws s3  --region us-east-1 cp $i s3://marthlab.rufus/ASC.out/$Out/ ; done
for i in *generator.V2.overlap.asembly.hash.fastq.*; do aws s3  --region us-east-1 cp $i s3://marthlab.rufus/ASC.out/$Out/ ; done

echo "done with everything "



