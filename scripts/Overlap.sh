#!/bin/bash
humanRef=$1
File=$2
FinalCoverage=$3
NameStub=$4.V2
HashList=$5
HashSize=$6
Threads=$7
 
SampleJhash=$8
ParentsJhash=$9

humanRefBwa=${10}
refHash=${11}

echo " you gave
File=$2
FinalCoverage=$3
NameStub=$4.V2
HashList=$5
HashSize=$6
Threads=$7
"

echo "final coveage is $FinalCoverage"

echo "@@@@@@@@@@@@@__IN_OVERLAP__@@@@@@@@@@@@@@@"
echo "human ref in Overlap is $humanRef"
echo "bwa human ref in Overlap is $humanRefBwa"
echo "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@"
mkdir ./TempOverlap/
mkdir ./Intermediates/
echo "Overlaping $File"

RDIR=/scratch/ucgd/lustre/u0991464/RUFUS.simulation.test/testStricterOverlap/RUFUS

OverlapHash=$RDIR/bin/Overlap
OverlapRebion2=$RDIR/bin/OverlapRegion
ReplaceQwithDinFASTQD=$RDIR/bin/ReplaceQwithDinFASTQD
ConvertFASTqD=$RDIR/bin/ConvertFASTqD.to.FASTQ
AnnotateOverlap=$RDIR/bin/AnnotateOverlap
#gkno=$RDIR/bin/gkno_launcher/gkno
bwa=$RDIR/bin/externals/bwa/src/bwa_project/bwa
samtools=$RDIR/bin/externals/samtools/samtools
RUFUSinterpret=$RDIR/bin/RUFUS.interpret
CheckHash=$RDIR/cloud/CheckJellyHashList.sh
OverlapSam=$RDIR/bin/OverlapSam
JellyFish=$RDIR/src/externals/jellyfish-2.2.5/bin/jellyfish


#if [ -s $NameStub.overlap.hashcount.fastq ]
#then 
#	echo "Skipping Overlap"
#else

if [ -s ./TempOverlap/$NameStub.sam.fastqd ]  
then 
	echo "skipping sam assemble"
else
	$bwa mem $humanRefBwa "$File" | $samtools sort -T $File -O bam - > $File.bam
	$samtools index $File.bam 
	#$gkno bwa-se -ps human  -q $File -id $File -s $File -o $File.bam -p ILLUMINA
	$OverlapSam <( $samtools view $File.bam ) .95 25 1 ./TempOverlap/$NameStub.sam $NameStub 1 $Threads
fi

if [ -s ./TempOverlap/$NameStub.1.fastqd ]
then 	
	echo "skipping first overlap"
else
	time $OverlapHash ./TempOverlap/$NameStub.sam.fastqd .98 50 1 FP 24 1 ./TempOverlap/$NameStub.1 1 $Threads #> $File.overlap.out
fi 
if [ -s ./TempOverlap/$NameStub.2.fastqd ]
then 
	echo "skipping second overlap"
else
	time $OverlapHash ./TempOverlap/$NameStub.1.fastqd .98 50 2 FP 15 1 ./TempOverlap/$NameStub.2 1 $Threads #>>  $File.overlap.out
fi
$ReplaceQwithDinFASTQD ./TempOverlap/$NameStub.2.fastqd > ./TempOverlap/$NameStub.3.fastqd
if [ -s ./TempOverlap/$NameStub.4.fastqd ]
then 
	echo "skipping overlap 4"
else 
	time $OverlapRebion2 ./TempOverlap/$NameStub.3.fastqd .98 50 2 ./TempOverlap/$NameStub.4 $NameStub 1 $Threads > /dev/null #>>  $File.overlap.out
fi
if [ -s ./TempOverlap/$NameStub.5.fastqd ]
then 
	echo "skipping ovelrap 5"
else
	time $OverlapRebion2 ./TempOverlap/$NameStub.4.fastqd .98 35 $FinalCoverage  ./TempOverlap/$NameStub.5 $NameStub 1 $Threads > /dev/null #>>  $File.overlap.out
fi 

if [ -s ./$NameStub.overlap.hashcount.fastq ]
then 
	echo "skipping final overlap work"
else

	$ReplaceQwithDinFASTQD ./TempOverlap/$NameStub.5.fastqd > ./$NameStub.overlap.fastqd
	$ConvertFASTqD ./$NameStub.overlap.fastqd > ./$NameStub.overlap.fastq

	echo "$AnnotateOverlap $HashList ./$NameStub.overlap.fastq TempOverlap/$NameStub.overlap.asembly.hash.fastq > ./$NameStub.overlap.hashcount.fastq"              
	      $AnnotateOverlap $HashList ./$NameStub.overlap.fastq TempOverlap/$NameStub.overlap.asembly.hash.fastq > ./$NameStub.overlap.hashcount.fastq
fi


if [ -s ./$NameStub.overlap.hashcount.fastq.bam ]
then 
	echo "skipping contig alignment" 
else
        $bwa mem -Y -E 0,0 -O 6,6  -d 500 -w 500 -L 0,0 $humanRefBwa ./$NameStub.overlap.hashcount.fastq | $samtools sort -T $File -O bam - > ./$NameStub.overlap.hashcount.fastq.bam
	$samtools index ./$NameStub.overlap.hashcount.fastq.bam
fi 


#############################################################################################################
if [ -e ./Intermediates/$NameStub.overlap.asembly.hash.fastq.ref.fastq ]
then 
	echo "skipping pull reference sequecnes"
else
    echo " $RDIR/bin/externals/bedtools2/src/bedtools2_project/bin/fastaFromBed -bed <( $RDIR/bin/externals/bedtools2/src/bedtools2_project/bin/bamToBed -i ./$NameStub.overlap.hashcount.fastq.bam) -fi $humanRef -fo ./Intermediates/$NameStub.overlap.asembly.hash.fastq.ref.fastq"
	$RDIR/bin/externals/bedtools2/src/bedtools2_project/bin/fastaFromBed -bed <( $RDIR/bin/externals/bedtools2/src/bedtools2_project/bin/bamToBed -i ./$NameStub.overlap.hashcount.fastq.bam) -fi $humanRef -fo ./Intermediates/$NameStub.overlap.asembly.hash.fastq.ref.fastq 
fi 

if [ -e ./Intermediates/$NameStub.overlap.hashcount.fastq.Jhash ]
then 
	echo "skipping var hash generationr"
else
	echo "$JellyFish count -m $HashSize -s 1G -t 20 -o ./Intermediates/$NameStub.overlap.hashcount.fastq.Jhash ./$NameStub.overlap.hashcount.fastq"
	$JellyFish count -m $HashSize -s 1G -t 20 -o ./Intermediates/$NameStub.overlap.hashcount.fastq.Jhash ./$NameStub.overlap.hashcount.fastq
	echo "$JellyFish dump  -c ./Intermediates/$NameStub.overlap.hashcount.fastq.Jhash > ./Intermediates/$NameStub.overlap.hashcount.fastq.Jhash.tab"
	$JellyFish dump  -c ./Intermediates/$NameStub.overlap.hashcount.fastq.Jhash > ./Intermediates/$NameStub.overlap.hashcount.fastq.Jhash.tab
fi 

if [ -e ././Intermediates/$NameStub.overlap.asembly.hash.fastq.ref.fastq.Jhash ] 
then 
	echo "skipping ref hash generation"
else
	echo " $JellyFish count -m $HashSize -s 1G -t 20 -o ././Intermediates/$NameStub.overlap.asembly.hash.fastq.ref.fastq.Jhash ././Intermediates/$NameStub.overlap.asembly.hash.fastq.ref.fastq"
	$JellyFish count -m $HashSize -s 1G -t 20 -o ././Intermediates/$NameStub.overlap.asembly.hash.fastq.ref.fastq.Jhash ././Intermediates/$NameStub.overlap.asembly.hash.fastq.ref.fastq
	echo "$JellyFish -c  ././Intermediates/$NameStub.overlap.asembly.hash.fastq.ref.fastq.Jhash > ././Intermediates/$NameStub.overlap.asembly.hash.fastq.ref.fastq.Jhash.tab"
	$JellyFish dump -c  ././Intermediates/$NameStub.overlap.asembly.hash.fastq.ref.fastq.Jhash > ././Intermediates/$NameStub.overlap.asembly.hash.fastq.ref.fastq.Jhash.tab
fi 
if [ -s Intermediates/$NameStub.overlap.asembly.hash.fastq.sample ]
then
        echo "skipping  Intermediates/$NameStub.overlap.asembly.hash.fastq.sample file already exitst"
else
	echo "starting hash lookup"
        bash $CheckHash $SampleJhash ./Intermediates/$NameStub.overlap.hashcount.fastq.Jhash.tab 0 > Intermediates/$NameStub.overlap.asembly.hash.fastq.sample
	echo "done with hash lookup"
fi
for parent in $ParentsJhash
        do
            if [ -s Intermediates/$NameStub.overlap.asembly.hash.fastq.$parent ]
            then
                echo "skiping Intermediates/$NameStub.overlap.asembly.hash.fastq.$parent already exists"
            else
                echo "-$parent-"
                echo " bash $CheckHash $parent ./Intermediates/$NameStub.overlap.hashcount.fastq.Jhash.tab 0 > Intermediates/$NameStub.overlap.asembly.hash.fastq.$parent "
                bash $CheckHash $parent ./Intermediates/$NameStub.overlap.hashcount.fastq.Jhash.tab 0 > Intermediates/$NameStub.overlap.asembly.hash.fastq.$parent
            fi
done



if [ -s Intermediates/$NameStub.overlap.asembly.hash.fastq.Ref.sample ]	
then 
	echo "skipping Intermediates/$NameStub.overlap.asembly.hash.fastq.Ref.sample"
else
	bash $CheckHash $SampleJhash  ././Intermediates/$NameStub.overlap.asembly.hash.fastq.ref.fastq.Jhash.tab 0 > Intermediates/$NameStub.overlap.asembly.hash.fastq.Ref.sample
fi 
for parent in $ParentsJhash
do
	if [ -s Intermediates/$NameStub.overlap.asembly.hash.fastq.Ref.$parent ]
        then
            echo "skipping Intermediates/$NameStub.overlap.asembly.hash.fastq.Ref.$parent already exitst"
        else
            echo "-$parent-"
            echo " bash $CheckHash $parent ././Intermediates/$NameStub.overlap.asembly.hash.fastq.ref.fastq.Jhash.tab 0 > Intermediates/$NameStub.overlap.asembly.hash.fastq.Ref.$parent"
            bash $CheckHash $parent ././Intermediates/$NameStub.overlap.asembly.hash.fastq.ref.fastq.Jhash.tab 0 > Intermediates/$NameStub.overlap.asembly.hash.fastq.Ref.$parent
       fi
done

parentCRString=""
c="-c"
cr="-cR"
space=" "


######################## BUILDING UP parent c and cR string ##############################
for parent in $ParentsJhash;
do
    parentCRString="$parentCRString -c Intermediates/$NameStub.overlap.asembly.hash.fastq.$parent -cR Intermediates/$NameStub.overlap.asembly.hash.fastq.Ref.$parent "
    echo "building up parentCR string in Overlap.sh... "
    echo "$parentCRString"
done

echo "final parent String is  $parentCRString"
##########################################################################################

if [ -s ./Intermediates/$NameStub.ref.RepRefHash ]
then
        echo "exclude already exists"
else
	echo "bash $CheckHash $refHash ./Intermediates/$NameStub.overlap.hashcount.fastq.Jhash.tab 1 > Intermediates/$NameStub.ref.RepRefHash"
	bash $CheckHash $refHash ./Intermediates/$NameStub.overlap.hashcount.fastq.Jhash.tab 1 > Intermediates/$NameStub.ref.RepRefHash
fi
wait

mkfifo check 

##samtools view ./$NameStub.overlap.hashcount.fastq.bam | $RUFUSinterpret -mod Intermediates/$NameStub.overlap.asembly.hash.fastq.sample -mQ 8 -r $humanRef -hf $HashList -o  ./$NameStub.overlap.hashcount.fastq.bam -m 1000000 -c $NameStub.overlap.asembly.hash.fastq.p1 -c $NameStub.overlap.asembly.hash.fastq.p2 -cR $NameStub.overlap.asembly.hash.fastq.Ref.p1 -cR $NameStub.overlap.asembly.hash.fastq.Ref.p2 -sR Intermediates/$NameStub.overlap.asembly.hash.fastq.Ref.sample -s Intermediates/$NameStub.overlap.asembly.hash.fastq.sample -e ./Intermediates/$NameStub.ref.RepRefHash 

$samtools view ./$NameStub.overlap.hashcount.fastq.bam | $RUFUSinterpret -mod Intermediates/$NameStub.overlap.asembly.hash.fastq.sample -mQ 0 -r $humanRef -hf $HashList -o  ./$NameStub.overlap.hashcount.fastq.bam -m 1000000 $(echo $parentCRString) -sR Intermediates/$NameStub.overlap.asembly.hash.fastq.Ref.sample -s Intermediates/$NameStub.overlap.asembly.hash.fastq.sample -e ./Intermediates/$NameStub.ref.RepRefHash 



#grep ^# ./$NameStub.overlap.hashcount.fastq.bam.vcf> ./$NameStub.overlap.hashcount.fastq.bam.vcf.sorted.vcf
#grep -v  ^# ./$NameStub.overlap.hashcount.fastq.bam.vcf | sort -k1,1 -k2,2n >> ./$NameStub.overlap.hashcount.fastq.bam.vcf.sorted.vcf
#bash $RDIR/scripts/VilterAutosomeOnly ./$NameStub.overlap.hashcount.fastq.bam.vcf.sorted.vcf > ./$NameStub.overlap.hashcount.fastq.bam.vcf.sorted.vcf.autosome.vcf
#
#$RDIR/bin/tabix/bgzip -f ./$NameStub.overlap.hashcount.fastq.bam.vcf.sorted.vcf
#$RDIR/bin/tabix/tabix ./$NameStub.overlap.hashcount.fastq.bam.vcf.sorted.vcf.gz
