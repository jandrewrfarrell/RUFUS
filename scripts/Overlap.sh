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
echo "Overlaping $File"

RDIR=/uufs/chpc.utah.edu/common/home/u0401321/finalRUFUStest/RUFUS

OverlapHash=$RDIR/bin/Overlap
OverlapRebion2=$RDIR/bin/OverlapRegion
ReplaceQwithDinFASTQD=$RDIR/bin/ReplaceQwithDinFASTQD
ConvertFASTqD=$RDIR/bin/ConvertFASTqD.to.FASTQ
AnnotateOverlap=$RDIR/bin/AnnotateOverlap
#gkno=$RDIR/bin/gkno_launcher/gkno
bwa=$RDIR/bin/externals/bwa/src/bwa_project/bwa
samtools=$RDIR/bin/externals/samtools/src/SAMTOOLS_PROJECT/samtools
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
	$OverlapSam <( $samtools view $File.bam ) .95 50 3 ./TempOverlap/$NameStub.sam $NameStub 1 $Threads
fi

if [ -s ./TempOverlap/$NameStub.1.fastqd ]
then 	
	echo "skipping first overlap"
else
	time $OverlapHash ./TempOverlap/$NameStub.sam.fastqd .98 50 2 FP 24 1 ./TempOverlap/$NameStub.1 1 $Threads #> $File.overlap.out
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
	time $OverlapRebion2 ./TempOverlap/$NameStub.3.fastqd .95 30 0 ./TempOverlap/$NameStub.4 $NameStub 1 $Threads > /dev/null #>>  $File.overlap.out
fi
if [ -s ./TempOverlap/$NameStub.5.fastqd ]
then 
	echo "skipping ovelrap 5"
else
	time $OverlapRebion2 ./TempOverlap/$NameStub.4.fastqd .95 30 $FinalCoverage  ./TempOverlap/$NameStub.5 $NameStub 1 $Threads > /dev/null #>>  $File.overlap.out
fi 

if [ -s ./$NameStub.overlap.hashcount.fastq ]
then 
	echo "skipping final overlap work"
else

	$ReplaceQwithDinFASTQD ./TempOverlap/$NameStub.5.fastqd > ./$NameStub.overlap.fastqd
	$ConvertFASTqD ./$NameStub.overlap.fastqd > ./$NameStub.overlap.fastq

	echo "$AnnotateOverlap $HashList ./$NameStub.overlap.fastq $NameStub.overlap.asembly.hash.fastq > ./$NameStub.overlap.hashcount.fastq"              
	      $AnnotateOverlap $HashList ./$NameStub.overlap.fastq $NameStub.overlap.asembly.hash.fastq > ./$NameStub.overlap.hashcount.fastq
fi


if [ -s ./$NameStub.overlap.hashcount.fastq.bam ]
then 
	echo "skipping contig alignment" 
else
        $bwa mem -Y -E 0,0 -O 6,6  -d 500 -w 500 -L 0,0 $humanRefBwa ./$NameStub.overlap.hashcount.fastq | $samtools sort -T $File -O bam - > ./$NameStub.overlap.hashcount.fastq.bam
	#$bwa mem "$humanRefBwa" ./"$NameStub".overlap.hashcount.fastq | samtools sort -T $File -O bam - > ./"$NameStub".overlap.hashcount.fastq.bam 	
	#$gkno bwa-se -ps human  -q ./$NameStub.overlap.hashcount.fastq -id ./$NameStub.overlap.hashcount.fastq -s ./$NameStub.overlap.hashcount.fastq -o ./$NameStub.overlap.hashcount.fastq.bam -p ILLUMINA
fi 


#############################################################################################################
if [ -e $NameStub.overlap.asembly.hash.fastq.ref.fastq ]
then 
	echo "skipping pull reference sequecnes"
else
    echo " $RDIR/bin/externals/bedtools2/src/bedtools2_project/bin/fastaFromBed -bed <( $RDIR/bin/externals/bedtools2/src/bedtools2_project/bin/bamToBed -i ./$NameStub.overlap.hashcount.fastq.bam) -fi $humanRef -fo $NameStub.overlap.asembly.hash.fastq.ref.fastq"
	$RDIR/bin/externals/bedtools2/src/bedtools2_project/bin/fastaFromBed -bed <( $RDIR/bin/externals/bedtools2/src/bedtools2_project/bin/bamToBed -i ./$NameStub.overlap.hashcount.fastq.bam) -fi $humanRef -fo $NameStub.overlap.asembly.hash.fastq.ref.fastq 
fi 

if [ -e ./$NameStub.overlap.hashcount.fastq.Jhash.tab ]
then 
	echo "skipping var hash generationr"
else
	echo "$JellyFish count -m $HashSize -s 1G -t 20 -o ./$NameStub.overlap.hashcount.fastq.Jhash ./$NameStub.overlap.hashcount.fastq"
	$JellyFish count -m $HashSize -s 1G -t 20 -o ./$NameStub.overlap.hashcount.fastq.Jhash ./$NameStub.overlap.hashcount.fastq
	echo "$JellyFish dump  -c ./$NameStub.overlap.hashcount.fastq.Jhash > ./$NameStub.overlap.hashcount.fastq.Jhash.tab"
	$JellyFish dump  -c ./$NameStub.overlap.hashcount.fastq.Jhash > ./$NameStub.overlap.hashcount.fastq.Jhash.tab
fi 

if [ -e ./$NameStub.overlap.asembly.hash.fastq.ref.fastq.Jhash.tab ] 
then 
	echo "skipping ref hash generation"
else
	echo " $JellyFish count -m $HashSize -s 1G -t 20 -o ./$NameStub.overlap.asembly.hash.fastq.ref.fastq.Jhash ./$NameStub.overlap.asembly.hash.fastq.ref.fastq"
	$JellyFish count -m $HashSize -s 1G -t 20 -o ./$NameStub.overlap.asembly.hash.fastq.ref.fastq.Jhash ./$NameStub.overlap.asembly.hash.fastq.ref.fastq
	echo "$JellyFish -c  ./$NameStub.overlap.asembly.hash.fastq.ref.fastq.Jhash > ./$NameStub.overlap.asembly.hash.fastq.ref.fastq.Jhash.tab"
	$JellyFish dump -c  ./$NameStub.overlap.asembly.hash.fastq.ref.fastq.Jhash > ./$NameStub.overlap.asembly.hash.fastq.ref.fastq.Jhash.tab
fi 
if [ -s $NameStub.overlap.asembly.hash.fastq.sample ]
then
        echo "skipping  $NameStub.overlap.asembly.hash.fastq.sample file already exitst"
else
	echo "starting hash lookup"
        bash $CheckHash $SampleJhash ./$NameStub.overlap.hashcount.fastq.Jhash.tab 0 > $NameStub.overlap.asembly.hash.fastq.sample
	echo "done with hash lookup"
fi
for parent in $ParentsJhash
        do
            if [ -s $NameStub.overlap.asembly.hash.fastq.$parent ]
            then
                echo "skiping $NameStub.overlap.asembly.hash.fastq.$parent already exists"
            else
                echo "-$parent-"
                echo " bash $CheckHash $parent ./$NameStub.overlap.hashcount.fastq.Jhash.tab 0 > "$NameStub".overlap.asembly.hash.fastq."$parent" "
                bash $CheckHash $parent ./$NameStub.overlap.hashcount.fastq.Jhash.tab 0 > $NameStub".overlap.asembly.hash.fastq."$parent
            fi
done



if [ -s $NameStub.overlap.asembly.hash.fastq.Ref.sample ]	
then 
	echo "skipping $NameStub.overlap.asembly.hash.fastq.Ref.sample"
else
	bash $CheckHash $SampleJhash  ./$NameStub.overlap.asembly.hash.fastq.ref.fastq.Jhash.tab 0 > $NameStub.overlap.asembly.hash.fastq.Ref.sample
fi 
for parent in $ParentsJhash
        do
            if [ -s $NameStub.overlap.asembly.hash.fastq.Ref.$parent ]
            then
                echo "skipping $NameStub.overlap.asembly.hash.fastq.Ref.$parent already exitst"
            else

                echo "-$parent-"
                echo " bash $CheckHash $parent ./$NameStub.overlap.asembly.hash.fastq.ref.fastq.Jhash.tab 0 > $NameStub.overlap.asembly.hash.fastq.Ref.$parent"
                   bash $CheckHash $parent ./$NameStub.overlap.asembly.hash.fastq.ref.fastq.Jhash.tab 0 > $NameStub.overlap.asembly.hash.fastq.Ref.$parent
           fi
        done

parentCRString=""
c="-c"
cr="-cR"
space=" "


######################## BUILDING UP parent c and cR string ##############################
for parent in $ParentsJhash;
do
    parentCRString="$parentCRString -c $NameStub.overlap.asembly.hash.fastq.$parent -cR $NameStub.overlap.asembly.hash.fastq.Ref.$parent "
    echo "building up parentCR string in Overlap.sh... "
    echo "$parentCRString"
done

echo "final parent String is  $parentCRString"
##########################################################################################

if [ -s ./$NameStub.ref.exclude ]
then
        echo "exclude already exists"
else
        bash $CheckHash $RefHash ./$NameStub.overlap.hashcount.fastq.Jhash.tab 1 > ./$NameStub.ref.RepRefHash
fi
wait

mkfifo check 
$samtools index ./$NameStub.overlap.hashcount.fastq.bam

##samtools view ./$NameStub.overlap.hashcount.fastq.bam | $RUFUSinterpret -mod $NameStub.overlap.asembly.hash.fastq.sample -mQ 8 -r $humanRef -hf $HashList -o  ./$NameStub.overlap.hashcount.fastq.bam -m 1000000 -c $NameStub.overlap.asembly.hash.fastq.p1 -c $NameStub.overlap.asembly.hash.fastq.p2 -cR $NameStub.overlap.asembly.hash.fastq.Ref.p1 -cR $NameStub.overlap.asembly.hash.fastq.Ref.p2 -sR $NameStub.overlap.asembly.hash.fastq.Ref.sample -s $NameStub.overlap.asembly.hash.fastq.sample -e ./$NameStub.ref.RepRefHash 

$samtools view ./$NameStub.overlap.hashcount.fastq.bam | $RUFUSinterpret -mod $NameStub.overlap.asembly.hash.fastq.sample -mQ 8 -r $humanRef -hf $HashList -o  ./$NameStub.overlap.hashcount.fastq.bam -m 1000000 $(echo $parentCRString) -sR $NameStub.overlap.asembly.hash.fastq.Ref.sample -s $NameStub.overlap.asembly.hash.fastq.sample -e ./$NameStub.ref.RepRefHash 



grep ^# ./$NameStub.overlap.hashcount.fastq.bam.vcf> ./$NameStub.overlap.hashcount.fastq.bam.vcf.sorted.vcf
grep -v  ^# ./$NameStub.overlap.hashcount.fastq.bam.vcf | sort -k1,1 -k2,2n >> ./$NameStub.overlap.hashcount.fastq.bam.vcf.sorted.vcf
bash $RDIR/scripts/VilterAutosomeOnly ./$NameStub.overlap.hashcount.fastq.bam.vcf.sorted.vcf > ./$NameStub.overlap.hashcount.fastq.bam.vcf.sorted.vcf.autosome.vcf

$RDIR/bin/tabix/bgzip -f ./$NameStub.overlap.hashcount.fastq.bam.vcf.sorted.vcf
$RDIR/bin/tabix/tabix ./$NameStub.overlap.hashcount.fastq.bam.vcf.sorted.vcf.gz
