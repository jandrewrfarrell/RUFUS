#!/bin/bash
File=$1
FinalCoverage=$2
NameStub=$3.V2
HashList=$4
HashSize=$5
Threads=$6
 
SampleJhash=$7
Parent1Jhash=$8
Parent2Jhash=$9
Parent3Jhash=$10

echo " you gave
File=$1
FinalCoverage=$2
NameStub=$3.V2
HashList=$4
HashSize=$5
Threads=$6

SampleJhash=$7
Parent1Jhash=$8
Parent2Jhash=$9
Parent3Jhash=$10
"

mkdir ./TempOverlap/
echo "Overlaping $File"

RDIR=/scratch/ucgd/lustre/u0991464/Projects/CEPH.1kg.cut0.5.v5/RUFUS

OverlapHash=$RDIR/bin/Overlap
OverlapRebion2=$RDIR/bin/OverlapRegion
ReplaceQwithDinFASTQD=$RDIR/bin/ReplaceQwithDinFASTQD
ConvertFASTqD=$RDIR/bin/ConvertFASTqD.to.FASTQ
AnnotateOverlap=$RDIR/bin/AnnotateOverlap
#gkno=$RDIR/bin/gkno_launcher/gkno
bwa=$RDIR/bin/bwa/bwa
#samtools=$RDIR/bin/gkno_launcher/tools/samtools/samtools
samtools=$RDIR/bin/samtools-1.6/samtools
RUFUSinterpret=$RDIR/bin/RUFUS.interpret
humanRef=/scratch/ucgd/lustre/u0991464/build_37_version_3/human_reference_v37_decoys
CheckHash=$RDIR/cloud/CheckJellyHashList.sh
OverlapSam=$RDIR/bin/OverlapSam
JellyFish=$RDIR//src/externals/jellyfish-2.2.5/bin/jellyfish



#if [ -s $NameStub.overlap.hashcount.fastq ]
#then 
#	echo "Skipping Overlap"
#else

if [ -s ./TempOverlap/$NameStub.sam.fastqd ]  
then 
	echo "skipping sam assemble"
else
	$bwa mem $humanRef $File | samtools sort -T $File -O bam - > $File.bam 
	$samtools index $File.bam 
	#$gkno bwa-se -ps human  -q $File -id $File -s $File -o $File.bam -p ILLUMINA
	$OverlapSam <( $samtools view $File.bam ) .95 75 3 ./TempOverlap/$NameStub.sam $NameStub 1 $Threads
fi

if [ -s ./TempOverlap/$NameStub.1.fastqd ]
then 	
	echo "skipping first overlap"
else
	time $OverlapHash ./TempOverlap/$NameStub.sam.fastqd .98 75 2 FP 24 1 ./TempOverlap/$NameStub.1 1 $Threads #> $File.overlap.out
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
	time $OverlapRebion2 ./TempOverlap/$NameStub.3.fastqd .97 50  0 ./TempOverlap/$NameStub.4 $NameStub 1 $Threads #>>  $File.overlap.out
fi
if [ -s ./TempOverlap/$NameStub.5.fastqd ]
then 
	echo "skipping ovelrap 5"
else
	time $OverlapRebion2 ./TempOverlap/$NameStub.4.fastqd .97 50  $FinalCoverage  ./TempOverlap/$NameStub.5 $NameStub 1 $Threads #>>  $File.overlap.out
fi 

if [ -s ./$NameStub.overlap.hashcount.fastq ]
then 
	echo "skipping final overlap work"
else

	$ReplaceQwithDinFASTQD ./TempOverlap/$NameStub.5.fastqd > ./$NameStub.overlap.fastqd
	$ConvertFASTqD ./$NameStub.overlap.fastqd > ./$NameStub.overlap.fastq
	$AnnotateOverlap $HashList ./$NameStub.overlap.fastq $NameStub.overlap.asembly.hash.fastq > ./$NameStub.overlap.hashcount.fastq 
fi


if [ -s ./$NameStub.overlap.hashcount.fastq.bam ]
then 
	echo "skipping contig alignment" 
else
	
	$bwa mem -Y -E 9,9 -O 4,4  -d 500 -w 500 -L 10,10     $humanRef ./$NameStub.overlap.hashcount.fastq | samtools sort -T $File -O bam - > ./$NameStub.overlap.hashcount.fastq.bam 	
	#$bwa mem   -Y  $humanRef ./$NameStub.overlap.hashcount.fastq | samtools sort -T $File -O bam - > ./$NameStub.overlap.hashcount.fastq.bam
	#$gkno bwa-se -ps human  -q ./$NameStub.overlap.hashcount.fastq -id ./$NameStub.overlap.hashcount.fastq -s ./$NameStub.overlap.hashcount.fastq -o ./$NameStub.overlap.hashcount.fastq.bam -p ILLUMINA
fi 

if [ -e $NameStub.overlap.asembly.hash.fastq.ref.fastq ]
then 
	echo "skipping pull reference sequecnes"
else 
	$RDIR/bin/bedtools2/bin/fastaFromBed -bed <( $RDIR/bin/bedtools2/bin/bamToBed -i ./$NameStub.overlap.hashcount.fastq.bam) -fi $RDIR/bin/gkno_launcher/resources/homo_sapiens/current/human_reference_v37.fa -fo $NameStub.overlap.asembly.hash.fastq.ref.fastq 
fi 

if [ -e ./$NameStub.overlap.hashcount.fastq.Jhash.tab ]
then 
	echo "skipping var hash generationr"
else
	$JellyFish count -C -m $HashSize -s 1G -t 20 -o ./$NameStub.overlap.hashcount.fastq.Jhash ./$NameStub.overlap.hashcount.fastq
	$JellyFish dump -c ./$NameStub.overlap.hashcount.fastq.Jhash > ./$NameStub.overlap.hashcount.fastq.Jhash.tab
fi 

if [ -e ./$NameStub.overlap.asembly.hash.fastq.ref.fastq.Jhash.tab ] 
then 
	echo "skipping ref hash generation"
else
	$JellyFish count -C -m $HashSize -s 1G -t 20 -o ./$NameStub.overlap.asembly.hash.fastq.ref.fastq.Jhash ./$NameStub.overlap.asembly.hash.fastq.ref.fastq
	$JellyFish dump -c ./$NameStub.overlap.asembly.hash.fastq.ref.fastq.Jhash > ./$NameStub.overlap.asembly.hash.fastq.ref.fastq.Jhash.tab
fi 
if [ -s $NameStub.overlap.asembly.hash.fastq.p1 ]
then
        echo "skipping hash lookup"
else
	echo "stargin hash lookup"
        bash $CheckHash $SampleJhash ./$NameStub.overlap.hashcount.fastq.Jhash.tab 0 > $NameStub.overlap.asembly.hash.fastq.sample
        bash $CheckHash $Parent1Jhash ./$NameStub.overlap.hashcount.fastq.Jhash.tab 0 > $NameStub.overlap.asembly.hash.fastq.p1
        bash $CheckHash $Parent2Jhash ./$NameStub.overlap.hashcount.fastq.Jhash.tab 0 > $NameStub.overlap.asembly.hash.fastq.p2
	echo "done with hash lookup"
fi


if [ -s $NameStub.overlap.asembly.hash.fastq.Ref.sample ]	
then 
	echo "skipping $NameStub.overlap.asembly.hash.fastq.Ref.sample"
else
	bash $CheckHash $SampleJhash  ./$NameStub.overlap.asembly.hash.fastq.ref.fastq.Jhash.tab 0 > $NameStub.overlap.asembly.hash.fastq.Ref.sample
	bash $CheckHash $Parent1Jhash ./$NameStub.overlap.asembly.hash.fastq.ref.fastq.Jhash.tab 0 > $NameStub.overlap.asembly.hash.fastq.Ref.p1
	bash $CheckHash $Parent2Jhash ./$NameStub.overlap.asembly.hash.fastq.ref.fastq.Jhash.tab 0 > $NameStub.overlap.asembly.hash.fastq.Ref.p2
	 
fi 

wait

mkfifo check 
$samtools index ./$NameStub.overlap.hashcount.fastq.bam
$samtools view ./$NameStub.overlap.hashcount.fastq.bam | $RUFUSinterpret -mod $NameStub.overlap.asembly.hash.fastq.sample -mQ 8 -r $humanRef.fa -hf $HashList -o  ./$NameStub.overlap.hashcount.fastq.bam -m 1000000 -c $NameStub.overlap.asembly.hash.fastq.p1 -c $NameStub.overlap.asembly.hash.fastq.p2 -cR $NameStub.overlap.asembly.hash.fastq.Ref.p1 -cR $NameStub.overlap.asembly.hash.fastq.Ref.p2 -sR $NameStub.overlap.asembly.hash.fastq.Ref.sample -s $NameStub.overlap.asembly.hash.fastq.sample 

grep ^# ./$NameStub.overlap.hashcount.fastq.bam.vcf> ./$NameStub.overlap.hashcount.fastq.bam.vcf.sorted.vcf
grep -v  ^# ./$NameStub.overlap.hashcount.fastq.bam.vcf | sort -k1,1 -k2,2n >> ./$NameStub.overlap.hashcount.fastq.bam.vcf.sorted.vcf
$RDIR/bin/tabix/bgzip -f ./$NameStub.overlap.hashcount.fastq.bam.vcf.sorted.vcf
$RDIR/bin/tabix/tabix ./$NameStub.overlap.hashcount.fastq.bam.vcf.sorted.vcf.gz
