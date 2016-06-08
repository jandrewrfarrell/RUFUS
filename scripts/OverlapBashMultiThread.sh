#!/bin/bash
File=$1
FinalCoverage=$2
NameStub=$3.V2
HashList=$4
Threads=$5

SampleJhash=$6
Parent1Jhash=$7
Parent2Jhash=$8
Parent3Jhash=$9

echo " you gave
File=$File
FinalCoverage=$FinalCoverage
NameStub=$NameStub
HashList=$HashList
Threads=$Threads
"

mkdir ./TempOverlap/
echo "Overlaping $File"

RDIR=/home/ubuntu/work/RUFUS

OverlapHash=$RDIR/bin/Overlap
OverlapRebion2=$RDIR/bin/OverlapRegion
ReplaceQwithDinFASTQD=$RDIR/bin/ReplaceQwithDinFASTQD
ConvertFASTqD=$RDIR/bin/ConvertFASTqD.to.FASTQ
AnnotateOverlap=$RDIR/bin/AnnotateOverlap
gkno=$RDIR/bin/gkno_launcher/gkno
samtools=$RDIR/bin/gkno_launcher/tools/samtools/samtools
RUFUSinterpret=$RDIR/bin/RUFUS.interpret
humanRef=$RDIR/bin/gkno_launcher/resources/homo_sapiens/build_37_version_3/human_reference_v37_decoys.fa
CheckHash=$RDIR/cloud/CheckJellyHashList.sh


time $OverlapHash $File .98 60 0 FP 24 1 ./TempOverlap/$NameStub.1 1 $Threads #> $File.overlap.out
time $OverlapHash ./TempOverlap/$NameStub.1.fastqd .98 50 2 FP 15 1 ./TempOverlap/$NameStub.2 1 $Threads #>>  $File.overlap.out
$ReplaceQwithDinFASTQD ./TempOverlap/$NameStub.2.fastqd > ./TempOverlap/$NameStub.3.fastqd
time $OverlapRebion2 ./TempOverlap/$NameStub.3.fastqd .95 30 0 ./TempOverlap/$NameStub.4 $NameStub 1 $Threads #>>  $File.overlap.out
time $OverlapRebion2 ./TempOverlap/$NameStub.4.fastqd .95 30 $FinalCoverage  ./TempOverlap/$NameStub.5 $NameStub 1 $Threads #>>  $File.overlap.out
$ReplaceQwithDinFASTQD ./TempOverlap/$NameStub.5.fastqd > ./$NameStub.overlap.fastqd
$ConvertFASTqD ./$NameStub.overlap.fastqd > ./$NameStub.overlap.fastq
$AnnotateOverlap $HashList ./$NameStub.overlap.fastq $NameStub.overlap.asembly.hash.fastq > ./$NameStub.overlap.hashcount.fastq 
bash $CheckHash $SampleJhash $NameStub.overlap.asembly.hash.fastq 0 > $NameStub.overlap.asembly.hash.fastq.sample
bash $CheckHash $Parent1Jhash $NameStub.overlap.asembly.hash.fastq 0 > $NameStub.overlap.asembly.hash.fastq.p1
bash $CheckHash $Parent2Jhash $NameStub.overlap.asembly.hash.fastq 0 > $NameStub.overlap.asembly.hash.fastq.p2
bash $CheckHash $Parent3Jhash $NameStub.overlap.asembly.hash.fastq 0 > $NameStub.overlap.asembly.hash.fastq.p3



$gkno bwa-se -ps human  -q ./$NameStub.overlap.hashcount.fastq -id ./$NameStub.overlap.hashcount.fastq -s ./$NameStub.overlap.hashcount.fastq -o ./$NameStub.overlap.hashcount.fastq.bam -p ILLUMINA

mkfifo check 
$samtools view ./$NameStub.overlap.hashcount.fastq.bam | $RUFUSinterpret -r $humanRef -hf $HashList -o  ./$NameStub.overlap.hashcount.fastq.bam -m 100000000 -c $NameStub.overlap.asembly.hash.fastq.p1 -c $NameStub.overlap.asembly.hash.fastq.p2 -c $NameStub.overlap.asembly.hash.fastq.p3

