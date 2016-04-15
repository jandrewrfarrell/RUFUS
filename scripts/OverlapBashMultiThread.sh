#!/bin/bash
File=$1
FinalCoverage=$2
NameStub=$3.V2
HashList=$4
Threads=$5

echo " you gave
File=$File
FinalCoverage=$FinalCoverage
NameStub=$NameStub
Threads=$Threads
"
mkdir ./TempOverlap/
echo "Overlaping $File"

RDIR=/uufs/chpc.utah.edu/common/home/u0991464/d1/home/farrelac/bin/RUFUS_git/

OverlapHash=$RDIR/bin/Overlap
OverlapRebion2=$RDIR/bin/OverlapRegion
ReplaceQwithDinFASTQD=$RDIR/bin/ReplaceQwithDinFASTQD
ConvertFASTqD=$RDIR/bin/ConvertFASTqD
AnnotateOverlap=$RDIR/bin/AnnotateOverlap
gkno=$RDIR/src/externals/gkno_launcher/gkno
samtools=$RDIR/src/externals/gkno_launcher/tools/samtools/samtools
RUFUS.interpret=$RDIR/bin/RUFUS.interpret
humanRef=$RDIR/src/externals/gkno_launcher/resources/homo_sapiens/build_37_version_3/human_reference_v37_decoys.fa



time $OverlapHash $File .98 60 0 FP 24 1 ./TempOverlap/$NameStub.1 1 $Threads #> $File.overlap.out
time $OverlapHash ./TempOverlap/$NameStub.1.fastqd .98 50 2 FP 15 1 ./TempOverlap/$NameStub.2 1 $Threads #>>  $File.overlap.out
$ReplaceQwithDinFASTQD ./TempOverlap/$NameStub.2.fastqd > ./TempOverlap/$NameStub.3.fastqd
time $OverlapRebion2 ./TempOverlap/$NameStub.3.fastqd .95 30 0 ./TempOverlap/$NameStub.4 $NameStub 1 $Threads #>>  $File.overlap.out
time $OverlapRebion2 ./TempOverlap/$NameStub.4.fastqd .95 30 $FinalCoverage  ./TempOverlap/$NameStub.5 $NameStub 1 $Threads #>>  $File.overlap.out
$ReplaceQwithDinFASTQD ./TempOverlap/$NameStub.5.fastqd > ./$NameStub.overlap.fastqd
$ConvertFASTqD ./$NameStub.overlap.fastqd > ./$NameStub.overlap.fastq
$AnnotateOverlap $HashList ./$NameStub.overlap.fastq > ./$NameStub.overlap.hashcount.fastq
$gkno bwa-se -r ~/d1/home/farrelac/references/current/human_reference_v37_decoys  -q ./$NameStub.overlap.hashcount.fastq -id ./$NameStub.overlap.hashcount.fastq -s ./$NameStub.overlap.hashcount.fastq -o ./$NameStub.overlap.hashcount.fastq.bam -p ILLUMINA

$samtools view ./$NameStub.overlap.hashcount.fastq.bam | $RUFUS.interpret -r $humanRef -hf $HashList -o  ./$NameStub.overlap.hashcount.fastq.bam -m 100000000

