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

OverlapHash=~/d1/home/farrelac/bin/RUFUSAwesome/Overlap19
OverlapRebion2=~/d1/home/farrelac/bin/RUFUSAwesome/OverlapRegion2

time $OverlapHash $File .98 60 0 FP 24 1 ./TempOverlap/$NameStub.1 1 $Threads #> $File.overlap.out  
time $OverlapHash ./TempOverlap/$NameStub.1.fastqd .98 50 2 FP 15 1 ./TempOverlap/$NameStub.2 1 $Threads #>>  $File.overlap.out
~/d1/home/farrelac/bin/RUFUSAwesome/ReplaceQwithDinFASTQD ./TempOverlap/$NameStub.2.fastqd > ./TempOverlap/$NameStub.3.fastqd 
time $OverlapRebion2 ./TempOverlap/$NameStub.3.fastqd .95 30 0 ./TempOverlap/$NameStub.4 $NameStub 1 $Threads #>>  $File.overlap.out
time $OverlapRebion2 ./TempOverlap/$NameStub.4.fastqd .95 30 $FinalCoverage  ./TempOverlap/$NameStub.5 $NameStub 1 $Threads #>>  $File.overlap.out
~/d1/home/farrelac/bin/RUFUSAwesome/ReplaceQwithDinFASTQD ./TempOverlap/$NameStub.5.fastqd > ./$NameStub.overlap.fastqd
~/d1/home/farrelac/bin/RUFUSAwesome/ConvertFASTqD.to.FASTQ ./$NameStub.overlap.fastqd > ./$NameStub.overlap.fastq
~/d1/home/farrelac/bin/RUFUSAwesome/AnnotateOverlap $HashList ./$NameStub.overlap.fastq > ./$NameStub.overlap.hashcount.fastq
~/bin/gkno_launcher/gkno bwa-se -r ~/d1/home/farrelac/references/current/human_reference_v37_decoys  -q ./$NameStub.overlap.hashcount.fastq -id ./$NameStub.overlap.hashcount.fastq -s ./$NameStub.overlap.hashcount.fastq -o ./$NameStub.overlap.hashcount.fastq.bam -p ILLUMINA

~/bin/samtools-1.2/samtools view ./$NameStub.overlap.hashcount.fastq.bam | ~/d1/home/farrelac/bin/RUFUSAwesome/RUFUS.interpret -r ~/d1/data/project_fasta/human_g1k_v37_decoy_phix.fasta -hf $HashList -o  ./$NameStub.overlap.hashcount.fastq.bam.V3  -m 100000000
