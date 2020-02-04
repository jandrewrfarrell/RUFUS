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

CDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"


RDIR=$CDIR/../

OverlapHash=$RDIR/bin/Overlap
OverlapRebion2=$RDIR/bin/OverlapRegion
ReplaceQwithDinFASTQD=$RDIR/bin/ReplaceQwithDinFASTQD
ConvertFASTqD=$RDIR/bin/ConvertFASTqD.to.FASTQ
AnnotateOverlap=$RDIR/bin/AnnotateOverlap
#gkno=$RDIR/bin/gkno_launcher/gkno
bwa=$RDIR/bin/externals/bwa/src/bwa_project/bwa
RUFUSinterpret=$RDIR/bin/RUFUS.interpret
CheckHash=$RDIR/scripts/CheckJellyHashList.sh
OverlapSam=$RDIR/bin/OverlapSam
JellyFish=$RDIR/bin/externals/jellyfish/src/jellyfish_project/bin/jellyfish
MOBList=$RDIR/resources/primate_non-LTR_Retrotransposon.fasta 

$bwa mem -t $Threads -Y -E 0,0 -O 6,6  -d 500 -w 500 -L 0,0 $MOBList ./$NameStub.overlap.hashcount.fastq | samtools sort -T $File -O sam - > ./TempOverlap/$NameStub.overlap.hashcount.fastq.MOB.sam
        samtools index ./$NameStub.overlap.hashcount.fastq.bam

#if [ -s $NameStub.overlap.hashcount.fastq ]
#then 
#	echo "Skipping Overlap"
#else

parentCRString=""
c="-c"
cr="-cR"
space=" "


######################## BUILDING UP parent c and cR string ##############################
for parent in $ParentsJhash;
do
    parentCRString="$parentCRString -c ./Intermediates/$NameStub.overlap.asembly.hash.fastq.$parent -cR ./Intermediates/$NameStub.overlap.asembly.hash.fastq.Ref.$parent "
done

#echo "final parent String is  $parentCRString"
##########################################################################################


mkfifo check 
samtools index ./$NameStub.overlap.hashcount.fastq.bam

echo "$RUFUSinterpret -mob ./TempOverlap/$NameStub.overlap.hashcount.fastq.MOB.sam -mod Intermediates/$NameStub.overlap.asembly.hash.fastq.sample -mQ 8 -r $humanRef -hf $HashList -o  ./$NameStub.overlap.hashcount.fastq.bam -m 1000000 $(echo $parentCRString) -sR Intermediates/$NameStub.overlap.asembly.hash.fastq.Ref.sample -s Intermediates/$NameStub.overlap.asembly.hash.fastq.sample -e ./Intermediates/$NameStub.ref.RepRefHash "

samtools view ./$NameStub.overlap.hashcount.fastq.bam | $RUFUSinterpret -mob ./TempOverlap/$NameStub.overlap.hashcount.fastq.MOB.sam -mod Intermediates/$NameStub.overlap.asembly.hash.fastq.sample -mQ 8 -r $humanRef -hf $HashList -o  ./$NameStub.overlap.hashcount.fastq.bam -m 100000000 $(echo $parentCRString) -sR Intermediates/$NameStub.overlap.asembly.hash.fastq.Ref.sample -s Intermediates/$NameStub.overlap.asembly.hash.fastq.sample -e ./Intermediates/$NameStub.ref.RepRefHash 


