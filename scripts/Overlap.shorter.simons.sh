#!/bin/bash

set -e 

humanRef=$1
File=$2
FinalCoverage=$3
NameStub=$4.V2
HashList=$5
HashSize=$6
Threads=$7
MaxAlleleSize=$8
speed=$9

SampleJhash=${10}
ParentsJhash=${11}

humanRefBwa=${12}
refHash=${13}
MaxCov=100000
echo " you gave
File=$2
FinalCoverage=$3
NameStub=$4.V2
HashList=$5
HashSize=$6
Threads=$7
"

echo "final coveage is $FinalCoverage"


echo "RUNNING THIS ONE"
echo "@@@@@@@@@@@@@__IN_OVERLAP__@@@@@@@@@@@@@@@"
echo "human ref in Overlap is $humanRef"
echo "bwa human ref in Overlap is $humanRefBwa"
echo "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@"
echo "Overlaping $File"

CDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"


RDIR=$CDIR/../

OverlapHash=$RDIR/bin/Overlap
OverlapRebion2=$RDIR/bin/OverlapRegion
OverlapRegionSmall=$RDIR/bin/OverlapRegion.small
ReplaceQwithDinFASTQD=$RDIR/bin/ReplaceQwithDinFASTQD
ConvertFASTqD=$RDIR/bin/ConvertFASTqD.to.FASTQ
AnnotateOverlap=$RDIR/bin/AnnotateOverlap
bwa=$RDIR/bin/externals/bwa/src/bwa_project/bwa
RUFUSinterpret=$RDIR/bin/RUFUS.interpret
CheckHash=$RDIR/scripts/CheckJellyHashList.sh
OverlapSam=$RDIR/bin/OverlapSam
JellyFish=$RDIR/bin/externals/jellyfish/src/jellyfish_project/bin/jellyfish
MOBList=$RDIR/resources/primate_non-LTR_Retrotransposon.fasta

#if [ -s $NameStub.overlap.hashcount.fastq ]
#then 
#	echo "Skipping Overlap"
#else


if [ -s ./$File.bam ] 
then 
	echo "skipping align"
else
	$bwa mem -t $Threads $humanRefBwa "$File" | samtools sort -T $File -O bam - > $File.bam
	samtools index $File.bam 
fi

if [ $( samtools view $File.bam| head | wc -l | awk '{print $1}') -eq "0" ]; then
        echo "ERROR: BWA failed on $File .  Either the files are exactly the same of something went wrong in previous step" 
        exit 100
fi

if [ "$speed" == "veryfast" ]
then
	echo "running very fast assembly"; 
	if [ -s ./TempOverlap/$NameStub.sam.fastqd ]
	then
	        echo "skipping sam assemble"
	else
	
	        $OverlapSam <( samtools view  -F 3328 $File.bam  ) .95 20 3 ./TempOverlap/$NameStub.sam $NameStub 1 $HashList $Threads
	fi 
	if [ -s ./TempOverlap/$NameStub.final.fastqd ]
	then 
		echo "skipping second assemble"
	else 
		$OverlapHash ./TempOverlap/$NameStub.sam.fastqd .99 75 $FinalCoverage $NameStub 15 1 ./TempOverlap/$NameStub.final 1 $Threads 
	fi
	
	if [ -s ./$NameStub.overlap.hashcount.fastq ]
	then
	        echo "skipping final overlap work"
	else
	
	        $ReplaceQwithDinFASTQD ./TempOverlap/$NameStub.final.fastqd > ./$NameStub.overlap.fastqd
	        $ConvertFASTqD ./$NameStub.overlap.fastqd > ./$NameStub.overlap.fastq
	
	        echo "$AnnotateOverlap $HashList ./$NameStub.overlap.fastq TempOverlap/$NameStub.overlap.asembly.hash.fastq > ./$NameStub.overlap.hashcount.fastq"              
	              $AnnotateOverlap $HashList ./$NameStub.overlap.fastq TempOverlap/$NameStub.overlap.asembly.hash.fastq > ./$NameStub.overlap.hashcount.fastq
	fi	
else
	echo "Running full assembly"; 
	
	if [ -s ./TempOverlap/$NameStub.sam.fastqd ]  
	then 
		echo "skipping sam assemble"
	else
	     
		#$OverlapSam <( samtools view  -F 3328 $File.bam | awk '$9 > 100 || $9 < -100 || $9==0' ) .95 20 1 ./TempOverlap/$NameStub.sam $NameStub 1 $HashList $Threads
		$OverlapSam <( samtools view  -F 3328 $File.bam ) .95 20 1 ./TempOverlap/$NameStub.sam $NameStub 1 $HashList $Threads
	fi
	
	#if [ $( wc -l ./TempOverlap/$NameStub.sam.fastqd | awk '{print $1}') -eq "0" ]; then
        if [ $( head ./TempOverlap/$NameStub.sam.fastqd | wc -l | awk '{print $1}') -eq "0" ]; then
		echo "ERROR Assembly produce output for ./TempOverlap/$NameStub.sam.fastqd"
		exit 100
	fi
	
	if [ -s ./TempOverlap/$NameStub.1.fastqd ]
	then 	
		echo "skipping first overlap"
	else
		$OverlapHash ./TempOverlap/$NameStub.sam.fastqd .98 100 1 FP 20 1 ./TempOverlap/$NameStub.1 0 $Threads #> $File.overlap.out
	fi
	
	if [ $( head ./TempOverlap/$NameStub.1.fastqd | wc -l | awk '{print $1}') -eq "0" ]; then
	        echo "ERROR Assembly produce output for ./TempOverlap/$NameStub.1.fastqd"
	        exit 100
	fi
	
	if [ -s ./TempOverlap/$NameStub.2.fastqd ]
	then 
		echo "skipping second overlap"
	else
		$OverlapHash ./TempOverlap/$NameStub.1.fastqd .98 75 2 FP 20 1 ./TempOverlap/$NameStub.2 1 $Threads #>>  $File.overlap.out
	fi
	
	if [ $( head ./TempOverlap/$NameStub.2.fastqd | wc -l  | awk '{print $1}') -eq "0" ]; then
	        echo "ERROR Assembly produce output for ./TempOverlap/$NameStub.2.fastqd"
	        exit 100
	fi
	
	if [ -s ./TempOverlap/$NameStub.3.fastqd ]
	then
	        echo "skipping third overlap"
	else
	        $OverlapHash ./TempOverlap/$NameStub.2.fastqd .98 50 2 $NameStub 20 1 ./TempOverlap/$NameStub.3 1 $Threads #>>  $File.overlap.out
	fi
	if [ $( head ./TempOverlap/$NameStub.3.fastqd | wc -l | awk '{print $1}') -eq "0" ]; then
	        echo "ERROR Assembly produce output for ./TempOverlap/$NameStub.3.fastqd"
	        exit 100
	fi
	
	if [ -s ./TempOverlap/$NameStub.4.fastqd ]
	then
	        echo "skipping fourth overlap"
	else
	        $OverlapRebion2 ./TempOverlap/$NameStub.3.fastqd .98 50 $FinalCoverage  ./TempOverlap/$NameStub.4 $NameStub 1 $Threads 
		# $OverlapHash ./TempOverlap/$NameStub.3.fastqd .98 25 $FinalCoverage $NameStub 15 1 ./TempOverlap/$NameStub.4 1 $Threads #>>  $File.overlap.out
	fi
	if [ $( head ./TempOverlap/$NameStub.4.fastqd | wc -l | awk '{print $1}') -eq "0" ]; then 
	        echo "ERROR Assembly produce output for ./TempOverlap/$NameStub.4.fastqd"
	        exit 100
	fi
 

	if [ -s ./$NameStub.overlap.hashcount.fastq ]
	then 
		echo "skipping final overlap work"
	else
	
		$ReplaceQwithDinFASTQD ./TempOverlap/$NameStub.4.fastqd > ./$NameStub.overlap.fastqd
		$ConvertFASTqD ./$NameStub.overlap.fastqd > ./$NameStub.overlap.fastq
	
		echo "$AnnotateOverlap $HashList ./$NameStub.overlap.fastq TempOverlap/$NameStub.overlap.asembly.hash.fastq > ./$NameStub.overlap.hashcount.fastq"              
		      $AnnotateOverlap $HashList ./$NameStub.overlap.fastq TempOverlap/$NameStub.overlap.asembly.hash.fastq > ./$NameStub.overlap.hashcount.fastq
	fi
fi

if [ $( head ./$NameStub.overlap.hashcount.fastq | wc -l | awk '{print $1}') -eq "0" ]; then 
        echo "ERROR Assembly produce output for ./$NameStub.overlap.hashcount.fastq"
        exit 100
fi

if [ -s ./$NameStub.overlap.hashcount.fastq.bam ]
then 
	echo "skipping contig alignment" 
else
#        $bwa mem -t $Threads -Y -E 0,0 -O 6,6  -d 500 -w 500 -L 2,2 $humanRefBwa ./$NameStub.overlap.hashcount.fastq | samtools sort -T $File -O bam - > ./$NameStub.overlap.hashcount.fastq.bam
#	$bwa mem -t $Threads -Y -E 0,0 -O 6,6  -L 2,2 $humanRefBwa ./$NameStub.overlap.hashcount.fastq | samtools sort -T $File -O bam - > ./$NameStub.overlap.hashcount.fastq.bam
         $bwa mem -t $Threads -Y  $humanRefBwa ./$NameStub.overlap.hashcount.fastq | samtools sort -T $File -O bam - > ./$NameStub.overlap.hashcount.fastq.bam
	 samtools index ./$NameStub.overlap.hashcount.fastq.bam
fi

if [ $( samtools view ./$NameStub.overlap.hashcount.fastq.bam | head | wc -l | awk '{print $1}') -eq "0" ]; then
        echo "ERROR: BWA failed on ./$NameStub.overlap.hashcount.fastq.bam .  Either the files are exactly the same of something went wrong in previous step" 
        exit 100
fi

echo "string hash lookup"
#############################################################################################################



parentCRString=""


######################## BUILDING UP parent c and cR string ##############################
for parent in $ParentsJhash;
do
    parentCRString="$parentCRString -c ./Intermediates/$NameStub.overlap.asembly.hash.fastq.$parent -cR ./Intermediates/$NameStub.overlap.asembly.hash.fastq.Ref.$parent "
done

#echo "final parent String is  $parentCRString"
##########################################################################################
echo "starting overlap index"
samtools index ./$NameStub.overlap.hashcount.fastq.bam
echo "done with overlap index" 
echo ""
echo "" 
echo ""
dumbFix=$(awk '{split($1, a, ".V2"); print a[1]}' <<< $NameStub)
echo "samtools view ./$NameStub.overlap.hashcount.fastq.bam | grep -v chrUn  | $RUFUSinterpret -mob ./Intermediates/$NameStub.overlap.hashcount.fastq.MOB.sam -mod $dumbFix.Jhash.histo.7.7.dist -mQ 1 -r $humanRef -hf $HashList -o  ./$NameStub.overlap.hashcount.fastq.bam -m $MaxAlleleSize $(echo $parentCRString) -sR Intermediates/$NameStub.overlap.asembly.hash.fastq.Ref.sample -s Intermediates/$NameStub.overlap.asembly.hash.fastq.sample -e ./Intermediates/$NameStub.ref.RepRefHash" 
samtools view ./$NameStub.overlap.hashcount.fastq.bam | grep -v chrUn  | $RUFUSinterpret -mob ./Intermediates/$NameStub.overlap.hashcount.fastq.MOB.sam -mod $dumbFix.Jhash.histo.7.7.dist -mQ 1 -r $humanRef -hf $HashList -o  ./$NameStub.overlap.hashcount.fastq.bam -m $MaxAlleleSize $(echo $parentCRString) -sR Intermediates/$NameStub.overlap.asembly.hash.fastq.Ref.sample -s Intermediates/$NameStub.overlap.asembly.hash.fastq.sample -e ./Intermediates/$NameStub.ref.RepRefHash 


