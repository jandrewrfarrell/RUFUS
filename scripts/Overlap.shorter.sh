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
if [ ! -d "./TempOverlap/" ]
then 
	mkdir ./TempOverlap/
else
	echo "TempOverlap already present"
fi
if [ ! -d "./Intermediates/" ]
then
	mkdir ./Intermediates/
else
	echo "Intermediates already present"
fi
echo "Overlaping $File"

CDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"


RDIR=$CDIR/../

AddSA=$RDIR/scripts/AddSAtoReadSame.pl
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
		#$OverlapSam <( samtools view  -F 3328 $File.bam | awk '$9 > 150 || $9 < -150 '  ) .99 25 3 ./TempOverlap/$NameStub.sam $NameStub 1 $HashList $Threads
	        $OverlapSam <( samtools view  -F 3328 $File.bam  ) .99 25 3 ./TempOverlap/$NameStub.sam $NameStub 1 $HashList $Threads
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
		echo "$OverlapHash ./TempOverlap/$NameStub.sam.fastqd .98 100 1 FP 20 1 ./TempOverlap/$NameStub.1 0 $Threads"
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
#	$bwa mem -t $Threads -Y -E 0,0 -O 6,6 -d 500 -w 500  -L 2,2 $humanRefBwa ./$NameStub.overlap.hashcount.fastq | samtools sort -T $File -O bam - > ./$NameStub.overlap.hashcount.fastq.bam
        $bwa mem -t $Threads -Y  $humanRefBwa ./$NameStub.overlap.hashcount.fastq | samtools sort -T $File -O bam - > ./$NameStub.overlap.hashcount.fastq.bam
	samtools index ./$NameStub.overlap.hashcount.fastq.bam
fi

if [ $( samtools view ./$NameStub.overlap.hashcount.fastq.bam | head | wc -l | awk '{print $1}') -eq "0" ]; then
        echo "ERROR: BWA failed on ./$NameStub.overlap.hashcount.fastq.bam .  Either the files are exactly the same of something went wrong in previous step" 
        exit 100
fi

echo "string hash lookup"
#############################################################################################################
echo "staring MOB check"
if [ -s ./Intermediates/$NameStub.overlap.hashcount.fastq.MOB.sam ]
then
	echo "skipping MOB alignemnt check "
else 
	$bwa mem -t $Threads -Y -E 0,0 -O 6,6  -d 500 -w 500 -L 0,0 $MOBList ./$NameStub.overlap.hashcount.fastq | samtools sort -T $File -O sam - > ./Intermediates/$NameStub.overlap.hashcount.fastq.MOB.sam
fi 


echo "starting reference pull "
if [ -e ./Intermediates/$NameStub.overlap.asembly.hash.fastq.ref.fastq ]
then 
	echo "skipping pull reference sequecnes"
else

	bedtools getfasta -bed <( bedtools bamtobed -i ./$NameStub.overlap.hashcount.fastq.bam |  awk '{s=$2-100; if (s<0) {print $1 "\t" 0  "\t" $3+100} else {print $1 "\t" s  "\t" $3+100}}'  ) -fi $humanRef -fo ./Intermediates/$NameStub.overlap.asembly.hash.fastq.ref.fastq 

fi 

echo "starting var hash generatrion" 
if [ -e ./Intermediates/$NameStub.overlap.hashcount.fastq.Jhash ]
then 
	echo "skipping var hash generationr"
else
	echo "$JellyFish count -m $HashSize -s 1G -t 20 -o ./Intermediates/$NameStub.overlap.hashcount.fastq.Jhash ./$NameStub.overlap.hashcount.fastq"
	$JellyFish count -m $HashSize -s 1G -t 1 -o ./Intermediates/$NameStub.overlap.hashcount.fastq.Jhash ./$NameStub.overlap.hashcount.fastq
	echo "$JellyFish dump  -c ./Intermediates/$NameStub.overlap.hashcount.fastq.Jhash > ./Intermediates/$NameStub.overlap.hashcount.fastq.Jhash.tab"
	$JellyFish dump  -c ./Intermediates/$NameStub.overlap.hashcount.fastq.Jhash > ./Intermediates/$NameStub.overlap.hashcount.fastq.Jhash.tab
fi 

echo "starting ref hash generation"
if [ -s ./Intermediates/$NameStub.overlap.asembly.hash.fastq.ref.fastq.Jhash ] 
then 
	echo "skipping ref hash generation"
else
	$JellyFish count -m $HashSize -s 1G -t 1 -o ./Intermediates/$NameStub.overlap.asembly.hash.fastq.ref.fastq.Jhash ./Intermediates/$NameStub.overlap.asembly.hash.fastq.ref.fastq
	$JellyFish dump -c  ./Intermediates/$NameStub.overlap.asembly.hash.fastq.ref.fastq.Jhash > ./Intermediates/$NameStub.overlap.asembly.hash.fastq.ref.fastq.Jhash.tab
fi
 
 echo "pull hashes from sample" 
if [ -s Intermediates/$NameStub.overlap.asembly.hash.fastq.sample ]
then
        echo "skipping  Intermediates/$NameStub.overlap.asembly.hash.fastq.sample file already exitst"
else
	echo "starting hash lookup this one"
        bash $CheckHash $SampleJhash ./Intermediates/$NameStub.overlap.hashcount.fastq.Jhash.tab 0 $MaxCov> Intermediates/$NameStub.overlap.asembly.hash.fastq.sample & 
	echo "done with hash lookup"
fi

echo "pull hashes from controls" 
for parent in $ParentsJhash
        do
            if [ -s Intermediates/$NameStub.overlap.asembly.hash.fastq.$parent ]
            then
                echo "skiping Intermediates/$NameStub.overlap.asembly.hash.fastq.$parent already exists"
            else
		echo "pulling  Intermediates/$NameStub".overlap.asembly.hash.fastq."$parent"
                bash $CheckHash $parent ./Intermediates/$NameStub.overlap.hashcount.fastq.Jhash.tab 0 $MaxCov> Intermediates/$NameStub".overlap.asembly.hash.fastq."$parent &
            fi
done

wait

if [ -s Intermediates/$NameStub.overlap.asembly.hash.fastq.Ref.sample ]	
then 
	echo "skipping Intermediates/$NameStub.overlap.asembly.hash.fastq.Ref.sample"
else
	bash $CheckHash $SampleJhash  ././Intermediates/$NameStub.overlap.asembly.hash.fastq.ref.fastq.Jhash.tab 0 $MaxCov> Intermediates/$NameStub.overlap.asembly.hash.fastq.Ref.sample&
fi

for parent in $ParentsJhash
do
    if [ -s ./Intermediates/$NameStub.overlap.asembly.hash.fastq.Ref.$parent ]
    then
        echo "skipping $NameStub.overlap.asembly.hash.fastq.Ref.$parent already exitst"
    else
	
        #echo "-$parent-"
        echo "  bash $CheckHash $parent ./Intermediates/$NameStub.overlap.asembly.hash.fastq.ref.fastq.Jhash.tab 0 $MaxCov> ./Intermediates/$NameStub.overlap.asembly.hash.fastq.Ref.$parent"
        bash $CheckHash $parent ./Intermediates/$NameStub.overlap.asembly.hash.fastq.ref.fastq.Jhash.tab 0 $MaxCov> ./Intermediates/$NameStub.overlap.asembly.hash.fastq.Ref.$parent &
    	echo "uncomment this"
    fi
done
wait
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
echo "here "
if [ -s ./Intermediates/$NameStub.ref.RepRefHash ]
then
        echo "Exclude already exists"
else
	if [ -z $refHash ]
	then 
		echo "refhash not provided, skipping"
		touch  Intermediates/$NameStub.ref.RepRefHash
	else
		
		echo "this one" 
		echo "bash $CheckHash $refHash ./Intermediates/$NameStub.overlap.hashcount.fastq.Jhash.tab 0 $MaxCov > Intermediates/$NameStub.ref.RepRefHash"
		bash $CheckHash $refHash ./Intermediates/$NameStub.overlap.hashcount.fastq.Jhash.tab 0 $MaxCov> Intermediates/$NameStub.ref.RepRefHash
		echo "outa this"
	fi
fi
wait

echo "starting overlap index"
samtools index ./$NameStub.overlap.hashcount.fastq.bam
echo "done with overlap index" 
echo ""
echo "" 
echo ""
dumbFix=$(awk '{split($1, a, ".V2"); print a[1]}' <<< $NameStub)
echo "$RUFUSinterpret -mob ./Intermediates/$NameStub.overlap.hashcount.fastq.MOB.sam -mod $dumbFix.Jhash.histo.7.7.dist -mQ 20 -r $humanRef -hf $HashList -o  ./$NameStub.overlap.hashcount.fastq.bam -m $MaxAlleleSize $(echo $parentCRString) -sR Intermediates/$NameStub.overlap.asembly.hash.fastq.Ref.sample -s Intermediates/$NameStub.overlap.asembly.hash.fastq.sample -e ./Intermediates/$NameStub.ref.RepRefHash"

samtools view ./$NameStub.overlap.hashcount.fastq.bam | perl $AddSA | grep -v chrUn  | $RUFUSinterpret -mob ./Intermediates/$NameStub.overlap.hashcount.fastq.MOB.sam -mod $dumbFix.Jhash.histo.7.7.dist -mQ 10 -r $humanRef -hf $HashList -o  ./$NameStub.overlap.hashcount.fastq.bam -m $MaxAlleleSize $(echo $parentCRString) -sR Intermediates/$NameStub.overlap.asembly.hash.fastq.Ref.sample -s Intermediates/$NameStub.overlap.asembly.hash.fastq.sample -e ./Intermediates/$NameStub.ref.RepRefHash 


