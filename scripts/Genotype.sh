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
RefHash=${11}

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

RDIR=/scratch/ucgd/lustre/u0991464/Projects/CEPH.1kg.cut0.5.v5/RUFUS

OverlapHash=$RDIR/bin/Overlap
OverlapRebion2=$RDIR/bin/OverlapRegion
ReplaceQwithDinFASTQD=$RDIR/bin/ReplaceQwithDinFASTQD
ConvertFASTqD=$RDIR/bin/ConvertFASTqD.to.FASTQ
AnnotateOverlap=$RDIR/bin/AnnotateOverlap
#gkno=$RDIR/bin/gkno_launcher/gkno
bwa=$RDIR/bin/bwa/bwa
samtools=$RDIR/bin/samtools-1.6/samtools
RUFUSinterpret=$RDIR/bin/RUFUS.interpret.onlytwoParents
CheckHash=$RDIR/cloud/CheckJellyHashList.sh
OverlapSam=$RDIR/bin/OverlapSam
JellyFish=$RDIR//src/externals/jellyfish-2.2.5/bin/jellyfish



#############################################################################################################
if [ -s $NameStub.overlap.asembly.hash.fastq.ref.fastq ]
then 
	echo "skipping pull reference sequecnes"
else
	echo " $RDIR/bin/bedtools2/bin/fastaFromBed -bed <( $RDIR/bin/bedtools2/bin/bamToBed -i ./$NameStub.overlap.hashcount.fastq.bam) -fi $humanRef -fo $NameStub.overlap.asembly.hash.fastq.ref.fastq"
	$RDIR/bin/bedtools2/bin/fastaFromBed -bed <( $RDIR/bin/bedtools2/bin/bamToBed -i ./$NameStub.overlap.hashcount.fastq.bam) -fi $humanRef -fo $NameStub.overlap.asembly.hash.fastq.ref.fastq 
fi 

if [ -s $NameStub.overlap.hashcount.fastq.Jhash.tab ]
then 
	echo "skipping var hash generationr"
else
	echo "$RDIR/bin/jellyfish/bin/jellyfish count -m $HashSize -s 1G -t 20 -o ./$NameStub.overlap.hashcount.fastq.Jhash ./$NameStub.overlap.hashcount.fastq"
	$RDIR/bin/jellyfish/bin/jellyfish count -m $HashSize -s 1G -t 20 -o ./$NameStub.overlap.hashcount.fastq.Jhash ./$NameStub.overlap.hashcount.fastq
	echo "$RDIR/bin/jellyfish/bin/jellyfish dump  -c ./$NameStub.overlap.hashcount.fastq.Jhash > ./$NameStub.overlap.hashcount.fastq.Jhash.tab"
	$RDIR/bin/jellyfish/bin/jellyfish dump  -c ./$NameStub.overlap.hashcount.fastq.Jhash > $NameStub.overlap.hashcount.fastq.Jhash.tab
fi 

if [ -s ./$NameStub.overlap.asembly.hash.fastq.ref.fastq.Jhash.tab ] 
then 
	echo "skipping ref hash generation"
else
	echo " $RDIR/bin/jellyfish/bin/jellyfish count -m $HashSize -s 1G -t 20 -o ./$NameStub.overlap.asembly.hash.fastq.ref.fastq.Jhash ./$NameStub.overlap.asembly.hash.fastq.ref.fastq"
	$RDIR/bin/jellyfish/bin/jellyfish count -m $HashSize -s 1G -t 20 -o ./$NameStub.overlap.asembly.hash.fastq.ref.fastq.Jhash ./$NameStub.overlap.asembly.hash.fastq.ref.fastq
	echo "$RDIR/bin/jellyfish/bin/jellyfish -c  ./$NameStub.overlap.asembly.hash.fastq.ref.fastq.Jhash > ./$NameStub.overlap.asembly.hash.fastq.ref.fastq.Jhash.tab"
	$RDIR/bin/jellyfish/bin/jellyfish dump -c  ./$NameStub.overlap.asembly.hash.fastq.ref.fastq.Jhash > ./$NameStub.overlap.asembly.hash.fastq.ref.fastq.Jhash.tab
fi 
if [ -s $NameStub.overlap.asembly.hash.fastq.sample ]
then
        echo "skipping  $NameStub.overlap.asembly.hash.fastq.sample file already exitst"
else
	echo "starting hash lookup"
        bash $CheckHash $SampleJhash ./$NameStub.overlap.hashcount.fastq.Jhash.tab 0 > $NameStub.overlap.asembly.hash.fastq.sample &
	echo "done with hash lookup"
fi
for parent in $ParentsJhash
        do
            if [ -s $NameStub".overlap.asembly.hash.fastq."$parent ]
            then
                echo "skiping $NameStub.overlap.asembly.hash.fastq.$parent already exists"
            else
                echo "-$parent-"
                echo " bash $CheckHash $parent ./$NameStub.overlap.hashcount.fastq.Jhash.tab 0 > "$NameStub".overlap.asembly.hash.fastq."$parent" "
                bash $CheckHash $parent ./$NameStub.overlap.hashcount.fastq.Jhash.tab 0 > $NameStub".overlap.asembly.hash.fastq."$parent &
            fi
done



if [ -s $NameStub.overlap.asembly.hash.fastq.Ref.sample ]	
then 
	echo "skipping $NameStub.overlap.asembly.hash.fastq.Ref.sample"
else
	bash $CheckHash $SampleJhash  ./$NameStub.overlap.asembly.hash.fastq.ref.fastq.Jhash.tab 0 > $NameStub.overlap.asembly.hash.fastq.Ref.sample &
fi 
for parent in $ParentsJhash
        do
            if [ -s $NameStub.overlap.asembly.hash.fastq.Ref.$parent ]
            then
                echo "skipping $NameStub.overlap.asembly.hash.fastq.Ref.$parent already exitst"
            else

                echo "-$parent-"
                echo " bash $CheckHash $parent ./$NameStub.overlap.asembly.hash.fastq.ref.fastq.Jhash.tab 0 > $NameStub.overlap.asembly.hash.fastq.Ref.$parent"
                   bash $CheckHash $parent ./$NameStub.overlap.asembly.hash.fastq.ref.fastq.Jhash.tab 0 > $NameStub.overlap.asembly.hash.fastq.Ref.$parent &
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
    parentCRString="$parentCRString -c $NameStub.overlap.asembly.hash.fastq.$parent -cR $NameStub.overlap.asembly.hash.fastq.Ref.$parent "
    echo "building up parentCR string in Overlap.sh... "
    echo "$parentCRString"
done

echo "final parent String is " "$parentCRString"
##########################################################################################

if [ -s ./$NameStub.ref.RepRefHash ]
then
        echo "exclude already exists"
else
        bash $CheckHash $RefHash ./$NameStub.overlap.hashcount.fastq.Jhash.tab 1 > ./$NameStub.ref.RepRefHash
fi




wait

mkfifo check 
samtools index ./$NameStub.overlap.hashcount.fastq.bam


samtools view ./$NameStub.overlap.hashcount.fastq.bam | $RUFUSinterpret -mod $NameStub.overlap.asembly.hash.fastq.sample -mQ 8 -r $humanRef -hf $HashList -o  ./$NameStub.overlap.hashcount.fastq.bam -m 1000000 $(echo $parentCRString) -sR $NameStub.overlap.asembly.hash.fastq.Ref.sample -s $NameStub.overlap.asembly.hash.fastq.sample -e ./$NameStub.ref.RepRefHash 



grep ^# ./$NameStub.overlap.hashcount.fastq.bam.vcf> ./$NameStub.overlap.hashcount.fastq.bam.vcf.sorted.vcf
grep -v  ^# ./$NameStub.overlap.hashcount.fastq.bam.vcf | sort -k1,1 -k2,2n >> ./$NameStub.overlap.hashcount.fastq.bam.vcf.sorted.vcf
bash $RDIR/scripts/VilterAutosomeOnly ./$NameStub.overlap.hashcount.fastq.bam.vcf.sorted.vcf > ./$NameStub.overlap.hashcount.fastq.bam.vcf.sorted.vcf.autosome.vcf

$RDIR/bin/tabix/bgzip -f ./$NameStub.overlap.hashcount.fastq.bam.vcf.sorted.vcf
$RDIR/bin/tabix/tabix ./$NameStub.overlap.hashcount.fastq.bam.vcf.sorted.vcf.gz
