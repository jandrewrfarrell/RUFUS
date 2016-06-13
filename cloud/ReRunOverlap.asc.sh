date
Parent1Generator=$1
Parent2Generator=$2
SiblingGenerator=$3
ProbandGenerator=$4
K=$5
Threads=$6
Out=$7

echo "You gave
Parent1Generator=$1
Parent2Generator=$2
SiblingGenerator=$2
ProbandGenerator=$4
K=$5
Threads=$6
Out=$7
"

if [ -z "$Out" ]
then
        echo "out file not specified"
        exit
fi

RDIR=/home/ubuntu/work/RUFUS
RUFUSmodel=$RDIR/bin/ModelDist
RUFUSbuild=$RDIR/bin/RUFUS.Build
RUFUSfilter=$RDIR/bin/RUFUS.Filter
DeDupDump=$RDIR/scripts/HumanDedup.grenrator.tenplate
PullSampleHashes=$RDIR/cloud/CheckJellyHashList.sh
RUFUS1kgFilter=$RDIR/bin/RUFUS.1kg.filter
RunJelly=$RDIR/cloud/RunJellyForRUFUS


MutantMinCov=4

aws s3 --region us-east-1 cp s3://marthlab.rufus/ASC.scripts/$Parent1Generator ./
aws s3 --region us-east-1 cp s3://marthlab.rufus/ASC.scripts/$Parent2Generator ./
aws s3 --region us-east-1 cp s3://marthlab.rufus/ASC.scripts/$SiblingGenerator ./
aws s3 --region us-east-1 cp s3://marthlab.rufus/ASC.scripts/$ProbandGenerator ./


aws s3 --region us-east-1 cp s3://marthlab.rufus/ASC.out/$Out/$ProbandGenerator.Mutations.fastq ./
aws s3 --region us-east-1 cp s3://marthlab.rufus/ASC.out/$Out/$ProbandGenerator.V2.overlap.hashcount.fastq.bam ./
aws s3 --region us-east-1 cp s3://marthlab.rufus/ASC.out/$Out/$ProbandGenerator.k$K_c$MutantMinCov.HashList ./
aws s3 --region us-east-1 cp s3://marthlab.rufus/ASC.out/$Out/$ProbandGenerator.V2.overlap.asembly.hash.fastq.p1 ./
aws s3 --region us-east-1 cp s3://marthlab.rufus/ASC.out/$Out/$ProbandGenerator.V2.overlap.asembly.hash.fastq.p2 ./
aws s3 --region us-east-1 cp s3://marthlab.rufus/ASC.out/$Out/$ProbandGenerator.V2.overlap.asembly.hash.fastq.p3 ./



if [ -e $ProbandGenerator.V2.overlap.hashcount.fastq.bam.vcf ]
then 
	echo "skipping overlap"
else

File=$ProbandGenerator.Mutations.fastq
FinalCoverage=5
NameStub=$ProbandGenerator
HashList=$ProbandGenerator.k$MutantMinCov.HashList


RUFUSinterpret=$RDIR/bin/RUFUS.interpret
humanRef=$RDIR/bin/gkno_launcher/resources/homo_sapiens/build_37_version_3/human_reference_v37_decoys.fa
CheckHash=$RDIR/cloud/CheckJellyHashList.sh

$samtools view ./$NameStub.overlap.hashcount.fastq.bam | $RUFUSinterpret -r $humanRef -hf $HashList -o  ./$NameStub.overlap.hashcount.fastq.bam -m 100000000 -c $NameStub.overlap.asembly.hash.fastq.p1 -c $NameStub.overlap.asembly.hash.fastq.p2 -c $NameStub.overlap.asembly.hash.fastq.p3	

fi

MutantMinCov=5
aws s3 --region us-east-1 cp s3://marthlab.rufus/ASC.out/$Out/$SiblingGenerator.V2.overlap.hashcount.fastq ./
aws s3 --region us-east-1 cp s3://marthlab.rufus/ASC.out/$Out/$SiblingGenerator.V2.overlap.hashcount.fastq.bam ./
aws s3 --region us-east-1 cp s3://marthlab.rufus/ASC.out/$Out/$SiblingGenerator.k$K_c$MutantMinCov.HashList ./
aws s3 --region us-east-1 cp s3://marthlab.rufus/ASC.out/$Out/$SiblingGenerator.V2.overlap.asembly.hash.fastq.p1 ./
aws s3 --region us-east-1 cp s3://marthlab.rufus/ASC.out/$Out/$SiblingGenerator.V2.overlap.asembly.hash.fastq.p2 ./
aws s3 --region us-east-1 cp s3://marthlab.rufus/ASC.out/$Out/$SiblingGenerator.V2.overlap.asembly.hash.fastq.p3 ./

if [ -e $SiblingGenerator.V2.overlap.hashcount.fastq.bam.vcf ]
then
	echo "skipping overlap"
else

File=$SiblingGenerator.Mutations.fastq
FinalCoverage=5
NameStub=$SiblingGenerator
HashList=$SiblingGenerator.k$MutantMinCov.HashList


RUFUSinterpret=$RDIR/bin/RUFUS.interpret
humanRef=$RDIR/bin/gkno_launcher/resources/homo_sapiens/build_37_version_3/human_reference_v37_decoys.fa
CheckHash=$RDIR/cloud/CheckJellyHashList.sh

$samtools view ./$NameStub.overlap.hashcount.fastq.bam | $RUFUSinterpret -r $humanRef -hf $HashList -o  ./$NameStub.overlap.hashcount.fastq.bam -m 100000000 -c $NameStub.overlap.asembly.hash.fastq.p1 -c $NameStub.overlap.asembly.hash.fastq.p2 -c $NameStub.overlap.asembly.hash.fastq.p3




fi

for i in *vcf*; do aws s3  --region us-east-1 cp $i s3://marthlab.rufus/ASC.out/$Out/ ; done 
for i in *bam; do aws s3  --region us-east-1 cp $i s3://marthlab.rufus/ASC.out/$Out/ ; done
for i in *generator.V2.overlap.asembly.hash.fastq.*; do aws s3  --region us-east-1 cp $i s3://marthlab.rufus/ASC.out/$Out/ ; done

echo "done with everything "



