if [ -e ./../runRufus.sh ]
then 
	echo "starting test"
else
	echo "could not find ./../runRufus.sh, are you in the RUFUS/runTest/ directory?"
	exit
fi 

./../runRufus.sh -s Child.mate1.fastq Child.mate2.fastq -q1 Child.mate1.fastq -q2 Child.mate2.fastq -c Mother.mate1.fastq Mother.mate2.fastq -c Father.mate1.fastq Father.mate2.fastq -k 25 -t 40  -r $PWD/../resources/references/small_test_human_reference_v37_decoys.fa $1
