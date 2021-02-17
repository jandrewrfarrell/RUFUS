if [ -e ./../runRufus.sh ]
then 
	echo "starting test"
else
	echo "could not find ./../runRufus.sh, are you in the RUFUS/runTest/ directory?"
	exit
fi 

./../runRufus.sh -s $PWD/Child.bam -c $PWD/Mother.bam -c $PWD/Father.bam -k 25 -t 40  -r $PWD/../resources/references/small_test_human_reference_v37_decoys.fa $1
