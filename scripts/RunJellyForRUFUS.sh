#!/bin/sh
set -e
GEN=$1
K=$2
T=$3
L=$4

CDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
RDIR=$CDIR/../


JELLYFISH="$RDIR/bin/externals/jellyfish/src/jellyfish_project/bin/jellyfish"
SORT="$RDIR/scripts/sort"

if [ -e "$GEN.Jhash" ]
then
	echo "Skipping jelly, $GEN.Jhash alreads exists"
else
	echo "Running jellyfish for $GEN"
	if [ -e $GEN.Jhash.temp ]; then 
		rm $GEN.Jhash.temp
	fi
	mkfifo $GEN.Jhash.temp
	if [ -e $GEN.fq ]; then 
		rm $GEN.fq
	fi
	mkfifo $GEN.fq
	bash $GEN | $RDIR/bin/PassThroughSamCheck $GEN.Jelly.chr > $GEN.fq &
	$JELLYFISH count --disk -m $K -L $L -s 8G -t $T -o $GEN.Jhash -C $GEN.fq
	rm $GEN.Jhash.temp
	rm $GEN.fq

	wait
fi

if [ ! -s  $GEN.Jhash.histo ]; then 
	$JELLYFISH histo -f -o $GEN.Jhash.histo $GEN.Jhash
fi
if [ $(awk '$2 > 0' $GEN.Jhash.histo | wc -l ) -eq "0" ]; then  
	echo "ERROR: jellyfish failed on the file $GEN"
	exit 100 
fi


exit
