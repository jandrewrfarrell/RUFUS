#!/bin/sh

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
    mkfifo $GEN.Jhash.temp
    mkfifo $GEN.fq
    bash $GEN | $RDIR/bin/PassThroughSamCheck $GEN.Jelly.chr > $GEN.fq &
    /usr/bin/time -v $JELLYFISH count --disk -m $K -L $L -s 8G -t $T -o $GEN.Jhash -C $GEN.fq
    /usr/bin/time -v $JELLYFISH histo -f -o $GEN.Jhash.histo $GEN.Jhash 
    rm $GEN.Jhash.temp
    rm $GEN.fq

	wait
fi
wait
exit
