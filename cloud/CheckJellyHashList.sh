#!/bin/sh

RDIR=/uufs/chpc.utah.edu/common/home/u0401321/finalRUFUStest/RUFUS
JellyFish=$RDIR/src/externals/jellyfish-2.2.5/bin/jellyfish
Jhash=$1
HashList=$2
MinCov=$3

$JellyFish query -s <(cat $HashList | awk '{print ">"$1"\n"$1}') $Jhash | awk -v var=$MinCov ' $2 >= var '  
