#!/bin/sh

RDIR=/uufs/chpc.utah.edu/common/home/u0991464/d1/home/farrelac/RUFUS/
Jhash=$1
HashList=$2
MinCov=$3
$RDIR/bin/jellyfish/bin/jellyfish query -s <(cat $HashList  | awk '{print ">"$1"\n"$1}') $Jhash |  awk -v var=$MinCov ' $2 >= var '   
