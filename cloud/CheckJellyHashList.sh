#!/bin/sh

RDIR=/scratch/ucgd/lustre/u0991464/Projects/EIEE.testRefactor/RUFUS
Jhash=$1
HashList=$2
MinCov=$3
$RDIR/bin/jellyfish/bin/jellyfish query -s <(cat $HashList | awk '{print ">"$1"\n"$1}') $Jhash | awk -v var=$MinCov ' $2 >= var '  
