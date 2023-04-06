#!/bin/sh

CDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
RDIR=$CDIR/../

JellyFish=$RDIR/bin/externals/jellyfish/src/jellyfish_project/bin/jellyfish
Jhash=$1
HashList=$2
MinCov=$3
MaxCov=$4

timeout 1h $JellyFish query -s <(cat $HashList | awk '{print ">"$1"\n"$1}') $Jhash | awk -v var=$MinCov ' $2 >= var ' | awk -v var=$MaxCov ' $2 <= var '
#cat $HashList | awk '{print $1}' | $JellyFish query -i $Jhash | awk -v var=$MinCov ' $2 >= var ' | awk -v var=$MaxCov ' $2 <= var '
