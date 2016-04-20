Jhash=$1
HashList=$2
../RUFUS/bin/jellyfish/bin/jellyfish query -s <(cat $HashList  | awk '{print ">"$1"\n"$1}') $Jhash | awk ' $2 >= 7 '  
