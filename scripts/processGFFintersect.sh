echo "" > $1.genes
for i in $(awk '$7 == "gene"' $1 | awk '{print $13}' ); do
	echo $i; 
	gene=$(awk '{split($1, a, ";Name="); split(a[2], b, ";"); print b[1]}' <<< $i);
	echo "gene = $gene";
	id=$(awk '{split($1, a, "ID=gene:"); split(a[2], b, ";"); print b[1]}' <<< $i);
	echo "id = $id"; 
	
	#grep $id  intersect.out 

	for j in $(grep  $id  $1 |  awk '$7 == "mRNA"' | awk '{print $13}'); do 
		ts=$(awk '{split($1, a, "ID=transcript:"); split(a[2], b, ";"); print b[1]}' <<< $j);
		grep $ts $1 | awk '$7 != "mRNA"'; 
		for i in $(grep $ts $1 | awk '$7 != "mRNA"'); do echo $gene; done | sort | uniq >> $1.genes
	done 
done 

