RDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

sample=$1
k=$2
Threads=$3
Min=$4
_arg_cramref=$5

    sampleFileName=$(basename "$1")
    echo "file name is" "$sampleFileName"
    sampleExtension="${sampleFileName##*.}"
    echo "file extension name is" "$sampleExtension"

    if  [[ "$sampleExtension" != "cram" ]] && [[ "$sampleExtension" != "bam" ]]  && [[ "$sampleExtension" != "generator" ]]
    then
        echo "The control bam/generator file" "$sample" " was not provided, or does not exist; killing run with non-zero exit status"
        kill -9 $$
    elif [[ "$sampleExtension" == "bam" ]]
    then
            sampleGenerator="$sampleFileName".generator
            ParentGenerators+=("$sampleGenerator")
            echo "samtools view -F 3328 $sample" > "$sampleGenerator"
            echo "You provided the control bam file" "$sample"
    elif [[ "$sampleExtension" == "cram" ]]
    then
            sampleGenerator="$sampleFileName".generator
            ParentGenerators+=("$sampleGenerator")
	    if [ "$_arg_cramref" == "" ]
            then
                echo "ERROR cram reference not provided for cram input";
                 kill -9 $$
            fi
            echo "samtools view -F 3328 -T $_arg_cramref $sample" > "$sampleGenerator"
            echo "You provided the control cram file" "$sample"    
    elif [[ "$sampleExtension" = "generator" ]]
    then
        sampleGenerator="$sampleFileName"
        ParentGenerators+=("$sampleGenerator")
        echo "You provided the control bam file" "$sample"
    fi

RunJelly=$RDIR/RunJellyForRUFUS.sh

bash $RunJelly $sampleGenerator $k $(echo $Threads -2 | bc) $Min 
