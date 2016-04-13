
perl -p -i -e "s/RDIR=.*\n/RDIR=$( echo $(pwd)| perl -p -i -e "s/\//\\\\\//g")/g" scripts/RunRUFUS.1000G.withDupRemoce.sh 
perl -p -i -e "s/RDIR=.*\n/RDIR=$( echo $(pwd)| perl -p -i -e "s/\//\\\\\//g")/g"  scripts/OverlapBashMultiThread.sh 


g++ src/AnnotateOverlap.cpp -o bin/AnnotateOverlap -std=gnu++0x
g++ src/ConvertFASTqD.to.FASTQ.cpp -o bin/ConvertFASTqD.to.FASTQ 
g++ src/ModelDist.cpp -o bin/ModelDist -fopenmp 
g++ src/Overlap.cpp -o bin/Overlap -std=gnu++0x -fopenmp 
g++ src/OverlapRegion.cpp -o bin/OverlapRegion -std=gnu++0x -fopenmp 
g++ src/ReplaceQwithDinFASTQD.cpp -o bin/ReplaceQwithDinFASTQD
g++ ./src/RUFUS.Filter.cpp -o ./bin/RUFUS.Filter -std=gnu++0x -fopenmp 
g++ src/RUFUS.Build.cpp -o bin/RUFUS.Build 
g++ ./src/RUFUS.interpret.cpp ./src/include/* -o ./bin/RUFUS.interpret -std=gnu++0x



cd src/externals/
tar -xvf jellyfish-2.2.5.tar.gz 
cd jellyfish-2.2.5
mkdir bin
./configure --prefix=/uufs/chpc.utah.edu/common/home/u0991464/d1/home/farrelac/bin/RUFUS_git/RUFUS/src/externals/jellyfish-2.2.5/bin
make
make install 
cd ../../../
