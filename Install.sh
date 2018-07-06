#install script for rufus 
echo "update script paths"
perl -p -i -e "s/RDIR=.*\n/RDIR=$( echo $(pwd)| perl -p -i -e "s/\//\\\\\//g")\n/g" scripts/RunRUFUS.1000G.withDupRemoce.sh
perl -p -i -e "s/RDIR=.*\n/RDIR=$( echo $(pwd)| perl -p -i -e "s/\//\\\\\//g")\n/g"  scripts/OverlapBashMultiThread.sh
perl -p -i -e "s/RDIR=.*\n/RDIR=$( echo $(pwd)| perl -p -i -e "s/\//\\\\\//g")\n/g"  scripts/RunJellyForRUFUS
perl -p -i -e "s/RDIR=.*\n/RDIR=$( echo $(pwd)| perl -p -i -e "s/\//\\\\\//g")\n/g"  scripts/HumanDedup.grenrator.tenplate
perl -p -i -e "s/RDIR=.*\n/RDIR=$( echo $(pwd)| perl -p -i -e "s/\//\\\\\//g")\n/g"  scripts/RunTumor.sh
perl -p -i -e "s/RDIR=.*\n/RDIR=$( echo $(pwd)| perl -p -i -e "s/\//\\\\\//g")\n/g"  cloud/RunJellyForRUFUS
perl -p -i -e "s/RDIR=.*\n/RDIR=$( echo $(pwd)| perl -p -i -e "s/\//\\\\\//g")\n/g"  cloud/CheckJellyHashList.sh
perl -p -i -e "s/RDIR=.*\n/RDIR=$( echo $(pwd)| perl -p -i -e "s/\//\\\\\//g")\n/g"  scripts/RunRUFUS.Trio.sh
perl -p -i -e "s/RDIR=.*\n/RDIR=$( echo $(pwd)| perl -p -i -e "s/\//\\\\\//g")\n/g"  scripts/Overlap.sh
perl -p -i -e "s/RDIR=.*\n/RDIR=$( echo $(pwd)| perl -p -i -e "s/\//\\\\\//g")\n/g"  scripts/OverlapBashMultiThread.trio.sh
perl -p -i -e "s/RDIR=.*\n/RDIR=$( echo $(pwd)| perl -p -i -e "s/\//\\\\\//g")\n/g"  scripts/RunRUFUS.1000G.sh 
perl -p -i -e "s/RDIR=.*\n/RDIR=$( echo $(pwd)| perl -p -i -e "s/\//\\\\\//g")\n/g"  scripts/RunRUFUS.sarcoma.sh
perl -p -i -e "s/RDIR=.*\n/RDIR=$( echo $(pwd)| perl -p -i -e "s/\//\\\\\//g")\n/g"  cloud/RunJellyForRUFUS.UGP
perl -p -i -e "s/RDIR=.*\n/RDIR=$( echo $(pwd)| perl -p -i -e "s/\//\\\\\//g")\n/g"  runRufus.sh



RUFUS_DIR=$(pwd)

echo "bulding rufus executables"
g++ -g src/AnnotateOverlap.cpp src/Util.cpp -o bin/AnnotateOverlap -std=gnu++11 -O3
g++ -g src/ConvertFASTqD.to.FASTQ.cpp src/Util.cpp -o bin/ConvertFASTqD.to.FASTQ -std=gnu++11 -O3
g++ -g src/ModelDist.cpp src/Util.cpp -o bin/ModelDist  -std=gnu++11 -fopenmp 
g++ -g src/Overlap.cpp src/Util.cpp -o bin/Overlap -std=gnu++11 -fopenmp -O3
g++ -g src/OverlapRegion.cpp src/Util.cpp -o bin/OverlapRegion -std=gnu++11 -fopenmp -O3
g++ -g src/OverlapSam.cpp src/Util.cpp -o bin/OverlapSam -std=gnu++11 -fopenmp -O3
g++ -g src/ReplaceQwithDinFASTQD.cpp src/Util.cpp -o bin/ReplaceQwithDinFASTQD -std=gnu++11  -O3
g++ -g ./src/RUFUS.Filter.cpp src/Util.cpp -o ./bin/RUFUS.Filter -std=gnu++11 -fopenmp -O3
g++ -g src/RUFUS.Build.cpp src/Util.cpp -o bin/RUFUS.Build -std=gnu++11 -fopenmp -O3
g++ -g ./src/RUFUS.interpret.cpp src/Util.cpp  ./src/include/* -o ./bin/RUFUS.interpret -std=gnu++11 -O3  -fpermissive
g++ -g ./src/RUFUS.1kg.filter.cpp src/Util.cpp -o ./bin/RUFUS.1kg.filter -std=gnu++11 -fopenmp -O3

echo "bulding external programs"
cd src/externals/


if [ -e ./jellyfish-2.2.5/bin/jellyfish ]
then
        echo "jellyfish already installed: skipping"
else

        tar -xvf jellyfish-2.2.5.tar.gz
        cd jellyfish-2.2.5
        mkdir bin
        ./configure --prefix=$RUFUS_DIR/bin/jellyfish
        make
        make install
        cd ..
fi

cd $RUFUS_DIR/bin/

#Make RufAlu
echo "!!!!!!!!!!!!!!!!! MAKING RUFALU IN BIN DIR !!!!!!!!!!!!!!!!!!!!!!"

if [ ! -d RufAlu ]; then
    git clone https://github.com/WilliamRichards2017/RufAlu.git
    cd RufAlu
    mkdir -p externals/external
    touch externals/external/CMakeLists.txt
    cd bin
    cmake ..
    make
    cd ../../
else 
    cd RufAlu/bin
    cmake ..
    make
    cd ../../
fi

echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"


#Make BWA

if [ ! -d bwa ]; then
    git clone https://github.com/lh3/bwa.git 
fi
    cd bwa
    make 
    cd ../
    
    #Make BedTools
if [ ! -d bedtools2 ]; then
    wget https://github.com/arq5x/bedtools2/releases/download/v2.25.0/bedtools-2.25.0.tar.gz
    tar -zxvf bedtools-2.25.0.tar.gz
fi

cd bedtools2
make
cd ../

#Make Samtools

if [ ! -d samtools-1.6 ]; then
    wget https://github.com/samtools/samtools/releases/download/1.6/samtools-1.6.tar.bz2
    tar -xvf samtools-1.6.tar.bz2
fi

cd samtools-1.6/
./configure --prefix=$(pwd)
make 
make install 
cd ../


cd ../

cd cloud
rm -rf jellyfish-MODIFIED-merge
tar -xvf jellyfish-2.2.5.tar.gz
mv jellyfish-2.2.5/ jellyfish-MODIFIED-merge/
cp merge_files.cc jellyfish-MODIFIED-merge/jellyfish/merge_files.cc 
cd jellyfish-MODIFIED-merge/
./configure --prefix=$(pwd)
make 
make install
cd ../
cd ../

