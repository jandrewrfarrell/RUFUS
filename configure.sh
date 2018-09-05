#!/bin/sh

echo "update script paths"
perl -p -i -e "s/RDIR=.*\n/RDIR=$( echo $(pwd)| perl -p -i -e "s/\//\\\\\//g")\n/g"  scripts/RunJellyForRUFUS
perl -p -i -e "s/RDIR=.*\n/RDIR=$( echo $(pwd)| perl -p -i -e "s/\//\\\\\//g")\n/g"  scripts/HumanDedup.grenrator.tenplate
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

cd $RUFUS_DIR/bin/externals
if [ ! -d "samtools" ]; then
wget https://github.com/samtools/samtools/releases/download/1.9/samtools-1.9.tar.bz2; 
tar xvjf samtools-1.9.tar.bz2;
mv samtools-1.9 samtools
fi
cd samtools; ./configure --prefix=${PWD}; make; make install; cd ../;

cd $RUFUS_DIR/src/externals
echo "current pwd is ${PWD}"

#if [ ! -d "jellyfish-2.2.5" ]; then
#wget https://github.com/gmarcais/Jellyfish/releases/download/v2.2.5/jellyfish-2.2.5.tar.gz
#tar -xvf jellyfish-2.2.5.tar.gz;
#fi
#cd jellyfish-2.2.5; ./configure --prefix=$PWD/bin/jellyfish; make; make install; cd ../;

if [ -e $RUFUS_DIR/src/externals/jellyfish-2.2.5/bin/jellyfish ]
then
        echo "jellyfish already installed: skipping"
else
    if [ ! -f $RUFUS_DIR/src/externals/jellyfish-2.2.5.tar.gz ]; then
	wget https://github.com/gmarcais/Jellyfish/releases/download/v2.2.5/jellyfish-2.2.5.tar.gz
    fi
    tar -xvf jellyfish-2.2.5.tar.gz
    cd jellyfish-2.2.5
    mkdir bin
    ./configure --prefix=$RUFUS_DIR/src/externals/jellyfish-2.2.5
    make
    make install
fi


cd $RUFUS_DIR/cloud

rm -rf jellyfish-MODIFIED-merg*
if [ ! -f jellyfish-2.2.5.tar.gz ]; then
wget https://github.com/gmarcais/Jellyfish/releases/download/v2.2.5/jellyfish-2.2.5.tar.gz
fi
tar -xvf jellyfish-2.2.5.tar.gz
mv jellyfish-2.2.5/ jellyfish-MODIFIED-merge/
cp merge_files.cc jellyfish-MODIFIED-merge/jellyfish/merge_files.cc 
cd jellyfish-MODIFIED-merge/
./configure --prefix=$RUFUS_DIR/cloud/jellyfish-MODIFIED_merge
make 
make install
cd ../
cd ../



echo "PWD FOR jelly-Modified-merge is $PWD"

#if [ ! -d "jellyfish-MODIFIED-merge" ]; then
#rm -rf jellyfish-MODIFIED-merge
#mv externals/jellyfish-2.2.5.tar.gz .
#tar -xvf jellyfish-2.2.5.tar.gz 
#mv jellyfish-2.2.5/ jellyfish-MODIFIED-merge/
#cp merge_files.cc jellyfish-MODIFIED-merge/jellyfish/merge_files.cc 
#fi
#cd jellyfish-MODIFIED-merge/
#make 
#make install
#./configure --prefix=$PWD/bin/jellyfish
#cd ../../
#echo "final pwd is $PWD"
