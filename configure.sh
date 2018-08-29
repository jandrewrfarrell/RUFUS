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


#if [ ! -s "samtools-1.9.tar.bz2" ]; then
#wget https://github.com/samtools/samtools/releases/download/1.9/samtools-1.9.tar.bz2; 
#tar xvjf samtools-1.9.tar.bz2;
#fi
#cd samtools-1.9; ./configure --prefix=${PWD}; make; make install; cd ../;

echo "current pwd is ${PWD}"

cd bin/externals/

if [ ! -d "jellyfish-2.2.5" ]; then
wget https://github.com/gmarcais/Jellyfish/releases/download/v2.2.5/jellyfish-2.2.5.tar.gz
tar -xzf jellyfish-2.2.5.tar.gz;
fi
cd jellyfish-2.2.5; ./configure --prefix=${PWD}; make; make install; cd ../;

#echo "PWD is ${PWD}"

if [ ! -d "jellyfish-MODIFIED-merge" ]; then
cp -r jellyfish-2.2.5/ jellyfish-MODIFIED-merge/
cp ../merge_files.cc jellyfish-MODIFIED-merge/jellyfish/merge_files.cc 
fi
cd jellyfish-MODIFIED-merge/
make 
make install
./configure --prefix=$(pwd)
cd ../../../
echo "final pwd is ${PWD}"
