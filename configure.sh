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

if [ ! -d "bin" ]; then
mkdir bin
fi
if [ ! -d "bin/externals" ]; then
cd bin && mkdir externals && cd externals && mkdir external && cd external && touch CMakeLists.txt
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


