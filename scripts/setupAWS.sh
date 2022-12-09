 sudo apt-get -y install gcc
 sudo apt-get -y install zlib
 sudo apt-get -y update
 sudo apt-get -y install build-essential
 sudo apt-get -y install zlib1g-dev
 sudo apt-get -y install libncurses5-dev
 sudo apt-get -y install bzip2
 sudo apt-get -y install libbzip2
 sudo apt-get -y install libbz2-dev
 sudo apt-get -y install liblzma-dev
 sudo apt-get -y install libcurl
 sudo apt -y install libcurl4-gnutls-dev
 sudo apt -y install libcurl4-nss-dev
 sudo apt -y install libcurl4-openssl-dev
 sudo apt-get -y install bedtools
 sudo apt-get -y install tabix
 sudo apt-get -y install samtools
 sudo apt-get -y install cmake

 git clone https://github.com/jandrewrfarrell/RUFUS.git
 cd RUFUS/
 git checkout dev
 git pull
 mkdir bin
 cd bin
 cmake ../
 make 
