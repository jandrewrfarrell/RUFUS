RUFUS
=====

K-mer based variant detection. 

Andrew Farrell PhD               
Research Associate          
Department of Human Genetics              
USTAR Center for Genetic Discovery   
Eccles Institute of Human Genetics   
University of Utah School of Medicine​
15 North 2030 East, Room 7140         
Salt Lake City, UT 84112-5330        
Email: JAndrewRFarrell@gmail.com         
http://marthlab.org/

This project is still under development and its not considered stable but you are free to use and any feedback is welcome. 

## Compiler Requirements

**RUFUS requires the use of the gcc/4.9.2 compiler.**  If you are not using the gcc/4.9.2 compiler, RUFUS will not build and install properly.  Please make sure that are using the gcc/4.9.2 compiler before you proceeed.  If you are loading your compiler from a module, cmake may not run the correct compiler path.  If this is the case, try replacing the
```
cmake ../
```
command with
```
cmake ../ -DCMAKE_C_COMPILER=$(which gcc) -DCMAKE_CXX_COMPILER=$(which g++)
```
to directly provide cmake with the path to your gcc compiler.

## Ubuntu dependencies
In order for RUFUS to run on a fresh Ubuntu build, all of the following packages must be installed:

**General**
```
sudo apt-get update
sudo apt-get install python
sudo apt-get install cmake
sudo apt-get install wget
```

**GCC-4.9** (c/c++ compiler)
```
sudo apt-get install build-essential
sudo add-apt-repository ppa:ubuntu-toolchain-r/test
sudo apt-get install g++-4.9
```

**zlib** (file compression library)
```
sudo apt-get install zlib1g-dev
sudo apt-get install libbz2-dev
```

**bzlib** (bz2 file compression library)
```
sudo apt-get install libbz2-dev
sudo apt-get install liblzma-dev
```

**bc** (floating point precision library)
```
sudo apt-get install bc
```

**Curse** (terminal control library)
```
sudo apt-get install libncurses5-dev
```

## Installing RUFUS

 **1) Download**
```
git clone https://github.com/jandrewrfarrell/RUFUS.git   
```

**2) Configure**
```
cd RUFUS
./configure.sh
```

__Note:__ ./configure must be run from the RUFUS root directory

**3) Build**
```
cd bin
cmake ..
make
```

If you get errors during the installation contact me at at JAndrewRFarrell@gmail.com or submit an issue.  NOTE, to date this only works on linux machines.


## Testing RUFUS

To make sure that RUFUS was successfully built, we provide users with a test run script to run RUFUS on a small test set of data with a small test reference, and default parameters.  To test RUFUS, simply run
```
cd testRun
bash runTest.sh
```

__NOTE:__ Make sure that runTest.sh is called directly from the testRun directory, or the testRun script will not be able to find the appropriate resources.

All data for this run is contained in the resources dir, and nothing needs to be provided by the user to test RUFUS.

At the end of a successfull test run, you should see a file named

```
testRun/Child.bam.generator.V2.overlap.hashcount.fastq.bam.vcf
```

This file should contain a single varient call.  The call should look exactly as follows: 

```
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  ./Child.bam     Mother.bam      Father.bam
5:177630000     12896   X-DeNovo        T       G       25      PASS    RN=NODE_Child.bam.generator.V2_0_L273_D22:10:12::MH0;MQ=60;cigar=273M;SB=0.454545;CVT=X;HD=-1_-1_-1_-1_-1_19_-1_19_19_-1_-1_-1_-1_-1_20_20_19_-1_-1_18_-1_18_-1_-1_18_-1_-1_;AO=19;VT=X GT:DP:RO:AO     0/1:39:20:19    0/0:23:23:0     0/0:23:23:0
```

If you were unable to reproduce this call, something went wront with the RUFUS install, and you should get the test run to work before proceeding further.  If you are unable to reproduce a successfull test run, please contact me at JAndrewRFarrell@gmail.scom 



## Running 

RUFUS is primarily used to find mutations unique to a  proband sample, that are not found in the control samples

Usage: 
```
./runRufus.sh [-s|--subject <arg>] [-c|--controls][<controls-1>] ... [<controls-n>]  [-t|--threads <arg>] [-k|--kmersize <arg>] [-r|--ref <arg>] [-m|--min <arg>] 
  [-h|--help]
 ```

 
 ```
-s,--subject: bam file containing the subject of interest (REQUIRED)

-c, --controls: bam files containing the control subjects (REQUIRED)

-t,--threads: number of threads to use (REQUIRED) (min 3)

-k,--kmersize: size of Kmer to use (REQUIRED)

-r,--ref: file path to the desired reference file to create VCF (REQUIRED)

-m,--min: overwrites the minimum k-mer count to call variant (OPTIONAL, Do not provide a min unless you are sure what you want)

-h,--help: HELP!!!!!!!!!!!!!!!
```


The command line should look something like this:

```
bash runRufus.sh --subject Child.bam --controls Mother.bam Father.bam  --kmersize 25 --threads 40 --ref human_reference_v37_decoys.fa
```

or 

```
bash runRufus.sh -s Child.bam -c Mother.bam Father.bam -k 25 -t 40 -r human_reference_v37_decoys.fa
```

**NOTE** runRUFUS must be provided with atleast three threads, as it assigns threads to processes as a proportion the total number of threads provided.  It does not matter if you have fewer than three cores, but you must simply provided at least three threads.

The flags can be provided any any order.

RUFUS can take any number control files (Must provide atleast one). 

Provide all of your control files after the [-c|--controls] flag

For Example:

```
bash runRufus.sh -s tumorT1.bam -c tumorT0.bam -k 25 -t 40 -r human_reference_v37_decoys.fa
```

or

```
bash runRufus.sh -s Proband.bam -c Mother.bam Father.bam Sibling1.bam Sibling2.bam -k 25 -t 40 -r human_reference_v37_decoys.fa
```



We recommend a kmer size of 25, 40 threads, and to NOT provide RUFUS with the optional --min parameter

## Providing a reference file.

After RUFUS has identified reads containing mutant kmers, the reads must be aligned to a reference fasta file.  Any fasta file can be used as a reference, as long as the fasta file has been indexed for BWA.  If a fasta has been indexed by bwa, there will be reference files with the following extensions: pac, .ann, .abm, .bwt, sa.  In order to prepare a reference fasta for bwa, simply type:

```
bwa index -a bwtsw reference.fa
samtools faidx reference.fa
```

This will produce the BWA index files, and the fasta file index respectively.  Make sure that the bwa index files and the fasta index file are in the same directory as reference.fa








