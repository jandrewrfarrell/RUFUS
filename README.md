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

## Download and install
```
git clone https://github.com/jandrewrfarrell/RUFUS.git   
cd RUFUS
bash Install.sh
```
This should install everything you need to use RUFUS.  If you get errors during the installation contact me at at JAndrewRFarrell@gmail.com

## Running 

RUFUS can be run on a proband sample, with any number of control samples.

Usage: 
./runRufus.sh [-s|--subject <arg>] [-r|--ref <arg>] [-t|--threads <arg>] [-k|--kmersize <arg>] [-m|--min <arg>] 
 
[-c|--controls][<controls-1>] ... [<controls-n>] [-h|--help]
 
 
-c, --controls: bam files containing the control subjects

-s,--subject: bam file containing the subject of interest (no default)

-r,--ref: file path to the desired reference file (no default)

-t,--threads: number of threads to use (no default)

-k,--kmersize: size of Kmer to use (no default)

-m,--min: overwrites the minimum k-mer count to call variant (Do not provide a min unless you are sure what you want)

-h,--help: HELP!!!!!!!!!!!!!!!
 


The command line should look something like this:

```
bash runRufus.sh --subject child.bam --controls Mother.bam Father.bam  --kmersize 25 --threads 40 --ref human_reference_v37_decoys.fa
```

or 

```
bash runRufus.sh -s child.bam -c Mother.bam Father.bam -k 25 -t 40 -r human_reference_v37_decoys.fa
```

We recommend a kmer size of 25, 40 threads, and to not provide RUFUS with a static minimum cutoff value




