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
This should install everything you need to use RUFUS.  If you get errors during the gkno installation contact Alistair Ward at award@genetics.utah.edu

## Running 

we have created a number of scripts to run rufus for typical experiments, the most stable one at the moment is 
 
```
RunRUFUS.Trio.sh
```

RunRUFUS.Trio takes 6 arguments, the command line should look like this 

```
bash RunRUFUS.Trio.sh Father.generator Mother.generator Child.generator 25 40 Child.generator.outstub
```
where 25 is the kmer size to be used and 40 is the number of threads to use.  You can chage those to what ever you want. 

RUFUS works on generator files.  You will need to make a generator for each one of your input bam files.  A generator is essentially a bash script that dumps sam formatted reads to the terminal.  It is ususally a basch script with a single line but it can be more complicated.  This lets us take in numerous data formats and lets you filter your input to remove duplicate reads ect.  An example of the generator file we usually use is below.  This is a single line file that will output only the primary alignments and ignore all duplicate reads.

```
samtools view -F 3328 yourbamfile.bam
```

