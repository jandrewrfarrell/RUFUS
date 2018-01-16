xy/*By ANDREW FARRELL 
 * the Marth Lab BC
 * 6.3 = 6.2 + testing sectioning the file 
 * now accepts nfastqd as input
 * 7 = hash based searching method for matches
 */
#include <algorithm>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <bitset>
#include <time.h>
#include<sys/stat.h>
#include<sys/types.h>

#include "Util.h"

using namespace std;

struct classcomp {
  bool operator() (const char& lhs, const char& rhs) const
  {return lhs<rhs;}
};


int main (int argc, char *argv[])
{
	int SearchHash = 30;
	int ACT = 0;
	float MinPercent ;
	int MinOverlap ;
	int MinCoverage ;			
	//cout << "you gave "<< argc << " Arguments"  << endl;
	if (argc != 2)
	{cout << "ERROR, wrong numbe of arguemnts\nCall is: FASTQD " << endl; return 0;}
		
	
	ifstream fastq;
	fastq.open (argv[1]);
	if ( fastq.is_open())
	{}	//cout << "##File Opend\n";
	else
	{
		cout << "Error, ParentHashFile could not be opened";
		return 0;
	}
	
	int lines = -1;

	string L1;
	string L2;
	string L3;
	string L4;
	string L5;
	string L6;
	int Rejects = 0;
	
	//read in fastq
	//int readcount = 0;
		int counter =0;
		//cout << "ATTENTION - Fastq+depth input detected, reading in FASTQD file \n";
		while (getline(fastq, L1))
                {
			counter++;
			if (counter % 100 == 1)
                        {
                        //        cout << "Read in " << counter << " lines and rejected " << Rejects << " reads\r";
                        }

                        getline(fastq, L2);
                        getline(fastq, L3);
                        getline(fastq, L4);
			getline(fastq, L5);
			getline(fastq, L6);
			

			cout << L1 << endl;
                        cout << L2 << endl;
			cout << L3 << endl;
			cout << L4 << endl;
                }

}
	
