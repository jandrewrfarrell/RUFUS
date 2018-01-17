
/*This vertion imcorporates both copy number and mutation detection in 1
 *  * it needs to be run in two sptes, first the build, then the filter
 *   * it is split up to allow distribution to a cluster */

#include <unistd.h>
#include <ios>

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unordered_map>
#include <bitset>
#include <time.h>
#include <sys/resource.h>
#include <stack>

#include "Util.h"

using namespace std;

int main (int argc, char *argv[])
{
	cout << "Call is PreBuiltMutHash Mutant.fq  firstpassfile hashsize MinQ HashCountThreshold window threads " << endl;
	double vm, rss, MAXvm, MAXrss;
	MAXvm = 0;
	MAXrss = 0;
	Util::process_mem_usage(vm, rss, MAXvm, MAXrss);
   	cout << "VM: " << vm << "; RSS: " << rss << endl;	
	
	
	int BufferSize = 1000;

 cout << "Paramaters are:\n  PreBuiltMutHash = " << argv[1] << "\n  Mutant.fq = " << argv[2] << "\n  out stub = " << argv[3] <<"\n  HashSize = " << argv[4] << "\n  MinQ = " << argv[5] << "\n  HashCountThreshold = " << argv[6]  << "\n  Window = " << argv[7] << "\n  Threads = " << argv[8] 
<< endl;
		// Read in file passed to the program on the command line
	
	string temp = argv[4];
        int HashSize = atoi(temp.c_str());	
	 temp = argv[5];
	int MinQ = atoi(temp.c_str());
	temp = argv[6];
	int HashCountThreshold = atoi(temp.c_str());
	temp = argv[7];
	int Window = atoi(temp.c_str());
	temp = argv[8];
	int Threads = atoi(temp.c_str());
	ifstream MutHashFile;		
	MutHashFile.open (argv[1]);
	if ( MutHashFile.is_open())
	{ cout << "Parent File open - " << argv[1] << endl;}	//cout << "##File Opend\n";
	else
	{
		cout << "Error, ParentHashFile could not be opened";
		return 0;
	}
	
	string filename = argv[2];
	ifstream MutFile;
 	if (filename == "stdin")
        {
                cout << "MutFile is STDIN" << endl;
                MutFile.open ("/dev/stdin");
        }
        else
        {
                cout << "MutFile is " << argv[2] << endl;
                MutFile.open (argv[2]);
        }
        if ( MutFile.is_open())
        {      cout << "##File Opend\n";}
        else
        {
                cout << "Error, MutFile could not be opened";
                return 0;
        }

	ofstream MutOutFile;
	string FirstPassFile = argv[3]; 
	FirstPassFile += ".Mutations.fastq";
	MutOutFile.open (FirstPassFile.c_str());
	if (MutOutFile.is_open())
	{}
	else
	{	
		cout << "ERROR, Output file could not be opened -" << argv[3] << endl;
		return 0;
	}

	
	string line;
	unordered_map  <unsigned long, int> Mutations;

	//std::unordered_map<std::bitmap, bool> ParentHashes;
	cout << "Reading in pre-built hash talbe\n";

	int lines = 0;

	string L1;
	string L2;
	string L3;
	string L4;
 
	unsigned long LongHash;

	bool notdone = true;
	
	cout << "starting " << endl;
	cout << "  Reading in MutHashFile" << endl;
	while (getline(MutHashFile, L1))
	{
		//cout << L1 << endl;
		vector<string>  temp;
		temp = Util::Split(L1, '\t');
		cout << temp[0] << " and " << temp[1] << endl;
		if (temp.size() == 2 && temp[0].size() == HashSize)
		{
		  unsigned long b = Util::HashToLong(temp[0].c_str());
		int c = atoi(temp[1].c_str()); 
		//cout << temp[0] << "\t" << temp[1] << endl << b << "\t" << LongToHash(b, 18) <<  endl<<endl;
		Mutations.insert(pair<unsigned long, int > (b, c));
		//break;
		}
	}
	MutHashFile.close();
	
	cout << "\nDone Hash Files\n";
	cout << "   Mutations Hash size is " << (int) Mutations.size() << endl;
	 
	
	int who = RUSAGE_SELF;
        struct rusage usage;
        int b = getrusage(RUSAGE_SELF, &usage);
        cout << "I am using " << usage.ru_maxrss << endl;

	cout << "VM: " << vm << "; RSS: " << rss <<"; maxVM: " << MAXvm << "; maxRSS: " << MAXrss << endl;
//	return 0;
	
	cout << "Starting Search " << endl;
	clock_t St,Et;
        St = clock();

	int found = 0;
	lines = 0;
	while (getline(MutFile, L1))
	{
		lines++;
		getline(MutFile, L2);
                getline(MutFile, L3);
                getline(MutFile, L4);
		
		vector<string> Buffer;
		Buffer.push_back(L1);
	 	Buffer.push_back(L2);
	 	Buffer.push_back(L3);
	 	Buffer.push_back(L4);
		 if (lines % 10000 > 1 && (lines % (10000+Threads) < Threads))
                 {
                        Et = clock();
                        float Dt = ((double) (Et - St)) * CLOCKS_PER_SEC;
                 	cout << "Read in " << lines << " lines: Found " << found << " Reads per sec = " << (float)lines/(float)Dt << " \r";
                 }
		int stackCount=0;
		while (Buffer.size()<Threads*4 && getline(MutFile, L1) )
		{
			lines++;
			getline(MutFile, L2);
                	getline(MutFile, L3);
                	getline(MutFile, L4);
			
			Buffer.push_back(L1);
			Buffer.push_back(L2);
                	Buffer.push_back(L3);
                	Buffer.push_back(L4);
		//	 cout << "adding  \n   " << L1 << "\n   " << L2 << "\n   " << L3 << "\n   " << L4 << endl;
		}
		//cout << " Buffer has " << Buffer.size() << " lines\n";
		#pragma omp parallel for shared(MutOutFile) num_threads(Threads)
		for (int BuffCount = 0; BuffCount < Buffer.size(); BuffCount+=4)
		{
			string B1, B2, B3, B4;
			vector<int> positions; 
			#pragma omp critical (update)
			{
				B1 = Buffer[BuffCount];
				B2 = Buffer[BuffCount +1];
				B3 = Buffer[BuffCount +2];
				B4 = Buffer[BuffCount +3];
			}
				
			int rejected = 0;
			int MutHashesFound = 0;
			cout << "working on " << B1 << "\n" << B2 << endl;
			for (int i=0; i<B2.length()-HashSize; i++)
			{
				cout << B2.c_str()[i] << "\t"; 
				string hash = B2.substr (i,HashSize);
				cout << hash << "\t"; 
				string Qhash = B4.substr (i,HashSize);	
				bool good = true;
				for (int j=0; j<HashSize; j++)
                        	{
	                               	int B = hash.c_str()[j];
        	                        if ( B == 78)
                	                {
                        	                good = false;
						cout << "N skipping frin " << i << " to " << i+j << endl;
						i = i+j;
                                	        break;
                                	}
					else
					{
						B = Qhash.c_str()[j];
						if (((int)B - 33) < MinQ)
						{
							good = false;
							cout << "Low q = " << (int)B - 33 << " skipping frin " << i << " to " << i+j << endl;
							i = i+j;
							break;
						}
					}
                        	}	
		 		if (good)
                        	{
				  unsigned long LongHash = Util::HashToLong(hash);
					
					if (Mutations.count(LongHash) > 0)
	                                {
	                                	MutHashesFound++;
						positions.push_back(i);
						cout << "found " << LongHash << "\t" << Util::LongToHash(LongHash, 25) << "\t" << Mutations[LongHash]; 
					}
					else if (Mutations.count(Util::HashToLong(Util::RevComp(hash))) > 0)
	                        	{
				
						MutHashesFound++;
                                                positions.push_back(i);
						cout << "found " << LongHash << "\t" << Util::LongToHash(Util::HashToLong(Util::RevComp(hash)), 25) << "\t" << Mutations[Util::HashToLong(Util::RevComp(hash))];	
					}
				}
	                        else
	                        {
					rejected++;
	                        }
				cout << endl;
			}
			cout << endl;	
			if (MutHashesFound >= HashCountThreshold and rejected < (B2.length()/2))
			{
				/*cout << B2 << endl;
				for (int i = 0; i<positions.size(); i++)
				{cout << "-" << positions[i]	;}
				cout << endl;*/
				int MaxCounter = -1;
                        	for (int i = 0; i<positions.size()-1; i++)
                        	{
                                	int counter = 1;
                                	for(int j = i+1; j<positions.size(); j++)
                                	{
                                	        if (positions[j] - positions[i] <= Window)
                                	        {counter++;}
						//else{break;}
                        	        }	
					//cout << "-"<< counter;
                	                if (counter > MaxCounter)
        	                        {MaxCounter = counter;}
	                        }
				//cout << endl;
				if (MaxCounter >= HashCountThreshold )
				{ 
					cout << "   kept - " << MaxCounter << endl;
					#pragma omp critical(MutWrite)
					{MutOutFile << B1 << ":MH" << MutHashesFound << endl << B2 << endl << B3 << endl << B4 << endl;}
					found ++;
				}
				//else
				//{cout << "    rejected " << endl;}
				
			}
		}
	}
	cout << "\nDone\n";
                     
	MutFile.close();
	MutOutFile.close();
	
	cout << "\nreally dont\n";
}
	
