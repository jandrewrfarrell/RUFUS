/*By ANDREW FARRELL
 * RUFUS.CheckHashFilter.cpp
 * TODO: describe funciton of file
 */

#include <bitset>
#include <fstream>
#include <ios>
#include <iostream>
#include <math.h>
#include <sstream>
#include <stack>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <sys/resource.h>
#include <time.h>
#include <unistd.h>
#include <unordered_map>
#include <vector>

#include "Util.h"

using namespace std;

int main(int argc, char *argv[]) 
{
	cout << "Call is PreBuiltMutHash Mutant.Mate1.fq Mutant.Mate2.fq firstpassfile hashsize MinQ HashCountThreshold threads " << endl;
	double vm, rss, MAXvm, MAXrss;
	MAXvm = 0;
	MAXrss = 0;
	Util::process_mem_usage(vm, rss, MAXvm, MAXrss);
	cout << "VM: " << vm << "; RSS: " << rss << endl;

	int BufferSize = 240;

	cout << "Paramaters are:\n	PreBuiltMutHash = " << argv[1]
			 << "\n	Mutant.mate1.fq = " << argv[2]
			 << "\n	Mutant.mate2.fq = " << argv[3] 
			 << "\n	out stub = " << argv[4]
			 << "\n	HashSize = " << argv[5] 
			 << "\n	MinQ = " << argv[6]
			 << "\n	HashCountThreshold = " << argv[7] 
			 << "\n	Threads = " << argv[8] << endl;
	// Read in file passed to the program on the command line

	string temp = argv[5];
	int HashSize = atoi(temp.c_str());
	temp = argv[6];
	int MinQ = atoi(temp.c_str());
	temp = argv[7];
	int HashCountThreshold = atoi(temp.c_str());
	temp = argv[8];
	int Threads = atoi(temp.c_str());
	ifstream MutHashFile;
	MutHashFile.open(argv[1]);
	if (MutHashFile.is_open()) {
		cout << "Parent File open - " << argv[1] << endl;
	}	// cout << "##File Opend\n";
	else {
		cout << "Error, ParentHashFile could not be opened";
		return 0;
	}

	ifstream MutFileM1;
	cout << "MutFile.mate1 is " << argv[2] << endl;
 
	MutFileM1.open(argv[2]);
 
	cout << "here " << endl;	
 
	if (MutFileM1.is_open()) {
		cout << "##File Opend\n";
	} else {
		cout << "Error, MutFile could not be opened";
		return 0;
	}

	ifstream MutFileM2;
	cout << "MutFile.mate2 is " << argv[3] << endl;
	MutFileM2.open(argv[3]);
	
	if (MutFileM2.is_open()) {
		cout << "##File Opend\n";
	} else {
		cout << "Error, MutFile could not be opened";
		return 0; 
	}


	ofstream MutOutFileM1;
	string FirstPassFile = argv[4];
	FirstPassFile += ".Mutations.Mate1.fastq";
	MutOutFileM1.open(FirstPassFile.c_str());
	if (MutOutFileM1.is_open()) {
	} else {
		cout << "ERROR, Output file could not be opened -" << argv[4] << endl;
		return 0;
	}

	ofstream MutOutFileM2;
	FirstPassFile = argv[4];
	FirstPassFile += ".Mutations.Mate2.fastq";
	MutOutFileM2.open(FirstPassFile.c_str());
	if (MutOutFileM2.is_open()) {
	} else {
		cout << "ERROR, Output file could not be opened -" << argv[4] << endl;
		return 0;
	}

	string line;
	unordered_map<unsigned long, int> Mutations;
	cout << "Reading in pre-built hash talbe\n";
	int lines = 0;
	string L1;
	unsigned long LongHash;
	bool notdone = true;
	cout << "starting " << endl;
	cout << "	Reading in MutHashFile" << endl;

	while (getline(MutHashFile, L1)) {
		vector<string> temp;
		temp = Util::Split(L1, '\t');

		if (temp.size() == 2) {
			unsigned long b = Util::HashToLong(temp[0]);
			unsigned long revb = Util::HashToLong(Util::RevComp(temp[0]));
			Mutations.insert(pair<unsigned long, int>(b, 0));
			Mutations.insert(pair<unsigned long, int>(revb, 0));
		} else if (temp.size() == 4) {
			unsigned long b = Util::HashToLong(temp[3]);
			unsigned long revb = Util::HashToLong(Util::RevComp(temp[3]));
			Mutations.insert(pair<unsigned long, int>(b, 0));
			Mutations.insert(pair<unsigned long, int>(revb, 0));
		}
		if (temp.size() == 1) {
			temp = Util::Split(L1, ' ');
			unsigned long b = Util::HashToLong(temp[0]);
			unsigned long revb = Util::HashToLong(Util::RevComp(temp[0]));
			Mutations.insert(pair<unsigned long, int>(b, 0));
			Mutations.insert(pair<unsigned long, int>(revb, 0));
		}
	}

	MutHashFile.close();
	cout << "\nDone Hash Files\n";
	cout << "	 Mutations Hash size is " << (int)Mutations.size() << endl;
	int who = RUSAGE_SELF;
	struct rusage usage;
	int b = getrusage(RUSAGE_SELF, &usage);
	cout << "I am using " << usage.ru_maxrss << endl;

	cout << "VM: " << vm << "; RSS: " << rss << "; maxVM: " << MAXvm << "; maxRSS: " << MAXrss << endl;
	cout << "Starting Search " << endl;
	clock_t St, Et;
	St = clock();
	int found = 0;
	lines = 0;
	string BufferMate1[2400];
	string BufferMate2[2400];

	while (getline(MutFileM1, L1)) 
	{
		lines++;
		BufferMate1[0] = L1;
		getline(MutFileM1, BufferMate1[1]);
		getline(MutFileM1, BufferMate1[2]);
		getline(MutFileM1, BufferMate1[3]);

		getline(MutFileM2, BufferMate2[0]);
		getline(MutFileM2, BufferMate2[1]);
		getline(MutFileM2, BufferMate2[2]);
		getline(MutFileM2, BufferMate2[3]);




		if (lines % 10000 > 1 && (lines % (10000 + Threads) < Threads)) {
			Et = clock();
			float Dt = ((double)(Et - St)) * CLOCKS_PER_SEC;
			cout << "Read in " << lines * (BufferSize/4) << " lines: Found " << found
					 << " Reads per sec = " << (float)lines / (float)Dt << " \r";
		}

		int pos = 4;

		while (getline(MutFileM1, BufferMate1[pos])) 
		{
			getline(MutFileM2, BufferMate2[pos]);
			pos++;
			if (pos >= BufferSize) {
				break;
			}
		}

		#pragma omp parallel for shared(MutOutFileM1, MutOutFileM2) num_threads(Threads)
		for (int BuffCount = 0; BuffCount < pos; BuffCount += 4) 
		{
			int MutHashesFound = 0;
			int streak = 0;
			int start = 0;

			for (int i = start; i < BufferMate1[BuffCount + 1].length()-1 ; i++) 
			{
			 	if (((int)BufferMate1[BuffCount + 3].c_str()[i] - 33) < MinQ || (int)BufferMate1[BuffCount + 1].c_str()[i] == 78) 
				{
					streak = 0;
					start = i ;//+ HashSize;
				}
				else
					streak++;

				if (streak >= HashSize ) 
				{
					if (Mutations.count(Util::HashToLong(	BufferMate1[BuffCount + 1].substr(i-HashSize+1, HashSize))) > 0) 
					{
						MutHashesFound++;
					}
				}
			}
			if (MutHashesFound >= HashCountThreshold )
                        {
                                #pragma omp critical(MutWrite)
                                {
                                        MutOutFileM1 << BufferMate1[BuffCount] << endl //<< ":MH" << MutHashesFound << endl
                                        << BufferMate1[BuffCount + 1] << endl
                                        << BufferMate1[BuffCount + 2] << endl
                                        << BufferMate1[BuffCount + 3] << endl;
                                        found++;
                                        MutOutFileM2 << BufferMate2[BuffCount] << endl //<< ":MH" << MutHashesFoundM2 << endl
                                        << BufferMate2[BuffCount + 1] << endl
                                        << BufferMate2[BuffCount + 2] << endl
                                        << BufferMate2[BuffCount + 3] << endl;
                                }
                        }
			else
			{
				int MutHashesFoundM2 = 0;
				int streakM2 = 0;
				int startM2 = 0;
	
				for (int i = startM2; i < BufferMate2[BuffCount + 1].length()-1 ; i++) 
				{
					if (((int)BufferMate2[BuffCount + 3].c_str()[i] - 33) < MinQ || (int)BufferMate2[BuffCount + 1].c_str()[i] == 78) 
					{
						streakM2 = 0;
						startM2 = i ;//+ HashSize;
	
					}
					else
						streakM2++;
	
					if (streakM2 >= HashSize )
					{
						if (Mutations.count(Util::HashToLong(	BufferMate2[BuffCount + 1].substr(i-HashSize+1, HashSize))) > 0) {
							MutHashesFoundM2++;
						}
					}
				}	
				if (MutHashesFoundM2 >= HashCountThreshold )
                        	{
                                	#pragma omp critical(MutWrite)
                                	{
                                	        MutOutFileM1 << BufferMate1[BuffCount] << endl //<< ":MH" << MutHashesFound << endl
                                	        << BufferMate1[BuffCount + 1] << endl
                                	        << BufferMate1[BuffCount + 2] << endl
                                	        << BufferMate1[BuffCount + 3] << endl;
                                	        found++;
                                	        MutOutFileM2 << BufferMate2[BuffCount] << endl //<< ":MH" << MutHashesFoundM2 << endl
                                	        << BufferMate2[BuffCount + 1] << endl
                                	        << BufferMate2[BuffCount + 2] << endl
                                	        << BufferMate2[BuffCount + 3] << endl;
                                	}
                        	}
			}
			/*if (MutHashesFound >= HashCountThreshold || MutHashesFoundM2 >= HashCountThreshold )
			{
				#pragma omp critical(MutWrite)
				{
					MutOutFileM1 << BufferMate1[BuffCount] << endl //<< ":MH" << MutHashesFound << endl
					<< BufferMate1[BuffCount + 1] << endl
					<< BufferMate1[BuffCount + 2] << endl
					<< BufferMate1[BuffCount + 3] << endl;
					found++;
					MutOutFileM2 << BufferMate2[BuffCount] << endl //<< ":MH" << MutHashesFoundM2 << endl
					<< BufferMate2[BuffCount + 1] << endl
					<< BufferMate2[BuffCount + 2] << endl
					<< BufferMate2[BuffCount + 3] << endl;
				}
			}*/

		}
	}
	MutFileM1.close();
	MutFileM2.close();
	MutOutFileM1.close();
	MutOutFileM2.close(); 
	cout << "\nDone running RUFUS.Filter.cpp\n";
}
