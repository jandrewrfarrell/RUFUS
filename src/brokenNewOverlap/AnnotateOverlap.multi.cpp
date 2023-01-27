/*This vertion imcorporates both copy number and mutation detection in 1
 *	* it needs to be run in two sptes, first the build, then the filter
 *	 * it is split up to allow distribution to a cluster */

#include <bitset>
#include <fstream>
#include <unistd.h>
#include <ios>
#include <iostream>
#include <math.h>
#include <sstream>
#include <stack>
#include <string>
#include <stdio.h>
#include <stdlib.h>
#include <sys/resource.h>
#include <time.h>
#include <unordered_map>
#include <vector>

#include "Util.h"

using namespace std;
class fasta{
	public: 

	string name; 
	string seq; 
	string sep; 
	string qual;

};

int main(int argc, char *argv[]) {
	int HashSize;
	double vm, rss, MAXvm, MAXrss;
	MAXvm = 0;
	MAXrss = 0;
	Util::process_mem_usage(vm, rss, MAXvm, MAXrss);
	int BufferSize = 1000;
	ifstream MutHashFile;
	MutHashFile.open(argv[1]);
	string filename = argv[2];
	ifstream MutFile;

	if (filename == "stdin") {
		MutFile.open("/dev/stdin");
	} else {
		MutFile.open(argv[2]);
	}

	string line;
	unordered_map<unsigned long, unsigned int> Mutations;
	ofstream HashOut;
	HashOut.open(argv[3]);
	string L1;
	string L2;
	string L3;
	string L4;
	unsigned long LongHash;
	bool notdone = true;

	while (getline(MutHashFile, L1)) {
		vector<string> temp;
		temp = Util::Split(L1, '\t');

		if (temp.size() == 2) {
			unsigned long b = Util::HashToLong(temp[0].c_str());
			Mutations.insert(
					 pair<unsigned long, unsigned int>(b, atoi(temp[1].c_str())));
			HashSize = temp[0].length();
		} else if (temp.size() == 4) {
			unsigned long b = Util::HashToLong(temp[3].c_str());
			Mutations.insert(
					 pair<unsigned long, unsigned int>(b, atoi(temp[2].c_str())));
			HashSize = temp[3].length();
		} else if (temp.size() == 1) {
			temp = Util::Split(L1, ' ');
			unsigned long b = Util::HashToLong(temp[0].c_str());
			Mutations.insert(
					 pair<unsigned long, unsigned int>(b, atoi(temp[1].c_str())));
			HashSize = temp[0].length();
		}
	}
	MutHashFile.close();
	int who = RUSAGE_SELF;
	struct rusage usage;
	int b = getrusage(RUSAGE_SELF, &usage);
	clock_t St, Et;
	St = clock();
	int found = 0;
	int size = 10000; 
	 
	while (getline(MutFile, L1)) {
		vector <fasta> buffer; 
		getline(MutFile, L2);
		getline(MutFile, L3);
		getline(MutFile, L4);
		
		fasta temp; 
		temp.name = L1;
		temp.seq = L2;
		temp.sep = L3;
		temp.qual = L4;
		buffer.push_back(temp); 
		int count = 1;
		while (count < size && getline(MutFile, L1)) 
		{
			getline(MutFile, L2);
			getline(MutFile, L3);
			getline(MutFile, L4);

			fasta temp;
			temp.name = L1;
			temp.seq = L2;
			temp.sep = L3;
			temp.qual = L4;
			buffer.push_back(temp);
			count++; 
		}
		#pragma omp parallel for  num_threads(50) 
		for (int b = 0; b < buffer.size(); b++)
		{
			fasta entry = buffer[b];
			unordered_map<unsigned long, unsigned int> LMutations;
			LMutations = Mutations; 
			vector<int> positions;
			vector<int> HashPos;
			for (int i = 0; i < entry.seq.size(); i++) {
				HashPos.push_back(0);
			}
			int rejected = 0;
			int MutHashesFound = 0;
			for (int i = 0; i < entry.seq.length() - HashSize; i++) {
				string hash = entry.seq.substr(i, HashSize);
				string Qhash = entry.qual.substr(i, HashSize);
				bool good = true;

				for (int j = 0; j < HashSize; j++) {
					int B = hash.c_str()[j];
					if (B == 78) {
						good = false;
						break;
					}
		
					
					int C = Qhash.c_str()[j];
					if ((int)C -33 < 3)
						good = false; 
				}
				if (good) {
					unsigned long LongHash = Util::HashToLong(hash);
					if (LMutations.count(LongHash) > 0) {
						positions.push_back(i);
	
						for (int j = i; j < i + HashSize; j++) {
							HashPos[j]++;
						}
					} else if (LMutations.count(Util::HashToLong(Util::RevComp(hash))) > 0) {
						positions.push_back(i);
	
						for (int j = i; j < i + HashSize; j++) {
							HashPos[j]++;
						}
					}
				}
			}
			
			{
				cout << entry.name << ":MH" << MutHashesFound << endl << entry.seq << endl << entry.sep << endl;

				if (HashPos[0] < 93)
					cout << char(HashPos[0] + 33);
				else
					cout << char(126);
				for (int i = 1; i < HashPos.size(); i++) {
					if (HashPos[i] < 93)
						cout << char(HashPos[i] + 33);
					else
						cout << char(126);
				}

				for (int i = 0; i < buffer[b].seq.size() - HashSize; i++) {
					string hash = buffer[b].seq.substr(i, HashSize);
					string rev = Util::RevComp(hash);
					if (hash < rev)
						HashOut << hash << " 1" << endl;
					else
						HashOut << rev << " 1" << endl;
				}
				cout << endl;
			}
		}
	}
	HashOut.close();
	MutFile.close();
}
