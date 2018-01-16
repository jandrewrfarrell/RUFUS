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
//	cout << "Call is PreBuiltMutHash Overlap.fq HashSize" << endl;
	int HashSize; 
	double vm, rss, MAXvm, MAXrss;
	MAXvm = 0;
	MAXrss = 0;
	Util::process_mem_usage(vm, rss, MAXvm, MAXrss);
	
	
	int BufferSize = 1000;
	
	ifstream MutHashFile;		
	MutHashFile.open (argv[1]);
	string filename = argv[2];
	ifstream MutFile;
 	if (filename == "stdin")
        {
                MutFile.open ("/dev/stdin");
        }
        else
        {
                MutFile.open (argv[2]);
        }
	string line;
	unordered_map  <unsigned long, unsigned int> Mutations;

	ofstream HashOut; 
	HashOut.open (argv[3]); 
	 

	int lines = 0;

	string L1;
	string L2;
	string L3;
	string L4;
 
	unsigned long LongHash;

	bool notdone = true;
	
	while (getline(MutHashFile, L1))
	{
		vector<string>  temp;
		temp = Util::Split(L1, '\t');
		//cout << "size = " << temp.size() << endl; 
		if (temp.size() == 2)
		{
		  unsigned long b = Util::HashToLong(temp[0].c_str());
			Mutations.insert(pair<unsigned long, unsigned int> (b, atoi(temp[1].c_str())));
			HashSize = temp[0].length();
		//	cout << "2 " << b << " = " << temp[0].c_str() << " - " << temp[1].c_str();	
		}
		else if(temp.size() ==4)
		{
			//cout << "here" << endl; 	
		  unsigned long b = Util::HashToLong(temp[3].c_str());
                        Mutations.insert(pair<unsigned long, unsigned int> (b, atoi(temp[2].c_str())));
                        HashSize = temp[3].length();
		//	cout << "4 " << b << " = " << temp[3].c_str() << " - " << temp[2].c_str(); 
		}
		else if (temp.size() == 1)
                {
		  temp = Util::Split(L1, ' ');
		  unsigned long b = Util::HashToLong(temp[0].c_str());
                        Mutations.insert(pair<unsigned long, unsigned int> (b, atoi(temp[1].c_str())));
                        HashSize = temp[0].length();
                 //       cout << "1 " << b << " = " << temp[0].c_str() << " - " << temp[1].c_str() << endl;
                }
	}
	MutHashFile.close();
	//cout << "Hash size = " << HashSize << endl;
	
	 
	
	int who = RUSAGE_SELF;
        struct rusage usage;
        int b = getrusage(RUSAGE_SELF, &usage);
//        cout << "I am using " << usage.ru_maxrss << endl;

//	cout << "VM: " << vm << "; RSS: " << rss <<"; maxVM: " << MAXvm << "; maxRSS: " << MAXrss << endl;
//	return 0;
	
//	cout << "Starting Search " << endl;
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
		
		vector<int> positions; 
		vector<int> HashPos;

		
		//cout << L2 << endl;
               	for (int i = 0; i < L2.size(); i++)
               	{HashPos.push_back(0);}
                //cout << "yay" << endl;

		int rejected = 0;
		int MutHashesFound = 0;
		for (int i=0; i<L2.length()-HashSize; i++)
		{
			string hash = L2.substr (i,HashSize);
			string Qhash = L4.substr (i,HashSize);	
                	

			bool good = true;
                        for (int j=0; j<HashSize; j++)
                        {
                        	int B = hash.c_str()[j];
                              	if ( B == 78)
                                {
                                        good = false;
                                        break;
                                }
                        }
                        if (good)
                        {
				
			  unsigned long LongHash = Util::HashToLong(hash);
				if (Mutations.count(LongHash) > 0)
	                	       {
					positions.push_back(i);
					for (int j = i; j<i+HashSize; j++)
					{ 
						HashPos[j]++;
						//if (HashPos[j] < Mutations[LongHash])
						//	HashPos[j] = Mutations[LongHash];	
					}
				}
				else if (Mutations.count(Util::HashToLong(Util::RevComp(hash))) > 0)
	               		{
                        	        positions.push_back(i);
					for (int j = i; j<i+HashSize; j++)
					{
						HashPos[j]++; 
						//if( HashPos[j] < Mutations[HashToLong(RevComp(hash))])
						//	HashPos[j] = Mutations[HashToLong(RevComp(hash))];
					}
				}
					
			}
		}
				
		//cout << L1 << ":MH" << MutHashesFound << endl << L2 << endl << L3 << endl << L4 << endl << L3 << endl;
		//cout << HashPos[0];
		//for (int i =1; i<HashPos.size(); i++) 
               // 	cout << " "<< HashPos[i];
                //cout << endl;
		//
	 	cout << L1 << ":MH" << MutHashesFound << endl << L2 << endl << L3 << endl ;
                if (HashPos[0] < 93)
			cout << char(HashPos[0]+33);
		else	
			cout << char(126);
                for (int i =1; i<HashPos.size(); i++)
              	{
			 if (HashPos[i] < 93)
                        	cout << char(HashPos[i]+33);
                	else
                        	cout << char(126); 
		}
		
		for (int i = 0; i < L2.size() - HashSize; i++)
		{
			string hash = L2.substr(i, HashSize); 
			string rev = Util::RevComp(hash); 
			if (hash < rev)
				HashOut << hash << " 1"<< endl;
			else
				HashOut << rev << " 1" << endl;	
		}
		
		cout << endl;
	}
	HashOut.close(); 
                  
	MutFile.close();
	
}
	
