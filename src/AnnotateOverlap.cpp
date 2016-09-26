
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

using namespace std;
int totalAdded =0;
int totalDeleted =0;

/////////////////////////
bool fncomp (char lhs, char rhs) {return lhs<rhs;}

struct classcomp {
  bool operator() (const char& lhs, const char& rhs) const
  {return lhs<rhs;}
};
///////////////////////////

const vector<string> Split(const string& line, const char delim) {
    vector<string> tokens;
    stringstream lineStream(line);
    string token;
    while ( getline(lineStream, token, delim) )
        tokens.push_back(token);
    return tokens;
}

unsigned long HashToLong (string hash)
{
	//cout << "booya" << hash << endl;
	bitset<64> HashBits;
	//#pragma omp parellel for
	for(int i=0; i<hash.length();i++)
	{
		if (hash.c_str()[i] == 'A')
		{
			HashBits[i*2] = 0;
			HashBits[i*2+1] = 0;
		}
		else  if (hash.c_str()[i] == 'C')
                {
                        HashBits[i*2] = 0;
                        HashBits[i*2+1] = 1;
                }
		else  if (hash.c_str()[i] == 'G')
                {
                        HashBits[i*2] = 1;
                        HashBits[i*2+1] = 0;
                }
		else  if (hash.c_str()[i] == 'T')
                {
                        HashBits[i*2] = 1;
                        HashBits[i*2+1] = 1;
                }
		else
		{
			cout << "ERROR, invalid character - " << hash.c_str()[i] << endl;
		}
	}
	//cout <<  HashBits.to_ulong() << "-" << endl;
	return HashBits.to_ulong();
}
/*
unsigned long HashToLong (string hash)
{
        //cout << "booya" << hash << endl;
        bitset<64> HashBits;
        int bitcout = 0;
        for(int i=0; i<hash.length();i++)
        {
                if (hash.c_str()[i] == 'A')
                {
                        HashBits[bitcout] = 0;
                        bitcout++;
                        HashBits[bitcout] = 0;
                        bitcout++;
                }
                else  if (hash.c_str()[i] == 'C')
                {
                        HashBits[bitcout] = 0;
                        bitcout++;
                        HashBits[bitcout] = 1;
                        bitcout++;
                }
                else  if (hash.c_str()[i] == 'G')
                {
                        HashBits[bitcout] = 1;
                        bitcout++;
                        HashBits[bitcout] = 0;
                        bitcout++;
                }
                else  if (hash.c_str()[i] == 'T')
                {
                        HashBits[bitcout] = 1;
                        bitcout++;
                        HashBits[bitcout] = 1;
                        bitcout++;
                }
                else
                {
                        cout << "ERROR, invalid character - " << hash.c_str()[i] << endl;
                }
        }
        //cout <<  HashBits.to_ulong() << "-" << endl;
        unsigned long booya;
	booya = HashBits.to_ulong();
	return booya;
}*/
string LongToHash (unsigned long LongHash, int HashSize)
{
        string value = "";
        //cout << LongHash << endl;
        bitset<64> test (LongHash);
        for (int i = 1; i < HashSize*2; i+=2)
        {
                //cout << "i-1=" << test[i-1] << " i=" << test[i] << " - ";
                if (test[i-1] == 0)
                {
                        if (test[i] == 0)
                        {
                                //cout << "A" << endl;
                                value = value + "A";
                        }
                        else
                        {
                                //cout << "C" << endl;
                                 value = value + "C";
                        }
                }
                else
                {
                        if (test[i] == 0)
                        {
                                //cout << "G" << endl;
                                 value = value + "G";
                        }
                        else
                        {
                                //cout << "T" << endl;
                                 value = value + "T";
                        }

                }


        }
//      for (int i = 0; i < HashSize*2; i++)
//      {
//              cout << i << " - " << test[i] << endl;
//      }
        return value;
//      char test [33];
//      itoa (4,test,2);
//      cout << test << endl;
//      cint >> test;
}
string RevComp (string Sequence)
{
        string NewString = "";
        //cout << "Start - " << Sequence << "\n";
        for(int i = Sequence.length()-1; i>=0; i+= -1)
        {
                char C = Sequence.c_str()[i];
        //      cout << C << endl;
                if (C == 'A')
                        NewString += 'T';
                else if (C == 'C')
                         NewString += 'G';
                else if (C == 'G')
                         NewString += 'C';
                else if (C == 'T')
                         NewString += 'A';
                else if (C == 'N')
                         NewString += 'N';
		else
                        cout << "ERROR IN RevComp - " << C << endl;


        }
        //cout << "end\n";
        return NewString;
}
void ProcessStack( stack <string>& Pstack, unordered_map  <unsigned long,  unsigned short >& Phash,unordered_map  <unsigned long,  unsigned short >& Mhash,  int MinMutDepth, int MaxParDepth, int side, unordered_map  <unsigned long, bool >& Mutations, unordered_map  <unsigned long, bool >& CopyNumber, int Par1x, int Mut1x) 
{
	//if (side == 0) cout << "Processing Parent " << endl;else cout << "Processing Mutant" << endl;
//	int size = Pstack.size();
//	 int count = 0;
	while (!Pstack.empty())
 //	#pragma omp parellel for
//	for ( int bam = 0; bam < size; bam++) 
	{
//		count++;
        	string L1 = Pstack.top();  
                Pstack.pop();
                vector<string> stuff;
                stuff = Split(L1, '\t');
                int count = atoi(stuff[1].c_str());
                unsigned short C;
                if (count > 65534){C = 65535;}else{C = count;}\

		//cout << "line = " << L1 << endl << "       " << stuff[0] << "\t" << C << endl;
                unsigned long  ParentLongHash = HashToLong(stuff[0]);        
		if (side == 0)
		{	
			if (Mhash.count(ParentLongHash) > 0)
                        {
                               	unsigned short M = Mhash[ParentLongHash];     			
		//		cout << "line = " << L1 << endl << "    P =" << stuff[0] << "\t" << C << endl <<"     M =" << stuff[0] << "\t" << M << endl;
				if ( M >= MinMutDepth && C <= MaxParDepth)
                                {
		//			 cout << "     Find mutation\n";
					Mutations[ParentLongHash] = 1;
                                        Mhash.erase(ParentLongHash);totalDeleted++;
				}
                                else
                                {
					float ParFold = (float) C / (float) Par1x;
                               	 	float MutFold = (float) M / (float) Mut1x;

                                	if ((fabs(ParFold - MutFold) > 1.9) && (C >= ( Par1x / 2 )) && (M >= ( Mut1x / 2 )) && (C <= ( Par1x * 20 )) && (M <= ( Mut1x * 20 )))
                                	{
		//				      cout << "     Find copyNum\n";
						CopyNumber[ParentLongHash] = 1;
                                        	Mhash.erase(ParentLongHash);totalDeleted++;
                                	}
					else 
					{Mhash.erase(ParentLongHash);totalDeleted++;}
				}
                        }
                        else
                        {
                                Phash.insert (pair<unsigned long,  unsigned short > (ParentLongHash,C));totalAdded++;
                        }
		}
		else if (side == 1)
		 {
                        if (Phash.count(ParentLongHash) > 0)
                        {
                               	unsigned short P = Phash[ParentLongHash];
		//		cout << "line = " << L1 << endl << "    M =" << stuff[0] << "\t" << C << endl <<"     P =" << stuff[0] << "\t" << P << endl;
				if ((int)C >= MinMutDepth && (int)P <= MaxParDepth)
                                {
		//			cout << "     Find mutation\n";
					Mutations[ParentLongHash] = 1;
                                        Phash.erase(ParentLongHash);totalDeleted++;
                                }
				else
                                {       
                                        float ParFold = (float) P / (float) Par1x;
                                        float MutFold = (float) C / (float) Mut1x;

                                        if ((fabs(ParFold - MutFold) > 1.9) && ((int)P >= ( Par1x / 2 )) && ((int)C >= ( Mut1x / 2 )) && ((int)P <= ( Par1x * 20 )) && ((int)C <= ( Mut1x * 20 )))
                                        {	
		//				cout << "     Find copyNum\n";
                                                CopyNumber[ParentLongHash] = 1;
                                                Phash.erase(ParentLongHash);totalDeleted++;
                                        }
                                	else{Phash.erase(ParentLongHash);totalDeleted++;}
                        	}
			}
                        else
                        {
                                Mhash.insert (pair<unsigned long,  unsigned short > (ParentLongHash,C));totalAdded++;
                        }
                }
	}
	//cout << size << " - " << count << endl;
return;
}
void process_mem_usage(double& vm_usage, double& resident_set, double& MAXvm, double& MAXrss)
{
   using std::ios_base;
   using std::ifstream;
   using std::string;

   vm_usage     = 0.0;
   resident_set = 0.0;

   // 'file' stat seems to give the most reliable results
   //
   ifstream stat_stream("/proc/self/stat",ios_base::in);

   // dummy vars for leading entries in stat that we don't care about
   //
   string pid, comm, state, ppid, pgrp, session, tty_nr;
   string tpgid, flags, minflt, cminflt, majflt, cmajflt;
   string utime, stime, cutime, cstime, priority, nice;
   string O, itrealvalue, starttime;

   // the two fields we want
   //
   unsigned long vsize;
   long rss;

   stat_stream >> pid >> comm >> state >> ppid >> pgrp >> session >> tty_nr
               >> tpgid >> flags >> minflt >> cminflt >> majflt >> cmajflt
               >> utime >> stime >> cutime >> cstime >> priority >> nice
               >> O >> itrealvalue >> starttime >> vsize >> rss; // don't care about the rest

   stat_stream.close();

   long page_size_kb = sysconf(_SC_PAGE_SIZE) / 1024; // in case x86-64 is configured to use 2MB pages
   vm_usage     = vsize / 1024.0;
   resident_set = rss * page_size_kb;
	if (vm_usage > MAXvm){MAXvm = vm_usage;}
	if (resident_set > MAXrss){MAXrss = resident_set;}
}
	
int main (int argc, char *argv[])
{
//	cout << "Call is PreBuiltMutHash Overlap.fq HashSize" << endl;
	int HashSize; 
	double vm, rss, MAXvm, MAXrss;
	MAXvm = 0;
	MAXrss = 0;
	process_mem_usage(vm, rss, MAXvm, MAXrss);
	
	
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
		temp = Split(L1, '\t');
		//cout << "size = " << temp.size() << endl; 
		if (temp.size() == 2)
		{
			unsigned long b = HashToLong(temp[0].c_str());
			Mutations.insert(pair<unsigned long, unsigned int> (b, atoi(temp[1].c_str())));
			HashSize = temp[0].length();
		//	cout << "2 " << b << " = " << temp[0].c_str() << " - " << temp[1].c_str();	
		}
		else if(temp.size() ==4)
		{
			//cout << "here" << endl; 	
			unsigned long b = HashToLong(temp[3].c_str());
                        Mutations.insert(pair<unsigned long, unsigned int> (b, atoi(temp[2].c_str())));
                        HashSize = temp[3].length();
		//	cout << "4 " << b << " = " << temp[3].c_str() << " - " << temp[2].c_str(); 
		}
		else if (temp.size() == 1)
                {
			temp = Split(L1, ' ');
                        unsigned long b = HashToLong(temp[0].c_str());
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
				
		       		unsigned long LongHash = HashToLong(hash);
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
				else if (Mutations.count(HashToLong(RevComp(hash))) > 0)
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
			string rev = RevComp(hash); 
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
	
