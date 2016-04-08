
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
	// 			1		2		3	     4	     5      6          7	    8       9
	cout << "Call is PreBuiltMutHash Mutant_mate1.fq Mutant_mate2.fq OutStub hashsize MinQ HashCountThreshold window threads " << endl;
	double vm, rss, MAXvm, MAXrss;
	MAXvm = 0;
	MAXrss = 0;
	process_mem_usage(vm, rss, MAXvm, MAXrss);
   	cout << "VM: " << vm << "; RSS: " << rss << endl;	
	
	
	int BufferSize = 1000;

 cout << "Paramaters are:\n  PreBuiltMutHash = " << argv[1] << "\n  Mutant.fq = " << argv[2] << "\n  out stub = " << argv[3] <<"\n  HashSize = " << argv[4] << "\n  MinQ = " << argv[5] << "\n  HashCountThreshold = " << argv[6]  << "\n  Window = " << argv[7] << "\n  Threads = " << argv[8] 
<< endl;
		// Read in file passed to the program on the command line
	
	string temp = argv[5];
        int HashSize = atoi(temp.c_str());	
	 temp = argv[6];
	int MinQ = atoi(temp.c_str());
	temp = argv[7];
	int HashCountThreshold = atoi(temp.c_str());
	temp = argv[8];
	int Window = atoi(temp.c_str());
	temp = argv[9];
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
	ifstream Mate1;
       	cout << "Mate1 is " << argv[2] << endl;
        Mate1.open (argv[2]);
        if ( Mate1.is_open())
        {      cout << "##File Opend\n";}
        else
        {
                cout << "Error, Mate1 could not be opened";
                return 0;
        }
	filename = argv[3];
        ifstream Mate2;
        cout << "Mate2 is " << argv[3] << endl;
        Mate1.open (argv[3]);
        if ( Mate1.is_open())
        {      cout << "##File Opend\n";}
        else
        {       
                cout << "Error, Mate1 could not be opened";
                return 0;
        }
	
	ofstream Mate1out;
	string FirstPassFile = argv[4]; 
	FirstPassFile += ".Mate1.filtered.fastq";
	Mate1out.open (FirstPassFile.c_str());
	if (Mate1out.is_open())
	{}
	else
	{	
		cout << "ERROR, Output file could not be opened -" << argv[3] << endl;
		return 0;
	}
	
	ofstream Mate2out;
        FirstPassFile = argv[4];
        FirstPassFile += ".Mate2.filtered.fastq";
        Mate2out.open (FirstPassFile.c_str());
        if (Mate2out.is_open())
        {}
        else
        {
                cout << "ERROR, Output file could not be opened -" << argv[3] << endl;
                return 0;
        }
	
	string line;
	unordered_map  <unsigned long, bool> Mutations;

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
		temp = Split(L1, ' ');
		unsigned long b = atol(temp[0].c_str());
		//cout << temp[0] << "\t" << temp[1] << endl << b << "\t" << LongToHash(b, 18) <<  endl<<endl;
		Mutations.insert(pair<unsigned long, bool> (b, 0));
		//break;
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
	while (getline(Mate1, L1))
	{
		lines++;
		getline(Mate1, L2);
                getline(Mate1, L3);
                getline(Mate1, L4);
		
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
		while (Buffer.size()<Threads*4 && getline(Mate1, L1) )
		{
			lines++;
			getline(Mate1, L2);
                	getline(Mate1, L3);
                	getline(Mate1, L4);
			
			Buffer.push_back(L1);
			Buffer.push_back(L2);
                	Buffer.push_back(L3);
                	Buffer.push_back(L4);
		//	 cout << "adding  \n   " << L1 << "\n   " << L2 << "\n   " << L3 << "\n   " << L4 << endl;
		}
		//cout << " Buffer has " << Buffer.size() << " lines\n";
		#pragma omp parallel for shared(Mate1out) num_threads(Threads)
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
			for (int i=0; i<B2.length()-HashSize; i++)
			{
				string hash = B2.substr (i,HashSize);
				string Qhash = B4.substr (i,HashSize);	
				bool good = true;
				for (int j=0; j<HashSize; j++)
                        	{
	                               	int B = hash.c_str()[j];
        	                        if ( B == 78)
                	                {
                        	                good = false;
                                	        break;
                                	}
					B = Qhash.c_str()[j];
					if (((int)B - 33) < MinQ)
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
	                                	MutHashesFound++;
						positions.push_back(i);
					}
					else if (Mutations.count(HashToLong(RevComp(hash))) > 0)
	                        	{
						MutHashesFound++;
                                                positions.push_back(i);
					}
				}
	                        else
	                        {
					rejected++;
	                        }
			}
				
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
				//	cout << "   kept - " << MaxCounter << endl;
					#pragma omp critical(MutWrite)
					{Mate1out << B1 << ":MH" << MutHashesFound << endl << B2 << endl << B3 << endl << B4 << endl;}
					found ++;
				}
				//else
				//{cout << "    rejected " << endl;}
				
			}
		}
	}
	cout << "\nDone\n";
                     
	Mate1.close();
	Mate1out.close();
	
	cout << "\nreally dont\n";
}
	
