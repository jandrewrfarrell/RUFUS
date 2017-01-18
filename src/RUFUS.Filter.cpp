
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
	bitset<64> HashBits;
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
	return HashBits.to_ulong();
}

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
string RevComp (string& Sequence)
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
	cout << "Call is PreBuiltMutHash Mutant.fq  firstpassfile hashsize MinQ HashCountThreshold window threads " << endl;
	double vm, rss, MAXvm, MAXrss;
	MAXvm = 0;
	MAXrss = 0;
	process_mem_usage(vm, rss, MAXvm, MAXrss);
   	cout << "VM: " << vm << "; RSS: " << rss << endl;	
	
	
	int BufferSize = 240;

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
		temp = Split(L1, '\t');
		if (temp.size() == 2)
		{
			unsigned long b = HashToLong(temp[0]);
			unsigned long revb = HashToLong(RevComp(temp[0]));
			//int c = atoi(temp[1].c_str()); 
			//cout << temp[0] << "\t" << temp[1] << endl << b << "\t" << LongToHash(b, 18) <<  endl<<endl;
			Mutations.insert(pair<unsigned long, int > (b, 0));
			Mutations.insert(pair<unsigned long, int > (revb, 0));
			//break;
		}
		else if (temp.size() == 4)
		{
			unsigned long b = HashToLong(temp[3]);
                	unsigned long revb = HashToLong(RevComp(temp[3]));
			//int c = atoi(temp[1].c_str()); 
                	//cout << temp[0] << "\t" << temp[1] << endl << b << "\t" << LongToHash(b, 18) <<  endl<<endl;
                	Mutations.insert(pair<unsigned long, int > (b, 0));
			Mutations.insert(pair<unsigned long, int > (revb, 0));
                	//break;
		}
		if (temp.size() == 1)
                {
		temp = Split(L1, ' ');
                unsigned long b = HashToLong(temp[0]);
                unsigned long revb = HashToLong(RevComp(temp[0]));
                //int c = atoi(temp[1].c_str()); 
                //cout << temp[0] << "\t" << temp[1] << endl << b << "\t" << LongToHash(b, 18) <<  endl<<endl;
                Mutations.insert(pair<unsigned long, int > (b, 0));
                Mutations.insert(pair<unsigned long, int > (revb, 0));
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
	string Buffer [2400];
	while (getline(MutFile, L1))
	{
		lines++;

                Buffer[0]=L1;
		getline(MutFile, Buffer[1]);
                getline(MutFile, Buffer[2]);
                getline(MutFile, Buffer[3]);
		
		 if (lines % 10000 > 1 && (lines % (10000+Threads) < Threads))
                 {
                        Et = clock();
                        float Dt = ((double) (Et - St)) * CLOCKS_PER_SEC;
                 	cout << "Read in " << lines << " lines: Found " << found << " Reads per sec = " << (float)lines/(float)Dt << " \r";
                 }
		int pos = 4; 
		while (getline(MutFile, Buffer[pos]) )
		{
			pos++;
			if (pos >=BufferSize)
				break;
			
		}
		#pragma omp parallel for shared(MutOutFile) num_threads(Threads)
		for (int BuffCount = 0; BuffCount < pos ; /*BuffCount < Buffer.size();*/ BuffCount+=4)
		{
			vector<int> positions; 
			int rejected = 0;
			int MutHashesFound = 0;
			bool good = true;
			int streak = 0; 
			int start = 0;  
			for (int i=0; i<Buffer[BuffCount+3].length()-HashSize+1; i++)
			{
				if ((int)Buffer[BuffCount+3].c_str()[i]-33<MinQ or (int)Buffer[BuffCount+1].c_str()[i] == 78)
				{
					streak = 0;
					start = i+1;
					rejected++; 
				}
				else
				{
					streak++; 
					if (streak >= HashSize-1)
						break;
					
				}	
			} 
			
			for (int i=start; i<Buffer[BuffCount+1].length()-HashSize; i++)
			{
				if (((int)Buffer[BuffCount+3].c_str()[i+HashSize-1] - 33) < MinQ or (int)Buffer[BuffCount+1].c_str()[i+HashSize-1] == 78) 
				{
					streak = 0;
					start = i+HashSize;
					for (int j=start; j<Buffer[BuffCount+3].length(); j++)
                        		{
                                		if ((int)Buffer[BuffCount+3].c_str()[j]-33 < MinQ or (int)Buffer[BuffCount+1].c_str()[j] == 78)
                                		{
                                		        streak = 0;
                                		        start = j+1;
                                			rejected++; 
						}
                                		else
                                		{
                                		        streak++;
                                        		if (streak >= HashSize)
                                        		        break;
                                		}
                        		}
					rejected++; 
					i = start;
			
				}
		 		if (streak  >= HashSize-1)
                        	{
				
					//unsigned long LongHash = HashToLong(Buffer[BuffCount+1].substr (i,HashSize));
					if (Mutations.count(HashToLong(Buffer[BuffCount+1].substr (i,HashSize))) > 0)
	                                {
						MutHashesFound++;
						positions.push_back(i);
					}
				}
			}
			//cout << MutHashesFound << endl;
			if (MutHashesFound >= HashCountThreshold and rejected < (Buffer[BuffCount+1].length()/2))
			{
			//cout << "yay keeping" << endl;
				int MaxCounter = -1;
                        	for (int i = 0; i<positions.size()-1; i++)
                        	{
                                	int counter = 1;
                                	for(int j = i+1; j<positions.size(); j++)
                                	{
                                	        if (positions[j] - positions[i] <= Window)
                                	        {counter++;}
                        	        }	
                	                if (counter > MaxCounter)
        	                        {MaxCounter = counter;}
	                        }
				if (MaxCounter >= HashCountThreshold )
				{
					#pragma omp critical(MutWrite)
					{MutOutFile << Buffer[BuffCount] << ":MH" << MutHashesFound << endl << Buffer[BuffCount+1] << endl << Buffer[BuffCount+2] << endl << Buffer[BuffCount+3] << endl;}
					found ++;
				}
				
			}
		}
	}
	cout << "\nDone\n";
                     
	MutFile.close();
	MutOutFile.close();
	
	cout << "\nreally dont\n";
}
	
