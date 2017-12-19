//this vertion pairs with ModeDist%
#include <unistd.h>
#include <ios>
#include <omp.h>

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
#include <sys/resource.h>
#include <stack>
#include <sys/time.h>

using namespace std;
bool WriteTrace = true;
ifstream ParentHash [100];
////Call is BitHashCompare Parent Mutant firstpassfile hashsize

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
string LongToHash (unsigned long LongHash, int HashSize)
{
    string value = "";
    bitset<64> test (LongHash);
    for (int i = 1; i < HashSize*2; i+=2)
    {
        if (test[i-1] == 0)
        {
            if (test[i] == 0){value = value + "A";}
            else{value = value + "C";}
        }
        else
        {
            if (test[i] == 0){value = value + "G";}
            else{value = value + "T";}
        }
    }
    return value;
}
string RevComp (string Sequence)
{
    string NewString = "";
    for(int i = Sequence.length()-1; i>=0; i+= -1)
    {
        char C = Sequence.c_str()[i];
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


    	ifstream MutHashFile;
    	string file = argv[1];
	MutHashFile.open (file);
    	if ( MutHashFile.is_open())
    	{}
    	else
    	{
    	    cout << "Error, MutHashFile could not be opened";
    	    return 0;
    	}
	
	cout << file << " is open " << endl; 	
	
	string MutLine = ""; 
	getline(MutHashFile, MutLine);
	getline(MutHashFile, MutLine);
	getline(MutHashFile, MutLine);
	       getline(MutHashFile, MutLine);
	               getline(MutHashFile, MutLine);
		               getline(MutHashFile, MutLine);
			       int counter =0; 
    while (getline(MutHashFile, MutLine) && counter < 13)
    {
	    counter ++; 
	    vector <string> temp; 
	    temp=Split(MutLine, '\t'); 
	    double num = 0.0;
	    for (int i = 3; i<temp.size(); i++){
		    num+=atof(temp[i].c_str());}	
	double value = num/(num+atof(temp[1].c_str()));
	//cout << temp[0] << "\t" << value << endl;  
	cout << value << endl;
	}
}

