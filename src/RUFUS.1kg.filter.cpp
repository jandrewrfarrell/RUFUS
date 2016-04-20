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
#include <unordered_map>
#include <bitset>
#include <time.h>
#include <sys/resource.h>
#include <stack>
#include <sys/time.h>
#include <fcntl.h>
#include <sys/types.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <errno.h>



using namespace std;
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
            //cout << "ERROR, invalid character - " << hash.c_str()[i] << endl;
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
		//      1	     2	           3                4        5            6          7        8
     cout << "Call is Parent_hash Mutant_Hash FilterLiestOutFile hashsize MaxParDepth MinMutDepth MaxDepth Threads" << endl;

    double vm, rss, MAXvm, MAXrss;
    MAXvm = 0;
    MAXrss = 0;
    process_mem_usage(vm, rss, MAXvm, MAXrss);

    string temp = argv[4];
    int HashSize = atoi(temp.c_str());
     temp = argv[5];
    int MaxParDepth = atoi(temp.c_str());

     temp = argv[6];
    int MinMutDepth = atoi(temp.c_str());
    temp = argv[7];
    int TooHigh = atoi(temp.c_str());
    temp = argv[8];
    int Threads = atoi(temp.c_str());

    int BufferSize = 10000 ;

 cout << "Paramaters are:\n  Parent_Hash = " << argv[1] <<
        "\n  Mutant_Hash = " << argv[2] <<
        "\n  FilterLiestOutFile = " << argv[3] <<
        "\n  HashSize = " << argv[4] <<
        "\n  Max Parent Depth = " << argv[5] <<
        "\n  Min Mut Depth = " << argv[6] <<
        "\n  Maximum count = " << argv[7] <<
        "\n  Threads = " << argv[8] << endl;

        // Read in file passed to the program on the command line


    /*ifstream ParentFile;
    ParentFile.open (argv[1]);
    if ( ParentFile.is_open())
    { cout << "Parent File open - " << argv[1] << endl;}    //cout << "##File Opend\n";
    else
    {
        cout << "Error, ParentFile could not be opened";
        return 0;
    }*/


    ifstream MutHashFile;
    MutHashFile.open (argv[2]);
    if ( MutHashFile.is_open())
    {}
    else
    {
        cout << "Error, MutHashFile could not be opened";
        return 0;
    }

    ofstream HashTable;
    string HashTablePath = argv[3];
    HashTable.open (HashTablePath.c_str());
    if (HashTable.is_open())
    {}
    else
    {
        cout << "Error HashTableFilter file couldnt be opened - " << HashTablePath << endl;
        return 0;
    }



    string line;
    //unordered_map  <unsigned long, vector<unsigned char> > Phash;
    unordered_map  <unsigned long, int > good;
    //std::unordered_map<std::bitmap, bool> ParentHashes;

    unsigned long int lines = 0;

    string L1;
    string L2;
    string L3;
    string L4;
    string ML1, PL1;
 

    clock_t St,Et;
    struct timeval start, end;
    //St = clock();
    gettimeofday(&start, NULL);
    bool notdone = true;
    int r= 0;
    cout << "Reading in Mutant Hash Table" << endl;
    while (notdone == true)
    {
        lines++;
        //if (lines == 100)
        //{break;}
        vector<string> Mstack;
        while (Mstack.size() < BufferSize && getline(MutHashFile, ML1) )
        {
            Mstack.push_back(ML1);
        }


        for ( int bam = 0; bam < Mstack.size(); bam++)
        {
            string B1 =  Mstack[bam];
            vector<string> stuff;
            stuff = Split(B1, '\t');
            int Bcount = atoi(stuff[1].c_str());
            //unsigned char Bcount;
            //if (Bcount > 250){C = 251;}else{C = Bcount;}
            //if (Bcount > 10000){C = 0;}
            unsigned long  MutantLongHash = HashToLong(stuff[0]);
	    //if (Bcount > MinMutDepth and Bcount < TooHigh )  
	    {good.insert(pair<unsigned long, int>(MutantLongHash, Bcount));}
            
        }




        if (Mstack.size() == 0 )
        {notdone = false;}
        if ((lines*BufferSize) % 1000 == 0 )
        {
            gettimeofday(&end, NULL);
            float Dt = end.tv_sec - start.tv_sec;
            process_mem_usage(vm, rss, MAXvm, MAXrss);
                cout << "Read in " << BufferSize * lines << " lines " <<
                        "; VM: " << vm <<
                        "; RSS: " << rss <<
                        " good = " << good.size() <<
                        ", TotalTime= " << Dt <<
                        " , second per line = " << Dt/(lines*BufferSize) <<"\r";
        }
    //    if (lines *BufferSize > 2000000){break;}
    }
    cout << "Done reading in mut table\nStarting filter\n";
    lines =0;


    #pragma omp parallel for num_threads(Threads) 
    for (int p = 0; p < Threads; p++)
    {
    	struct timeval start, end;
	gettimeofday(&start, NULL);
    	
	long int fd;
	struct stat sb;
    	char* fName = argv[1]; 
    	char *data;
    	char *lineMap;
    	// map the file
    	fd = open(fName, O_RDONLY);
    	fstat(fd, &sb);
    	long int pageSize;
    	pageSize = sysconf(_SC_PAGE_SIZE);

    	long int off = 0;
   
    	// get lines
    	long int count = 0;
    	char *fileptr = NULL;
    	string foo = "";
    	unordered_map <char, bool >  chars;
    	chars.insert(pair<char, bool>('A', true));
    	chars.insert(pair<char, bool>('C', true));
    	chars.insert(pair<char, bool>('G', true));
    	chars.insert(pair<char, bool>('T', true));
	
 	long int ChunkSize = sb.st_size/pageSize/Threads; 
	off =  ChunkSize * p * pageSize; 
	cout << "total File size = " << sb.st_size << endl; 
	cout << "Chunc Size = " << ChunkSize<< endl;
	cout << "thread " << p +1 << endl;
	cout << "I go from " << off << "    to    " <<  sb.st_size/pageSize/Threads  * (p+1) << endl; 
	long pages = 0;
    	while( off  <=   ChunkSize  * pageSize * (p+1)) 
    	{
		pages++;
		//	cout << "off = "<< off << endl;
		fileptr = (char*)mmap64(NULL, pageSize, PROT_READ, MAP_PRIVATE | MAP_POPULATE, fd, off);
		data = fileptr;
		for (int i = 0; i < pageSize; i++)
        	{
			if ( chars.count(data[i])>0)
			{
				foo+=data[i];
			}
			else
			{
				if (foo.length() == HashSize)
				{
					unsigned long ParHash = HashToLong(foo);
                        		int Count;
                        		Count =	 good.count(ParHash);
			//		cout << foo.length() << " - " << foo << endl;
                        		if (Count > 0)
                        		{
						gettimeofday(&end, NULL);
						float Dt = end.tv_sec - start.tv_sec;
						#pragma omp critical (boom)
                        			{cout << "   FoundHash -" << foo << "- in thread " << p << " " << foo << " - " << ParHash << " - % done =" <<  ((float) pages/(float) ChunkSize ) * 100.0  << ", Pages/sec = " << ((float) pages * (float) Threads)/Dt <<endl;}
                        		        //#pragma omp cricial(a)
						{good[ParHash] = -1;}
                        		}
				}
				foo = "";
			}
        	}
	 	munmap(fileptr, pageSize);
  		off += pageSize;
    	}
	
    }

    cout << "\nDone reading Parent File\n";

    //ParentFile.close();
    int who = RUSAGE_SELF;
    struct rusage usage;
    int b = getrusage(RUSAGE_SELF, &usage);
    cout << "I am using " << usage.ru_maxrss << endl;


    unordered_map  <unsigned long, int >::iterator pos;
    for (pos = good.begin(); pos != good.end(); ++pos)
    {
	if (pos->second > MinMutDepth and pos->second < TooHigh)
        {HashTable << pos->first << "\t" << "0" << "\t" << pos->second << "\t" << LongToHash(pos->first, HashSize) << endl;}
    }

    HashTable.close();

cout << "\nreally dont\n";
}

