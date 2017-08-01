
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
#include <map>
#include <ctime>
#include <cmath>
#include <fcntl.h>
#include <sys/types.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <errno.h>

#define NUMINTS  (1000)
#define FILESIZE (NUMINTS * sizeof(int))

using namespace std;


int HashSize = 25; 
/////////////////////////
const vector<string> Split(const string& line, const char delim) {
    vector<string> tokens;
    stringstream lineStream(line);
    string token;
    while ( getline(lineStream, token, delim) )
        tokens.push_back(token);
    return tokens;
}

int checkPage( char *data, string hash, long int pageSize, string line)
{
        //cout << "checking page" << "with size = " << pageSize<< endl;
        bool firstNew = false;
        for (int i = 0; i <  pageSize; i++)
        {
        //      cout << i << endl;
        //      cout << data[i] << endl;
                if (data[i] == '\n')
                {
        //              cout << "n found";
                        if (firstNew != true)
                                firstNew = true;
                        else
                        {
                                vector<string> stuff;
                                stuff = Split(line, ' ');
                                string PageHash = stuff[0];
        //                      cout << "PageHash " << endl;
                                if (hash == PageHash)
                                {
          //                              cout << "found a hash " << hash  << " - " << PageHash;
                                        return atoi(stuff[1].c_str());
                                }
                                line = "";
                        }
                }
                else if (firstNew == true)
                {
        //              cout << "first Newline found" << endl;
                        line += data[i];
                }

        }
        return 0;
}

void ProcessPage(  char *data, string& PageFirstHash, string& PageLastHash, long int pageSize)
{
        string line = "";
	bool firstNew = false;\
        for (int i = 0; i < pageSize; i++)
        {
                if (data[i] == '\n')
                {
		
                        if (firstNew == true)
                                break;
                        else
                                firstNew = true;
                }
                else if (firstNew == true)
                        line += data[i];
        }

        vector<string> stuff;
        stuff = Split(line, ' ');
        PageFirstHash = stuff[0];
        firstNew = false;
        line = "";
        for (int i = pageSize-1; i > 0; i+=-1)
        {
                if ( data[i] == '\n')
                {
                        if (firstNew == true)
                                break;
                        else
                                firstNew = true;
                }
                else if(firstNew == true)
                        line = data[i] + line;
        }
        stuff = Split(line, ' ');
        PageLastHash = stuff[0];

}



int search(long int& fd, string hash, char* fileptr)
{
        //cout << "searching for " << hash << endl;
        char *data;
        struct stat sb;
        fstat(fd, &sb);

        long int pageSize;
        pageSize = sysconf(_SC_PAGE_SIZE);
        long int NumPages = sb.st_size/pageSize;
        //cout << "Number of pages = " << NumPages << endl;
 //       char *fileptr = NULL;

        long int off = 0;
        long int firstPos;
        long int lastPos;
        firstPos = 0;
        fileptr = (char*)mmap64(NULL, pageSize, PROT_READ, MAP_PRIVATE | MAP_POPULATE, fd, firstPos*pageSize);
        data =  fileptr;
        // above should get me first page
        string FirstPageFirstHash;
        string FirstPageLastHash;
        ProcessPage(data, FirstPageFirstHash, FirstPageLastHash, pageSize);
        if (munmap(fileptr, pageSize) == -1) {
        	perror("Error un-mmapping the file");
        }
	//cout << FirstPageFirstHash << endl << FirstPageLastHash << endl;
        //quck check to see if on first page
        if (hash >= FirstPageFirstHash and hash <= FirstPageLastHash)
        {
              //cout << "found on first page" << endl;
         	fileptr = (char*)mmap64(NULL, pageSize, PROT_READ, MAP_PRIVATE | MAP_POPULATE, fd, firstPos*pageSize);
                int val =  checkPage(data, hash, pageSize, "");
		if (munmap(fileptr, pageSize) == -1) {
                	perror("Error un-mmapping the file");
        	}
		return val; 
        }
	if (hash < FirstPageFirstHash)
	{
		//cout << "HASH NOT IN FILE " << hash << endl;
		return 0; 
	}
        fileptr = (char*)mmap64(NULL, pageSize, PROT_READ, MAP_PRIVATE | MAP_POPULATE, fd, pageSize*(NumPages-1));
        data = fileptr;
        lastPos = NumPages-1;
        //the above should get me the last two pages, we take two to ensure the last pages isnt just one character or something like that

        string LastPageFirstHash;
        string LastPageLastHash;
        ProcessPage(data, LastPageFirstHash, LastPageLastHash, pageSize);
        if (munmap(fileptr, pageSize) == -1) {
                perror("Error un-mmapping the file");
        }
	//cout << LastPageFirstHash << endl << LastPageLastHash << endl;
        //quck check to see if on last page
        if (hash >= LastPageFirstHash and hash <= LastPageLastHash)
        {
              	fileptr = (char*)mmap64(NULL, pageSize, PROT_READ, MAP_PRIVATE | MAP_POPULATE, fd, pageSize*(NumPages-1));
		//cout << "found on last page" << endl;
                int val =  checkPage(data, hash, pageSize, "");
        	if (munmap(fileptr, pageSize) == -1) {
                	perror("Error un-mmapping the file");
        	}
		return val; 
	}
	if (hash > LastPageLastHash)
	{
	        //cout << "HASH NOT IN FILE " << hash << endl;
                return 0;
        }
        //start the search
        int counter = 0;
        while (true)
        {
              //cout << "ON LOOP " << counter << endl << endl;
                counter++;
                long int currentPage = lastPos - ((lastPos-firstPos)/2);
              //cout << "checking page " << currentPage << " last = " << lastPos << " and first = " << firstPos << endl;;
                if (currentPage == lastPos or currentPage == firstPos or lastPos - firstPos <= 3)
                {
                        string extra = "";
                 //     cout << "\nenvoked this" << endl;
                        fileptr = (char*)mmap64(NULL, pageSize*5, PROT_READ, MAP_PRIVATE | MAP_POPULATE, fd, pageSize * (firstPos-1));
                        data =  fileptr;
                     // cout << "made it here" << endl;
                        int val = checkPage(data, hash, pageSize*5, extra);
			if (munmap(fileptr, pageSize*5) == -1) {
                		perror("Error un-mmapping the file");
        		}
			return val; 
                }
		//cout << " fileptr = (char*)mmap64(NULL, " << pageSize*2 <<", PROT_READ, MAP_PRIVATE | MAP_POPULATE, fd, " << pageSize<<" * " << currentPage <<");"<< endl;
                fileptr = (char*)mmap64(NULL, pageSize, PROT_READ, MAP_PRIVATE | MAP_POPULATE, fd, pageSize*currentPage);
		data =  fileptr;
                string CurrentPageFirstHash;
                string CurrentPageLastHash;
		ProcessPage(data, CurrentPageFirstHash, CurrentPageLastHash, pageSize);
		if (munmap(fileptr, pageSize) == -1) {
                	perror("Error un-mmapping the file");
        	}
		//cout << " with " << CurrentPageFirstHash << " and " << CurrentPageLastHash << endl;
                if (hash >= CurrentPageFirstHash and hash <= CurrentPageLastHash)
                {
			fileptr = (char*)mmap64(NULL, pageSize, PROT_READ, MAP_PRIVATE | MAP_POPULATE, fd, pageSize*currentPage);
                        int val = checkPage(data, hash, pageSize, "");
			if (munmap(fileptr, pageSize) == -1) {
                        perror("Error un-mmapping the file");
                	}
			return val;
                }
                else
                {
                        if (hash < CurrentPageFirstHash)
                        {
                  //            cout << "hash " << hash << " is greater than " << CurrentPageFirstHash << " looking above" << endl;
                                lastPos = currentPage;
                                LastPageFirstHash = CurrentPageFirstHash;
                                LastPageLastHash = CurrentPageLastHash;
                        }
                        else if (hash > CurrentPageLastHash)
                        {
                    //          cout << "hash \n" << hash << " is less than \n" << CurrentPageLastHash << " looking below" << endl;
                                firstPos = currentPage;
                                FirstPageFirstHash = CurrentPageFirstHash;
                                FirstPageLastHash = CurrentPageLastHash;
                        }
                }

        }
	close(fd); 

}

bool fncomp (char lhs, char rhs) {return lhs<rhs;}

struct classcomp {
  bool operator() (const char& lhs, const char& rhs) const
  {return lhs<rhs;}
};
///////////////////////////
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
//	ifstream testthis [100]; 
//	testthis[0].open("./test.txt"); 
//	string boom2; 
//	while (getline(testthis[0], boom2))
//	{
//		cout << boom2 << endl;
//	}
//	return 0;  
	//************************************************
	//my arg parser

	string helptext;	
	helptext = \
"\
RUFUS.interpret: converts RUFUS aligned contigs into a VCF \n\
By Andrew Farrell\n\
   The Marth Lab\n\
\n\ 
options:\
  -h [ --help ]	  Print help message\n\
  -hf arg		Path to HashFile from RUFUS.build\n\
  -o  arg		Output stub\n\
  -c  arg		Path to sorted.tab file for the parent sample\n\
  -hs arg		Hash size (default = 25)\n\
";
	
	string RefFile = ""; 
	string HashListFile = "" ; 	
	string outStub= "";
	int  Reff1kgArgPos = -1; 
	 for(int i = 1; i< argc; i++)
	{
	//	cout << i << " = " << argv[i] << endl; 
	}
	//cout <<"****************************************************************************************" << endl;
	for(int i = 1; i< argc; i++)
	{
		string p = argv[i];
	//	cout << i << " = " << argv[i]<< endl;
		if( p == "-h")
		{
			//print help 
			cout << helptext << endl;
			return 0; 
		}
		else if (p == "-o")
		{
			outStub =  argv[i+1];
			i++;
		}
		else if (p == "-hf")
		{
			HashListFile =  argv[i+1];
			i++;
		} 
		else if (p == "-c")
                {
        //                cout << "Par Hash = " << argv[i+1] << endl;
                        Reff1kgArgPos = i+1;
                        i=i+1;
                }
		else if (p == "-hs")
                {
          //              cout << "Hash Size= " << argv[i+1] << endl;
                        HashSize = atoi(argv[i+1]);
                        i+=1;
                }
		else
		{
			cout << "ERROR: unkown command line paramater -" <<  argv[i] << "-"<< endl;
			return 0; 
		}
		
	}
	//check values 
	if (Reff1kgArgPos == -1)
	{
		cout << "ERROR 1kg Reference required" << endl;
		return 0; 
	}
	if (HashListFile == "")
	{
		cout << "Error HashList required" << endl;
		return 0; 
	}
	if (outStub == "")
	{
		cout << "ERROR out file stub required " << endl;
		return -1; 
		
	}
	//***********************************************
	

      	ifstream HashList;
      	HashList.open (HashListFile);
      	if ( HashList.is_open())
      	{	}// cout << "HashList Open " << HashListFile << endl;}   //cout << "##File Opend\n";
      	else
      	{
      		cout << "Error, HashList could not be opened";
      		return 0;
      	}
	ofstream Outfile; 
        Outfile.open(outStub);
	if (Outfile.is_open())
	{ cout <<  "Out file open " << outStub << endl;}
	else
        {
                cout << "Error, out could not be opened";
                return 0;
        }

	string line; 
	 
	long int Reader1kg;
        char* fName = argv[Reff1kgArgPos];
        Reader1kg =  open(fName, O_RDONLY);
	while ( getline(HashList, line))
	{
		vector<string> temp = Split(line, '\t');
		string hash = "";
		string hashcount = "";
                if (temp.size() ==2)
		{
			hash  = temp[0];
			hashcount = temp[1];
		
			//cout << temp[0]<<endl;
		}
		else 
		{
			hash = temp[3]; 
			hashcount = temp[2]; 
		}
                string rev = RevComp(hash);
                char *fileptr = NULL;
                if (hash < rev)
               	{
			int count = search(Reader1kg, hash, fileptr);
			cout << hash << " " << count << endl; 
			if (count == 0)
				Outfile << hash << "\t" << hashcount << endl;
			
		}
                else
                {	
			int count = search(Reader1kg, rev, fileptr); 
			cout << hash << " " << count << endl;
			if (count == 0)
				 Outfile << hash << "\t" << hashcount << endl;
		}

                if (munmap(fileptr, FILESIZE) == -1) {
                	perror("Error un-mmapping the file");
                }
       }
       close(Reader1kg);
	Outfile.close();	
}
	
