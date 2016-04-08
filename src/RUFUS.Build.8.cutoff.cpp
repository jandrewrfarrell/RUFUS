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
    double z= 0;
    double one = 1;
    cout << "0/1 = " << z/1 << endl;

	if ("AAA" < "AAC")
		cout << "yup" << endl;
	else
		cout << "nope" << endl;

    double MaxCount = 65534;
    double vm, rss, MAXvm, MAXrss;
    MAXvm = 0;
    MAXrss = 0;
    process_mem_usage(vm, rss, MAXvm, MAXrss);
//////////////////////lets parse some command line paramaters////////////////////////////
        string helptext;
        helptext = \
"\
RUFUS.interpret: converts RUFUS aligned contigs into a VCF \n\
By Andrew Farrell\n\
   The Marth Lab\n\
\n\
options:\
  -h [ --help ]          Print help message\n\
  -c  arg                Control or parent file pair, Hash file then dist file separeted by space\n\
  -s  arg                mutant or subject file paie, Hash file then dist file separeted by space \n\
  -o  arg                Output stub\n\
  -hs arg		 Hash size used in kmer counting\n\
  -mS arg		 minimum coverage in the subject\n\
  -mC arg 		 maximum coverage in the control samples to consider 0\n\
  -max arg		 Maximum coverage to consider\n\
  -t  arg 		 Number of threds to use\n\
";


	vector <string> ParentHashFilePaths; 
        string MutHashFilePath = "" ;
        string OutFileStub = "";
	int HashSize = -1;
	double MinCovSubject = -1.0; //probability of a mutation; 
	double MaxCovControl = -1.0; //probablitliy of a copy nubmer event; 
	double Pvalue = -1.0; //Pvlaue threshold
	int threads = 1;
	long MaxCoverage = 100000000;  
		
        for(int i = 1; i< argc; i++)
        {
                string p = argv[i];
                cout << i << " = " << argv[i]<< endl;
                if( p == "-h")
                {
                        //print help
                        cout << helptext << endl;
                        return 0;
                }
                else if (p == "-c")
                {
			cout << "Par Hash = " << argv[i+1] << endl;
                        ParentHashFilePaths.push_back(argv[i+1]);
                        i=i+1;
                }
                else if (p == "-s")
                {
			cout << "Sub Hash = " << argv[i+1] << endl;
			cout << "Sub Dist = " << argv[i+2] << endl;
                        MutHashFilePath = argv[i+1];
			i+=1;
                }
                else if (p == "-o")
                {
                        OutFileStub = argv[i+1];
                        i++;
                }
     		else if (p == "-hs")
                {
                       	HashSize =  atoi(argv[i+1]);
                        cout << "HS = " << HashSize << endl;
			i++;
                }
		else if (p == "-mS")
                {
                        MinCovSubject =  atof(argv[i+1]);
                 	cout << "MinCovSubject = " << MinCovSubject << endl; 
		       	i++;
                }
		else if (p == "-mC")
                {
                        MaxCovControl =  atof(argv[i+1]);
                	cout << "MaxCovControl = " << MaxCovControl << endl;
		        i++;
                }
		else if (p == "-max")
		{
			MaxCoverage = atoi(argv[i+1]);
			i++; 
		}
		else if (p == "-t")
                {
                        threads =  atoi(argv[i+1]);
                        i++;
                }
                else
                {
                        cout << "ERROR: unkown command line paramater -" <<  argv[i] << "-"<< endl;
                        return 0;
                }

        }	
	//check all the inputs :(
	if (ParentHashFilePaths.size() == 0 )
	{
		cout << "Error in control file inputs, either empty or number of file not the same for both" << endl;
		return 0; 
	}
	if (MutHashFilePath == "" )
	{
		cout << "Error subject file required, curnt = " << MutHashFilePath << endl; 
		return 0; 
	}	
	if (OutFileStub == "")
	{
		cout << "not out file provided, using Subject Hash path as prefix" << endl;
		OutFileStub = MutHashFilePath; 
	}
	if (HashSize == -1)
	{
		cout << "Error Hash Size required" << endl;
		return 0; 
	}
	if (MinCovSubject <0)
	{
		cout << "Error Minimum coverage in the subject required and must be greather than 0" << endl;
		return 0; 
	}
	if (MaxCovControl <0)
	{
		cout << "Error Maximum Coverage in the control must be set and must be greaterthan or equal to 0" << endl;
		return 0; 
	}
/////////////////////////////////donzo with command line paramaters 

	string LastLine = ""; 
	for (int c = 0; c < HashSize; c++)
	{
		LastLine += "T"; 
	}
	LastLine += "	0"; 
    //long BufferSize = 1000000 ;

	//in files 
//	ifstream ParentHash [100]; 
	
	for(int pi = 0; pi<ParentHashFilePaths.size(); pi++)
	{
		ParentHash[pi].open(ParentHashFilePaths[pi].c_str()); 
		if ( ParentHash[pi].is_open())
    		{ cout << "Parent Hash File open - " << ParentHashFilePaths[pi] << endl;}    //cout << "##File Opend\n";
    		else
    		{
        		cout << "Error, ParentHashFile could not be opened" << endl << ParentHashFilePaths[pi] << endl;
        		return 0;
    		}
		

	}


    	ifstream MutHashFile;
    	MutHashFile.open (MutHashFilePath.c_str());
    	if ( MutHashFile.is_open())
    	{}
    	else
    	{
    	    cout << "Error, MutHashFile could not be opened";
    	    return 0;
    	}
	
	
	
	
	//out files 
	ofstream MutMutHashTable;
    	string MutHashTablePath = OutFileStub ;
    	MutMutHashTable.open (MutHashTablePath.c_str());
    	if (MutMutHashTable.is_open())
    	{}else{cout << "Error MutHashTableFilter file couldnt be opened - " << MutHashTablePath << endl; return 0;}
	
    	ofstream Trace;
    	string TracePath = OutFileStub;
    	TracePath += ".Buld.Log";
    	Trace.open (TracePath.c_str());
    	if (Trace.is_open())
    	{}else{cout << "Error MutHashTableFilter file couldnt be opened - " << TracePath << endl; return 0;}
	
	
//////////////lets get this started



string line; 

////////////////////////////////////////////////////////////////////////////////////////////////////////


    long MutCount = 0; 
    long lines = 0;

    string ML1, PL1;
 

    clock_t St,Et;
    struct timeval start, end;
    gettimeofday(&start, NULL);
    bool notdone = true;
    int r= 0;
    cout << "starting " << endl;
    string MutLine; 
    //prime loop, grab first hash from each parent file
	long int HetCount = 0; 
    vector<string> ParLines; 
    //cout << "Priming the parent read" << endl;
    for (int pi = 0; pi<ParentHashFilePaths.size(); pi++)
    {
//	cout << "File " << pi << " = " << ParentHashFilePaths[pi] << endl;
	string line; 
	getline(ParentHash[pi], line); 
	ParLines.push_back(line); 
//	cout << "first line = " << line << endl; 
    }
	
//    for (int i =0; i < ParLines.size(); i++)
//    {
//	cout << ParLines[i] << endl;
//    }
    while (getline(MutHashFile, MutLine))
    {
        lines++;
 	//if (lines > 1000000)
        //{
        //        cout <<  MutLine << endl;
        //        return 0;
        //}

//	cout << "Mut Line is " << MutLine << endl;
	vector<string> Mline = Split(MutLine, '\t');
	int Acount =  atoi(Mline[1].c_str());
//	cout << "Acount = " << Acount << endl; 
	if (Acount >= MinCovSubject and Acount <= MaxCoverage)
	{
//		cout << "MutPassedCoverage" << endl;
//		cout << "reading parents" << endl;
		for (int pi = 0; pi<ParentHashFilePaths.size(); pi++)
    		{
//			cout << "Parent File " << pi<< endl;

			vector<string> Pline = Split(ParLines[pi], '\t');	
			while (Pline[0]<Mline[0])
			{
				string last = ParLines[pi]; 
				getline(ParentHash[pi], ParLines[pi]);
				if(ParLines[pi].size() == 0)
 				{
					cout << "Parent " << ParentHashFilePaths[pi] << " died at " << last << endl;
         				ParLines[pi] = LastLine;
 				}	
				
				Pline = Split(ParLines[pi], '\t');
			} //will need to add here check to see if at eof
		}

		stringstream Bcounts; 
		int TotalParentDepth = 0; 
		for (int pi = 0; pi < ParentHashFilePaths.size(); pi++)
        	{
        	    vector<string> Pline = Split(ParLines[pi], '\t');
  //      	    cout << "Par values are " << Pline[0] << "and " << Pline[1] << endl;
		    if (Pline[0]==Mline[0]) 
        	    {
			int Bcount = atoi(Pline[1].c_str());
			TotalParentDepth+=Bcount; 
			
			Bcounts << Bcount << "-"; 
			stringstream ss;
        	        ss << Pline[0] << '\t' << Bcount << '\t' << Acount << '\t'; 
		    
		    }
        	}
//		cout << "TotalParentDepth = " << TotalParentDepth << endl;
		///this is the simplest form, needs more thought
		if (Acount >= MinCovSubject and TotalParentDepth <= MaxCovControl)
		{
//			cout << "Found One Sub = " << Acount << " Control = " << MaxCovControl << endl;
        	        #pragma omp critical(h)
        	        { MutMutHashTable << HashToLong(Mline[0].c_str()) << "\t" << TotalParentDepth << "\t" << Acount << "\t" << Mline[0] << endl; HetCount++;}
		}
				 
	}
        if ((lines) % 10000 == 0 )
        {
            gettimeofday(&end, NULL);
            double Dt = end.tv_sec - start.tv_sec;
            process_mem_usage(vm, rss, MAXvm, MAXrss);
            cout << "Read in " << lines << " lines " << "; VM: " << vm << "; RSS: " << rss << " muts = "  << HetCount << " , lines per second = " << 1.0/(Dt/(lines)) <<"\r";
        }
	//if (WriteTrace){Trace << TraceBuffer.str();}
    }
	
    cout << "\nDone reading Parent File\n";

    int who = RUSAGE_SELF;
    struct rusage usage;
    int b = getrusage(RUSAGE_SELF, &usage);
    cout << "I am using " << usage.ru_maxrss << endl;

    MutMutHashTable.close();
    Trace.close();

cout << "\nreally dont\n";
}

