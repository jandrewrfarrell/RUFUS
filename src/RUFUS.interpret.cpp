
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
#include "externals/fastahack/Fasta.h"
#include <fcntl.h>
#include <sys/types.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <errno.h>

#define NUMINTS  (1000)
#define FILESIZE (NUMINTS * sizeof(int))

using namespace std;

vector <unordered_map <unsigned long int, int> > ParentHashes; 
unordered_map  <unsigned long int, int> MutantHashes; 

FastaReference Reff;
int HashSize = 25; 
int totalDeleted; 
int totalAdded;
int MaxVarentSize = 1000; 
ofstream VCFOutFile;
ofstream BEDOutFile;
ofstream BEDBigStuff;
ofstream BEDNotHandled;
ofstream Invertions; 
ofstream Translocations; 
ofstream Unaligned; 
map <string, int> Hash; 
/////////////////////////
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
	//	      cout << "n found";
			if (firstNew != true)
				firstNew = true;
			else
			{
				vector<string> stuff;
				stuff = Split(line, '\t');
				string PageHash = stuff[0];
	//		      cout << "PageHash " << endl;
				if (hash == PageHash)
				{
	  //			      cout << "found a hash " << hash  << " - " << PageHash;
					return atoi(stuff[1].c_str());
				}
				line = "";
			}
		}
		else if (firstNew == true)
		{
	//	      cout << "first Newline found" << endl;
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
	stuff = Split(line, '\t');
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
	stuff = Split(line, '\t');
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
	 //     cout << "found on first page" << endl;
	 	fileptr = (char*)mmap64(NULL, pageSize, PROT_READ, MAP_PRIVATE | MAP_POPULATE, fd, firstPos*pageSize);
		int val =  checkPage(data, hash, pageSize, "");
		if (munmap(fileptr, pageSize) == -1) {
			perror("Error un-mmapping the file");
		}
		return val; 
	}
	if (hash < FirstPageFirstHash)
	{
		cout << "HASH NOT IN FILE " << hash << endl;
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
	//	cout << "found on last page" << endl;
		int val =  checkPage(data, hash, pageSize, "");
		if (munmap(fileptr, pageSize) == -1) {
			perror("Error un-mmapping the file");
		}
		return val; 
	}
	if (hash > LastPageLastHash)
	{
		cout << "HASH NOT IN FILE " << hash << endl;
		return 0;
	}
	//start the search
	int counter = 0;
	while (true)
	{
	 //     cout << "ON LOOP " << counter << endl << endl;
		counter++;
		long int currentPage = lastPos - ((lastPos-firstPos)/2);
	   //   cout << "checking page " << currentPage << " last = " << lastPos << " and first = " << firstPos << endl;;
		if (currentPage == lastPos or currentPage == firstPos or lastPos - firstPos < 3)
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
	//	cout << " with " << CurrentPageFirstHash << " and " << CurrentPageLastHash << endl;
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
	  //		    cout << "hash " << hash << " is greater than " << CurrentPageFirstHash << " looking above" << endl;
				lastPos = currentPage;
				LastPageFirstHash = CurrentPageFirstHash;
				LastPageLastHash = CurrentPageLastHash;
			}
			else if (hash > CurrentPageLastHash)
			{
	    //		  cout << "hash \n" << hash << " is less than \n" << CurrentPageLastHash << " looking below" << endl;
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
int GetReadOrientation(int flag)
{
	// returns 0 = forward, 1 = reverse
	if (flag == 0 or flag == 2048)
		return 0;
	else if (flag == 16 or flag == 2064)
		return 1;
	else
		cout <<  "error flag not recognised " << flag << endl;

	return -1;
}
string getHash(string seq, int j, int HashSize)
{
	int bases = 0;
	string NewHash = "";
	while (bases<HashSize and j<seq.size())
	{
		if (seq.c_str()[j] == 'A' or seq.c_str()[j] == 'C' or seq.c_str()[j] == 'G' or seq.c_str()[j] == 'T')
		{
			NewHash+=seq.c_str()[j];
			bases++;
		}
		j++;
	}
	return NewHash;
}

string RevComp (string Sequence)
{
	string NewString = "";
	//cout << "Start - " << Sequence << "\n";
	//cout << Sequence.length() << endl; 
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
		{
			cout << "ERROR IN RevComp - " << C << " " ;
			NewString += C; 
		}

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

class SamRead
{
	public:
	string name;
	int flag;
	string chr;
	int pos;
	int mapQual;
	string cigar;
	string seq;
	string qual;
	string RefSeq; 
	string originalSeq;
	string originalQual; 
	string cigarString;
	string strand;  
	vector<int> alignments; 
	vector<int> Positions; 
	vector<string> ChrPositions; 
	bool first; // = true;  
	bool combined; // = false; 

	void parse(string read);
	void getRefSeq(); 
	void processCigar();
	void parseMutations( char *argv[] );
	void processMultiAlignment(); 
	void write();
	void writeVertical();
	void writetofile(ofstream &out);
	void flipRead(); 
};

void SamRead::flipRead()
{
	cout <<"FLIPPING reads not on the same strand";
	write(); 
	string FlipSeq = ""  ;
	string FlipQual = ""   ; 
	string FlipRefSeq = "";
	string FlipCigarString = ""; 
	string FlipStrand = "";
	vector<int> FlipPos;
	vector<string> FlipChrPos; 
	for (int i = seq.size() -1;  i >=0; i--)
	{
	//	FlipSeq += seq.c_str()[i]; 
		FlipQual += qual.c_str()[i];
	//	FlipRefSeq += RefSeq.c_str()[i];
		FlipCigarString += cigarString.c_str()[i];
		FlipStrand += '-'; 
		FlipPos.push_back(Positions[i]);
		FlipChrPos.push_back(ChrPositions[i]);
	}
	FlipSeq = RevComp(seq); 
	FlipRefSeq = RevComp(RefSeq);

	seq = FlipSeq; 
	qual = FlipQual; 
	RefSeq = FlipRefSeq;
	cigarString = FlipCigarString; 
	strand = FlipStrand; 
	Positions = FlipPos; 
	ChrPositions = FlipChrPos; 
	
	write(); 
	
}

void SamRead::processMultiAlignment()
{
	//check if this is a mis-joined contig

}
void SamRead::write()
{
	cout << name << endl;
	cout << "   flag = " << flag << endl;
	cout << "   mapQual = " << mapQual << endl;
	cout << "   Strand = " << GetReadOrientation(flag) << endl;
	cout << "   Alignments = " << alignments.size() << endl;
	cout << "   chr " << chr << " - " << pos << " qual = " << mapQual <<  endl;
	cout << "   cigar  = " << cigar << endl;
	cout << "   Seq    = " << seq << endl;
	cout << "   Qual   = " << qual << endl; 
	cout << "   Cigar  = " << cigarString << endl; 
	cout << "   RefSeq = " << RefSeq << endl;
	cout << "   strand = " << strand << endl;
	cout << "   RefPositions: ";
	for (int i =0; i < Positions.size(); i++)
		cout << Positions[i] << " \t"; 
	cout << endl;
	cout << "   RefChromoso: "; 
	for (int i =0; i < ChrPositions.size(); i++)
	       	cout << ChrPositions[i] << " \t";
	cout << endl;
	
}
void SamRead::writetofile(ofstream &out)
{
	out << name << endl;
	out << "   flag = " << flag << endl;
	out << "   mapQual = " << mapQual << endl; 
       out << "   Strand = " << GetReadOrientation(flag) << endl;
	out << "   Alignments = " << alignments.size() << endl;
	out << "   chr " << chr << " - " << pos << " qual = " << mapQual <<  endl;
	out << "   cigar  = " << cigar << endl;
	out << "   Seq    = " << seq << endl;
	out << "   Qual   = " << qual << endl;
	out << "   Cigar  = " << cigarString << endl;
	out << "   RefSeq = " << RefSeq << endl;
}

void SamRead::writeVertical()
{
	cout << name << endl;
	cout << "   flag = " << flag << endl;
	cout << "   chr " << chr << " - " << pos << " qual = " << mapQual <<  endl;
	cout << "   cigar  = " << cigar << endl;
	cout << "   Alignments = " << alignments.size() << endl; 
	for (int i =0; i < seq.size(); i++)
		cout << seq.c_str()[i] << "\t" << qual.c_str()[i] << "\t" << cigarString.c_str()[i] << "\t" << RefSeq.c_str()[i] << "\t" << ChrPositions[i] << "\t" << Positions[i] << endl;
}
string compressVar(string line, int start, string& StructCall)
{      
	cout << "compressing var" << endl; 
	char current = line.c_str()[0];
	int currentCount = 1;
	string CV = ""; 
	for (int i = 1; i< line.size(); i++)
	{       
		cout << current << endl;
		if (line.c_str()[i] == current)
		{       
			currentCount++;
		}
		else    
		{       
			if (currentCount > 2)
			{       
				ostringstream convert;
				convert << currentCount; 
				CV+=convert.str();
				CV+=current;
				
				ostringstream convertEND; 
				int end = currentCount+start; 
				convertEND << end; 

				
				if (current == 'T')
				{
					cout << "YAAAY STRUCT" << endl; 
					StructCall = "SVTYPE=DUP;END="; 
					StructCall += convertEND.str(); 
					StructCall += ";SVLEN="; 
					StructCall += convert.str(); 
					StructCall += ";";
					cout << StructCall << endl; 
				}
			}
			else if (currentCount ==2)
			{       
				CV+=current;
				CV+=current;
			}
			else if (currentCount == 1)
			{       
				CV+=current;
			}
			else    
			{       
				cout << "ERROR in compress " << current << " " << currentCount << endl;
			}
			
			current = line.c_str()[i];
			currentCount = 1;
		}
	}
	if (currentCount > 2)
	{       
		ostringstream convert;
		convert << currentCount;
		CV+=convert.str();
		CV+=current;
		
		ostringstream convertEND;
		int end = currentCount+start;
		convertEND << end;


		if (current == 'T')
		{
			cout << "YAAAY STRUCT" << endl;
			StructCall = "SVTYPE=DUP:TANDEM;END=";
			StructCall += convertEND.str();
			StructCall += ";SVLEN=";
			StructCall += convert.str();
			StructCall += ";";
			cout << StructCall << endl;
		}
	}
	else if (currentCount ==2)
	{       
		CV+=current;
		CV+=current;
	}
	else if (currentCount == 1)
	{       
		CV+=current;
	}
	else
	{
		cout << "ERROR in compress " << current << " " << currentCount << endl;
	}
	
	return CV;
}
void SamRead::parseMutations( char *argv[])
{
	
	cout << "Parsing Mutations " << endl;
	write(); 
	string StructCall = "ST"; 	
	vector<bool> PeakMap; 
	int max = -1; 
	int maxSpot = -1;
	int i =0;
	while (i<qual.size())
	{
		int j;
		for ( j = i+1; j < qual.size()-1; j++)
		{
			if(qual.c_str()[j+1] < qual.c_str()[j])
			{
				max =  qual.c_str()[j];
				while (qual.c_str()[j+1] < qual.c_str()[j] && j < qual.size()-1)
				{
					j++;
				}
				break; 	
			}		
		}
		int k; 
		for ( k = i; k <=j; k++)
		{
			if (qual.c_str()[k] == max)
				PeakMap.push_back(1);
			else	
				PeakMap.push_back(0);
		}
		i = i +k; 
	}	
	
	//cout << "testhing hash search" << endl;
	vector <string> hashes;
	vector <bool> varHash; 
	for (int i = 0; i <  seq.size() - HashSize; i++)
	{
		string newHash = "";
		newHash += seq.c_str()[i];
		int count = 0; 
		if ((cigarString.c_str()[i] != 'D' and cigarString.c_str()[i] != 'R' and cigarString.c_str()[i] != 'H'))
		{
			for (int j = 1; j<seq.size() - i and count < HashSize-1 ; j++)
			{
				
				if (cigarString.c_str()[i+j] != 'D' and cigarString.c_str()[i+j] != 'R'  and cigarString.c_str()[i+j] != 'H')
				{
					newHash += seq.c_str()[i+j]; 
					count++;
		//			cout << cigarString.c_str()[i+j] << " - " << i+j<< endl;
		//			cout << newHash << endl; 
				}
	//			else
		//		cout << "yay cigar was not ok - "   << cigarString.c_str()[i+j]<< endl;
			}
		}
		hashes.push_back(newHash);
		//cout << newHash << endl; 
		if (Hash.count(newHash) > 0 or Hash.count(RevComp(newHash)) > 0)
			varHash.push_back(true); 
		else
			varHash.push_back(false); 
	}
	//cout << "made hash list to check " << endl; 
	vector <vector<int>> parentCounts;
	//vector <long int> ParentHash;
	for(int pi = 0; pi<ParentHashes.size(); pi++)
	{
		vector <int> counts;
		for(int i = 0; i< hashes.size(); i++)
		{
			string hash =  hashes[i];
			//cout << "Hash = " << hash << endl;
			bool checkHash = true; 
			for (int j = 0; j < 25; j++)
			{
				if (!(hash[j] == 'A' or hash[j] == 'C' or hash[j] == 'G' or hash[j] == 'T'))
				{
					checkHash = false; 
					break; 
				}
			}
			if (checkHash)
			{
				//cout << "passed" << endl;
				unsigned long int LongHash = HashToLong(hash);
				//cout << "herhe " << LongHash << endl;
				if (ParentHashes[pi].count(LongHash) >0)
				{
					//cout << "exists" << endl;
					counts.push_back(ParentHashes[pi][LongHash]);
				}
				else
				{
					//cout << "dosnt exist" << endl;
					counts.push_back(0);
				}
			}
			else 
			{
				//cout << "didnt pass" << endl;
				counts.push_back(-1);
			} 
		}
		parentCounts.push_back(counts);
	}
	//cout << "here" << endl;
	for(int i =0; i < hashes.size(); i++)
	{
		cout << i << "\t" << hashes[i] << "\t" << varHash[i] << "\t" << PeakMap[i] << "\t" << (int) qual.c_str()[i]-33;
		for (int j = 0; j < parentCounts.size(); j++)
		{
				cout << "\t" <<  parentCounts[j][i];
		}
		cout << endl;
	}

	
	//cout << "Peak map = " << PeakMap.size() << endl;
	//cout << "qual string = " << qual.size() << endl;
	for (int i = 0; i < qual.size(); i++)
	{
		int qualcount =  qual.c_str()[i] - 33; 
	//	cout <<  qual.c_str()[i] << " - " << qualcount << " - " << PeakMap[i] << seq.c_str()[i]<< endl;
	}

	string reff = "";
	string alt = "";
	string varType = "";
	for(int i = 0; i<cigarString.size(); i++)
	{
		reff = ""; 
		alt = "";
		varType = ""; 
		if ((cigarString.c_str()[i] == 'X' or cigarString.c_str()[i] == 'I' or cigarString.c_str()[i] == 'D' or cigarString.c_str()[i] == 'T'/*  or cigarString.c_str()[i] == 'S' *or cigarString.c_str()[i] == 'H'*/) and RefSeq.c_str()[i] != 'N')
		{
			int size = -1; 
			int startPos = i; 
			bool  AnyBasesOver0  = false; 
			string Denovo = "inherited"; 
			//varType += cigarString.c_str()[i]; 
			if (qual.c_str()[i] > '!')
				AnyBasesOver0 = true;
			if (PeakMap[i] == 1)
				Denovo = "DeNovo";
			for(int j = 0; j< cigarString.size() - i; j++)
			{
				//if(cigarString.c_str()[i+j] == cigarString.c_str()[i])
				if(cigarString.c_str()[i+j] == 'X' or cigarString.c_str()[i+j] == 'D' or cigarString.c_str()[i+j] == 'I' or cigarString.c_str()[i+j] == 'T' /*or cigarString.c_str()[i+j] == 'S' or cigarString.c_str()[i+j] == 'H'*/)
				{
					size = j; 
					if (qual.c_str()[i+j] > '!')
						AnyBasesOver0 = true;
					if (PeakMap[i+j] == 1)
						Denovo = "DeNovo";
	
				}
				else //if (qual.c_str()[i+j] == '!')
					break; 
			}
			bool LowCov = false; 
			int low = i - 30; 
			if (low < 0){low = 0;}
			cout << "checking bases " << low << " to " << i+size+5 << endl;
			for(int j = low; j < i+size+5 and j < hashes.size(); j++)
			{
				bool AllPar0 = true; 
				for (int k = 0; k < parentCounts.size(); k++)
				{
					if (parentCounts[k][j] != 0 )
						AllPar0 = false;	
				}
				if ( AllPar0 )
				{
					varHash[j] = true; 
				}
				if (hashes[j].size() == HashSize and varHash[j] == false)
				{
					 cout << "base " << j ; 
					 
					for (int k = 0; k < parentCounts.size(); k++)
					{
						if (parentCounts[k][j] <= 5 and parentCounts[k][j] > 0)
						{
							cout << "\t" << parentCounts[k][j] ; 
							LowCov = true;
						}
						 
					}
					cout << endl;
				}
			}
			if (AnyBasesOver0)  //enabling this will only report varites covered by hashes 
			{
				if ( cigarString.c_str()[i] == 'I' or cigarString.c_str()[i] == 'D' or cigarString.c_str()[i] == 'T' /*or cigarString.c_str()[i] == 'S' or cigarString.c_str()[i] == 'H'*/)
				{
					//reff+=RefSeq.c_str()[i-1];
					//alt+=seq.c_str()[i-1];
					//startPos = i-1;
						
					for (int k = 1; i-k >= 0; k++)
					{
						if (ChrPositions[i-k] == "nope")
						{}
						else
						{
							reff+=RefSeq.c_str()[i-k]; 
							alt+=seq.c_str()[i-k];
							startPos = i-k;
							break;
						}
					}  
				}
				//if ( cigarString.c_str()[i] == 'I' )
				//	varType += "INS";
				//if (cigarString.c_str()[i] == 'D')
				//	varType += "DEL";  	
				//if (cigarString.c_str()[i] == 'X')
				//	varType += "SNP"; 
				for(int j = 0; j<= size; j++)
				{
					if (RefSeq.c_str()[i+j] == 'A' or RefSeq.c_str()[i+j] == 'C' or RefSeq.c_str()[i+j] == 'G' or RefSeq.c_str()[i+j] == 'T')
						reff+=RefSeq.c_str()[i+j]; 
					if (seq.c_str()[i+j] == 'A' or seq.c_str()[i+j] == 'C' or seq.c_str()[i+j] == 'G' or seq.c_str()[i+j] == 'T')
						alt+=seq.c_str()[i+j]; 
					varType += cigarString.c_str()[i+j]; 
				}
	
				//***********Build up hash depth nehborhod********************
				//cout << seq << endl;
				//for (int j = 0; j <i; j++)
				//{	cout << " ";}
				//cout << alt << endl;
					
			
				int lower = i-HashSize; 
				if (lower < 0){lower =0;}
				int upper = i+alt.length()+1; 
				vector <int> HashCounts; 
				vector <int> HashCountsOG; 
				for (int j = lower; j<upper; j++)
				{
				//	for (int k = 0; k<j; k++)
				//	{       cout << " ";}
					string hash = getHash(seq, j, HashSize);
					string Revhash = RevComp(hash);
				//	cout << hash << "\t";
					if (Hash.count(hash) > 0)
					{
						HashCounts.push_back(Hash[hash]);
						HashCountsOG.push_back(Hash[hash]);
				//		cout << Hash[hash] << endl;
					}
					else if (Hash.count(Revhash) > 0)
					{
						HashCounts.push_back(Hash[Revhash]);
						HashCountsOG.push_back(Hash[Revhash]);
				//	cout << Hash[Revhash] << endl;
					}
					else
					{
						HashCounts.push_back(-1); 
				//		cout << -1 << endl;
					}
				}
				
				//***********check that the alese are only baess**************
				bool good = true; 
				for (int j = 0; j<reff.size(); j++)
				{	if (reff.c_str()[j] != 'A' or reff.c_str()[j] != 'C' or reff.c_str()[j] != 'G' or reff.c_str()[j] != 'T'){good = false;}}
				for (int j = 0; j<alt.size(); j++)
				{	if (alt.c_str()[j] != 'A' or alt.c_str()[j] != 'C' or alt.c_str()[j] != 'G' or alt.c_str()[j] != 'T'){good = false;}}
				if (good = false)
				{
					cout <<"ERROR in SNP detect"<< endl;
					cout <<endl<< chr << "\t" << pos+i << "\t" << reff << "\t" << alt << endl;
					write();
				}
				 //***********check that the alese are only baess done**************
				
				//write();
				string CompressedVarType = compressVar(varType, Positions[startPos], StructCall); 
					
				cout <<  chr << "\t" << pos+i << "\t" << CompressedVarType /*"."*/ << "\t" << reff << "\t" << alt << "\t" << HashCountsOG.size() << "\t" << varType << "\t" << "." << "\t" << "." << "\t" << "." << endl;
				if (LowCov)
				{
					cout << "LOW COVERAGE" << endl;
					Denovo = "LowCov"; 
				}
				else
				   	cout << "GOOD COVERAGE" << endl;
				cout << ChrPositions[startPos] << "\t" <<Positions[startPos] << "\t" << CompressedVarType <<"-" <<Denovo /*"."*/  << "\t" << reff << "\t" << alt << "\t" << HashCountsOG.size() << "\t" << "." << "\t" << StructCall << "CVT=" << CompressedVarType << ";HD=";	
				VCFOutFile << ChrPositions[startPos] << "\t" <<Positions[startPos] << "\t" << CompressedVarType <<"-" <<Denovo /*"."*/  << "\t" << reff << "\t" << alt << "\t" << HashCountsOG.size() << "\t" << "." << "\t"  << StructCall <<";RN=" << name << ";MQ=" << mapQual << ";cigar=" << cigar << ";" << ";CVT=" << CompressedVarType << ";HD="; 
				for (int j = 0; j < HashCounts.size(); j++)
				{	VCFOutFile << HashCounts[j] << "_"; }
				if (HashCountsOG.size()>0)
				{
					std::sort (HashCountsOG.begin(), HashCountsOG.end());
					VCFOutFile << ";AO=" << HashCountsOG[HashCountsOG.size()/2];
				}
				else
				{
					 VCFOutFile << ";AO=" << "-1";
				} 
				VCFOutFile <<  ";VT=" <<  varType << "\t" << "GT:AK" << "\t" << "0/1:"<< "7" << endl; 
				BEDOutFile << chr << "\t" << pos+i << "\t" <<  pos+i+size << "\t" << chr << ":" << pos+i << ":" << (int)(reff.length() - alt.length()) << ":" << HashCountsOG.size() << endl;
				i+=size;
			} 
		}
		/*else if (cigarString.c_str()[i] == 'D'  and RefSeq.c_str()[i] != 'N')
		{
			cout << "DELETION" << endl;
			//write(); 
		      	int size = -1;
			bool AnyBasesOver0 = false; 
			for(int j = 0; j< cigarString.size() - i; j++)
			{
				if(cigarString.c_str()[i+j] == 'D' or cigarString.c_str()[i+j] == 'X' or cigarString.c_str()[i+j] == 'I' )//or cigarString.c_str()[i+j] == 'S' or cigarString.c_str()[i+j] == 'H')
				{
					size = j;
					 if (RefSeq.c_str()[i+j] == 'A' or RefSeq.c_str()[i+j] == 'C' or RefSeq.c_str()[i+j] == 'G' or RefSeq.c_str()[i+j] == 'T')
						reff+=RefSeq.c_str()[i+j];
					if (seq.c_str()[i+j] == 'A' or seq.c_str()[i+j] == 'C' or seq.c_str()[i+j] == 'G' or seq.c_str()[i+j] == 'T')
						alt+=seq.c_str()[i+j];	
					
					if (qual.c_str()[i+j] > '!')
						AnyBasesOver0 = true; 
				}
				else
					break;
			}
			cout << "Raw reff = " << reff << endl << "raw alt  = " << alt << endl;
			if (AnyBasesOver0 and size < MaxVarentSize)
			{
				cout << "size = " << size<< endl;
				reff = RefSeq.c_str()[i-1] + reff;
				alt = seq.c_str()[i-1] + alt;

				cout << "reff = " << reff << endl << "alt  = " << alt << endl;	
				//for(int j = 0; j<= size; j++)
				//{
				//	reff+=RefSeq.c_str()[i+j];
				//}

				 //***********Build up hash depth nehborhod********************
				//cout << seq << endl;
				//for (int j = 0; j <i; j++)
				//{       cout << " ";}
				//cout << alt << endl;
	
	
				int lower = i-HashSize;
				if (lower < 0){lower =0;}
				int upper = i+2;
				vector <int> HashCounts;
				vector <int> HashCountsOG;
				for (int j = lower; j<upper; j++)
			       	{
       				 //	for (int k = 0; k<j; k++)
       				 //	{       cout << " ";}
       				 	string hash = getHash(seq, j, HashSize);
       				 	string Revhash = RevComp(hash);
				//	cout << hash << "\t";
					if (Hash.count(hash) > 0)
					{
						HashCounts.push_back(Hash[hash]);
						HashCountsOG.push_back(Hash[hash]);
				//		cout << Hash[hash] << endl;
					}
					else if (Hash.count(Revhash) > 0)
					{
						HashCounts.push_back(Hash[Revhash]);
						HashCountsOG.push_back(Hash[Revhash]);
				//		cout << Hash[Revhash] << endl;
					}
					else
					{
						HashCounts.push_back(-1);
				//		cout << -1 << endl;
					}
				}


				//***********check that the alese are only baess**************
				bool good = true;
				for (int j = 0; j<reff.size(); j++)
				{       if (reff.c_str()[j] != 'A' or reff.c_str()[j] != 'C' or reff.c_str()[j] != 'G' or reff.c_str()[j] != 'T'){good = false;}}
				for (int j = 0; j<alt.size(); j++)
       			 	{	if (alt.c_str()[j] != 'A' or alt.c_str()[j] != 'C' or alt.c_str()[j] != 'G' or alt.c_str()[j] != 'T'){good = false;}}
       			 	if (good = false)
       			 	{
       				 	cout <<"ERROR in DELETION detect"<< endl;
       				 	cout <<endl<< chr << "\t" << pos+i << "\t" << reff << "\t" << alt << endl;
       				 	write();
       			 	}	
       			  	//***********check that the alese are only baess done**************
				//write(); 
				cout <<chr << "\t" << pos+i-1 << "\t" << "." << "\t" << reff << "\t" << alt << "\t" << HashCountsOG.size() << "\t" << "DEL" << "\t" << "." << "\t" << "." << "\t" << "." << endl;
				//  "#CHROM 	      POS	       ID	      REF	    ALT	     QUAL	   FILTER	     INFO	  FORMAT	Sample";
				
				VCFOutFile << ChrPositions[i] << "\t" << Positions[i-1] << "\t" << "." << "\t" << reff << "\t" << alt << "\t" << HashCountsOG.size() << "\t" << "DEL" << "\t" "HD=";;  
				for (int j = 0; j < HashCounts.size(); j++)
				{       VCFOutFile << HashCounts[j] << "_"; }	
				if (HashCountsOG.size()>0)
				{
					std::sort (HashCountsOG.begin(), HashCountsOG.end());
					VCFOutFile << ";AO=" << HashCountsOG[HashCountsOG.size()/2];
				}
				else
				{
					 VCFOutFile << ";AO=" << "-1";
				}
				VCFOutFile << "\t" << "." << "\t" << "."<< endl;

				BEDOutFile << ChrPositions[i] << "\t" << Positions[i-1] << "\t" <<  pos+i+size << endl;
       				i+=size;
			}
			else
				cout << "nobase :):" << endl; 
		}*/
		/*else if ((cigarString.c_str()[i] == 'I' )  and RefSeq.c_str()[i] != 'N')
		{
			cout << "Working on insertion" << endl;
			bool AnyBasesOver0 = false;	
			int size = -1;
			for(int j = 0; j< cigarString.size() - i; j++)
			{
			       	if(cigarString.c_str()[i+j] == 'X' or cigarString.c_str()[i+j] == 'D' or cigarString.c_str()[i+j] == 'I' or cigarString.c_str()[i+j] == 'S' or cigarString.c_str()[i+j] == 'H')
				{
					size = j;
					if (qual.c_str()[i+j] > '!')
						AnyBasesOver0 = true;
				}
				else
					break;
			}
			if (AnyBasesOver0 and size < MaxVarentSize)
			{
				reff += RefSeq.c_str()[i-1];
				alt += 	seq.c_str()[i-1];
		
       				for(int j = 0; j<= size; j++)
				{
					if (seq.c_str()[i+j] == 'A' or seq.c_str()[i+j] == 'G' or seq.c_str()[i+j] == 'C' or seq.c_str()[i+j] == 'T')
	       			 		alt+=seq.c_str()[i+j];
					//else
					//	alt+=RefSeq.c_str()[i+j];
	       		 	}
				for(int j = 0; j<= size; j++)
				{
					if (RefSeq.c_str()[i+j] == 'A' or RefSeq.c_str()[i+j] == 'G' or RefSeq.c_str()[i+j] == 'C' or RefSeq.c_str()[i+j] == 'T')
						reff+=seq.c_str()[i+j];
				}
				  //***********Build up hash depth nehborhod********************
				//cout << seq << endl;
				//for (int j = 0; j <i; j++)
				//{       cout << " ";}
				//cout << alt << endl;
	
	
				int lower = i-HashSize;
				if (lower < 0){lower =0;}
				int upper = i+size+2;
				vector <int> HashCounts;
				vector <int> HashCountsOG;
				for (int j = lower; j<upper; j++)
				{
				//	for (int k = 0; k<j; k++)
				//	{       cout << " ";}
					string hash = getHash(seq, j, HashSize);
					string Revhash = RevComp(hash);
				//	cout << hash << "\t";
					if (Hash.count(hash) > 0)
					{
						HashCounts.push_back(Hash[hash]);
						HashCountsOG.push_back(Hash[hash]);
			//			cout << Hash[hash] << endl;
			       	 	}
					else if (Hash.count(Revhash) > 0)
					{
						HashCounts.push_back(Hash[Revhash]);
						HashCountsOG.push_back(Hash[Revhash]);
			  //      		cout << Hash[Revhash] << endl;
					}
					else
					{
						HashCounts.push_back(-1);
			    //    		cout << -1 << endl;
					}
				}
				//***********check that the alese are only baess**************
				bool good = true;
				for (int j = 0; j<reff.size(); j++)
				{       if (reff.c_str()[j] != 'A' or reff.c_str()[j] != 'C' or reff.c_str()[j] != 'G' or reff.c_str()[j] != 'T'){good = false;}}
				for (int j = 0; j<alt.size(); j++)
				{	if (alt.c_str()[j] != 'A' or alt.c_str()[j] != 'C' or alt.c_str()[j] != 'G' or alt.c_str()[j] != 'T'){good = false;}}
				if (good = false)
				{
					cout <<"ERROR in INSERTION detect"<< endl;
					cout <<endl<< chr << "\t" << pos+i << "\t" << reff << "\t" << alt << endl;
					write();
				}
				//***********check that the alese are only baess done**************
	
	
 				//write(); 
				//cout <<  chr << "\t" << pos+i-1 << "\t" << "." << "\t" << reff << "\t" << alt << "\t" << "100" << "\t" << "INS" << "\t" << "." << "\t" << "." << "\t" << "." << endl;
				//	    "#CHROM	      POS	       ID		REF	       ALT	       QUAL	      FILTER		 INFO	       FORMAT	     Sample";
	
				VCFOutFile << ChrPositions[i] << "\t" << Positions[i-1] << "\t" << "." << "\t" << reff << "\t" << alt << "\t" << HashCountsOG.size() << "\t" << "INS" << "\t" "HD="; 
				for (int j = 0; j < HashCounts.size(); j++)
				{       VCFOutFile << HashCounts[j] << "_"; }
				if (HashCountsOG.size()>0)
				{
					std::sort (HashCountsOG.begin(), HashCountsOG.end());
					VCFOutFile << ";AO=" << HashCountsOG[HashCountsOG.size()/2];
				}
				else
				{
					 VCFOutFile << ";AO=" << "-1";
				}
				VCFOutFile << "\t" << "." << "\t" << "."<< endl;
				BEDOutFile << ChrPositions[i] << "\t" << Positions[i-1] << "\t" <<  Positions[i]+size << endl;
       			 	i+=size;	
			}			
		}
		/*else if (cigarString.c_str()[i] == 'S'  and RefSeq.c_str()[i] != 'N')
		{
			bool AnyBasesOver0 = false;
			int size = -1;
			for(int j = 0; j< cigarString.size() - i; j++)
			{
				if(cigarString.c_str()[i+j] == 'S' or cigarString.c_str()[i+j] == 'I' )
				{
					size = j;
					if (qual.c_str()[i+j] > '!')
						AnyBasesOver0 = true;
				}
				else
					break;
			}
			if (AnyBasesOver0)
			{
				reff += RefSeq.c_str()[i-1];
				alt +=  seq.c_str()[i-1];

				for(int j = 0; j<= size; j++)
				{
					alt+=seq.c_str()[i+j];
				}
				//***********check that the alese are only baess**************
				bool good = true;
				for (int j = 0; j<reff.size(); j++)
				{       if (reff.c_str()[j] != 'A' or reff.c_str()[j] != 'C' or reff.c_str()[j] != 'G' or reff.c_str()[j] != 'T'){good = false;}}
				for (int j = 0; j<alt.size(); j++)
				{	if (alt.c_str()[j] != 'A' or alt.c_str()[j] != 'C' or alt.c_str()[j] != 'G' or alt.c_str()[j] != 'T'){good = false;}}
				if (good = false)
				{
					cout <<"ERROR in INSERTION detect"<< endl;
					cout <<endl<< chr << "\t" << pos+i << "\t" << reff << "\t" << alt << endl;
					write();
				}
				//***********check that the alese are only baess done**************


				//write();
				//cout <<  chr << "\t" << pos+i-1 << "\t" << "." << "\t" << reff << "\t" << alt << "\t" << "100" << "\t" << "INS" << "\t" << "." << "\t" << "." << "\t" << "." << endl;
				//	  "#CHROM	   POS	       ID	      REF	    ALT	     QUAL	   FILTER	     INFO	  FORMAT	Sample";

				VCFOutFile << chr << "\t" << pos+i-1 << "\t" << "." << "\t" << reff << "\t" << alt << "\t" << "100" << "\t" << "INS-SofClipped" << "\t" << "." << "\t" << "." << "\t" << "." << endl;
				BEDOutFile << chr << "\t" << pos+i-1 << "\t" <<  pos+i+size << endl;
				i+=size;
			}
		}*/
	}
}
void SamRead::processCigar()
{
	string num = ""; 
	for (int i = 0; i < cigar.length(); i++)
	{
		if (cigar.c_str()[i] >= 48 and cigar.c_str()[i] <= 57)
			num = num + cigar.c_str()[i]; 
		else 
		{	
			int number = atoi(num.c_str()); 
			for(int j = 0; j < number; j++)
			{cigarString += cigar.c_str()[i];}
			num = ""; 
		}
	}
	
}

void SamRead::getRefSeq()
{
	originalSeq = seq; 
	originalQual = qual; 
	RefSeq = ""; 
	string NewSeq = "";
	string NewQual = "";
	string NewCigar = "";
	string NewStrand = ""; 
	vector<int> NewPositions; 
	vector<string> NewChromosome;
	int InsOffset = 0;
	
	if ( Reff.sequenceNameStartingWith(chr) == "") //come back to, need to check if chr is in reference 
	{
		cout << "ERROR chr " << chr << " not found\n"; 
		return; 
	}

   //correct star position of the read to account for Hard and soft clipped bases as we are counting those now
	for (int i = 0; i< cigarString.size(); i++)
	{
		if (cigarString.c_str()[i]== 'H')
		{}
		else
		{
			pos = pos-i;
			break;
		}
	}
	for (int i = 0; i< cigarString.size(); i++)
      	{
		if (cigarString.c_str()[i]== 'S')
		{}
		else
		{
			pos = pos-i;
			break;
		}
	}
	



	int Roffset = 0; 
	int Coffset = 0; 
	for (int i =0; i<cigarString.length(); i++)
	{
		if (cigarString.c_str()[i] == 'M')
		{
			//cout << "yay cigar = M" << endl;
			RefSeq += toupper(Reff.getSubSequence(chr, i+pos-1+Roffset, 1).c_str()[0]);
			NewSeq += seq.c_str()[i-Coffset]; 
			NewQual += qual.c_str()[i-Coffset];
			NewPositions.push_back(pos+i-InsOffset);
			NewChromosome.push_back(chr);
			if (   toupper(Reff.getSubSequence(chr,i+pos-1+Roffset, 1).c_str()[0]) == seq.c_str()[i-Coffset])
				NewCigar +='M'; 
			else
				NewCigar +='X';

		}
		else if(cigarString.c_str()[i] == 'I')
		{
			// cout << "yay cigar = I" << endl;
			RefSeq+='-'; 
			Roffset+= -1;
			NewSeq += seq.c_str()[i-Coffset];
			NewQual += qual.c_str()[i-Coffset];
			NewCigar += 'I';
			InsOffset++; 
			NewPositions.push_back(pos+i-InsOffset);
			NewChromosome.push_back(chr);
		}
		else if(cigarString.c_str()[i] == 'D')
		{
			// cout << "yay cigar = D" << endl;
			NewSeq += '-';
			NewQual += ' '; 
			Coffset++;
			RefSeq += toupper(Reff.getSubSequence(chr, i+pos-1+Roffset, 1).c_str()[0]);
			NewCigar += 'D';
		 	NewPositions.push_back(pos+i-InsOffset);
			NewChromosome.push_back(chr);
		}

		else if(cigarString.c_str()[i] == 'H')
		{
			RefSeq += 'H';
			NewSeq += 'H';
			NewQual += ' '; 
			Coffset++; 
			NewCigar+='H';
		//	Roffset+= -1;
			NewPositions.push_back(-1);
			NewChromosome.push_back("nope");
		}
		else if(cigarString.c_str()[i] == 'S')
		{
			//I have an error here, when a read starts with S its ofset wrong
			RefSeq += '-';
			NewSeq += seq.c_str()[i-Coffset];
			NewQual += qual.c_str()[i-Coffset]; 
			NewCigar+='S';
		//	Roffset += -1; 
			NewPositions.push_back(pos+i-InsOffset); //NewPositions.push_back(-1);
			NewChromosome.push_back(chr); //NewChromosome.push_back("nope");
		}
		//{ cout << "yay cigar = H" << endl;}
		else{cout << "well shit" << endl;}
		//cout << "yay" << endl;
	}
	seq = NewSeq;
	cigarString = NewCigar; 
	qual.clear();
	char lastQ = ' ';

	//cout << "qual = " << qual; 
	for (int i =0; i<NewQual.size(); i++)
	{       
	 	if (NewQual.c_str()[i] == ' ')
			qual += lastQ;
		else 
		{
			qual+=NewQual.c_str()[i];
			lastQ = NewQual.c_str()[i];
		}
	}
	
	//correct star position of the read to account for Hard and soft clipped bases as we are counting those now
/*	for (int i = 0; i< cigarString.size(); i++)
	{
		if (cigarString.c_str()[i]== 'H')
		{}
		else
		{
			pos = pos-i; 		
			break;
		}
	}
	for (int i = 0; i< cigarString.size(); i++)
      {
		if (cigarString.c_str()[i]== 'S')
		{}
		else
		{
			pos = pos-i;
			break;
		}
	}*/
	//cout << "qual = " << qual;
	//build up the string that will indidcate if the read has to be flipped
	for (int i = 0; i<qual.size(); i++)
	{
		strand+="+"; 
	}
	Positions = NewPositions;
	ChrPositions = NewChromosome;
}
void SamRead::parse(string read)
{
	vector <string> temp = Split(read, '\t');
	name = temp[0];
	flag = atoi(temp[1].c_str());
	chr = temp[2];
	pos = atoi(temp[3].c_str());
	mapQual = atoi(temp[4].c_str());
	cigar = temp[5];
	seq = temp[9];
	if (temp[10] == "*")
	{
		cout << "correcting missing qualtiy" << endl; 
		string newQual = ""; 
		for (int i = 0; i < seq.size(); i++)
		{
			newQual+='5'; 
		}
		temp[10] = newQual; 
	}
	qual = temp[10];
	alignments.clear(); 
	
	first = true; 
	combined = false; 
	string NewSeq = ""; 
	for (int i = 0; i < seq.size(); i++)
		NewSeq+=toupper(seq.c_str()[i]);
	
	seq = NewSeq; 
	processCigar(); 
}

SamRead BetterWay(vector<SamRead> reads)
{
	
	int A = 0; 
	int B = 1;
	//cout << "B = " << B << endl;
	cout << "Working on " << reads[A].name << ", size = " << reads[B].pos -reads[A].pos << endl;
	vector<vector<int>> AlignmentPos;
	vector<vector<string>> AlignmentChr;
	int ALastGoodRef = -1;
	string ALastGoodChr = "nope"; 
	int BLastGoodRef = -1; 
	string BLastGoodChr = "nope";


	vector <string> NewSeqs; 
	vector <string> NewQuals; 
	vector <string> NewRefs; 
	vector <string> NewCigars; 

	for(int i = 0; i<reads.size(); i++)
	{
		NewSeqs.push_back("");
		NewQuals.push_back(""); 
		NewRefs.push_back("");
		NewCigars.push_back("");
	}
	
	
	for (int i = 0; i < reads.size(); i++)
	{
		cout << "Pre Lining up reads " << i << endl;
		reads[i].write();
	}
	int Acount = 0; 
	int Bcount = 0;
	 if (B == 1 and GetReadOrientation(reads[A].flag) != GetReadOrientation(reads[B].flag))
	{
		cout <<"FLIPPING reads not on the same strand"; 
		reads[B].flipRead(); 
	}
	while (Acount < reads[A].cigarString.size() and Bcount < reads[B].cigarString.size())
	{
		vector<int> currentPos;
		vector<string> currentChr;
		//need to get all the reads lined up with the same number of bases, taking account of I's and D's that change length//
		if (reads[A].cigarString.c_str()[Acount] == 'D' and reads[B].cigarString.c_str()[Acount] != 'D')
		{
			currentPos.push_back(reads[A].Positions[Acount]);
			currentPos.push_back(-1);

			currentChr.push_back(reads[A].ChrPositions[Acount]);
			currentChr.push_back("nope");
		
			ALastGoodRef = reads[A].Positions[Acount];
			ALastGoodChr = reads[A].ChrPositions[Acount];
			
			NewSeqs[A]+=reads[A].seq.c_str()[Acount];
			NewQuals[A]+=reads[A].qual.c_str()[Acount];
			NewRefs[A]+=reads[A].RefSeq.c_str()[Acount];
			NewCigars[A]+=reads[A].cigarString.c_str()[Acount];

			Acount++; 

			NewSeqs[B]+='-';
			NewQuals[B]+='!';
			NewRefs[B]+="-";
			NewCigars[B]+='R';
		}
		else if (reads[A].cigarString.c_str()[Acount] != 'D' and reads[B].cigarString.c_str()[Acount] == 'D')
		{
			currentPos.push_back(-1);
			currentPos.push_back(reads[B].Positions[Bcount]);
			
			currentChr.push_back("nope");
			currentChr.push_back(reads[B].ChrPositions[Bcount]);

			BLastGoodRef = reads[B].Positions[Bcount];
			BLastGoodChr = reads[B].ChrPositions[Bcount];

			NewSeqs[B]+=reads[B].seq.c_str()[Bcount];
			NewQuals[B]+=reads[B].qual.c_str()[Bcount];
			NewRefs[B]+=reads[B].RefSeq.c_str()[Bcount];
			NewCigars[B]+=reads[B].cigarString.c_str()[Bcount];
			
			Bcount++;

			NewSeqs[A]+='-';
			NewQuals[A]+='!';
			NewRefs[A]+='-';
			NewCigars[A]+='R'; 
		}
		else
		{
			if (reads[A].cigarString.c_str()[Acount] == 'H' or reads[A].cigarString.c_str()[Acount] == 'S')
			{
				currentPos.push_back(-1); 
				currentChr.push_back("nope");

				NewSeqs[A]+=reads[A].seq.c_str()[Acount];
				NewQuals[A]+=reads[A].qual.c_str()[Acount];
				NewRefs[A]+=reads[A].RefSeq.c_str()[Acount];
				NewCigars[A]+=reads[A].cigarString.c_str()[Acount];
				
				Acount++; 
		
			}
			else if (reads[A].cigarString.c_str()[Acount] == 'M' or reads[A].cigarString.c_str()[Acount] == 'X' or reads[A].cigarString.c_str()[Acount] == 'D') 
			{
				currentPos.push_back(reads[A].Positions[Acount]);
				currentChr.push_back(reads[A].ChrPositions[Acount]);
				ALastGoodRef = reads[A].Positions[Acount]; 
				ALastGoodChr = reads[A].ChrPositions[Acount]; 
				
				NewSeqs[A]+=reads[A].seq.c_str()[Acount];
				NewQuals[A]+=reads[A].qual.c_str()[Acount];
				NewRefs[A]+=reads[A].RefSeq.c_str()[Acount];
				NewCigars[A]+=reads[A].cigarString.c_str()[Acount];
				
				Acount++; 
			}
			else if (reads[A].cigarString.c_str()[Acount] == 'I')
			{
				currentPos.push_back(ALastGoodRef); 
				currentChr.push_back(ALastGoodChr);

				NewSeqs[A]+=reads[A].seq.c_str()[Acount];
				NewQuals[A]+=reads[A].qual.c_str()[Acount];
				NewRefs[A]+=reads[A].RefSeq.c_str()[Acount];
				NewCigars[A]+=reads[A].cigarString.c_str()[Acount];
				
				Acount++; 
			}
			else 
				cout << "WTF, cigar = " << reads[A].cigarString.c_str()[Acount]; 
			
			
			
			if (reads[B].cigarString.c_str()[Bcount] == 'H' or  reads[B].cigarString.c_str()[Bcount] == 'S')
			{
				currentPos.push_back(-1);
				currentChr.push_back("nope");
				
				NewSeqs[B]+=reads[B].seq.c_str()[Bcount];
				NewQuals[B]+=reads[B].qual.c_str()[Bcount];
				NewRefs[B]+=reads[B].RefSeq.c_str()[Bcount];
				NewCigars[B]+=reads[B].cigarString.c_str()[Bcount];
			
				Bcount++; 
			}
			else if (reads[B].cigarString.c_str()[Bcount] == 'M' or reads[B].cigarString.c_str()[Bcount] == 'X' or reads[B].cigarString.c_str()[Bcount] == 'D')
			{
				currentPos.push_back(reads[B].Positions[Bcount]);
				currentChr.push_back(reads[B].ChrPositions[Bcount]);
				BLastGoodRef = reads[B].Positions[Bcount];
				BLastGoodChr = reads[B].ChrPositions[Bcount]; 
			
				NewSeqs[B]+=reads[B].seq.c_str()[Bcount];
				NewQuals[B]+=reads[B].qual.c_str()[Bcount];
				NewRefs[B]+=reads[B].RefSeq.c_str()[Bcount];
				NewCigars[B]+=reads[B].cigarString.c_str()[Bcount];
				
				Bcount++; 
			}
			else if (reads[B].cigarString.c_str()[Bcount] == 'I')
			{
				currentPos.push_back(BLastGoodRef);
				currentChr.push_back(BLastGoodChr);
	
				NewSeqs[B]+=reads[B].seq.c_str()[Bcount];
				NewQuals[B]+=reads[B].qual.c_str()[Bcount];
				NewRefs[B]+=reads[B].RefSeq.c_str()[Bcount];
				NewCigars[B]+=reads[B].cigarString.c_str()[Bcount];
	
				Bcount++; 
			}
			else
				cout << "WTF, cigar = " << reads[B].cigarString.c_str()[Bcount];  
		
		}
		AlignmentPos.push_back(currentPos);
		AlignmentChr.push_back(currentChr);
	}

	//write  
	for (int i =0; i< reads.size(); i++)
	{
		reads[i].Positions.clear(); 
		reads[i].ChrPositions.clear(); 
		reads[i].seq = NewSeqs[i];
		reads[i].qual = NewQuals[i];
		reads[i].RefSeq = NewRefs[i]; 
		reads[i].cigarString = NewCigars[i];
		for (int j = 0; j<AlignmentPos.size(); j++)
		{reads[i].Positions.push_back(AlignmentPos[j][i]);}
		for (int j = 0; j<AlignmentPos.size(); j++)
		{reads[i].ChrPositions.push_back(AlignmentChr[j][i]);}
		
	}

	cout << "Post adjustment" << endl;
	reads[A].write(); 
	reads[B].write(); 

	bool deletion = true;
	string NewCigar = "";
	string NewSeq = "";
	string NewQual = "";
	string NewRef = "";
	vector<int>NewPos; 
	vector<string> NewChr; 
	
	char LastAlignedQ = ' ';
	int LastAlignedPos = -1; 
	string LastAlignedChr = "nope"; 


	//set LastAlignedPos to the first base with an aligned base
	bool notfound = true; 
	int base = 0; 
	while (notfound) 
	{
		for (int i = 0; i< reads.size(); i++) 
		{
			if (reads[i].Positions[base] > -1)
			{	
				LastAlignedPos = reads[i].Positions[base];
				LastAlignedChr = reads[i].ChrPositions[base];
				notfound = false;
				break; 
			}
		}
		if (notfound)
			base++;
	}
	

	cout << "************INTO**********" << endl;
	if (B == 1 and GetReadOrientation(reads[A].flag) == GetReadOrientation(reads[B].flag)) //if they are on the same strand, and there are only two split reads  
	{
		for (int i =base; i < reads[A].seq.size(); i++)
		{
			if (reads[A].Positions[i] > -1) //if this base is aligned in A 
			{
				if (reads[A].Positions[i] - LastAlignedPos > 1) //indicates a deletion
				{
					if (LastAlignedQ == '!' or reads[A].qual.c_str()[i] == '!' )
					{
						cout << "well fuck this shit A" << endl;
						return reads[A]; 
					}
					//if(reads[A].ChrPositions[i] == LastAlignedChr and abs(reads[A].Positions[i] -LastAlignedPos ) < MaxVarentSize ) //must be on the same chromosome
					if(reads[A].chr == reads[B].chr and abs(reads[A].Positions[i] -LastAlignedPos ) < MaxVarentSize ) //must be on the same chromosome
					{
						cout << "wel this dosnt make any sense" << endl; //reads are in order in the bam so A should always be downstream of B, theus the deletion shoould be detected in B
						BEDBigStuff << reads[A].chr << "\t" << LastAlignedPos << "\t" << reads[A].Positions[i] << "\t" << "Deletion" << endl;
						for (int j = LastAlignedPos; j<reads[A].Positions[i]-1; j++)
						{		
							NewCigar += 'D';
							NewSeq += '-';
							NewQual += LastAlignedQ;
							NewRef+= toupper(Reff.getSubSequence(reads[A].chr, i, 1).c_str()[0]);
							NewPos.push_back(j);
							NewChr.push_back(reads[A].ChrPositions[i]);
						}
					}
					else
					{
						//if( reads[A].ChrPositions[i] == LastAlignedChr and abs(reads[A].Positions[i] -LastAlignedPos ) >= MaxVarentSize )
						if( reads[A].chr == reads[B].chr and abs(reads[A].Positions[i] -LastAlignedPos ) >= MaxVarentSize )
						{
							Translocations << "too big " << abs(reads[A].Positions[i] -LastAlignedPos ) << endl;
							reads[A].writetofile(Translocations);
							reads[B].writetofile(Translocations);
							Translocations << endl << endl;	
							
						}
						else
						{
							if ((reads[A].chr == "hs37d5" and reads[B].chr != "hs37d5" ) or  (reads[A].chr != "hs37d5" and reads[B].chr == "hs37d5" ))
							{
								Translocations << "mobil elemnt " << endl;
								reads[A].writetofile(Translocations); 
								reads[B].writetofile(Translocations); 
								Translocations << endl << endl;
							}
							else
							{
								Translocations << "we got a translocation" << endl;
								reads[A].writetofile(Translocations);
								reads[B].writetofile(Translocations);
								Translocations << endl << endl;

							}	
						}
					}
				}
				else if (reads[A].Positions[i] -LastAlignedPos < 0 and abs(reads[A].Positions[i] -LastAlignedPos ) < MaxVarentSize) // indicates a possible insertion or tandem duplication
				{
					 if (LastAlignedQ == '!'  or reads[A].qual.c_str()[i] == '!')
					{
						cout << "well fuck this shit B" << endl;
						return reads[A];
					}
					cout << "this could be one A, last = " << LastAlignedPos << " Current = " << reads[A].Positions[i] << " chr = " << LastAlignedChr << " and " << reads[A].ChrPositions[i] << endl;
//					if (reads[A].ChrPositions[i] == LastAlignedChr )
					if (reads[A].chr == reads[B].chr )
					{
						cout << "This is an insertion A at base " << i <<  endl;
						BEDBigStuff << reads[A].chr << "\t" << reads[A].Positions[i] << "\t" << LastAlignedPos << "\tTandemDup" << endl;
						cout << "tadem dup" << endl;
						int j = 0; 
						for( j = i; j<reads[A].seq.size() and reads[A].Positions[j] < LastAlignedPos ; j++)
						{
							NewCigar += 'T'; //'I';
							NewSeq += reads[A].seq.c_str()[j];
							NewQual += reads[A].qual.c_str()[j];
							NewRef+= '-';
							NewPos.push_back(reads[A].Positions[i]);
							NewChr.push_back(reads[A].ChrPositions[i]);
						}
						i=j;
						//find the last base that was aligned, htere can be novel insertion stuff so you cant jsut take the last base
						int k; 
						for (k =reads[A].Positions.size()-1; k >=0; k+= -1)
						{
							if (reads[A].Positions[k]+1 > 1)
								break;
						}
						//cout << "j= " << reads[A].Positions[k]+1 << " < " <<  LastAlignedPos << " - " << endl;
						for( j = reads[A].Positions[k]+1; j < LastAlignedPos; j++)
						{
							NewCigar += 'T'; //'I';
							NewSeq += toupper(Reff.getSubSequence(reads[A].chr, j, 1).c_str()[0]);
							NewQual += '!';
							NewRef+= '-';
							NewPos.push_back(j);
							NewChr.push_back(reads[A].chr);
						}
						cout << "yaya finished" << endl;
						 
					}
					else
					{
						if ((reads[A].chr == "hs37d5" and reads[B].chr != "hs37d5" ) or  (reads[A].chr != "hs37d5" and reads[B].chr == "hs37d5" ))
						{
							Translocations << "mobil elemnt " << endl;
							reads[A].writetofile(Translocations);
							reads[B].writetofile(Translocations);
							Translocations << endl << endl;
						}
						Translocations << "we got a translocation" << endl;
						reads[A].writetofile(Translocations);
						reads[B].writetofile(Translocations);
						Translocations << endl << endl;
					}
	
				}
				else if( reads[A].Positions[i] -LastAlignedPos < 0  and abs(reads[A].Positions[i] -LastAlignedPos ) >= MaxVarentSize )
				{
					 if (LastAlignedQ == '!'  or reads[A].qual.c_str()[i] == '!')
					{
						cout << "well fuck this shit C" << endl;

						return reads[A];
					}
					Translocations << "too big " << abs(reads[A].Positions[i] -LastAlignedPos ) << endl;
					reads[A].writetofile(Translocations);
					reads[B].writetofile(Translocations);
					Translocations << endl << endl;
				}

				
				if (i<reads[A].cigarString.size())
				{
					NewCigar += reads[A].cigarString.c_str()[i];
					NewSeq += reads[A].seq.c_str()[i];
					NewQual += reads[A].qual.c_str()[i];
					NewRef+= reads[A].RefSeq.c_str()[i];
					NewPos.push_back(reads[A].Positions[i]);
					NewChr.push_back(reads[A].ChrPositions[i]);
		       	 		LastAlignedQ = reads[A].qual.c_str()[i];
					LastAlignedPos = reads[A].Positions[i];
					LastAlignedChr = reads[A].ChrPositions[i];
					//cout <<	 "A " << reads[A].seq.c_str()[i] << " " << reads[A].cigarString.c_str()[i] << " " << reads[A].RefSeq.c_str()[i] <<" " << reads[A].ChrPositions[i] << " " << reads[A].Positions[i] << NewSeq << endl;
				}	
			}
			else if (reads[B].Positions[i] > -1)
			{
				if (reads[B].Positions[i] - LastAlignedPos > 1)
				{
					 if (LastAlignedQ == '!'  or reads[B].qual.c_str()[i] == '!')
					{
						
						cout << "well fuck this shit D" << endl;

						return reads[A];
					}
					//if(reads[B].ChrPositions[i] == LastAlignedChr  and abs(reads[B].Positions[i] - LastAlignedPos ) < MaxVarentSize )
					if(reads[B].chr == reads[A].chr  and abs(reads[B].Positions[i] - LastAlignedPos ) < MaxVarentSize )
					{
						cout << "striahgtup deletion, size = " <<  abs(reads[B].Positions[i] -LastAlignedPos ) << " at base " << i << " from Position " << LastAlignedPos << " to " << reads[B].Positions[i] << endl; 
						BEDBigStuff << reads[B].chr << "\t" << LastAlignedPos << "\t" << reads[B].Positions[i]  << "\t" << "Deletion" << endl;
//						cout << "Inserting reff sequence from " << LastAlignedPos+1 << " to " << reads[B].Positions[i] << " = " << reads[B].Positions[i]-LastAlignedPos << endl;
						for (int j = LastAlignedPos; j< reads[B].Positions[i]-1; j++)
			       	 		{
							//cout << j << " - " << j - LastAlignedPos<< endl;
			       				NewCigar += 'D';
			       				NewSeq += '-';
			       		 		NewQual += LastAlignedQ;
			       		 		NewRef+= toupper(Reff.getSubSequence(reads[B].chr, j, 1).c_str()[0]);
							NewPos.push_back(j);
							NewChr.push_back(reads[B].ChrPositions[i]);

	//						char tmp = toupper(Reff.getSubSequence(reads[B].chr, j, 1).c_str()[0]); 
	  //		      			cout << "R " << '-' << " " << 'D' << " " << tmp << " " << reads[B].chr << " " << j << " " << NewSeq << endl;
	//						cout << "R " << '-' << " " << 'D' << " " << tmp << " " << reads[B].chr << " " << j << " " << NewRef << endl << endl;
						}
					}
					else
					{
						// if( reads[A].ChrPositions[i] == LastAlignedChr and abs(reads[B].Positions[i] -LastAlignedPos ) >= MaxVarentSize )
						if( reads[A].chr == reads[B].chr and abs(reads[B].Positions[i] -LastAlignedPos ) >= MaxVarentSize )
						{
							Translocations << "too big " << abs(reads[A].Positions[i] -LastAlignedPos ) << endl;
							reads[A].writetofile(Translocations);
							reads[B].writetofile(Translocations);
							Translocations << endl << endl;
						}
						else
						{
							if ((reads[A].chr == "hs37d5" and reads[B].chr != "hs37d5" ) or  (reads[A].chr != "hs37d5" and reads[B].chr == "hs37d5" ))
							{
								Translocations << "mobil elemnt " << endl;
								reads[A].writetofile(Translocations);
								reads[B].writetofile(Translocations);
								Translocations << endl << endl;
							}
							else	
							{
								Translocations << "we got a translocation" << endl;
								reads[A].writetofile(Translocations);
								reads[B].writetofile(Translocations);
								Translocations << endl << endl;
							}
						}
					}

				}
				else if (reads[B].Positions[i] -LastAlignedPos < 0 and abs(reads[B].Positions[i] - LastAlignedPos ) < MaxVarentSize) // indicates a possible insertion or tandem duplication
				{
					 if (LastAlignedQ == '!' or reads[B].qual.c_str()[i] == '!' )
					{
						cout << "well fuck this shit E" << endl;

						return reads[A];
					}
					 // cout << "this could be one B, last = " << LastAlignedPos << " Current = " << reads[B].Positions[i] << " chr = " << LastAlignedChr << " and " << reads[B].ChrPositions[i] << endl;
					//if (reads[B].ChrPositions[i] == LastAlignedChr  
					if (reads[B].chr == reads[A].chr  )
					{
						cout << "This is an insertion B at base " << i << endl;

						BEDBigStuff << reads[B].chr << "\t" << reads[B].Positions[i] << "\t" << LastAlignedPos << "\tTandemDup" << endl;
					//	cout << "tadem dup" << endl;
						int j = 0; 
						for(j = i; j < reads[B].seq.size() and reads[B].Positions[j] <= LastAlignedPos; j++)
						{
							NewCigar += 'T';
							NewSeq += reads[B].seq.c_str()[j];
							NewQual += reads[B].qual.c_str()[j];
							NewRef+= '-';
							NewPos.push_back(reads[B].Positions[i]);
							NewChr.push_back(reads[B].ChrPositions[i]);
						}
						i=j;
						//need to find last base that was alinged, insetion can mess this up so you can just take the last base pos 
						int k;
						for (k =reads[B].Positions.size()-1; k >=0; k+= -1)
						{
							if (reads[B].Positions[k]+1 > 1)
								break;
						}
						//cout << "j= " << reads[B].Positions[k]+1 << " < " <<  LastAlignedPos << " - " << endl;
						for( j = reads[B].Positions[k]+1; j < LastAlignedPos; j++) 
						{
							NewCigar += 'T';
							NewSeq += toupper(Reff.getSubSequence(reads[B].chr, j, 1).c_str()[0]);
							NewQual += '!';
							NewRef+= '-';
							NewPos.push_back(j);
							NewChr.push_back(reads[B].chr);
						}
						
					}
					else
					{
						cout << "we got a translocation" << endl;
						if ((reads[A].chr == "hs37d5" and reads[B].chr != "hs37d5" ) or  (reads[A].chr != "hs37d5" and reads[B].chr == "hs37d5" ))
						{
							cout << "mobil elemnt " << endl;
							//reads[A].write();
							//reads[B].write();
						}
					}

				}
				else if( reads[B].Positions[i] -LastAlignedPos < 0  and abs(reads[B].Positions[i] -LastAlignedPos ) >= MaxVarentSize )
				{
					 if (LastAlignedQ == '!' or reads[B].qual.c_str()[i] == '!')
					{
						cout << "well fuck this shit F" << endl;
						return reads[A];
					}
					Translocations << "too big " << abs(reads[A].Positions[i] -LastAlignedPos ) << endl;
					reads[A].writetofile(Translocations);
					reads[B].writetofile(Translocations);
					Translocations << endl << endl;
				}


				if(i<reads[B].cigarString.size())
				{
					NewCigar += 'M'; //reads[B].cigarString.c_str()[i];
					NewSeq += reads[B].seq.c_str()[i];
					NewQual	+= reads[B].qual.c_str()[i];
					NewRef+= reads[B].RefSeq.c_str()[i];
					NewPos.push_back(reads[B].Positions[i]);
					NewChr.push_back(reads[B].ChrPositions[i]);

					LastAlignedQ = reads[B].qual.c_str()[i];
					LastAlignedPos = reads[B].Positions[i];
					LastAlignedChr = reads[B].ChrPositions[i];
					//	 cout << "B " << reads[B].seq.c_str()[i] << " " << reads[B].cigarString.c_str()[i] << " " << reads[B].RefSeq.c_str()[i] <<" " << reads[B].ChrPositions[i] << " " << reads[B].Positions[i] << NewSeq << endl;
				}
			}
			else
			{
				if (reads[A].cigarString.c_str()[i] == 'S')
				{
					NewCigar += reads[A].cigarString.c_str()[i];
					NewSeq +=reads[A].seq.c_str()[i];
					NewQual += reads[A].qual.c_str()[i];
					NewRef+= reads[A].RefSeq.c_str()[i];
					NewPos.push_back(reads[A].Positions[i]);
					NewChr.push_back(reads[A].ChrPositions[i]);

				}
				else if (reads[B].cigarString.c_str()[i] == 'S')
				{
					NewCigar += reads[B].cigarString.c_str()[i];
					NewSeq +=reads[B].seq.c_str()[i];
					NewQual += reads[B].qual.c_str()[i];
					NewRef+= reads[B].RefSeq.c_str()[i];
					NewPos.push_back(reads[B].Positions[i]);
					NewChr.push_back(reads[B].ChrPositions[i]);

				}
				else
					cout << "no base works :(" << endl;
				
			}
		}

		//fix inernal S bases
		int First = -1; 
		int Last = -1; 
		string NewNewCigar = ""; 
		for (int i =0; i < NewCigar.size(); i++)
		{
			if(NewCigar.c_str()[i] != 'S' and NewCigar.c_str()[i] != 'H')
			{
				First = i; 
				break; 
			}
		}
		for (int i = NewCigar.size() -1; i>=0; i--)
		{
			if(NewCigar.c_str()[i] != 'S' and NewCigar.c_str()[i] != 'H')
			{
				Last = i; 
				break; 
			}
		}

		for (int i = 0; i<NewCigar.size(); i++)
		{
			if (i > First and i < Last)
			{
				if (NewCigar.c_str()[i] == 'S' or NewCigar.c_str()[i] == 'H')
					NewNewCigar+= 'I';
				else 
					NewNewCigar+=NewCigar.c_str()[i]; 
			}
			else 
				NewNewCigar+=NewCigar.c_str()[i];
					
		}
		NewCigar = NewNewCigar; 
		//cout << NewSeq << endl << NewQual << endl << NewCigar << endl << NewRef << endl;
		reads[A].first = true;
//		cout <<"redoing alignemnt number" << endl;
		int temp = reads[A].alignments[0];
		reads[A].alignments.clear();
		reads[A].alignments.push_back(temp);
//		cout << "alignemtns should be " << reads[A].alignments.size() << endl;

       		reads[A].cigarString = NewCigar;
		reads[A].seq = NewSeq;
		reads[A].qual = NewQual;
		reads[A].RefSeq = NewRef;
		reads[A].Positions.clear();
		reads[A].ChrPositions.clear();
		reads[A].Positions = NewPos; 
		reads[A].ChrPositions = NewChr; 
		reads[A].combined = true;
		reads[B].combined = true;
		
		//reads[A].write(); 		
		//{
		//reads[A].write(); 
	       	//reads[A].writeVertical();
		//}		
	


		cout <<	 "Done with " << reads[A].name << endl;
	}
	else
	{
		//this needs to be the invertion stuff 
		if( reads[A].chr == reads[B].chr )//and abs(reads[A].Positions[i] -LastAlignedPos ) <= MaxVarentSize )
		{
			Translocations << "INVERSION "  << endl;
			reads[A].writetofile(Translocations);
			reads[B].writetofile(Translocations);
			Translocations << endl << endl;
			string Acig = "";
			string Bcig = "";	
			for (int i = 0; i < reads[A].seq.size(); i++)
			{
				char Ab = reads[A].cigarString.c_str()[i]; 
				char Bb = reads[B].cigarString.c_str()[i]; 		
				if ((reads[A].cigarString.c_str()[i] == 'M' or reads[A].cigarString.c_str()[i] == 'X') and (reads[B].cigarString.c_str()[i] == 'S' or reads[B].cigarString.c_str()[i] == 'H'))
					Bb = 'U';
				if ((reads[B].cigarString.c_str()[i] == 'M' or reads[B].cigarString.c_str()[i] == 'X') and (reads[A].cigarString.c_str()[i] == 'S' or reads[A].cigarString.c_str()[i] == 'H'))
                              		Ab = 'U'; 
				Acig += Ab; 
				Bcig += Bb; 
				   	
			}
			reads[A].cigarString = Acig; 
			reads[B].cigarString = Bcig; 
			cout << "invertion adjust string"; 
			reads[A].write(); 
			reads[B].write(); 	
		}
		else if ((reads[A].chr == "hs37d5" and reads[B].chr != "hs37d5" ) or  (reads[A].chr != "hs37d5" and reads[B].chr == "hs37d5" ))
		{
			Translocations << "mobil elemnt " << endl;
			reads[A].writetofile(Translocations);
			reads[B].writetofile(Translocations);
			Translocations << endl << endl;
		}
		else
		{
			Translocations << "we got a translocation and invertion" << endl;
			reads[A].writetofile(Translocations);
			reads[B].writetofile(Translocations);
			Translocations << endl << endl;
			
		}

			
		cout << "*********SKIPPING************\ndifference strands" << endl;
		BEDNotHandled << "Different strands" << endl;
		BEDNotHandled << reads[A].chr << "\t" << reads[A].pos << "\t" << reads[A].pos+reads[A].seq.size() << "\t" << reads[A].name << "\t" << reads[A].cigar << endl;
		reads[A].writetofile(BEDNotHandled);
		BEDNotHandled << reads[B].chr << "\t" << reads[B].pos << "\t" << reads[B].pos+reads[B].seq.size() << "\t" << reads[B].name << "\t" << reads[B].cigar << endl;	
		reads[B].writetofile(BEDNotHandled);
		
		BEDNotHandled << endl << endl;


		Invertions << reads[A].chr << "\t" << reads[A].pos << "\t" << reads[B].pos << "\t" << reads[B].pos- reads[A].pos << endl;

	}
	reads[A].write(); 
	return reads[A];
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
  -sam  arg		Path to input SAM file, omit for stdin\n\
  -r  arg		Path to reference file \n\
  -hf arg		Path to HashFile from RUFUS.build\n\
  -o  arg		Output stub\n\
  -m  arg		Maximum varient size: default 1Mb\n\
			(Sorry it has to be a num, no 1kb, must be 1000\n\
  -c  arg		Path to sorted.tab file for the parent sample\n\
  -s  arg 		Path to sorted.tab file for the subject sample\n\
  -mQ arg		Minimum map quality to consider varients in\n\
";
	
	string MutHashFilePath = "" ;
	MaxVarentSize = 1000000;
	string RefFile = ""; 
	string HashListFile = "" ; 	
	string samFile = "stdin"; 
	string outStub= "";
	int MinMapQual = 0; 
	for(int i = 1; i< argc; i++)
	{
		cout << i << " = " << argv[i] << endl; 
	}
	cout <<"****************************************************************************************" << endl;
	vector <int> ParentHashFilePaths; 
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
		else if (p == "-r")
		{
			RefFile = argv[i+1];
			i=i+1;
			 cout << "YAAAY added RefFile = " << RefFile << endl;	
		}
		else if (p == "-sam")
		{
			samFile =  argv[i+1];
			i++;
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
		else if (p == "-m")
		{
			MaxVarentSize =  atoi(argv[i+1]);
			i++;
			cout << "YAAAY added MaxVarSize = " << MaxVarentSize << endl;
		}
		else if (p == "-c")
		{
			cout << "Par Hash = " << argv[i+1] << endl;
			ParentHashFilePaths.push_back(i+1);
			i=i+1;
		}
		else if (p == "-s")
		{
			cout << "Sub Hash = " << argv[i+1] << endl;
			MutHashFilePath = argv[i+1];
			i+=1;
		}
		else if (p == "-mQ")
                {
                        cout << "Min Mapping Qualtiy = " << argv[i+1] << endl;
                        MinMapQual = atoi(argv[i+1]);
                        i+=1;
                }
		else
		{
			cout << "ERROR: unkown command line paramater -" <<  argv[i] << "-"<< endl;
			return 0; 
		}
		
	}
	//check values 
	if (RefFile == "")
	{
		cout << "ERROR Reference required" << endl;
		return 0; 
	}
	if (HashListFile == "")
	{
		cout << "Error HashList required" << endl;
		return 0; 
	}
	if (outStub == "")
	{
		if (samFile != "stdin")
			outStub = samFile; 
		else
		{
			cout << "ERROR out file stub required " << endl;
			return -1; 
		}
	}

	for (int i = 0; i < ParentHashFilePaths.size(); i++)
	{
		ifstream reader; 
		reader.open (argv[ParentHashFilePaths[i]]); 
		string line = "";
		unordered_map <unsigned long int, int> hl;  
		while (getline(reader, line))
		{
			vector <string> temp = Split(line, ' '); 
			unsigned long hash = HashToLong(temp[0]); 
			hl[hash] = atoi(temp[1].c_str());
			hash = HashToLong(RevComp(temp[0])); 
			hl[hash] = atoi(temp[1].c_str()); 
		}
		ParentHashes.push_back(hl); 
		reader.close(); 
	}
	cout << "check parent thing" << endl; 
	for(int i =0; i < ParentHashes.size(); i++)
	{
		cout << "sample " << i << endl;
	}
	
	ifstream reader;
	reader.open (MutHashFilePath);
	string line = "";
	while (getline(reader, line))
	{
		
		vector <string> temp = Split(line, ' ');
		unsigned long hash = HashToLong(temp[0]);
		MutantHashes[hash] = atoi(temp[1].c_str());
		hash = HashToLong(RevComp(temp[0]));
		MutantHashes[hash] = atoi(temp[1].c_str());
	}
	reader.close(); 
			
	
	//***********************************************
	//cout << "Call is Reference Contigs.fa OutStub HashList MaxVarientSize" << endl;
	double vm, rss, MAXvm, MAXrss;
	MAXvm = 0;
	MAXrss = 0;
	process_mem_usage(vm, rss, MAXvm, MAXrss);
   	cout << "VM: " << vm << "; RSS: " << rss << endl;	

	
	int BufferSize = 1000;

	Reff.open(RefFile);


      	ifstream HashList;
      	HashList.open (HashListFile);
      	if ( HashList.is_open())
      	{	 cout << "HashList Open " << HashListFile << endl;}   //cout << "##File Opend\n";
      	else
      	{
      		cout << "Error, HashList could not be opened";
      		return 0;
      	}
	line = "";
	getline(HashList, line);
	cout << "line = " << line << endl; 
	vector<string> temp = Split(line, '\t');
	cout << "split = " << temp[0] << " and " << temp[1] << endl;
	if (temp.size() ==4)
	{ 
		HashSize = temp[3].length(); 
		Hash.insert(pair<string, int>(temp[3], atoi(temp[2].c_str())));
	
		while ( getline(HashList, line))
		{
			vector<string> temp = Split(line, '\t');
			Hash.insert(pair<string, int>(temp[3], atoi(temp[2].c_str())));
			//cout << "added pair " << temp[3] << "\t" << temp[2] << endl;
		}
		HashList.close(); 
		cout << "done with HashList" << endl;
	}
	else if (temp.size() ==2)
	{
		HashSize = temp[0].length();
		Hash.insert(pair<string, int>(temp[0], atoi(temp[1].c_str())));

		while ( getline(HashList, line))
		{
			vector<string> temp = Split(line, '\t');
			Hash.insert(pair<string, int>(temp[0], atoi(temp[1].c_str())));
			//cout << "added pair " << temp[3] << "\t" << temp[2] << endl;
		}
		HashList.close();
		cout << "done with HashList" << endl;
	}
	else if (temp.size() ==1)
	{
		vector<string> temp = Split(line, ' ');
		HashSize = temp[0].length();
		Hash.insert(pair<string, int>(temp[0], atoi(temp[1].c_str())));

		while ( getline(HashList, line))
		{
			vector<string> temp = Split(line, ' ');
			Hash.insert(pair<string, int>(temp[0], atoi(temp[1].c_str())));
			//cout << "added pair " << temp[3] << "\t" << temp[2] << endl;
		}
		HashList.close();
		cout << "done with HashList" << endl;
	}
	//map<string, int>::iterator it;
	//for ( it = Hash.begin(); it != Hash.end(); it++ )
	//{
	//	cout << "-"<<it->first<<"-" << "\t" << it->second << endl;
	//}

	ifstream SamFile;
       	if (samFile == "stdin")
	{
		cout << "Sam File is STDIN" << endl;
		SamFile.open ("/dev/stdin");
	}
	else
	{
		cout << "Sam File is " << samFile << endl;
		SamFile.open (samFile);
	}
	if ( SamFile.is_open())
	{      cout << "Sam File Opend\n";}
	else
	{
		cout << "Error, SamFile could not be opened";
		return 0;
	}

	string boom = outStub;
	VCFOutFile.open(boom+ ".vcf"); 
	BEDOutFile.open(boom+ ".vcf.bed"); 
	BEDBigStuff.open(boom+ ".vcf.Big.bed");
	BEDNotHandled.open(boom+ ".vcf.NotHandled.bed");
	Invertions.open(boom+".vcf.invertions.bed");
	Translocations.open(boom+ ".vcf.Translocations");
	Unaligned.open(boom+"vcf.Unaligned");

	//write VCF header
	VCFOutFile << "##fileformat=VCFv4.1" << endl;
	VCFOutFile << "##fileDate=" << time(0) << endl;
	VCFOutFile << "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">" << endl;
	VCFOutFile << "##FORMAT=<ID=AK,Number=1,Type=Integer,Description=\"Alternte Kmer Count\">" << endl;
	VCFOutFile << "##INFO=<ID=AO,Number=A,Type=Integer,Description=\"Alternate allele observations, with partial observations recorded fractionally\">"<<endl;
	VCFOutFile << "##INFO=<ID=HD,Number=A,Type=String,Description=\"Hash counts for each k-mer overlapping the vareint, -1 indicates no info\">"<< endl;
	VCFOutFile << "##INFO=<ID=RN,Number=1,Tinype=String,Description=\"Name of contig that produced the call\">"<< endl;
	VCFOutFile << "##INFO=<ID=MQ,Number=1,Tinype=Integer,Description=\"Mapping quality of the contig that created the call\">"<< endl;
	VCFOutFile << "##INFO=<ID=cigar,Number=1,Tinype=String,Description=\"Cigar string for the contig that created the call\">"<< endl;
	VCFOutFile << "##INFO=<ID=CVT,Number=1,Tinype=String,Description=\"Compressed Varient Type\">"<< endl;
	VCFOutFile << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t";
	VCFOutFile << outStub << endl;

	int lines = 0;

	line="";
 
	unsigned long LongHash;

	cout << "Reading in Sam File"  << endl;
	map <string, int> Names; 
	vector <SamRead> reads; 
	int counter = 0; 
	while (getline(SamFile, line))
	{
		//cout << line << endl;
		counter ++; 
		SamRead read; 
		read.parse(line);
		read.getRefSeq();
		//if (read.mapQual > 0)
		{
			reads.push_back(read);
			if (counter%100 == 0)
				cout << "read " << counter << " entries " << char(13); 
		}
	}
	cout << endl;
	cout << "Read in " << reads.size() << " reads " << endl;
	
	cout << "procesing split reads" << endl;
	for (int i = 0; i < reads.size(); i++) 
	{
		cout << "processing read " << reads[i].name << endl;
		if (reads[i].alignments.size() == 0)
		{
			reads[i].alignments.push_back(i);
			for (int j = i+1; j < reads.size(); j++)
			{
				if (strcmp(reads[i].name.c_str(), reads[j].name.c_str()) == 0)
				{
					cout << "found mate " << reads[j].name << endl;
				//	if (reads[i].chr == reads[j].chr)
					{
						reads[i].alignments.push_back(j); 
						reads[j].first = false;
					}
				//	else
				//		cout << "skip that shit" << endl; 
				}
			}
		}
		for (int j = 0; j<reads[i].alignments.size(); j++)
			reads[reads[i].alignments[j]].alignments = reads[i].alignments; 

	}

	for (int i = 0; i < reads.size(); i++)
	{
		SamRead read = reads[i]; 
		cout << read.name  << endl;
		if (strcmp(read.chr.c_str(), "*") == 0)
		{
			read.writetofile(Unaligned);
			Unaligned << endl;
		}
		//if  (read.first)
		{	
			//if it looks like a simple split read alignment colaps it into a single read
			if (read.first and read.alignments.size() == 2 )
			{
				cout << "atempting colaps" << endl;
				cout << read.name << endl;
				vector<SamRead> R; 
				for(int j =0; j< read.alignments.size(); j++)
				{
					//these better be sorted by position
					//if (reads[read.alignments[j]].chr == read.chr)
					{
						R.push_back(reads[read.alignments[j]]);
						cout << reads[read.alignments[j]].name << endl;
					}
				}
				if (R.size() == 2 & /*R[0].chr == R[1].chr & */ R[0].mapQual > 0 & R[1].mapQual > 0)
				{
					read = BetterWay(R); 
				}
			}
			else if(read.first and read.alignments.size() >2)
			{
				BEDNotHandled << "too many alignments" << endl;
				BEDNotHandled << read.chr << "\t" << read.pos << "\t" << read.pos+read.seq.size() << "\t" << read.name << "\t" << read.cigar << endl;
				cout << "too many alignments" << endl;
				cout << read.chr << "\t" << read.pos << "\t" << read.pos+read.seq.size() << "\t" << read.name << "\t" << read.cigar << endl;
				for (int j = 0; j< read.alignments.size(); j++)
				{
					SamRead mate = reads[read.alignments[j]]; 
					cout << j << "\t" << mate.chr << "\t" << mate.pos << "\t" << mate.pos+mate.seq.size() << "\t" << mate.name << "\t" << mate.cigar << endl;
					 BEDNotHandled << mate.chr << "\t" << mate.pos << "\t" << mate.pos+mate.seq.size() << "\t" << mate.name << "\t" << mate.cigar << endl;
					mate.writetofile(BEDNotHandled);
				}
				BEDNotHandled << endl << endl;
			}
			//read.write(); 	
			//if (read.first && read.alignments.size() ==1 )
			// 	read.parseMutations();
			//else if (read.combined == false )
			if (read.mapQual > MinMapQual and read.alignments.size() <=2)
			{
				read.parseMutations(argv); 
				//BEDNotHandled << read.chr << "\t" << read.pos << "\t" << read.pos+read.seq.length() << "\t" << read.name << endl;
			}
		}
	}
	VCFOutFile.close(); 
	BEDOutFile.close();
	BEDBigStuff.close();
	BEDNotHandled.close();
	Invertions.close(); 
	cout << "\nreally dont\n";
}
	
