
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

vector <string> ParNames; 
vector <unordered_map <unsigned long int, int> > ParentHashes;
unordered_map  <unsigned long int, int> MutantHashes; 
unordered_map <unsigned long int, int> ExcludeHashes;
FastaReference Reff;
int CurrentSVeventID = 0; 
int MaxBND = 0; 
int HashSize = 25; 
int totalDeleted; 
int totalAdded;
int MaxVarentSize = 1000; 
int ParLowCovThreshold = 7; 
int SegThreshold = 10;
ofstream VCFOutFile;
ofstream BEDOutFile;
ofstream BEDBigStuff;
ofstream BEDNotHandled;
ofstream Invertions; 
ofstream Translocations; 
ofstream Translocationsbed; 
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
			//cout << "ERROR, invalid character - " << hash.c_str()[i] << endl;
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
	bool b[16];
	
       	for (int j = 0;  j < 16;  ++j){
		b [j] =  0 != (flag & (1 << j));	
	}
	//cout << "bits are " << b[0] << " " << b[1] << " " << b[2] << " " << b[3] << " " << b[4] << " " << b[5] << " " << b[6] << " " << b[7] << " " << b[8] << " " << b[9] << " " << b[10] << " " << b[11] << " " << b[12] << endl;
	//cout << "flag is " << flag << endl; 
	//cout << "orientation is " << b[4] << endl; 
	return b[4];
	/*if (flag == 0 or flag == 2048){
		cout << "0 forward" << endl;
		return 0;}
	else if (flag == 16 or flag == 2064){
		cout << "1 reverse" << endl; 
		return 1;}
	else
		cout <<  "error flag not recognised " << flag << endl;

	return -1;*/
}
string getHash(string seq, int j, int HashSize)
{
	int bases = 0;
	string NewHash = "";
	if (j <seq.size())
	{
		if (seq.c_str()[j] == 'A' or seq.c_str()[j] == 'C' or seq.c_str()[j] == 'G' or seq.c_str()[j] == 'T') //if the first base is not in a real seqeunce, return a blank hash
		{
			while (bases<HashSize and j<seq.size())
			{
				if (seq.c_str()[j] == 'A' or seq.c_str()[j] == 'C' or seq.c_str()[j] == 'G' or seq.c_str()[j] == 'T')
				{
					NewHash+=seq.c_str()[j];
					bases++;
				}
				j++;
			}
			//cout << "grabbed hash " << NewHash << endl;
		}
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
		else if (C == '-')
			NewString += "-"; 
		else
		{
			cout << "ERROR IN RevComp - " << C << " " << endl;
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
string ShittyGenotyper(int Alt, int Ref)
{
	if (Alt ==0 and Ref ==0)
		return ".";
	else if (Alt == 0 and Ref > 1)
		return "0/0";
	else if (Alt >0 and Ref ==0)
		return "1/1";
	else if ((double) Alt / ((double) Ref + (double) Alt)  >.85)
		return "1/1";
	else if ((double) Alt / ((double) Ref + (double) Alt)  <.15)
		return "0/0";
	else
		return "0/1";
}
class MobRead
{
	public:
	string name; 
	string cigarString; 
	int flag; 
	bool FlagBits[16]; 
	string chr;
	int pos;
	string cigar;
	int mapQual;
	string seq;
	string qual;
	void parse(string read);
	int AS = -1; 
	void processCigar(); 
	void write(); 
};
void MobRead::write() 
{
	cout << name << endl;
	cout << "   flag = " << flag << endl;
	cout << "   mapQual = " << mapQual << endl;
	cout << "   Strand = " << GetReadOrientation(flag) << endl;
	cout << "   chr " << chr << " - " << pos << " qual = " << mapQual <<  endl;
	cout << "   cigar  = " << cigar << endl;
	cout << "   Seq    = " << seq << endl;
	cout << "   Qual   = " << qual << endl;
	cout << "   Cigar  = " << cigarString << endl;

}
	

void MobRead::processCigar()
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
void MobRead::parse(string read)
{
	//cout << "Mob parsing " << read << endl;
	vector <string> temp = Split(read, '\t');
	name = temp[0];
	//cout << "name " << name;
	flag = atoi(temp[1].c_str());
	chr = temp[2];
	pos = atoi(temp[3].c_str());
	mapQual = atoi(temp[4].c_str());
	cigar = temp[5];
	seq = temp[9];
	qual = temp[10];
	
	for (int i = 11; i < temp.size(); i++)
	{
		if (temp[i].length() > 3)
		{
			if (temp[i].c_str()[0]== 'A' && temp[i].c_str()[1]== 'S' && temp[i].c_str()[2]== ':')
			{
				vector <string> temp2 = Split(temp[i], ':'); 
				AS = atoi(temp2[2].c_str()); 
			}
		}
	}

	for (int j = 0;  j < 16;  ++j){
		FlagBits [j] =  0 != (flag & (1 << j));
	}
	processCigar(); 
}
class SamRead
{
	public:
	
	bool parsed = false; 
	bool MobAligned = false;
	bool AllA = false; 
	int SVeventid = 0; 
	int BNDid =0; 
	string clipPattern = "";
	int isSplitRead = 0; 
	string MobContig = "none"; 
	int PolyA = 0; 
	string name;
	int flag;
	bool FlagBits[16];
	string chr;
	int pos;
	int mapQual;
	int AlignScore; 
	int MobAS = -1; 
	string cigar;
	string seq;
	string qual;
	string RefSeq; 
	string originalSeq;
	string originalQual; 
	string cigarString;
	string strand;  
	string phase = "none"; 
	float StrandBias; 
	string strands; 
	int forward; 
	int reverse; 
	bool UsedForBigVar; 
	vector<int> alignments; 
	vector<int> Positions; 
	vector<string> ChrPositions;
	int AlignmentSegments; 
	int AlignmentSegmentsCigar;
	vector<long> MutAltCounts; 
	vector<long> MutRefCounts;
	vector<long> MutContigCounts; 
	vector<long> MutHashListCounts; 

	vector<vector <long>> RefAltCounts; 
	vector<vector <long>> RefRefCounts; 
	vector<string> AltKmers;
	vector<string> RefKmers;	
	bool first; // = true;  
	bool combined; // = false; 
	vector<bool> PeakMap; 

	void createPeakMap(); 
	void parse(string read);
	void getRefSeq();
	void CountAlignmentSegments();
	void CountAlignmentSegmentsCigar();
	void processCigar();
	void parseInsertions( SamRead B); 
	void parseMutations( char *argv[], vector<SamRead>& reads );
	void GetModes(int pos, string alt, string reff, int &MutRefMode, int &MutAltMode, vector <int> &ParRefModes, vector <int> &ParAltModes,  vector <int> &HashCounts, vector <int> &HashCountsOG, int &PossibleVarKmer);
	//string ShittyGenotyper(int Alt, int Ref); 	
	int GetSupportingHashCount(int pos, string alt,  string reff);
	void processMultiAlignment(); 
	void write();
	void writeVertical();
	void writetofile(ofstream &out);
	void flipRead(); 
	void LookUpKmers();
	void CheckPhase(); 
	void FixTandemRef();
	int isPolyA(vector<SamRead> r);
	int CheckParentCov(int &mode); 
	bool StartsWithAlign(int &pos, string &insert); 
	bool EndsWithAlign(int &pos, string &insert); 
	bool StartsWithAlignAtPeak(int &pos, string &insert, int &Kdepth);
	bool EndsWithAlignAtPeak(int &pos, string &insert, int &Kdepth);
	void checkMob(unordered_map <string, MobRead> m); 
	int sigBreakPoint(); 
	int CountBasesAligned(int start); 
	int BreakPoint(); 
	string createStructGenotype(int pos); 
	string ClipPattern(); 
	string filterSV(); 
	bool CheckEndsAlign(); 
	int CheckBasesAligned();
	int SVCheckParentsForLowCov(int spot); 

	
	vector <string> hashes;
	vector <string> hashesRef;
	vector <bool> varHash;
	vector <bool> candidateHash; 
	vector <vector<int>> parentCounts;
	vector <vector<int>> parentCountsReference;
	vector <int> mutCounts;
	vector <int> mutCountsRef;
	void BuildUpHashCountTable(); 
	void GetQualityHashes(int &Mut, int &Pos, int spot);  
	string getClippedSequence(int pos, string type);
};
string SamRead::getClippedSequence(int pos, string type)
{
	int start = -1; 
	int end = -1;
	if (type == "mc")
	{
		start = pos; 
		end = seq.size(); 
	}
	else if( type == "cm")
	{
		start = 0; 
		end = pos; 
	}
	else 
	{
		cout<< "unrecognised type = " << type; 
		return ""; 
	}
	stringstream returnseq; 
	for(int i = start; i < end; i++)
	{
		returnseq<< seq.c_str()[i]; 
	}
	return returnseq.str(); 
}
string SamRead::filterSV()
{
	string Filter = "";

	if(StrandBias>=0)
	{
		if (StrandBias >0.99 || StrandBias < 0.01)
			Filter+="SB;"; 
	}
	if(AlignmentSegments > SegThreshold || AlignmentSegmentsCigar > SegThreshold)
	{
		Filter+="PA;"; 
	}
	return Filter; 
}

string SamRead::ClipPattern()
{
	//cout << "checking clip pattern"; 
	char last = 'a'; 
	string pattern = "";
	int count = 0;
	
	if (cigarString.c_str()[0] == 'H' || cigarString.c_str()[0] == 'S')
	{
		last = 'c';
	}
	else
	{
		last = 'm';
	}

	for (int i =1; i<cigarString.size(); i++)
	{
		if (cigarString.c_str()[i] == 'H' || cigarString.c_str()[i] == 'S')
		{
			if(last == 'c')
			{
				count++; 
			}
			else 
			{
				if( count > 10)
				{
					pattern +=last;
				}
				last = 'c';
				count = 1; 
			}
		}
		else
		{
			if (last == 'm')
			{
				count ++; 
			}
			else
			{
				if(count> 10)
				{
					pattern +=last;
				}
				last = 'm';
				count = 1; 
			}
		}
	}
	if(count> 10)
	{
		pattern +=last;
	}
	
	//write();
	//cout << "pattern = " << pattern << endl;
	return pattern; 
};
void SamRead::GetQualityHashes(int &Mut, int &Pos, int spot)
{
	cout << "struct GettingHashes" << endl;
	write();
	vector <int> MutAltCounts;
	vector <int> MutRefCounts;
	vector <int> a;
	a.clear();
	///build up mutant kmer depths for ALT 
	int start = spot - HashSize+1;
	if (start <0)
		start = 0;

	string hash  ="";
	string ref = "";
	string lastHash = ""; 
	int PossibleVarKmer = 0;
	cout << "starting loop, start = " << start << " spot = " << spot << " seq size = " << seq.size() << endl; 

	for(int i = start; i <= spot && i < seq.size()-HashSize; i++)
	{
		hash = seq.substr(i, HashSize);
		ref = RefSeq.substr(i, HashSize);
		unsigned long int LongHash = HashToLong(hash);
		if ( hash != ref and (ExcludeHashes[HashToLong(hash)]<1 &&  ExcludeHashes[HashToLong(RevComp(hash))]<1) and hash != lastHash)
		{
			if (Hash.count(hash) > 0)
			{
				MutAltCounts.push_back(Hash[hash]);
			}
			else if (Hash.count(RevComp(hash)) > 0)
			{
				MutAltCounts.push_back(Hash[RevComp(hash)]);
			}
	
		//	if ( hash != ref and (ExcludeHashes[HashToLong(hash)]<1 &&  ExcludeHashes[HashToLong(RevComp(hash))]<1) and hash != lastHash)
			{
	
				PossibleVarKmer++;
		//	      cout << "Different" << endl;
			}

		}
		lastHash = hash; 
	}
	Mut = MutAltCounts.size();
	Pos = PossibleVarKmer; 

}
bool CheckGenotypes(string genotypes)
{
	vector <string> temp = Split(genotypes, '\t');
	if (temp.size() < 1)
		return false; 
	for(int i =0; i < temp.size(); i++)
	{
		if (temp[i].c_str()[0] == '.')
			return false;
	}
	return true; 
}
int SamRead::SVCheckParentsForLowCov(int spot)
{
	int MinParCov = 1; 
	cout << "CheckingParLowCov" << endl; 
	write(); 
	vector <int> MutAltCounts; 
	vector<vector <int>> sParAltCounts; 
	vector <int> a;
	a.clear(); 
	for (int i =0; i<ParentHashes.size(); i++)
	{
		sParAltCounts.push_back(a); 
	
	}
	
	int start = spot - HashSize+1;
	if (start <0)
		start = 0; 
	string hash  ="";

	vector <int> streak; 
	for (int i =0; i<sParAltCounts.size(); i++)
	{
		streak.push_back(0);
	}
	for(int i = start; i < spot && spot < seq.size(); i++)
	{
		hash = seq.substr(i, HashSize);
		cout << "hash = " << hash << endl; 
		cout << i << "\t" << spot << "\t" << hash  ;
		unsigned long int LongHash = HashToLong(hash);
		if (ExcludeHashes[HashToLong(hash)]<1 && ExcludeHashes[HashToLong(RevComp(hash))]<1)	
		{	
			int numlow = 0; 
			for (int k =0; k<sParAltCounts.size(); k++)
			{
				if (ParentHashes[k].count(LongHash)>0 )
				{
					cout << "\t" << ParentHashes[k][LongHash]; 
					if (ParentHashes[k][LongHash] > 0 && ParentHashes[k][LongHash] <= ParLowCovThreshold )
					{
						numlow ++; 
					}
				}
				else
					cout << "\t" << "-1";
			}	
			if(numlow == 1)
			{
				
				cout << " keeping"; 
				for (int k =0; k<sParAltCounts.size(); k++)
				{
					if (ParentHashes[k].count(LongHash)>0 )
					{
						if (ParentHashes[k][LongHash] > MinParCov && ParentHashes[k][LongHash] <=ParLowCovThreshold)
						{
                                                        streak[k]++;
							cout << "\ts=" << streak[k]; 
							if (streak[k] >=3)
								sParAltCounts[k].push_back(ParentHashes[k][LongHash]);
						}
						else
							streak[k]=0; 
					}
					else
						streak[k]=0; 
				}
			}
			else
			{
				cout << " common reject";
				for (int k =0; k<sParAltCounts.size(); k++)
				{
					streak[k]=0; 
				}
			}
		}
		else
			cout << "common hash"; 

		cout << endl; 

	}
	int LowCovPar = 0; 
	cout << "summing low cov par" << endl; 
	for (int i =0; i<sParAltCounts.size(); i++)
	{
		cout << "parent[" << i << "] = " << sParAltCounts[i].size()<< endl ;
		if (sParAltCounts[i].size() >=1)
			LowCovPar ++;  
	}
	return LowCovPar; 
}
string SamRead::createStructGenotype(int spot)
{
	if (spot <= 0)
		return ""; 

	cout << "struct genotyper" << endl; 
	write(); 
	vector <int> MutAltCounts; 
	vector <int> MutRefCounts; 
	vector<vector <int>> sParAltCounts; 
	vector<vector <int>> sParRefCounts; 
	vector <int> a;
	a.clear(); 
	for (int i =0; i<ParentHashes.size(); i++)
	{
		sParAltCounts.push_back(a); 
		sParRefCounts.push_back(a);
	
	}
	cout << "RefSizes before anlying is added are " << endl;
	for (int i =0; i<sParRefCounts.size(); i++)
	{
		cout << sParRefCounts[i].size() << endl;
	}
	///build up mutant kmer depths for ALT 
	int start = spot - HashSize;
	if (start <0)
		start = 0; 
	string hash  =""; 
	cout << "hereTESTING " << endl; 
	for(int i = start; i < spot && i+HashSize < seq.size(); i++)
	{
		hash = seq.substr(i, HashSize);
		unsigned long int LongHash = HashToLong(hash);
		
		if (Hash.count(hash) > 0)//if (MutantHashes.count(LongHash))
		{
			MutAltCounts.push_back(Hash[hash]);
			cout << hash << "\t" << Hash[hash] << " i = " << i << " seq.size = " <<  seq.size() << " HashSize = " << HashSize ; 
			for (int i =0; i<sParAltCounts.size(); i++)
			{
				if (ParentHashes[i].count(LongHash)>0)
				{
					sParAltCounts[i].push_back(ParentHashes[i][LongHash]);
					cout << "\t" << ParentHashes[i][LongHash];
				}	
				else
				{
					cout << "\t-1";
				}
			}
		}
		else if (Hash.count(RevComp(hash)) > 0)
		{ 
			MutAltCounts.push_back(Hash[RevComp(hash)]);
			cout << hash << "\t" << Hash[RevComp(hash)] << " i = " << i << " seq.size = " <<  seq.size() << " HashSize = " << HashSize;
			for (int i =0; i<sParAltCounts.size(); i++)
			{       
				if (ParentHashes[i].count(LongHash)>0)
				{
					sParAltCounts[i].push_back(ParentHashes[i][LongHash]);
					cout << "\t" << ParentHashes[i][LongHash];
				}
				else
				{
					cout << "\t-1";
				}
			}
		}
		else 
		{	
			cout << hash << "\tnone";
		}
		
		cout << endl; 

	}
	
	cout << "RefSizes after anlying is added are " << sParRefCounts.size() << endl;
	for (int i =0; i<sParRefCounts.size(); i++)
	{
		cout << sParRefCounts[i].size() << endl;
	}



	///Build up mutant kmer depths for REF 
	cout << "chr = " << chr << " pos = " << pos << "+"<< spot <<"-"<<HashSize << " = " << pos+spot-HashSize << endl;   
	int pullstart = pos+spot-HashSize;
	int pullend = HashSize+HashSize; 
	if (pullstart < 0)
	{
		pullend += pullstart; 
		pullstart = 0;
	}
	string refs = Reff.getSubSequence(chr, pos+spot-HashSize, HashSize+HashSize);
	cout << "refs = " << refs << endl ; 
	cout << "hereTESTING" << endl; 
	for (int i = 0; i+HashSize < refs.size() && refs.size() > 0 ; i++)
	{
		hash = refs.substr(i, HashSize);
		cout << "hash = " << hash << " i = " << i << " refs.size = " << refs.size() << endl; 
		unsigned long int LongHash = HashToLong(hash);
		if (MutantHashes.count(LongHash) > 0)
		{
			MutRefCounts.push_back(MutantHashes[LongHash]); 
			cout << hash << "\t" << MutantHashes[LongHash]; 
		}
		else 
		{       
			cout << hash << "\t-1";
		}
		cout << "here" << endl; 
		for (int j =0; j<sParRefCounts.size(); j++)
		{
			if (ParentHashes[j].count(LongHash)>0)
			{
				sParRefCounts[j].push_back(ParentHashes[j][LongHash]);
				cout << "\t" << ParentHashes[j][LongHash]; 
			}
			else
			{
				cout << "\t-1";
			}
		}
		cout << endl; 
	}
	cout << "RefSizes are" << endl; 
	for (int i =0; i<sParRefCounts.size(); i++)
	{
		cout << sParRefCounts[i].size() << endl;
	}

	sort (MutAltCounts.begin(), MutAltCounts.end());
	sort (MutRefCounts.begin(), MutRefCounts.end());
	
	for (int i =0; i<sParRefCounts.size(); i++)
	{
		sort (sParRefCounts[i].begin(), sParRefCounts[i].end()); 
	}
	for (int i =0; i<sParAltCounts.size(); i++)
	{
		sort (sParAltCounts[i].begin(), sParAltCounts[i].end());
	}
	cout << "done sorting";	

	int MutAlt; 
	int MutRef; 
	vector <int> ParAlt; 
	vector <int> ParRef; 

	if (MutAltCounts.size()>0)
		MutAlt = MutAltCounts[0]; 
	else
		MutAlt = 0; 
	
	if (MutRefCounts.size()>0)
		MutRef = MutRefCounts[0];
	else    
		MutRef = 0;
	cout << "build Mut Alelles " << MutAlt << " " << MutRef << endl; 

	for (int i =0; i<sParAltCounts.size(); i++)
	{
		if (sParAltCounts[i].size()>0)
			ParAlt.push_back(sParAltCounts[i][0]);
		else
			ParAlt.push_back(0); 
	}

	for (int i =0; i<sParRefCounts.size(); i++)
	{
		if (sParRefCounts[i].size()>0)
			ParRef.push_back(sParRefCounts[i][0]);
		else    
			ParRef.push_back(0);
	} 

	cout << "counts are: \nMut: " << MutAlt << "\t" <<MutRef << endl;
	for (int i =0; i < ParAlt.size(); i++)
	{
		cout << ParAlt[i] << "\t" << ParRef[i] << endl;
	}
	cout << "done with counts" << endl;
	stringstream ss;
	ss << ShittyGenotyper(MutAlt, MutRef) << ":" << MutAlt+MutRef << ":" << MutRef << ":" << MutAlt; 
	for (int i =0; i<sParAltCounts.size(); i++)
	{
		ss <<"\t" << ShittyGenotyper(ParAlt[i], ParRef[i]) << ":" << ParAlt[i]+ParRef[i] << ":" << ParRef[i] << ":" << ParAlt[i];
	}
	cout << "returning " << ss << endl;
	return ss.str(); 
} 
int SamRead::BreakPoint()
{
	for (int i = 1; i<seq.size(); i++)
	{
		if ((cigarString.c_str()[i-1] == 'H' || cigarString.c_str()[i-1] == 'S') && (cigarString.c_str()[i] == 'M' || cigarString.c_str()[i] == 'X' || cigarString.c_str()[i] == 'D' || cigarString.c_str()[i] == 'I' ))
		{
			return i;
		}
		else if  ((cigarString.c_str()[i] == 'H' || cigarString.c_str()[i] == 'S') && (cigarString.c_str()[i-1] == 'M' || cigarString.c_str()[i-1] == 'X' || cigarString.c_str()[i-1] == 'D' || cigarString.c_str()[i-1] == 'I'))
		{
			return i;
		}
	}
	return -1;
}
bool BreakpointInUnalignedCenter(SamRead A, SamRead B)
{
	if (GetReadOrientation(A.flag) != GetReadOrientation(B.flag))
	{
		B.flipRead();
	}
	
	int delFixA = 0; 
	int delFixB = 0;  
	bool StartAlign = false; 
	bool EndAlign = false; 
	bool InUnAlign = false; 
	int centerPeak = 0; 
	for (int i = 0; i + delFixA < A.seq.size() && i + delFixB < B.seq.size() ; i++)
	{

		///////keep everything lined up /////
		while (A.seq.c_str()[i+delFixA] =='-')
		{
			cout << "fixing base" << endl; 
			delFixA++;
			
		}
		while(B.seq.c_str()[i+delFixB] == '-')
		{
			delFixB++;
		}
		/////////////////////////////////////
		if ( (A.cigarString.c_str()[i+delFixA] != 'H' && A.cigarString.c_str()[i+delFixA] != 'S') || (B.cigarString.c_str()[i+delFixB] != 'S' && B.cigarString.c_str()[i+delFixB] != 'H'))
		{
			if (StartAlign == false && EndAlign == false && InUnAlign == false)
				StartAlign = true;
			if (StartAlign == true && EndAlign == false && InUnAlign == true)
				EndAlign = true; 
		}
		else if((A.cigarString.c_str()[i+delFixA] == 'H' || A.cigarString.c_str()[i+delFixA] == 'S') && (B.cigarString.c_str()[i+delFixB] == 'S' || B.cigarString.c_str()[i+delFixB]== 'H'))
		{
			InUnAlign = true; 
			if (centerPeak == 0)
			{
				EndAlign = false; 
			}
			if (A.PeakMap[i+delFixA] == true || B.PeakMap[i+delFixB] == true)
			{
				centerPeak++;
			}
		}
		else 
		{
			cout << "not sure how I got here " << endl; 
		}
	}
	
	if (StartAlign == true && EndAlign == true && InUnAlign == true && centerPeak > 0)
	{
		cout << "BreakpointInUnalignedCenter true with " << centerPeak << " peaks" << endl; 
		return true; 
	}
	else
		return false; 


	
}
//bool BreakpointInUnalignedCenter(SamRead A, SamRead B)
//{
//	if (A.name == B.name)
//	{
//		if (A.clipPattern == "cm" && B.clipPattern == "mc" && GetReadOrientation(A.flag) == GetReadOrientation(B.flag))
//		{
//			int end = A.BreakPoint();
//			int start = B.BreakPoint();
//			for (int i = start; i <=end; i++)
//			{
//				if (A.PeakMap[i])
//					return true; 
//			}
//		}
//		if (B.clipPattern == "cm" && A.clipPattern == "mc" && GetReadOrientation(A.flag) == GetReadOrientation(B.flag))
//		{
//			int end =B.BreakPoint();
//			int start = A.BreakPoint();
//			for (int i = start; i <=end; i++)
//			{
//				if (A.PeakMap[i])
//					return true;
//			}
//		}
//		//////NEED TO ADD STUFF IF THEY ARE NOT ON THE SAME STRAND
//
//	}
//	return false; 
//}

int SamRead::CountBasesAligned(int start)
{
	int count = 0; 
	for (int i = start; i<seq.size(); i++)
	{
		if  (cigarString.c_str()[i] == 'M' || cigarString.c_str()[i] == 'X' || cigarString.c_str()[i] == 'D' || cigarString.c_str()[i] == 'I' )
		{
			count++; 
		}
		else
			return count; 
	}
	return count; 
}
int SamRead::sigBreakPoint()
{
	for (int i = 1; i<seq.size(); i++)
	{
		if ((cigarString.c_str()[i-1] == 'H' || cigarString.c_str()[i-1] == 'S') && (cigarString.c_str()[i] == 'M' || cigarString.c_str()[i] == 'X' || cigarString.c_str()[i] == 'D' || cigarString.c_str()[i] == 'I' ))
		{
			if (PeakMap[i-1] || PeakMap[i])
				return i; 
		}
		else if  ((cigarString.c_str()[i-1] == 'M' || cigarString.c_str()[i-1] == 'X' || cigarString.c_str()[i-1] == 'D' || cigarString.c_str()[i-1] == 'I') && (cigarString.c_str()[i] == 'H' || cigarString.c_str()[i] == 'S')  )
		{
			if (PeakMap[i-1] || PeakMap[i])
				return i;
		}
	}
	return -1; 

};
int SamRead::isPolyA(vector <SamRead> R)
{
	int MinPolyAsize = 10; 
	cout << "in poly a test" << endl; 
	cout << "size of array = " << R.size() << endl; 
	int start = -1; 
	int end = -1;
	char base = 'f'; 
	bool clipped = false;
	bool atpeak = false;
	cout << "Checking PolyAstatus" << endl; 
	write(); 
	for (int j = 0; j<R.size(); j++)
	{
		//R[j].write(); 
		//cout << R[j].seq.size() << " != " << seq.size() << endl;
		//if (R[j].seq.size() != seq.size())
		//{
		//	cout << "warning not all the same size" << endl;
		//	return -1; 
		//}
		if (GetReadOrientation(flag) != GetReadOrientation(R[j].flag))
		{
			R[j].flipRead(); 
		}
	}
	write(); 
	for (int j = 0; j<R.size(); j++){R[j].write();}
	cout << "got all the reads in the same orientation" << endl;  
	bool check = false; 
	
	int delFix = 0; 
	vector <int> fix; 
	for (int i =0; i<R.size(); i++)
	{
		int f = 0; 
		fix.push_back(f);
	}
	
	
	 
	for (int i = 0; i + delFix < seq.size(); i++)
	{

			int u =0; 
			////////////////////////////////////////////////
			//cout << "i = " << i << "\tdelFix= " << delFix << "\t" << seq.c_str()[i+delFix] << "\t" << base << "\t" << cigarString.c_str()[i+delFix] ;  
			for (int j =0; j < R.size(); j++)
			{
				cout << "\t" << fix[j] << "\t" << R[j].seq.c_str()[i+fix[j]] << "\t" << R[j].cigarString.c_str()[i+fix[j]]; 
			}
			cout << endl; 
			while (seq.c_str()[i+delFix] =='-')
			{
				cout << "fixing base" << endl; 
				delFix++;
				
			}
			for (int j =0; j < R.size(); j++)
			{
				while(R[j].seq.c_str()[i+fix[j]] == '-')
				{
					fix[j]++;
				}
			}
			///////////////////////////////////////////////
			check = false;
			for (int j =0; j < R.size(); j++)
			{
				if (R[j].mapQual > 0 && R[j].cigarString.c_str()[i+fix[j]] != 'S' && R[j].cigarString.c_str()[i+fix[j]] != 'H')
					check =true;
			}

		if ( base == 'f' && (seq.c_str()[i+delFix] == 'T' || seq.c_str()[i+delFix] == 'A') && ( cigarString.c_str()[i+delFix] == 'H' || cigarString.c_str()[i+delFix] == 'S') && check == false)
		{
			base = seq.c_str()[i+delFix]; 
			start = i+delFix;
			cout << "found a start " << base << " " << start  << endl; 
		}
		else if (base != 'f' && seq.c_str()[i+delFix] == base && ( cigarString.c_str()[i+delFix] == 'H' || cigarString.c_str()[i+delFix] == 'S') && check == false)
		{
			cout << "on a run " << start << " " << base << endl; 
		}
		else if (base != 'f' && (( seq.c_str()[i+delFix] != base || ( cigarString.c_str()[i+delFix] != 'H' && cigarString.c_str()[i+delFix] != 'S')) || check == true ))
		{
			end = i+delFix;
			int size = end - start; 
			cout << "ended run " << end << " " << seq.c_str()[i] << " " << end << " so size = "  << size << endl; 
			if (size> MinPolyAsize)
			{
				cout << "chekcing peak and clipped" << endl; 
				for(int j = start; j <=end; j++)
				{
					if ( PeakMap[j] == true)
					{	atpeak = true;}
					if ( cigarString.c_str()[j] == 'H' || cigarString.c_str()[j] == 'S')
					{	clipped = true;} 
				}
				cout << "done checking peak and clipped" << endl; 
			}
			cout << "clipped = " << clipped << " and atpeak = " << atpeak << endl; 
			if (clipped == true  && atpeak  == true)
			{
				cout << "found poly A with start " << start << "and end " << end << "with base " << base << endl; 
				if (clipPattern == "mc")
					return start;
				else if (clipPattern == "cm")
					return end; 

			}
			clipped = false; 
			atpeak = false; 
			base = 'f'; 
			start = -1;
			end = -1; 
		}
	}
	if (base != 'f' &&  seq.c_str()[seq.size()-1] == base && check == false )
	{
		end = seq.size()-1;
		int size = end - start;
		cout << "ended run " << end << " " << seq.c_str()[seq.size()-1] << " " << end << " so size = "  << size << endl;
		if (size> MinPolyAsize)
		{
			cout << "chekcing peak and clipped" << endl;
			for(int j = start; j <=end; j++)
			{
				if ( PeakMap[j] == true)
				{       atpeak = true;}
				if ( cigarString.c_str()[j] == 'H' || cigarString.c_str()[j] == 'S')
				{       clipped = true;}
			}
		}
	//		cout << "clipped = " << clipped << " and atpeak = " << atpeak << endl;
		if (clipped == true  && atpeak  == true)
		{
	//	cout << "found poly A with start " << start << "and end " << end << "with base " << base << endl;
			return start;
		}
	}

	return -1; 
}

void SamRead::checkMob(unordered_map <string, MobRead> m)
{
	//cout << "mob check of " << name << endl; 
	if (m.count(name) > 0)
	{
		MobAligned = true; 
		MobContig = m[name].chr; 
		MobAS = m[name].AS; 
	//	cout << "mob found " << MobContig << endl;
	}


}

void SamRead::BuildUpHashCountTable()
{
  /////////////////Building up varHash and hash lists ///////////// 
	//cout << "Building up varHash" << endl;
	for (int i = 0; i <  seq.size() - HashSize; i++)
	{
		string newHash = "";
		string newHashRef = "";
		newHash += seq.c_str()[i];
		newHashRef += RefSeq.c_str()[i];
		int count = 0;
		////can i replace this with get hash ?  
		if ((cigarString.c_str()[i] != 'D' and cigarString.c_str()[i] != 'R' and cigarString.c_str()[i] != 'H'))
		{
			for (int j = 1; j<seq.size() - i and count < HashSize-1 ; j++)
			{

				if (cigarString.c_str()[i+j] != 'D' and cigarString.c_str()[i+j] != 'R'  and cigarString.c_str()[i+j] != 'H')
				{
					newHash += seq.c_str()[i+j];
					newHashRef += RefSeq.c_str()[i+j];
					count++;
				}
			}
		}
		hashes.push_back(newHash);
		  hashesRef.push_back(newHashRef);
		if (Hash.count(newHash) > 0 or Hash.count(RevComp(newHash)) > 0)
			varHash.push_back(true);
		else
			varHash.push_back(false);
	}
	///////////////////////////////////////////////////


	///////////////////building up parent hash counts //////////////////
	//cout << "Bulding Par hash counts" << endl;
	//vector <long int> ParentHash;
	for(int pi = 0; pi<ParentHashes.size(); pi++)
	{
		vector <int> counts;
		vector <int> countsRef;
		for(int i = 0; i< hashes.size(); i++)
		{
			string hash =  hashes[i];
			string hashRef = hashesRef[i];
			bool checkHash = true;
			for (int j = 0; j < HashSize; j++)
			{
				if (!(hash[j] == 'A' or hash[j] == 'C' or hash[j] == 'G' or hash[j] == 'T'))
				{
					checkHash = false;
					break;
				}
			}
			if (checkHash)
			{
				unsigned long int LongHash = HashToLong(hash);
				if (ParentHashes[pi].count(LongHash) >0)
					counts.push_back(ParentHashes[pi][LongHash]);
				else
					counts.push_back(0);
				 unsigned long int LongHashRef = HashToLong(hashRef);
				if (ParentHashes[pi].count(LongHashRef) >0)
					countsRef.push_back(ParentHashes[pi][LongHashRef]);
				else
					countsRef.push_back(0);
			}
			else
			{
				counts.push_back(-1);
				countsRef.push_back(-1);
			}
		}
		parentCounts.push_back(counts);
		parentCountsReference.push_back(countsRef);
	}
	 //////////////////////////////////////////////////////////////////
			/////////////////////bulid Mut counts/////////////////////////////
	//cout << "bulding mut counts" << endl; 
	//cout << hashes.size() << endl;
	//cout << hashesRef.size() << endl; 
	for(int i = 0; i< hashes.size(); i++)
	{
	//	cout << i<< endl;
		string hash =  hashes[i];
	//	cout << "   hash = " << hash << endl;
		string hashRef = hashesRef[i];
	//	cout << "RefHash = " << hashRef << endl;
		bool checkHash = true;
		for (int j = 0; j < HashSize; j++)
		{
			if (!(hash[j] == 'A' or hash[j] == 'C' or hash[j] == 'G' or hash[j] == 'T'))
			{
				checkHash = false;
				break;
			}
		}
	//	cout << "check hash = " << checkHash << endl;
		if (checkHash)
		{
			unsigned long int LongHash = HashToLong(hash);
			if (MutantHashes.count(LongHash) >0)
			{       mutCounts.push_back(MutantHashes[LongHash]);}
			else
			{       mutCounts.push_back(0);}

			unsigned long int LongHashRef = HashToLong(hashRef);
			if (MutantHashes.count(LongHashRef) >0)
			{       mutCountsRef.push_back(MutantHashes[LongHashRef]);}
			else
			{       mutCountsRef.push_back(0);}
		}
		else
		{
			mutCounts.push_back(-1);
			mutCountsRef.push_back(-1);
		}
	}
	/////////////////////////////////////////////////////////////////////

	 ////////////////////write out vertical table/////////////////////////   
//	cout << "writing hashes out vert" << endl;
//	for(int i =0; i < hashes.size(); i++)
//	{
//		cout << i+pos << "\t" <<  i  << "\t" << hashes[i] << "\t" << varHash[i] << "\t" << PeakMap[i] << "\t" << (int) qual.c_str()[i]-33;
//		cout << "\t" << "MutVar-" << mutCounts[i]; 
//		for (int j = 0; j < parentCounts.size(); j++)
//		{
//				cout << "\t" <<  parentCounts[j][i];
//		}
//		cout << "\t" << "MutRef-" << mutCountsRef[i];
//		for (int j = 0; j < parentCountsReference.size(); j++)
//		{
//				cout << "\t" <<  parentCountsReference[j][i];
//		}
//		cout << endl;
//
//	}
	////////////////////////////////////////////////////////////////
}
int SamRead::GetSupportingHashCount(int pos, string alt,  string reff)
{
				int Count =0; 
				int lower = pos-HashSize;
				if (lower < 0){lower =0;}
				int upper = pos+alt.length()+reff.length();//-1;
	//			cout << pos<<" + "<< alt.length() << " + " << reff.length() << endl;
				if (upper > MutRefCounts.size()){
	//				cout << "this is going to break " << upper << " > " << MutRefCounts.size() << endl;
					upper = MutRefCounts.size();
				}

				for (int j = lower; j<upper; j++)
				{
					if (Hash.count(AltKmers[j]) > 0 and Hash[AltKmers[j]] > 0)
						Count++; 
					else if (Hash.count(RevComp(AltKmers[j])) > 0 and Hash[RevComp(AltKmers[j])])
						Count++; 
				}
	return Count; 
}
void SamRead::GetModes(int pos, string alt,  string reff, int &MutRefMode, int &MutAltMode, vector <int> &ParRefModes, vector <int> &ParAltModes,  vector <int> &HashCounts, vector <int> &HashCountsOG, int &PossibleVarKmer)
{
	int lower = pos-HashSize+1;
	if (lower < 0){lower =0;}
	int upper = pos+alt.length()+reff.length()-1;
	cout << pos<<" + "<< alt.length() << " + " << reff.length() << endl;
	if (upper > MutRefCounts.size()){
		cout << "this is going to break " << upper << " > " << MutRefCounts.size() << endl;
		upper = MutRefCounts.size();
	}

	//////////////chekcing allele frequencies /////////// 
	//vector <int> HashCountsOG;
	vector <int> varMutRefCounts;
	vector <int> varMutAltCounts;
	vector<vector <int>> varParRefCounts;
	vector<vector <int>> varParAltCounts;
	vector <int> temp;
	for(int pi = 0; pi<ParentHashes.size(); pi++){
		varParRefCounts.push_back(temp);
		varParAltCounts.push_back(temp);
	}
	string nonspecific = "nonspecific";
	//vector <float> freqs;
	//cout << "checking NonSpecic Kmers";
	string LastAtlKmer = "boomba";
	for (int j = lower; j<upper; j++)
	{
	//	cout <<  AltKmers[j] << endl << RefKmers[j] <<endl;
		if ( AltKmers[j] != RefKmers[j] and (ExcludeHashes[HashToLong(AltKmers[j])]<1 && ExcludeHashes[HashToLong(RevComp(AltKmers[j]))]<1) and AltKmers[j] != LastAtlKmer)
		{
			
			PossibleVarKmer++;
	//		cout << "Different" << endl;
		} 
	//	else
	//		cout << "SAME" << endl; 
		LastAtlKmer = AltKmers[j]; 
		if(MutRefCounts[j]>0 and AltKmers[j] != RefKmers[j] ) //and MutRefCounts[j]<200 and (ExcludeHashes[HashToLong(RefKmers[j])]<2 or ExcludeHashes[HashToLong(RevComp(RefKmers[j]))]<2)) //needs to be fixed, should be based on cov not cutoff of 200 
			varMutRefCounts.push_back(MutRefCounts[j]);
		if (MutAltCounts[j]>0 and AltKmers[j] != RefKmers[j] and MutAltCounts[j]<200 and (Hash.count(AltKmers[j]) > 0 or Hash.count(RevComp(AltKmers[j])) > 0)  and (ExcludeHashes[HashToLong(AltKmers[j])]<1 && ExcludeHashes[HashToLong(RevComp(AltKmers[j]))]<1)) //needs to be fixed, should be based on cov not cutoff of 200 
			varMutAltCounts.push_back(MutAltCounts[j]);
		
		for (int pi=0; pi < varParRefCounts.size(); pi++){
			if (RefRefCounts[pi][j] >0 and AltKmers[j] != RefKmers[j] )//and RefRefCounts[pi][j] <200 and (ExcludeHashes[HashToLong(AltKmers[j])]<2 or ExcludeHashes[HashToLong(RevComp(AltKmers[j]))]<2)) //needs to be fixed, should be based on cov not cutoff of 200 
				varParRefCounts[pi].push_back(RefRefCounts[pi][j]);
			if (RefAltCounts[pi][j] >0 and AltKmers[j] != RefKmers[j] and RefAltCounts[pi][j] < 200 and (Hash.count(AltKmers[j]) > 0 or Hash.count(RevComp(AltKmers[j])) > 0)   and (ExcludeHashes[HashToLong(RefKmers[j])]<1 && ExcludeHashes[HashToLong(RevComp(RefKmers[j]))]<1)) //needs to be fixed, should be based on cov not cutoff of 200 
				varParAltCounts[pi].push_back(RefAltCounts[pi][j]);

		}

		if (Hash.count(AltKmers[j]) > 0 and Hash[AltKmers[j]] > 0 and AltKmers[j] != RefKmers[j] )
			HashCountsOG.push_back(Hash[AltKmers[j]]);
		else if (Hash.count(RevComp(AltKmers[j])) > 0 and Hash[RevComp(AltKmers[j])] and AltKmers[j] != RefKmers[j] )
			HashCountsOG.push_back(Hash[RevComp(AltKmers[j])]);
		
					
		if (Hash[AltKmers[j]] > 0 and AltKmers[j] != RefKmers[j])
			HashCounts.push_back(Hash[AltKmers[j]]);
		else if (Hash[RevComp(AltKmers[j])] > 0 and AltKmers[j] != RefKmers[j])
			 HashCounts.push_back(Hash[RevComp(AltKmers[j])]);
		else
			HashCounts.push_back(-1);
	}
       // float freq = 0;
       // if (freqs.size() > 0){
       //	 for (int i =0; i<freqs.size(); i++){
       //		 freq+=freqs[i];
       //	 }
       //	 freq = freq/freqs.size();
       //` }
	//////////////////////////////////////////////

/*	cout << "<><><><><>MutRef<><><><><><>" << endl ;
	for (int s =0; s<varMutRefCounts.size(); s++){
		cout << varMutRefCounts[s] << " " ;
	}
	cout << endl << "<><><><><>MutAlt<><><><><><>" << endl;
	for (int s =0; s<varMutAltCounts.size(); s++)
	{
		cout << varMutAltCounts[s] << " " ;
	}
	cout << endl;
*/
	sort (varMutRefCounts.begin(), varMutRefCounts.end());
	sort (varMutAltCounts.begin(), varMutAltCounts.end());

/*	cout << "<><><><><>MutRefSorted<><><><><><>" << endl ;
	for (int s =0; s<varMutRefCounts.size(); s++)
	{
		cout << varMutRefCounts[s] << " " ;
	}
	cout << endl << "<><><><><>MutAltSorted<><><><><><>" << endl;
	for (int s =0; s<varMutAltCounts.size(); s++)
	{
		cout << varMutAltCounts[s] << " " ;
	
	}
	cout << endl;
*/
	for(int pi = 0; pi<varParRefCounts.size(); pi++)
	{
		 sort (varParRefCounts[pi].begin(), varParRefCounts[pi].end());
/*		cout << "<><><><><>Ref" << pi << "<><><><><<><>" << endl;
		for (int s =0; s<varParRefCounts[pi].size(); s++){
			cout << varParRefCounts[pi][s] << " " ;
		}
		cout << endl;
*/	}
	///////////////////////////////////////////////
	if (varMutRefCounts.size() >1)
		MutRefMode = varMutRefCounts[0];
		//MutRefMode = varMutRefCounts[(varMutRefCounts.size())/2]; //// switch this for line above to get the mode, right now were taking the min 
	else if (varMutRefCounts.size() ==1)
		MutRefMode = varMutRefCounts[0]; 
	else
		MutRefMode = 0;
	if (varMutAltCounts.size() >1)
		 MutAltMode = varMutAltCounts[0];
		//MutAltMode= varMutAltCounts[(varMutAltCounts.size()-2)/2]; //// switch this for line above to get the mode, right now were taking the min 
	else if (varMutAltCounts.size() ==1)
		MutAltMode = varMutAltCounts[0];
	else
		MutAltMode=0;
	for(int pi = 0; pi<varParRefCounts.size(); pi++)
	{
		if (varParRefCounts[pi].size() >1)
			ParRefModes.push_back(varParRefCounts[pi][0]);
			//ParRefModes.push_back(varParRefCounts[pi][((varParRefCounts[pi].size())/2)]);  //// switch this for line above to get the mode, right now were taking the min 
		else if (varParRefCounts[pi].size() ==1)
			ParRefModes.push_back(varParRefCounts[pi][0]);
		else
			ParRefModes.push_back( 0);
	}
	for(int pi =0; pi<varParAltCounts.size(); pi++)
	{
		if (varParAltCounts[pi].size()>1)
			ParAltModes.push_back(varParAltCounts[pi][0]);
			// ParAltModes.push_back(varParAltCounts[pi][((varParAltCounts[pi].size())/2)]); //// switch this for line above to get the mode, right now were taking the min 
		else if (varParAltCounts[pi].size()==1)
			ParAltModes.push_back(varParAltCounts[pi][0]);
		else
			ParAltModes.push_back( 0); 
	}
}
void SamRead::CountAlignmentSegmentsCigar()
{
	AlignmentSegmentsCigar = 0;
	char last = cigar.c_str()[0];
	for (int i =1; i < cigar.size(); i++)
	{
	  if (cigar.c_str()[i] == 'M' or cigar.c_str()[i] == 'S' or cigar.c_str()[i] == 'H' or cigar.c_str()[i] == 'D' or cigar.c_str()[i] == 'I')
	  {}
	  else if (last == 'M' or last == 'S' or last == 'H' or last == 'I' or last == 'D')
	  {
	    AlignmentSegmentsCigar++;
	  }
	  last = cigar.c_str()[i];
	}
	if (last == 'M' or last == 'S' or last == 'H' or last == 'I' or last == 'D')
	{
		AlignmentSegmentsCigar++;
	}
}
void SamRead::CountAlignmentSegments()
{
	AlignmentSegments = 0;
	char last = cigarString.c_str()[0]; 
	for (int i =1; i < cigarString.size(); i++)
	{
	  if (cigarString.c_str()[i] == 'M')
	  {}
	  else if (last == 'M')
	  {
	    AlignmentSegments++; 
	  }
	  last = cigarString.c_str()[i];
	}
	if (last == 'M')
	{
		AlignmentSegments++;
	}
}
int SamRead::CheckParentCov(int &mode)
{
	bool good = true;
	int lowC = 0;  
	vector<int> cov; 
	for (int pi = 0; pi < RefRefCounts.size(); pi++){
		for (int i = 0; i < RefRefCounts[pi].size(); i++){
			if (RefKmers[i] != ""){
				int ParRef = 0;
				int ParAlt = 0; 
				if (RefAltCounts[pi][i] > 0)
					ParAlt = RefAltCounts[pi][i]; 
				if (RefRefCounts[pi][i] > 0)
					ParRef = RefRefCounts[pi][i]; 
				cov.push_back(ParRef+ParAlt); 
				if (ParRef+ParAlt > 0 && ParRef+ParAlt < 10)
					lowC++; 
			}
		}
	}	
	if(cov.size()>1){
	
		sort (cov.begin(), cov.end());
		mode = cov[cov.size()/2];
	}
	else
		mode = -1; 

	return lowC; 
}

void SamRead::flipRead()
{
//	cout <<"FLIPPING reads not on the same strand";
//	write(); 
	string FlipSeq = ""  ;
	string FlipQual = ""   ; 
	string FlipRefSeq = "";
	string FlipCigarString = ""; 
	string FlipStrand = "";
	string FlipclipPattern=""; 
	vector<bool> FlipPeakMap ; 
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
		FlipPeakMap.push_back(PeakMap[i]);
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
	PeakMap=FlipPeakMap; 
	for (int i = clipPattern.size()-1; i >=0; i--)
	{
		FlipclipPattern += clipPattern.c_str()[i];
	}
	clipPattern = FlipclipPattern; 
//	write(); 
	
}

void SamRead::processMultiAlignment()
{
	//check if this is a mis-joined contig

}
void SamRead::write()
{
	cout << name << endl;
	cout << "   parsed= " << parsed << endl; 
	cout << "   flag = " << flag << endl;
	cout << "   mapQual = " << mapQual << endl;
	cout << "   Strand = " << GetReadOrientation(flag) << endl;
	cout << "   Alignments = " << alignments.size() << endl;
	cout << "   SVeventID = " << SVeventid << endl; 
	cout << "   MobAligned = " << MobAligned << endl; 
	cout << "   MobContig = " << MobContig << endl; 
	cout << "   MobAS = " << MobAS << endl; 
	cout << "   AllA = " << AllA << endl; 
	cout << "   isSplitRead = " << isSplitRead << endl; 
	cout << "   chr " << chr << " - " << pos << " qual = " << mapQual <<  endl;
	cout << "   clip  = " << clipPattern << endl;
	cout << "   cigar  = " << cigar << endl;
	cout << "   Seq    = " << seq << endl;
	cout << "   Qual   = " << qual << endl; 
	cout << "   Cigar  = " << cigarString << endl; 
	cout << "   RefSeq = " << RefSeq << endl;
	cout << "   strand = " << strand << endl;
	cout << "  PeakMap = "; 
	for (int i =0; i < PeakMap.size(); i++)
		{cout << PeakMap[i]; }
	cout << endl;
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
	out << "   AlignScore = " << AlignScore << endl; 
	out << "   clip  = " << clipPattern << endl;
	out << "   cigar  = " << cigar << endl;
	out << "   Seq    = " << seq << endl;
	out << "   Qual   = " << qual << endl;
	out << "   Cigar  = " << cigarString << endl;
	out << "   RefSeq = " << RefSeq << endl;
	out << "   PeakMap= "; 
	for (int i =0; i < PeakMap.size(); i++)
	{
		out <<  PeakMap[i] ;
	}
	out << endl; 
	out << "   PMSize = " << PeakMap.size() << endl;
}

void SamRead::writeVertical()
{
	cout << "ParentHashes size = " <<  ParentHashes.size() << "RefAltCounts size " << RefAltCounts.size() << endl;
	cout << name << endl;
	cout << "   flag = " << flag << endl;
	cout << "   chr " << chr << " - " << pos << " qual = " << mapQual <<  endl;
	cout << "   clip  = " << clipPattern << endl;
	cout << "   cigar  = " << cigar << endl;
	cout << "   Alignments = " << alignments.size() << endl; 
	cout << "   AlignScore = " << AlignScore << endl; 
	for (int i =0; i < seq.size(); i++){

		cout << seq.c_str()[i] << "\t" << qual.c_str()[i] << "\t" << cigarString.c_str()[i] << "\t" << RefSeq.c_str()[i] << "\t" << ChrPositions[i] << "\t" << Positions[i] << "\t" << MutAltCounts[i] << "\t" << MutRefCounts[i] << "\t" << MutHashListCounts[i] << "\t" ;
		cout << "\tParents";
		for (int pi=0; pi < RefAltCounts.size(); pi++){
			cout << "\t" << RefAltCounts[pi][i] << "\t" << RefRefCounts[pi][i];
		}
		cout << "\t" << RefKmers[i] << "\t" << AltKmers[i]; 
		cout<< endl;
	}
}
void SamRead::CheckPhase()
{
	cout << "Checking Phasing" << endl; 
	cout << "ParentHashes size = " <<  ParentHashes.size() << " RefAltCounts size " << RefAltCounts.size() << endl;
	cout << name << endl;
	cout << "   flag = " << flag << endl;
	cout << "   chr " << chr << " - " << pos << " qual = " << mapQual <<  endl;
	cout << "   cigar  = " << cigar << endl;
	cout << "   Alignments = " << alignments.size() << endl;
	cout << "   AlignScore = " << AlignScore << endl;
	write(); 
	vector <int> phased ;
	phased.push_back(0); 
	phased.push_back(0); 
	for (int i =0; i < seq.size(); i++){
		cout << seq.c_str()[i] << "\t" << qual.c_str()[i] << "\t" << cigarString.c_str()[i] << "\t" << RefSeq.c_str()[i] << "\t" << ChrPositions[i] << "\t" << Positions[i] << "\tMutContigCounts=" << MutContigCounts[i] << "\t" << MutAltCounts[i] << "\t" << MutRefCounts[i] << "\t" << MutHashListCounts[i] << "\t" ;
		//cout << "made it all the way here " << endl; 
		cout << "\tParents";
		for (int pi=0; pi < RefAltCounts.size(); pi++){
			cout << "\t" << RefAltCounts[pi][i] << "\t" << RefRefCounts[pi][i];
		}
		bool p = false; 
		if (RefAltCounts.size()>=2)
		{
			if (RefAltCounts[0][i] == 0 & RefAltCounts[1][i] > 3 && MutContigCounts[i] > 2	      && RefAltCounts[0][i] >= 0 && RefAltCounts[1][i] >= 0 && RefRefCounts[0][i] >= 0 && RefRefCounts[1][i] >=0)
			{	
				p=true;
				phased[1]++; 
			}	
			else if(RefAltCounts[0][i] > 3 & RefAltCounts[1][i] == 0 && MutContigCounts[i] > 2		 && RefAltCounts[0][i] >= 0 && RefAltCounts[1][i] >= 0 && RefRefCounts[0][i] >= 0 && RefRefCounts[1][i] >=0)
			{	
				p=true;
				phased[0]++;
			}
			else if(RefRefCounts[0][i] ==0 && RefRefCounts[1][i] > 3  && MutContigCounts[i] < -2	       && RefAltCounts[0][i] >= 0 && RefAltCounts[1][i] >= 0 && RefRefCounts[0][i] >= 0 && RefRefCounts[1][i] >=0)
			{	
				p=true; 
				phased[1]++;
			}
			else if(RefRefCounts[0][i] > 3 && RefRefCounts[1][i] == 0   && MutContigCounts[i] < -2		&& RefAltCounts[0][i] >= 0 && RefAltCounts[1][i] >= 0 && RefRefCounts[0][i] >= 0 && RefRefCounts[1][i] >=0)
			{	
				p=true; 
				phased[0]++;
			}
		}
		cout << "\t" << RefKmers[i] << "\t" << AltKmers[i];
		if (p)
		{
			cout << " FOUND ONE " << phased[0] << "-" << phased[1] << "\t"; 
		}
		cout<< endl;
	}
	if ( phased[0]>0 and phased[1] ==0)
	{
		cout << "PHASED CONTIG " << phased[0] << "-" << phased[1]  << endl; 
		

		 ostringstream convert;
		 convert << phased[0]; 

		phase = "PHASED-" + convert.str() + "-" + ParNames[0]; 
	}
	else if( phased[0]==0 and phased[1] >0)
	{

		cout << "PHASED CONTIG " << phased[0] << "-" << phased[1]  << endl; 
		ostringstream convert;
		convert << phased[1];
		phase = "PHASED-" + convert.str() + "-" + ParNames[1]; 
	}
	else if (phased[0]>0 and phased[1] >0)
	{
	
		cout << "Conflicting PHASED CONTIG " << phased[0] << "-" << phased[1]  << endl; 
		ostringstream convert;
		ostringstream convert2;
		convert << phased[1];
		convert2 << phased[0]; 

		phase = "ConflictingPHASED-" + convert.str() + "-" + convert2.str(); 
	}
	cout << "done checking phasing" << endl; 
}
string compressVar(string line, int start, string& StructCall)
{      
	//cout << "compressing var" << endl; 
	char current = line.c_str()[0];
	int currentCount = 1;
	string CV = ""; 
	for (int i = 1; i< line.size(); i++)
	{       
	//	cout << current << endl;
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

				
				if (current == 'Y')
				{
	//				cout << "YAAAY STRUCT" << endl; 
					StructCall = "SVTYPE=DUP;END="; 
					StructCall += convertEND.str(); 
					StructCall += ";SVLEN="; 
					StructCall += convert.str(); 
					StructCall += ";";
	//				cout << StructCall << endl; 
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
	//			cout << "ERROR in compress " << current << " " << currentCount << endl;
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


		if (current == 'Y')
		{
	//		cout << "YAAAY STRUCT" << endl;
			StructCall = "SVTYPE=DUP:TANDEM;END=";
			StructCall += convertEND.str();
			StructCall += ";SVLEN=";
			StructCall += convert.str();
			StructCall += ";";
	//		cout << StructCall << endl;
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
	//	cout << "ERROR in compress " << current << " " << currentCount << endl;
	}
	
	return CV;
}
char next(string& qual , int i)
{

	for(int j = i+1; j < qual.size(); j++)
	{
		if (qual[j] != qual[i] || qual[j] == '!')
			return qual[j]; 
	}
	return qual[i]; 

}
char last(string& qual, int i )
{
	for(int j = i-1; j >=0; j+= -1)
        {
                if (qual[j] != qual[i] || qual[j] == '!')
                        return qual[j];
        }
        return qual[i]; 

}
void SamRead::createPeakMap()
{
	vector<bool> tempPeakMap;
	for (int i =0; i< qual.size()-1; i++)
	{

		if (qual[i] <='!')
		{
			tempPeakMap.push_back(0);
		}
		else
		{
			if (qual[i] >= last(qual, i) && qual[i] >= next(qual, i))
				tempPeakMap.push_back(1);
			else
				tempPeakMap.push_back(0);

		}
	}
	tempPeakMap.push_back(0);

	// I hate one time corrections, but here on is to correct if ther is a del 
	for (int i =0; i< qual.size(); i++)
	{
		if (seq[i] == '-'){
			tempPeakMap[i] == tempPeakMap[i-1];
		}
	}

	PeakMap.clear();	
	PeakMap = tempPeakMap;
}
/*void SamRead::createPeakMap()
{
	vector<bool> tempPeakMap;
	int i =0; 
	int max = -1; 
	int maxSpot = -1; 
	int last = 0; 
	for (int i =0; i< qual.size(); i++)
	{

		if (qual[i] <='!')
		{
			tempPeakMap.push_back(0);
		}
		else
		{
			int j = i; 
			max = qual[j]; 
			while ( j < qual.size() and qual[j] > '!' )
			{
				if (max < qual[j])
				{max = qual[j];}
				j++;
			}
			j = j-1;
			for ( int k = i; k < qual.size() and k <= j ; k++)
			{
				if (qual[k]==max and cigarString[k] != 'H')
					tempPeakMap.push_back(1);
				else
					tempPeakMap.push_back(0);
			}
			i = j;
		}
	}

	// I hate one time corrections, but here on is to correct if ther is a del 
	for (int i =0; i< qual.size(); i++)
	{
		if (seq[i] == '-'){
			tempPeakMap[i] == tempPeakMap[i-1];
		}
	}

	PeakMap.clear();	
	PeakMap = tempPeakMap;
}*/
/*void SamRead::createPeakMap()
{
//	cout << "crateing PeakMap" << endl;
	vector<bool> tempPeakMap;
	int i =0; 
	int max = -1; 
	int maxSpot = -1; 
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
			       tempPeakMap.push_back(1);
			else
				tempPeakMap.push_back(0);
		}
		i = i +k;
	}
	for (int i = tempPeakMap.size()-1; i < qual.size(); i++)
	{
		tempPeakMap.push_back(0); 
	}
	PeakMap.clear();	
	PeakMap = tempPeakMap; 
//	for (int i =0; i< PeakMap.size(); i++)
//		cout << i << "-" << PeakMap[i] << endl;

}*/
void SamRead::parseMutations( char *argv[], vector<SamRead>& reads)
{
	cout << "In Parsing Mutations " << endl;
	BuildUpHashCountTable();
	createPeakMap();
	write(); 
	
	
	string StructCall = ""; 
	string reff = "";
	string alt = "";
	string varType = "";
	for(int i = 0; i<cigarString.size(); i++)
	{
		reff = ""; 
		alt = "";
		varType = "";
		//find the first variant base
		if ((cigarString.c_str()[i] == 'X' or cigarString.c_str()[i] == 'I' or cigarString.c_str()[i] == 'D' or cigarString.c_str()[i] == 'Y'/*  or cigarString.c_str()[i] == 'S' *or cigarString.c_str()[i] == 'H'*/) and RefSeq.c_str()[i] != 'N')
		{
			cout << "found a " << cigarString.c_str()[i] << endl;
			cout << "at pos " << i << " so pos " << i+pos << endl; 
			int size = -1; 
			int startPos = i; 
			bool  AnyBasesOver0  = false; 
			string Denovo = "inherited"; 
			
			if (qual.c_str()[i] > '!')
				AnyBasesOver0 = true;
			if (PeakMap[i] == 1)
				Denovo = "DeNovo";
			
			for(int j = 0; j< cigarString.size() - i; j++)
			{
				if(cigarString.c_str()[i+j] == 'X' or cigarString.c_str()[i+j] == 'D' or cigarString.c_str()[i+j] == 'I' or cigarString.c_str()[i+j] == 'Y' /*or cigarString.c_str()[i+j] == 'S' or cigarString.c_str()[i+j] == 'H'*/)
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
			cout << "size =" << size<< endl;
			
			if (AnyBasesOver0)  //enabling this will only report varites covered by hashes 
			{
				
				if ( cigarString.c_str()[i] == 'I' or cigarString.c_str()[i] == 'D' or cigarString.c_str()[i] == 'Y' /*or cigarString.c_str()[i] == 'S' or cigarString.c_str()[i] == 'H'*/)
				{
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
				
				/////////build the alleles and var type/////////
				for(int j = 0; j<= size; j++)
				{
					if (RefSeq.c_str()[i+j] == 'A' or RefSeq.c_str()[i+j] == 'C' or RefSeq.c_str()[i+j] == 'G' or RefSeq.c_str()[i+j] == 'T')
						reff+=RefSeq.c_str()[i+j]; 
					if (seq.c_str()[i+j] == 'A' or seq.c_str()[i+j] == 'C' or seq.c_str()[i+j] == 'G' or seq.c_str()[i+j] == 'T')
						alt+=seq.c_str()[i+j]; 
					varType += cigarString.c_str()[i+j]; 
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
				////////Starting Checks /////////////     
				int MutRefMode;
				int MutAltMode;
				vector <int> ParRefModes;
				vector <int> ParAltModes;
				vector <int> HashCounts;
			       	vector <int> HashCountsOG;
				int PossibleAltKmer=0; 
				GetModes(i, alt, reff, MutRefMode, MutAltMode, ParRefModes, ParAltModes, HashCounts, HashCountsOG, PossibleAltKmer);
				int SupportingHashes = GetSupportingHashCount(i, alt,  reff);
				string Genotype = ShittyGenotyper(MutAltMode, MutRefMode); 
				string CompressedVarType = compressVar(varType, Positions[startPos], StructCall); 
					cout <<  chr << "\t" << pos+i << "\t" << CompressedVarType /*"."*/ << "\t" << reff << "\t" << alt << "\t" << SupportingHashes << "\t" << varType << "\t" << "." << "\t" << "." << "\t" << "." << endl;
				////////////////generatre parent genotypes and check///////////////////////
				vector <string> ParGenotypes; 
				for (int p = 0; p< ParRefModes.size(); p++)
				{
					ParGenotypes.push_back(ShittyGenotyper(ParAltModes[p], ParRefModes[p]) );
				}
				
				cout << endl;
				////////////////check that parents have enough coverage////////////////////
				cout << "PAR LOW COV CHECK" << endl; 
				int NumLowCov = 0; 
				int low = i-HashSize-10; 
				if (low < 0)
				  low = 0; 

				for(int k = low ; k <= i+10 and k < hashes.size(); k++)
				{
	       				for (int j = 0; j < parentCounts.size(); j++)
					{
						int sum = 0; 
						if (hashesRef[k] == hashes[k])
						{sum = parentCountsReference[j][k]; cout <<hashesRef[k]<< "\t" << hashes[k] << "\tParSame-" << sum; }
						else
						{sum = parentCounts[j][k] + parentCountsReference[j][k]; cout <<hashesRef[k]<< "\t" << hashes[k] << "\tParDiff-" << sum;  }
						
						if (sum <= 7 and parentCounts[j][k] + parentCountsReference[j][k] > 2 )
						{
							NumLowCov++;
							 cout << "\tLOWCOV-" << NumLowCov ;
						}
					}
					cout << endl;

				}
				//////////////////////////////check if the parenst contain any of mut hashes//////////////////////////////////////////
				bool LowCov = false; 
				int lowCount = 0;
				low = i - HashSize ; 
				if (low < 0){low = 0;}
				cout << "checking bases " << low << " to " << i+size+5 << endl;
				vector<int> streak; 
				for (int k = 0; k < parentCounts.size(); k++)
				{
					streak.push_back(0); 
				}
				for(int j = low; j <= i+size and j < hashes.size(); j++)
				{

					if (hashesRef[j] != hashes[j])
					{
						 cout << "Checking Par Hash " << hashes[j] << "\t" << hashesRef[j]; 
					  	if ((ExcludeHashes[HashToLong(hashes[j])]<1 && ExcludeHashes[HashToLong(RevComp(hashes[j]))]<1))
						{
							for (int k = 0; k < parentCounts.size(); k++)
			       		  		{
								//cout << "Checking Par Hash " << hashes[j] << "\t" << parentCounts[k][j]  << "\t" << hashesRef[j] << "\t" << parentCountsReference[k][j];
								cout << "\t" << parentCounts[k][j] << "\t" << parentCountsReference[k][j]; 
								float varFreq = 1;
	
								if (parentCountsReference[k][j] > 0)
								{
								  varFreq = (double)parentCounts[k][j]/((double)parentCountsReference[k][j] + (double)parentCounts[k][j]); 
								}
								cout << "\tvarFreq=" << varFreq; 
								
								if (parentCounts[k][j] >= 1 and parentCounts[k][j] <= ParLowCovThreshold and varFreq > .02 )//and parentCountsReference[k][j]<150 ) //if (parentCounts[k][j] <= 5 and parentCounts[k][j] > 0 )
								{
									streak[k]++; 
									cout << " LC HASH FOUND streak = " << streak[k] ;
									if (streak[k] >=3)
									{
										cout << " LC STREAK FOUND streak = " <<  streak[k] ; 
			       			     				LowCov = true;
							    			lowCount++;
									}
								}
								else
								{
									streak[k] = 0;
									cout << "streak " << k << " reset to 0 - " << streak[k] ; 
								}
							}
					  	}
						else
						{ cout << "      HASH FOUND IN REF " << hashes[j]  << " - " << ExcludeHashes[HashToLong(hashes[j])]  << " - " << ExcludeHashes[HashToLong(RevComp(hashes[j]))] ;}
						cout << endl; 		
					}
					else
						cout << "hashes are same, skpping " << hashes[j] << "\t" << hashesRef[j];
				}

				cout << "done with LC check, LC = " << LowCov << " and LowCount = " << lowCount << endl; 
				///////////////////final filter check/////////////////////////////////////////
				string Filter = ".";
				string InfoFilter = ""; 
 				if (Genotype.find("1") == std::string::npos) {
					Denovo = "Mosaic";
				} 
				if (AlignmentSegments > SegThreshold or AlignmentSegmentsCigar > SegThreshold)
				{
				  	Denovo = "PoorAlignment"; 
					stringstream ss;
					ss << AlignmentSegments << "-" << AlignmentSegmentsCigar; 
					Denovo += ss.str();
					if (Filter == ".")
						Filter = "";  
					Filter+="PA";
					InfoFilter+="PA";
					InfoFilter+=ss.str();
					InfoFilter+=",";
					Filter+=";";
				}
				if (NumLowCov > 3)
				{
					Denovo = "ParLowCovRegion";
					cout << "ParLowCov " << NumLowCov << endl;
					if (Filter == ".")
						Filter = "";
					stringstream ss;
					ss << NumLowCov; 
					Filter+="PLC";
					InfoFilter+="PLC"; 
					InfoFilter+=ss.str();
					Filter+=";";  
					InfoFilter+=",";
				}
				//if (LowCov)
				if (lowCount >=2)
				{
					cout << "LOW COVERAGE" << endl;
					Denovo = "Inherited";
					stringstream ss;
					ss << lowCount;
					Denovo += ss.str(); 
					if (Filter == ".")
						Filter = "";
					Filter+="LCH";
					InfoFilter+="LCH";
					InfoFilter+=ss.str();
					Filter+=";";
					InfoFilter+=",";
				}
				else
				   	cout << "GOOD COVERAGE" << endl;
				if (StrandBias >= 0){

					if (StrandBias >0.99999 or StrandBias < 0.00001)
					{	
						Denovo = "StrandBias";
						stringstream ss;
						ss << StrandBias; 
						if (Filter == ".")
							Filter = "";
						Filter+="SB";
						InfoFilter+="SB";
						InfoFilter+=ss.str();
						Filter+=";";
						InfoFilter+=",";
					}
				}	 
				if (Denovo == "DeNovo" and Filter == ".")
					Filter = "PASS";
				if (InfoFilter=="")
					InfoFilter="PASS"; 


				//checking if the aligment is a split that shows a good event to remove it from consideratin for the less likely sv stuff
				cout << "ok lets see if this is a big event" << endl; 
				
				///if (Denovo == "DeNovo")
				{
					cout << "yup DeNovo" << endl; 
					if (isSplitRead > 0)
					{
						cout << "yup more than one alignment" << endl; 
						if ( varType.find("D") != std::string::npos  || varType.find("Y") != std::string::npos || varType.find("I") != std::string::npos )
						{
							cout << "YAAAAAY found a deletion, let me se these SV's" << endl;
							for (int w = 0; w < alignments.size(); w++)
							{
								reads[alignments[w]].SVeventid = -1; 
							}
						}
					}
				}
				for (int p = 0; p< ParRefModes.size(); p++)
				{
					//if (ParGenotypes[p].find("1") != std::string::npos) {
						 //Denovo = "PresentInParents"; }
				}
				/////////////////////////////////////////////////////// 
				 cout << "startpos = " << startPos << " chrsize = " << ChrPositions.size() << endl; 	
				 cout << ChrPositions[startPos] << "\t" << endl;
				 cout << Positions[startPos] << "\t"  << endl;
				 cout << CompressedVarType <<"-"  << endl;
				 cout << Denovo /*"."*/  << "\t"  << endl;
				 cout << reff << "\t"  << endl;
				 cout << alt << "\t"  << endl;
				 cout << SupportingHashes << "\t"  << endl;
				 cout << Filter << "\t"  << endl;
				 cout << StructCall  << endl;
				 cout <<"RN=" << name  << endl;
				 cout << ";MQ=" << mapQual  << endl;
				 cout << ";cigar=" << cigar  << endl;
				 cout << ";" << "CVT=" << CompressedVarType << ";HD=" << endl;

				double Score = ((double)SupportingHashes/(double)PossibleAltKmer) * 100.0; 
			////////////////////////Writing var out to file/////////////////////////
				cout       << ChrPositions[startPos] << "\t" <<Positions[startPos] << "\t" << CompressedVarType <<"-" <<Denovo /*"."*/  << "\t" << reff << "\t" << alt << "\t" << SupportingHashes << "\t" << Filter << "\t" << StructCall <<"FEX=" << InfoFilter << ";RN=" << name << ";MQ=" << mapQual << ";cigar=" << cigar << ";" << "CVT=" << CompressedVarType << ";HD="; 
				
				VCFOutFile << ChrPositions[startPos] << "\t" <<Positions[startPos] << "\t" << CompressedVarType <<"-" << Denovo /*"."*/  << "\t" << reff << "\t" << alt << "\t" << Score << /*SupportingHashes << "/" << PossibleAltKmer << */ "\t" << Filter << "\t" << "PH=" << phase << ";FEX=" << InfoFilter << ";FS=" << SupportingHashes << "/" << PossibleAltKmer << ";RN=" << name << ";MQ=" << mapQual << ";cigar=" << cigar << ";SB=" << StrandBias << ";" << "AS=" << AlignmentSegments << "-" << AlignmentSegmentsCigar << ";" << "CVT=" << CompressedVarType << ";HD="; 
				//VCFOutFile << ChrPositions[startPos] << "\t" <<Positions[startPos] << "\t" << CompressedVarType <<"-" << Denovo /*"."*/  << "\t" << reff << "\t" << alt << "\t" << Score << /*SupportingHashes << "/" << PossibleAltKmer << */ "\t" << Filter << "\t" << StructCall << "PH=" << phase << ";FEX=" << InfoFilter << ";FS=" << SupportingHashes << "/" << PossibleAltKmer << ";RN=" << name << ";MQ=" << mapQual << ";cigar=" << cigar << ";SB=" << StrandBias << ";" << "AS=" << AlignmentSegments << "-" << AlignmentSegmentsCigar << ";" << "CVT=" << CompressedVarType << ";HD="; 

				for (int j = 0; j < HashCounts.size(); j++) 
				{	
					cout       << HashCounts[j] << "_";
					VCFOutFile << HashCounts[j] << "_"; 
				}
				VCFOutFile << ";AO=" << MutAltMode;
				cout       << ";VT=" <<  varType << "\t" ;
				VCFOutFile <<  ";VT=" <<  varType << "\t" ;
				cout       << "GT:DP:RO:AO" << "\t" << Genotype << ":" << MutRefMode + MutAltMode << ":" << MutRefMode << ":" << MutAltMode ;
				VCFOutFile << "GT:DP:RO:AO" << "\t" << Genotype << ":" << MutRefMode + MutAltMode << ":" << MutRefMode << ":" << MutAltMode ; 
				int ParentMode; 
				int lowC = CheckParentCov(ParentMode); 
				//      cout << ":" << lowC << ":" << ParentMode << ":" << StrandBias;; 
				//VCFOutFile << ":" << lowC << ":" << ParentMode << ":" << StrandBias;;
				for (int p = 0; p< ParRefModes.size(); p++)
				{
					VCFOutFile << "\t" << ParGenotypes[p] << ":" << ParAltModes[p] + ParRefModes[p] << ":" << ParRefModes[p] << ":" << ParAltModes[p]; 
				}
				cout << endl; 
				VCFOutFile << endl;
					
				BEDOutFile << chr << "\t" << pos+i << "\t" <<  pos+i+size << "\t" << chr << ":" << pos+i << ":" << (int)(reff.length() - alt.length()) << ":" << SupportingHashes << endl;
				i+=size;
				

				cout << "\nModes\t" << ChrPositions[startPos] << "\t" << Positions[startPos] << "\t" << CompressedVarType <<"-" <<Denovo /*"."*/  << "\t" << reff << "\t" << alt <<"\t" << MutRefMode <<"\t" << MutAltMode;
				for(int pi = 0; pi<ParentHashes.size(); pi++){
					 cout << "\t" << ParRefModes[pi];
					 }
				cout << endl;
			} 
		}
	}
		cout << "Out parsingMutations" << endl;
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
	isSplitRead= 0; 
	for (int i = 0; i<cigarString.length(); i++)
	{
		if (cigarString.c_str()[i] == 'H' || cigarString.c_str()[i] == 'S')
			isSplitRead++; 
	}
	
}
void SamRead::FixTandemRef()
{
//	cout << "FOUND TANDEM" << endl;
//	write();
	//writeVertical();  
	string lastChr = "nope"; 
	int lastPos = -1; 
	string NewRef = ""; 
	for (int i = 0; i<seq.size(); i++){
//		cout << i << " " << seq[i] << " " <<  cigarString[i] << " " << RefSeq[i] << endl;
		if (cigarString[i]=='Y'){
			NewRef+= toupper(Reff.getSubSequence(lastChr, lastPos+1, 1).c_str()[0]); 
			lastPos++; 
		}
		else{
			NewRef+=RefSeq[i]; 
			lastChr = ChrPositions[i]; 
			lastPos = Positions[i];
		}
	}
	RefSeq = NewRef; 
//	cout << "done tandtem fix" << endl;
//	write();
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
		{
			if (NewCigar.c_str()[i] == 'D')
				qual += lastQ; //'!';
			else
				qual += '!';
		} 
		else 
		{
			qual+=NewQual.c_str()[i];
			lastQ = NewQual.c_str()[i];
		}
	}
	for (int i = 0; i<qual.size(); i++)
	{
		strand+="+"; 
	}
	Positions = NewPositions;
	ChrPositions = NewChromosome;
	//************************Lookup Kmer counts ***************************//
	LookUpKmers();
	vector <long> blank;
	
			

	//**********************************************************************//
	CountAlignmentSegments(); 	
	CountAlignmentSegmentsCigar();
//	cout << "After getRefSeq";
	//FullOutwriteVertical(); 
}

void SamRead::LookUpKmers() 
{
//	cout << "SeqSize = " << seq.size() << " RefSize = " << RefSeq.size() << endl;
	vector <long> blank;
	RefAltCounts.clear();
	RefRefCounts.clear();
	MutHashListCounts.clear();
	MutContigCounts.clear(); 
	MutAltCounts.clear();
	MutRefCounts.clear();
	RefKmers.clear();
	AltKmers.clear();
	for (int pi = 0; pi < ParentHashes.size(); pi++){
		RefAltCounts.push_back(blank);
		RefRefCounts.push_back(blank);
	}
	for (int j = 0; j<seq.size(); j++)
	{
		//hash table already has revcomp counted in it
		string hash = getHash(seq, j, HashSize);
		string Refhash = getHash(RefSeq, j, HashSize);
		RefKmers.push_back(Refhash); 
		AltKmers.push_back(hash);
		if (hash != ""){
			if (MutantHashes.count(HashToLong(hash)) > 0 ){
				if (hash == Refhash)
					MutContigCounts.push_back(MutantHashes[HashToLong(hash)]*-1);
				else
					MutContigCounts.push_back(MutantHashes[HashToLong(hash)]);
			}
			else
			{
				MutContigCounts.push_back(0);
			}


			if (hash == Refhash){
				MutAltCounts.push_back(0);
				for (int pi = 0; pi < ParentHashes.size(); pi++){
					RefAltCounts[pi].push_back(0);
				}
			}
			else{
				if (MutantHashes.count(HashToLong(hash)) > 0 ){
					MutAltCounts.push_back(MutantHashes[HashToLong(hash)]);
			       	} 
				else
					MutAltCounts.push_back(-1);
				
				for (int pi = 0; pi < ParentHashes.size(); pi++){
					if (ParentHashes[pi].count(HashToLong(hash)) > 0){
						RefAltCounts[pi].push_back(ParentHashes[pi][HashToLong(hash)]);
					}
					else
						RefAltCounts[pi].push_back(-1);
				}
			}

			if (Hash.count(hash) > 0){
					MutHashListCounts.push_back(Hash[hash]);
				}
				else
					MutHashListCounts.push_back(-1);


		}
		else{
			MutContigCounts.push_back(-3); 
			MutAltCounts.push_back(-3);
			MutHashListCounts.push_back(-3); 
			for (int pi = 0; pi < ParentHashes.size(); pi++){
				RefAltCounts[pi].push_back(-3);
			}
		}

		if(Refhash != ""){
			if (MutantHashes.count(HashToLong(Refhash)) > 0 )
			{
				MutRefCounts.push_back(MutantHashes[HashToLong(Refhash)]);
			}
			else
				MutRefCounts.push_back(-1);

			for (int pi = 0; pi < ParentHashes.size(); pi++){
				if (ParentHashes[pi].count(HashToLong(Refhash)) > 0){
						RefRefCounts[pi].push_back(ParentHashes[pi][HashToLong(Refhash)]);
				}
				else
					RefRefCounts[pi].push_back(-1);
			}

		}
		else{

			MutRefCounts.push_back(-3);
			for (int pi = 0; pi < ParentHashes.size(); pi++){
				RefRefCounts[pi].push_back(-3);
			}
		}

	}
	//cout << "donezo" << endl;
}
void SamRead::parse(string read)
{
  	
	cout << "parsing " << read << endl;
	vector <string> temp = Split(read, '\t');
	cout << "boom" << endl;
	name = temp[0];
	cout << "name " << name << endl;  
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
	
	UsedForBigVar = false;UsedForBigVar = false;  
	first = true; 
	combined = false; 
	string NewSeq = ""; 
	for (int i = 0; i < seq.size(); i++)
		NewSeq+=toupper(seq.c_str()[i]);

	seq = NewSeq; 
	
	processCigar();
	cout << "working on strand bias" << endl; 
	cout << temp[0] << endl; 
	vector <string> temp2 = Split(name, ':');
	if (temp2.size() >= 2){
		cout << "break it down " << endl;
		strands = temp2[1];
		cout << "strands = " << strands << endl;
		forward = 0;
		reverse = 0; 
		forward = atoi(temp2[1].c_str()); 
		reverse = atoi(temp2[2].c_str()); 
		cout << "forward = " << forward << endl; 
		cout << "reverse = " << reverse << endl; 
		if (forward+reverse ==0)
			StrandBias = 1; 
		else
			StrandBias = ((float)forward)/((float)forward+(float)reverse); 
		cout << "strand bias = " << StrandBias << endl; 
	}
	else{
		cout << "no strand data " << temp2.size() << endl; 
		strands = ""; 
		StrandBias = -1; 
		forward = -1; 
		reverse = -1;
	}
	AlignScore = 0;
	for (int i = 11; i< temp.size(); i++){
	  vector<string> astemp = Split(temp[i], ':');
	  if (astemp[0] == "AS"){
	    AlignScore = atoi(astemp[2].c_str());
	  }
	}
	cout << "getting flag bits " << endl;
	for (int j = 0;  j < 16;  ++j){
		FlagBits [j] =  0 != (flag & (1 << j));
	}
	cout << "Read Pared = " << FlagBits[0] << endl;
	cout << "read mapped in proper pair = " << FlagBits[1] << endl; 
	cout << "read unmapped = " << FlagBits[2] << endl; 
	cout << "mate unmapped = " << FlagBits[3] << endl; 
	cout << "read reverse strand = " << FlagBits[4] << endl; 
	cout << "mate referse strand = " << FlagBits[5] << endl; 
	cout << "first in pair =" << FlagBits[6] << endl; 
	cout << "second in pair =" << FlagBits[7] << endl; 
	cout << "not primary alignment =" << FlagBits[8] << endl; 
	cout << "read fails platform or vendor quality checks =" << FlagBits[9] << endl; 
	cout << "read is PCR or optical duplicate =" << FlagBits[10] << endl; 
	cout << "supplementary alignment =" << FlagBits[11] << endl; 
}
int findBreak(SamRead& read)
{
	char Afirst = read.cigarString.c_str()[0];


//	cout << "Afirst = " << Afirst  << endl;
//	cout << "starting A check " << endl;
	if (Afirst == 'H' or Afirst == 'S')
	{
//		cout << "forward" << endl;
		for (int i =0; i < read.seq.size(); i++)
		{
			if  (read.cigarString.c_str()[i] == 'H' or read.cigarString.c_str()[i] == 'S')
			{
				//keep going
//				 cout << i <<  " == " << read.cigarString.c_str()[i] << ' ' << read.seq.c_str()[i]<< endl;
			}
			else
			{
//				cout << i <<  " == " << read.cigarString.c_str()[i]  << ' ' << read.seq.c_str()[i] << endl;
//				cout << "fond break at " << i << endl;
				return i;
			}
		}
	}
	else
	{
//		cout << "reverse" << endl;
		for (int i = read.seq.size()-1; i >= 0; i += -1)
		{
			if  (read.cigarString.c_str()[i] == 'H' or read.cigarString.c_str()[i] == 'S')
			{
//				cout << i <<  " == " << read.cigarString.c_str()[i]  << ' ' << read.seq.c_str()[i] << endl;
				//keep going
			}
			else
			{
//				cout << i <<  " == " << read.cigarString.c_str()[i]  << ' ' << read.seq.c_str()[i] << endl;
//				cout << "fond break at " << i << endl;
				return i;
		       }
		}
	}


}
SamRead BetterWay(vector<SamRead> reads)
{
	cout << "In BetterWay" << endl;
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
	for (int i =0; i < base; i++)
	{

		NewCigar += reads[A].cigarString.c_str()[i];
					NewSeq +=reads[A].seq.c_str()[i];
					NewQual += reads[A].qual.c_str()[i];
					NewRef+= reads[A].RefSeq.c_str()[i];
					NewPos.push_back(reads[A].Positions[i]);
					NewChr.push_back(reads[A].ChrPositions[i]);
	}
	//corect qualites so everyone has the same ones, H will produce no quality 
	cout << "checking quals" << endl;
	string bestQual = reads[0].qual; 
	for (int i =0; i<reads.size(); i++)
	{	
		bool h = false; 
		cout << "read " << reads[i].name << endl;
		for (int j = 0; j < reads[i].RefSeq.size(); j ++)
		{
			if (reads[i].RefSeq.c_str()[j] == 'H')
			{
				h = true;
			}
		}
		if (h)
		{
			cout << "contains H" << reads[i].RefSeq << endl;
		}
		else
		{
			cout << "does not contain H " << reads[i].RefSeq << endl;
			bestQual = reads[i].qual; 
		}
		 
		
						
	}
	for (int i =0; i<reads.size(); i++)
	{
		reads[i].qual = bestQual; 
	}
		
	for (int i =0; i< reads.size(); i++)
	{
		reads[i].createPeakMap();
	}

	cout << "Post adjustment" << endl;
	reads[A].write();
	reads[B].write();
//	for (int i =0; i< reads[A].PeakMap.size(); i++)
//	      cout << i << "-" << reads[A].PeakMap[i] << "-" << reads[B].PeakMap[i] << endl;
	cout << "************INTO**********" << endl;
	if (B == 1 and GetReadOrientation(reads[A].flag) == GetReadOrientation(reads[B].flag)) //if they are on the same strand, and there are only two split reads  
	{
//	  	cout << "same strand and only split reads" << endl;
		for (int i =base; i < reads[A].seq.size(); i++)
		{
//			cout << "base " << i << " of " << reads[A].seq.size() << " or " << reads[B].seq.size() << endl;  
			if (reads[A].Positions[i] > -1) //if this base is aligned in A 
			{
//				 cout << A << " is aligned " << reads[A].seq.c_str()[i] << " " << reads[A].Positions[i] << " " <<  reads[A].chr << endl; 
				if (reads[A].Positions[i] - LastAlignedPos > 1) //indicates a deletion
				{
					if (LastAlignedQ == '!' or reads[A].qual.c_str()[i] == '!' )
					{
//						cout << "well fuck this shit A" << endl;
						//return reads[A]; 
					}
					//if(reads[A].ChrPositions[i] == LastAlignedChr and abs(reads[A].Positions[i] -LastAlignedPos ) < MaxVarentSize ) //must be on the same chromosome
					if(reads[A].chr == reads[B].chr and abs(reads[A].Positions[i] -LastAlignedPos ) <= (MaxVarentSize+1000) ) //must be on the same chromosome
					{
//						cout << "wel this dosnt make any sense" << endl; //reads are in order in the bam so A should always be downstream of B, theus the deletion shoould be detected in B
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
						if( reads[A].chr == reads[B].chr and abs(reads[A].Positions[i] -LastAlignedPos ) >= MaxVarentSize + 1000 )
						{
							int Abreak = findBreak(reads[A]);
							int Bbreak = findBreak(reads[B]);
//							cout << " if( "<<reads[A].PeakMap[Abreak]<<" == 1 or "<<reads[A].PeakMap[Abreak-1]<<" == 1 or " << reads[B].PeakMap[Bbreak]<<" == 1 or "<<reads[B].PeakMap[Bbreak -1]<<" == 1)";
							if( /*1==1 or*/ (reads[A].PeakMap[Abreak] == 1 or reads[A].PeakMap[Abreak-1] == 1) and ( reads[B].PeakMap[Bbreak] == 1 or reads[B].PeakMap[Bbreak -1] == 1) and Abreak > 0 and Bbreak > 0)
							{
								cout <<  "INVERSION written to file"  << endl;
								Translocations << "Too Big, Same strand and chr "<< abs(reads[A].Positions[i] -LastAlignedPos ) << endl;
								reads[A].writetofile(Translocations);
								reads[B].writetofile(Translocations);
								Translocations << endl << endl;
								Translocationsbed << reads[A].chr << "\t" << reads[A].Positions[Abreak]-200 << "\t" <<  reads[A].Positions[Abreak]+200 << endl <<  reads[B].chr << "\t" << reads[B].Positions[Bbreak]-200 << "\t" <<  reads[B].Positions[Bbreak]+200 << endl; 
						       }
							else
							{
								cout <<  "INVERSION skipped"  << endl;
							}
						//	if( /*1==1 or*/ reads[A].PeakMap[i] == 1 or reads[A].PeakMap[i-1] == 1 or reads[B].PeakMap[i] == 1 or reads[B].PeakMap[i -1] == 1)
						//	{
						//		Translocations << "TOO BIG " << abs(reads[A].Positions[i] -LastAlignedPos ) << endl;
						//		reads[A].writetofile(Translocations);
						//		reads[B].writetofile(Translocations);
						//		Translocations << endl << endl;	
						//	}
						}
						else
						{
							if ((reads[A].chr == "hs37d5" and reads[B].chr != "hs37d5" ) or  (reads[A].chr != "hs37d5" and reads[B].chr == "hs37d5" ))
							{
								int Abreak = findBreak(reads[A]);
								int Bbreak = findBreak(reads[B]);
//								cout << " if( "<<reads[A].PeakMap[Abreak]<<" == 1 or "<<reads[A].PeakMap[Abreak-1]<<" == 1 or " << reads[B].PeakMap[Bbreak]<<" == 1 or "<<reads[B].PeakMap[Bbreak -1]<<" == 1)";
								if( /*1==1 or*/ (reads[A].PeakMap[Abreak] == 1 or reads[A].PeakMap[Abreak-1] == 1) and ( reads[B].PeakMap[Bbreak] == 1 or reads[B].PeakMap[Bbreak -1] == 1) and Abreak > 0 and Bbreak > 0 )
								{
									cout <<  "INVERSION written to file"  << endl;
									Translocations << "Possible mob event "<< abs(reads[A].Positions[i] -LastAlignedPos ) << endl;
									reads[A].writetofile(Translocations);
									reads[B].writetofile(Translocations);
									Translocations << endl << endl;
									Translocationsbed << reads[A].chr << "\t" << reads[A].Positions[Abreak]-200 << "\t" <<  reads[A].Positions[Abreak]+200 << endl <<  reads[B].chr << "\t" << reads[B].Positions[Bbreak]-200 << "\t" <<  reads[B].Positions[Bbreak]+200 << endl;
								}
								else
								{
									cout <<  "INVERSION skipped"  << endl;
								}
						//		if( /*1==1 or*/ reads[A].PeakMap[i] == 1 or reads[A].PeakMap[i-1] == 1 or reads[B].PeakMap[i] == 1 or reads[B].PeakMap[i -1] == 1)
						//		{
						//			Translocations << "mobil elemnt " << endl;
						//			reads[A].writetofile(Translocations); 
						//			reads[B].writetofile(Translocations); 
						//			Translocations << endl << endl;
						//		}
							}
							else
							{
							int Abreak = findBreak(reads[A]);
							int Bbreak = findBreak(reads[B]);
//							cout << " if( "<<reads[A].PeakMap[Abreak]<<" == 1 or "<<reads[A].PeakMap[Abreak-1]<<" == 1 or " << reads[B].PeakMap[Bbreak]<<" == 1 or "<<reads[B].PeakMap[Bbreak -1]<<" == 1)";
							if( /*1==1 or*/ (reads[A].PeakMap[Abreak] == 1 or reads[A].PeakMap[Abreak-1] == 1) and ( reads[B].PeakMap[Bbreak] == 1 or reads[B].PeakMap[Bbreak -1] == 1) and Abreak > 0 and Bbreak > 0)
							{
								cout <<  "INVERSION written to file"  << endl;
								Translocations << "Translocataion, same strand "<< abs(reads[A].Positions[i] -LastAlignedPos ) << endl;
								reads[A].writetofile(Translocations);
								reads[B].writetofile(Translocations);
								Translocations << endl << endl;
								Translocationsbed << reads[A].chr << "\t" << reads[A].Positions[Abreak]-200 << "\t" <<  reads[A].Positions[Abreak]+200 << endl <<  reads[B].chr << "\t" << reads[B].Positions[Bbreak]-200 << "\t" <<  reads[B].Positions[Bbreak]+200 << endl; 
						       }
							else
							{
								cout <<  "INVERSION skipped"  << endl;
							}	
							//	if( /*1==1 or*/ reads[A].PeakMap[i] == 1 or reads[A].PeakMap[i-1] == 1 or reads[B].PeakMap[i] == 1 or reads[B].PeakMap[i -1] == 1)
							//	{
							//		Translocations << "we got a translocation" << endl;
							//		reads[A].writetofile(Translocations);
							//		reads[B].writetofile(Translocations);
							//		Translocations << endl << endl;
							//	}
							}	
						}
					}
				}
				else if (reads[A].Positions[i] -LastAlignedPos < 0 and abs(reads[A].Positions[i] -LastAlignedPos ) <= MaxVarentSize +1000) // indicates a possible insertion or tandem duplication
				{
					 if (LastAlignedQ == '!'  or reads[A].qual.c_str()[i] == '!')
					{
						cout << "well fuck this shit B" << endl;
						//return reads[A];
					}
//					cout << "this could be one A, last = " << LastAlignedPos << " Current = " << reads[A].Positions[i] << " chr = " << LastAlignedChr << " and " << reads[A].ChrPositions[i] << endl;
//					if (reads[A].ChrPositions[i] == LastAlignedChr )
					if (reads[A].chr == reads[B].chr )
					{
						cout << "This is an insertion A at base " << i <<  endl;
						BEDBigStuff << reads[A].chr << "\t" << reads[A].Positions[i] << "\t" << LastAlignedPos << "\tTandemDup" << endl;
						cout << "tadem dup" << endl;
						int j = 0; 
						for( j = i; j<reads[A].seq.size() and reads[A].Positions[j] < LastAlignedPos ; j++)
						{
							NewCigar += 'Y'; //'I';
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
							NewCigar += 'Y'; //'I';
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
							int Abreak = findBreak(reads[A]);
							int Bbreak = findBreak(reads[B]);
//							cout << " if( "<<reads[A].PeakMap[Abreak]<<" == 1 or "<<reads[A].PeakMap[Abreak-1]<<" == 1 or " << reads[B].PeakMap[Bbreak]<<" == 1 or "<<reads[B].PeakMap[Bbreak -1]<<" == 1)";
							if( /*1==1 or*/ (reads[A].PeakMap[Abreak] == 1 or reads[A].PeakMap[Abreak-1] == 1) and ( reads[B].PeakMap[Bbreak] == 1 or reads[B].PeakMap[Bbreak -1] == 1) and Abreak > 0 and Bbreak > 0)
							{
//								cout <<  "INVERSION written to file"  << endl;
								Translocations << "Possible mob event "<< abs(reads[A].Positions[i] -LastAlignedPos ) << endl;
								reads[A].writetofile(Translocations);
								reads[B].writetofile(Translocations);
								Translocations << endl << endl;
								Translocationsbed << reads[A].chr << "\t" << reads[A].Positions[Abreak]-200 << "\t" <<  reads[A].Positions[Abreak]+200 << endl <<  reads[B].chr << "\t" << reads[B].Positions[Bbreak]-200 << "\t" <<  reads[B].Positions[Bbreak]+200 << endl;
							}
							else
							{
//								cout <<  "INVERSION skipped"  << endl;
							}
							//if( /*1==1 or*/ reads[A].PeakMap[i] == 1 or reads[A].PeakMap[i-1] == 1 or reads[B].PeakMap[i] == 1 or reads[B].PeakMap[i -1] == 1)
							//{
							//	Translocations << "mobil elemnt " << endl;
							//	reads[A].writetofile(Translocations);
							//	reads[B].writetofile(Translocations);
							//	Translocations << endl << endl;
							//}
						}
				 		int Abreak = findBreak(reads[A]);
						int Bbreak = findBreak(reads[B]);
//					       	cout << " if( "<<reads[A].PeakMap[Abreak]<<" == 1 or "<<reads[A].PeakMap[Abreak-1]<<" == 1 or " << reads[B].PeakMap[Bbreak]<<" == 1 or "<<reads[B].PeakMap[Bbreak -1]<<" == 1)";
						if( /*1==1 or*/ (reads[A].PeakMap[Abreak] == 1 or reads[A].PeakMap[Abreak-1] == 1) and ( reads[B].PeakMap[Bbreak] == 1 or reads[B].PeakMap[Bbreak -1] == 1) and Abreak > 0 and Bbreak > 0)
						{
//							cout <<  "INVERSION written to file"  << endl;
							Translocations << "Translocation, same strand "<< abs(reads[A].Positions[i] -LastAlignedPos ) << endl;
							reads[A].writetofile(Translocations);
							reads[B].writetofile(Translocations);
							Translocations << endl << endl;
				 			Translocationsbed << reads[A].chr << "\t" << reads[A].Positions[Abreak]-200 << "\t" <<  reads[A].Positions[Abreak]+200 << endl <<  reads[B].chr << "\t" << reads[B].Positions[Bbreak]-200 << "\t" <<  reads[B].Positions[Bbreak]+200 << endl;
					       	}
						else
						{
//							cout <<  "INVERSION skipped"  << endl;
						}
				//		if( /*1==1 or*/ reads[A].PeakMap[i] == 1 or reads[A].PeakMap[i-1] == 1 or reads[B].PeakMap[i] == 1 or reads[B].PeakMap[i -1] == 1)
				//		{
				//			Translocations << "we got a translocation" << endl;
				//			reads[A].writetofile(Translocations);
				//			reads[B].writetofile(Translocations);
				//			Translocations << endl << endl;
				//		}
					}
	
				}
				else if( reads[A].Positions[i] -LastAlignedPos < 0  and abs(reads[A].Positions[i] -LastAlignedPos ) >= MaxVarentSize +1000 )
				{
					 if (LastAlignedQ == '!'  or reads[A].qual.c_str()[i] == '!')
					{
						cout << "well fuck this shit C" << endl;

					//	return reads[A];
					}
					int Abreak = findBreak(reads[A]);
					int Bbreak = findBreak(reads[B]);
//					cout << " if( "<<reads[A].PeakMap[Abreak]<<" == 1 or "<<reads[A].PeakMap[Abreak-1]<<" == 1 or " << reads[B].PeakMap[Bbreak]<<" == 1 or "<<reads[B].PeakMap[Bbreak -1]<<" == 1)";
					if( /*1==1 or*/ (reads[A].PeakMap[Abreak] == 1 or reads[A].PeakMap[Abreak-1] == 1) and ( reads[B].PeakMap[Bbreak] == 1 or reads[B].PeakMap[Bbreak -1] == 1) and Abreak > 0 and Bbreak > 0)
					{
//						cout <<  "INVERSION written to file"  << endl;
						if (reads[A].chr == reads[B].chr)
							Translocations << "TOO BIG 3 "<< abs(reads[A].Positions[i] -LastAlignedPos ) << endl;
						else 
							Translocations << "Translocation 3 "<< abs(reads[A].Positions[i] -LastAlignedPos ) << endl;
						reads[A].writetofile(Translocations);
						reads[B].writetofile(Translocations);
						Translocations << endl << endl;
						Translocationsbed << reads[A].chr << "\t" << reads[A].Positions[Abreak]-200 << "\t" <<  reads[A].Positions[Abreak]+200 << endl <<  reads[B].chr << "\t" << reads[B].Positions[Bbreak]-200 << "\t" <<  reads[B].Positions[Bbreak]+200 << endl; 
				       }
					else
					{
//						cout <<  "INVERSION skipped"  << endl;
					}
			//		if( /*1==1 or*/ reads[A].PeakMap[i] == 1 or reads[A].PeakMap[i-1] == 1 or reads[B].PeakMap[i] == 1 or reads[B].PeakMap[i -1] == 1)
			//		{
			//			Translocations << "TOO BIG " << abs(reads[A].Positions[i] -LastAlignedPos ) << endl;
			//			reads[A].writetofile(Translocations);
			//			reads[B].writetofile(Translocations);
			//			Translocations << endl << endl;
			//		}
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
//				 cout << B << " is aligned  " << reads[A].seq.c_str()[i] << " " << reads[A].Positions[i] << " " << reads[A].chr << endl;
				if (reads[B].Positions[i] - LastAlignedPos > 1)
				{
//					cout << B << " is aligned infront of A " << endl; 
					 if (LastAlignedQ == '!'  or reads[B].qual.c_str()[i] == '!')
					{
						
						cout << "well fuck this shit D" << endl;

						//return reads[A];
					}
					//if(reads[B].ChrPositions[i] == LastAlignedChr  and abs(reads[B].Positions[i] - LastAlignedPos ) < MaxVarentSize )
					cout << "yup this one " << "abs(" << reads[B].Positions[i]<< "  - " << LastAlignedPos <<" ) < " << MaxVarentSize << endl; 
					cout << abs(reads[B].Positions[i] - LastAlignedPos )  << endl; 
					if(reads[B].chr == reads[A].chr  and abs(reads[B].Positions[i] - LastAlignedPos ) <= MaxVarentSize +1000 )
					{
						cout << "striahgtup deletion, size = " <<  abs(reads[B].Positions[i] -LastAlignedPos ) << " at base " << i << " from Position " << LastAlignedPos << " to " << reads[B].Positions[i] << endl; 
						BEDBigStuff << reads[B].chr << "\t" << LastAlignedPos << "\t" << reads[B].Positions[i]  << "\t" << "Deletion" << endl;
						cout << "Inserting reff sequence from " << LastAlignedPos+1 << " to " << reads[B].Positions[i] << " = " << reads[B].Positions[i]-LastAlignedPos << endl;
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
						cout << "done inserting sequence" << endl; 
					}
					else
					{
						cout << "fond one way too big" << endl; 
						// if( reads[A].ChrPositions[i] == LastAlignedChr and abs(reads[B].Positions[i] -LastAlignedPos ) >= MaxVarentSize )
						if( reads[A].chr == reads[B].chr and abs(reads[B].Positions[i] -LastAlignedPos ) >= MaxVarentSize +1000 )
						{
							int Abreak = findBreak(reads[A]);
							int Bbreak = findBreak(reads[B]);
//							cout << " if( "<<reads[A].PeakMap[Abreak]<<" == 1 or "<<reads[A].PeakMap[Abreak-1]<<" == 1 or " << reads[B].PeakMap[Bbreak]<<" == 1 or "<<reads[B].PeakMap[Bbreak -1]<<" == 1)";
							if( /*1==1 or*/ (reads[A].PeakMap[Abreak] == 1 or reads[A].PeakMap[Abreak-1] == 1) and ( reads[B].PeakMap[Bbreak] == 1 or reads[B].PeakMap[Bbreak -1] == 1) and Abreak > 0 and Bbreak > 0)
							{
								cout <<  "INVERSION written to file"  << endl;
								Translocations << "TOO BIG 2 "<< abs(reads[A].Positions[i] -LastAlignedPos ) << endl;
								reads[A].writetofile(Translocations);
								reads[B].writetofile(Translocations);
								Translocations << endl << endl;
								Translocationsbed << reads[A].chr << "\t" << reads[A].Positions[Abreak]-200 << "\t" <<  reads[A].Positions[Abreak]+200 << endl <<  reads[B].chr << "\t" << reads[B].Positions[Bbreak]-200 << "\t" <<  reads[B].Positions[Bbreak]+200 << endl;
							}
							else
							{
								cout <<  "INVERSION skipped"  << endl;
							}
						//	if( /*1==1 or*/ reads[A].PeakMap[i] == 1 or reads[A].PeakMap[i-1] == 1 or reads[B].PeakMap[i] == 1 or reads[B].PeakMap[i -1] == 1)
						//	{
						//		Translocations << "TOO BIG " << abs(reads[A].Positions[i] -LastAlignedPos ) << endl;
						//		reads[A].writetofile(Translocations);
						//		reads[B].writetofile(Translocations);
						//		Translocations << endl << endl;
						//	}
						}
						else
						{
							if ((reads[A].chr == "hs37d5" and reads[B].chr != "hs37d5" ) or  (reads[A].chr != "hs37d5" and reads[B].chr == "hs37d5" ))
							{
								int Abreak = findBreak(reads[A]);
								int Bbreak = findBreak(reads[B]);
//								cout << " if( "<<reads[A].PeakMap[Abreak]<<" == 1 or "<<reads[A].PeakMap[Abreak-1]<<" == 1 or " << reads[B].PeakMap[Bbreak]<<" == 1 or "<<reads[B].PeakMap[Bbreak -1]<<" == 1)";
								if( /*1==1 or*/ (reads[A].PeakMap[Abreak] == 1 or reads[A].PeakMap[Abreak-1] == 1) and ( reads[B].PeakMap[Bbreak] == 1 or reads[B].PeakMap[Bbreak -1] == 1) and Abreak > 0 and Bbreak > 0)
								{
//									cout <<  "INVERSION written to file"  << endl;
									Translocations << "Possible mob event "<< abs(reads[A].Positions[i] -LastAlignedPos ) << endl;
									reads[A].writetofile(Translocations);
									reads[B].writetofile(Translocations);
									Translocations << endl << endl;
									Translocationsbed << reads[A].chr << "\t" << reads[A].Positions[Abreak]-200 << "\t" <<  reads[A].Positions[Abreak]+200 << endl <<  reads[B].chr << "\t" << reads[B].Positions[Bbreak]-200 << "\t" <<  reads[B].Positions[Bbreak]+200 << endl;
								}
								else
								{
//									cout <<  "INVERSION skipped"  << endl;
								}
							//	if( /*1==1 or*/ reads[A].PeakMap[i] == 1 or reads[A].PeakMap[i-1] == 1 or reads[B].PeakMap[i] == 1 or reads[B].PeakMap[i -1] == 1)
							//	{
							//		Translocations << "mobil elemnt " << endl;
							//		reads[A].writetofile(Translocations);
							//		reads[B].writetofile(Translocations);
							//		Translocations << endl << endl;
							//	}
							}
							else	
							{
								int Abreak = findBreak(reads[A]);
								int Bbreak = findBreak(reads[B]);
//								cout << " if( "<<reads[A].PeakMap[Abreak]<<" == 1 or "<<reads[A].PeakMap[Abreak-1]<<" == 1 or " << reads[B].PeakMap[Bbreak]<<" == 1 or "<<reads[B].PeakMap[Bbreak -1]<<" == 1)";
								if( /*1==1 or*/ (reads[A].PeakMap[Abreak] == 1 or reads[A].PeakMap[Abreak-1] == 1) and ( reads[B].PeakMap[Bbreak] == 1 or reads[B].PeakMap[Bbreak -1] == 1) and Abreak > 0 and Bbreak > 0)
								{
//									cout <<  "INVERSION written to file"  << endl;
									Translocations << "Translocation 2 "<< abs(reads[A].Positions[i] -LastAlignedPos ) << endl;
									reads[A].writetofile(Translocations);
									reads[B].writetofile(Translocations);
									Translocations << endl << endl;
									Translocationsbed << reads[A].chr << "\t" << reads[A].Positions[Abreak]-200 << "\t" <<  reads[A].Positions[Abreak]+200 << endl <<  reads[B].chr << "\t" << reads[B].Positions[Bbreak]-200 << "\t" <<  reads[B].Positions[Bbreak]+200 << endl; 
							       }
								else
								{
//									cout <<  "INVERSION skipped"  << endl;
								}

							
							//	if( /*1==1 or*/ reads[A].PeakMap[i] == 1 or reads[A].PeakMap[i-1] == 1 or reads[B].PeakMap[i] == 1 or reads[B].PeakMap[i -1] == 1)
							//	{
							//		Translocations << "we got a translocation" << endl;
							//		reads[A].writetofile(Translocations);
							//		reads[B].writetofile(Translocations);
							//		Translocations << endl << endl;
							//	}
							}
						}
					}

				}
				else if (reads[B].Positions[i] -LastAlignedPos < 0 and abs(reads[B].Positions[i] - LastAlignedPos ) < MaxVarentSize +1000 ) // indicates a possible insertion or tandem duplication
				{
					cout << "IN POSSIBLE TANDEM DUP BIT" << endl; 
					 if (LastAlignedQ == '!' or reads[B].qual.c_str()[i] == '!' )
					{
						cout << "well fuck this shit E" << endl;

					//	return reads[A];
					}
					 // cout << "this could be one B, last = " << LastAlignedPos << " Current = " << reads[B].Positions[i] << " chr = " << LastAlignedChr << " and " << reads[B].ChrPositions[i] << endl;
					//if (reads[B].ChrPositions[i] == LastAlignedChr  
					if (reads[B].chr == reads[A].chr  )
					{
//						cout << "This is an insertion B at base " << i << endl;

						BEDBigStuff << reads[B].chr << "\t" << reads[B].Positions[i] << "\t" << LastAlignedPos << "\tTandemDup" << endl;
					//	cout << "tadem dup" << endl;
						int j = 0; 
						for(j = i; j < reads[B].seq.size() and reads[B].Positions[j] <= LastAlignedPos; j++)
						{
							NewCigar += 'Y';
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
							NewCigar += 'Y';
							NewSeq += toupper(Reff.getSubSequence(reads[B].chr, j, 1).c_str()[0]);
							NewQual += '!';
							NewRef+= '-';
							NewPos.push_back(j);
							NewChr.push_back(reads[B].chr);
						}
						
					}
					else
					{
//						cout << "we got a translocation" << endl;
						if ((reads[A].chr == "hs37d5" and reads[B].chr != "hs37d5" ) or  (reads[A].chr != "hs37d5" and reads[B].chr == "hs37d5" ))
						{
//							cout << "mobil elemnt " << endl;
							//reads[A].write();
							//reads[B].write();
						}
					}

				}
				else if( reads[B].Positions[i] -LastAlignedPos < 0  and abs(reads[B].Positions[i] -LastAlignedPos ) >= MaxVarentSize +1000 )
				{
					 if (LastAlignedQ == '!' or reads[B].qual.c_str()[i] == '!')
					{
						cout << "well fuck this shit F" << endl;
						//return reads[A];
					}
				 	int Abreak = findBreak(reads[A]);
					int Bbreak = findBreak(reads[B]);



//		       		 	cout << " if( "<<reads[A].PeakMap[Abreak]<<" == 1 or "<<reads[A].PeakMap[Abreak-1]<<" == 1 or " << reads[B].PeakMap[Bbreak]<<" == 1 or "<<reads[B].PeakMap[Bbreak -1]<<" == 1)";
					if( /*1==1 or*/ (reads[A].PeakMap[Abreak] == 1 or reads[A].PeakMap[Abreak-1] == 1) and ( reads[B].PeakMap[Bbreak] == 1 or reads[B].PeakMap[Bbreak -1] == 1) and Abreak > 0 and Bbreak > 0)
					{
//						cout <<  "INVERSION written to file"  << endl;
						Translocations << "TOO BIG 1 "<< abs(reads[A].Positions[i] -LastAlignedPos ) << endl;
						reads[A].writetofile(Translocations);
						reads[B].writetofile(Translocations);
						Translocations << endl << endl;
	 					Translocationsbed << reads[A].chr << "\t" << reads[A].Positions[Abreak]-200 << "\t" <<  reads[A].Positions[Abreak]+200 << endl <<  reads[B].chr << "\t" << reads[B].Positions[Bbreak]-200 << "\t" <<  reads[B].Positions[Bbreak]+200 << endl;
			       		}
					else
					{
//				 		cout <<  "INVERSION skipped"  << endl;
					}

				//	if( /*1==1 or*/ reads[A].PeakMap[i] == 1 or reads[A].PeakMap[i-1] == 1 or reads[B].PeakMap[i] == 1 or reads[B].PeakMap[i -1] == 1)
				//       	{
				//		Translocations << "TOO BIG " << abs(reads[A].Positions[i] -LastAlignedPos ) << endl;
				//		reads[A].writetofile(Translocations);
				//		reads[B].writetofile(Translocations);
				//		Translocations << endl << endl;
				//	}
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
				//cout << "no base is aligned" << endl; 
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
	
		int UnalignedCount = 0;
		for (int i =0; i<NewCigar.size(); i++)
		{
			if (NewCigar.c_str()[i] == 'H' || NewCigar.c_str()[i] == 'S')
				UnalignedCount++;
		}
//		cout << "Unaligned Baises = " << UnalignedCount << " %= " << (double)UnalignedCount/(double)NewCigar.size() << endl;
		if (UnalignedCount < 150)
		{
			NewCigar = NewNewCigar; 
			//cout << NewSeq << endl << NewQual << endl << NewCigar << endl << NewRef << endl;
			reads[A].first = true;
//			cout <<"redoing alignemnt number" << endl;
			//int temp = reads[A].alignments[0];
			//reads[A].alignments.clear();
			//reads[A].alignments.push_back(temp);
//			cout << "alignemtns should be " << reads[A].alignments.size() << endl;
	
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
	
		}
		else
		{
//			cout << "Skipping" << endl;
		}
		cout <<	 "Done with " << reads[A].name << endl;
	}
	else //if (3==1)
	{
	  	cout << "here somehow" << endl; 
		//this needs to be the invertion stuff 
		//need to add here looing at every base 
		if( reads[A].chr == reads[B].chr )//and abs(reads[A].Positions[i] -LastAlignedPos ) <= MaxVarentSize )
		{
			//need to check each base and find the breakpoint 
			cout << "found possible Inversion " << endl;
			


			int Abreak = findBreak(reads[A]); 
			int Bbreak = findBreak(reads[B]);
			
			

//			cout << " if( "<<reads[A].PeakMap[Abreak]<<" == 1 or "<<reads[A].PeakMap[Abreak-1]<<" == 1 or " << reads[B].PeakMap[Bbreak]<<" == 1 or "<<reads[B].PeakMap[Bbreak -1]<<" == 1)";
			//if( /*1==1 or*/ reads[A].PeakMap[Abreak] == 1 or reads[A].PeakMap[Abreak-1] == 1 or reads[B].PeakMap[Bbreak] == 1 or reads[B].PeakMap[Bbreak -1] == 1) and Abreak > 0 and Bbreak > 0
			if( /*1==1 or*/ (reads[A].PeakMap[Abreak] == 1 or reads[A].PeakMap[Abreak-1] == 1) and ( reads[B].PeakMap[Bbreak] == 1 or reads[B].PeakMap[Bbreak -1] == 1) and Abreak > 0 and Bbreak > 0)
			{	
//				cout <<  "INVERSION written to file"  << endl;
				Translocations << "INVERSION"  << endl;
				reads[A].writetofile(Translocations);
				reads[B].writetofile(Translocations);
				Translocations << endl << endl;
				Translocationsbed << reads[A].chr << "\t" << reads[A].Positions[Abreak]-200 << "\t" <<  reads[A].Positions[Abreak]+200 << endl <<  reads[B].chr << "\t" << reads[B].Positions[Bbreak]-200 << "\t" <<  reads[B].Positions[Bbreak]+200 << endl;
			}
			else
			{
//				 cout <<  "INVERSION skipped"  << endl;
			}
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
//			cout << "invertion adjust string"; 
//			reads[A].write(); 
//			reads[B].write(); 	
		}
		else if ((reads[A].chr == "hs37d5" and reads[B].chr != "hs37d5" ) or  (reads[A].chr != "hs37d5" and reads[B].chr == "hs37d5" ))
		{
			 int Abreak = findBreak(reads[A]);
			int Bbreak = findBreak(reads[B]);



//			cout << " if( "<<reads[A].PeakMap[Abreak]<<" == 1 or "<<reads[A].PeakMap[Abreak-1]<<" == 1 or " << reads[B].PeakMap[Bbreak]<<" == 1 or "<<reads[B].PeakMap[Bbreak -1]<<" == 1)";
			//if( /*1==1 or*/ reads[A].PeakMap[Abreak] == 1 or reads[A].PeakMap[Abreak-1] == 1 or reads[B].PeakMap[Bbreak] == 1 or reads[B].PeakMap[Bbreak -1] == 1) and Abreak > 0 and Bbreak > 0
			 if( /*1==1 or*/ (reads[A].PeakMap[Abreak] == 1 or reads[A].PeakMap[Abreak-1] == 1) and ( reads[B].PeakMap[Bbreak] == 1 or reads[B].PeakMap[Bbreak -1] == 1) and Abreak > 0 and Bbreak > 0)
			{
				Translocations << "mobil elemnt inverted" << endl;
				reads[A].writetofile(Translocations);
				reads[B].writetofile(Translocations);
				Translocations << endl << endl;
				Translocationsbed << reads[A].chr << "\t" << reads[A].Positions[Abreak]-200 << "\t" <<  reads[A].Positions[Abreak]+200 << endl <<  reads[B].chr << "\t" << reads[B].Positions[Bbreak]-200 << "\t" <<  reads[B].Positions[Bbreak]+200 << endl;	
			}
		}
		else
		{
			 int Abreak = findBreak(reads[A]);
			int Bbreak = findBreak(reads[B]);



//			cout << " if( "<<reads[A].PeakMap[Abreak]<<" == 1 or "<<reads[A].PeakMap[Abreak-1]<<" == 1 or " << reads[B].PeakMap[Bbreak]<<" == 1 or "<<reads[B].PeakMap[Bbreak -1]<<" == 1)";
			//if( /*1==1 or*/ reads[A].PeakMap[Abreak] == 1 or reads[A].PeakMap[Abreak-1] == 1 or reads[B].PeakMap[Bbreak] == 1 or reads[B].PeakMap[Bbreak -1] == 1) and Abreak > 0 and Bbreak > 0
			 if( /*1==1 or*/ (reads[A].PeakMap[Abreak] == 1 or reads[A].PeakMap[Abreak-1] == 1) and ( reads[B].PeakMap[Bbreak] == 1 or reads[B].PeakMap[Bbreak -1] == 1) and Abreak > 0 and Bbreak > 0)
			{
				Translocations << "we got a translocation and invertion" << endl;
				reads[A].writetofile(Translocations);
				reads[B].writetofile(Translocations);
				Translocations << endl << endl;
				Translocationsbed << reads[A].chr << "\t" << reads[A].Positions[Abreak]-200 << "\t" <<  reads[A].Positions[Abreak]+200 << endl <<  reads[B].chr << "\t" << reads[B].Positions[Bbreak]-200 << "\t" <<  reads[B].Positions[Bbreak]+200 << endl;
			}	
		}

			
//		cout << "*********SKIPPING************\ndifference strands" << endl;
		BEDNotHandled << "Different strands" << endl;
		BEDNotHandled << reads[A].chr << "\t" << reads[A].pos << "\t" << reads[A].pos+reads[A].seq.size() << "\t" << reads[A].name << "\t" << reads[A].cigar << endl;
		reads[A].writetofile(BEDNotHandled);
		BEDNotHandled << reads[B].chr << "\t" << reads[B].pos << "\t" << reads[B].pos+reads[B].seq.size() << "\t" << reads[B].name << "\t" << reads[B].cigar << endl;	
		reads[B].writetofile(BEDNotHandled);
		
		BEDNotHandled << endl << endl;


		Invertions << reads[A].chr << "\t" << reads[A].pos << "\t" << reads[B].pos << "\t" << reads[B].pos- reads[A].pos << endl;

	}
	cout << "made it here" << endl; 
	//**********************Adding K-mer lookup stuff here *********************************
	
	//**************************************************************************************
	//reads[A].FixTandemRef();
	if (reads[B].phase != "none" && reads[A].phase == "none")
	{
		reads[A].phase == reads[B].phase;
	}
	cout << "starting kmer lookup" << endl; 
	reads[A].LookUpKmers();
	cout << "ReAdjustedKmers" <<endl;
	reads[A].write(); 
	//FullOutreads[A].writeVertical(); 
	return reads[A];
	cout << "Out Of SamReadBetterWay " << endl;
}
int SamRead::CheckBasesAligned()
{
  int longest = 0;
  int count = 0; 
 for ( int j = 0; j< cigarString.size(); j++)
 {
	if (cigarString.c_str()[j] != 'H' and cigarString.c_str()[j] != 'S')
	{count++; }
	else
	{
	  if (count > longest){longest = count;}
	  count = 0; 
	}
   
 }
 if (count > longest){longest = count;}
 return longest; 
}
bool SamRead::CheckEndsAlign()
{
  int StartAlign = 0;
  int j;
  for ( j = 10; j< cigarString.size(); j++)
  {
    if (cigarString.c_str()[j] != 'H' and cigarString.c_str()[j] != 'S')
    {StartAlign++;}
    else
    {break;}
  }
  int EndAlign = 0;
  int i;
  for ( i = cigarString.size()-10; i>=0; i--)
  {
    if (cigarString.c_str()[i] != 'H' and  cigarString.c_str()[i] != 'S')
    {EndAlign++;}
    else
    {break;}
  }

  if (StartAlign > 20 or EndAlign > 20)
  {
    return true; 
  }

return false; 
}


bool SamRead::StartsWithAlignAtPeak(int &pos, string &insert, int &Kdepth)
{
	//cout << "starst with " << endl;
  	///////////////////
  	int EndClip = 0;
  	int i;
  	int PeakBases = 0; 
  	for ( i = cigarString.size()-1; i>=0; i--)
  	{
  	  	if (cigarString.c_str()[i] == 'H' or cigarString.c_str()[i] == 'S')
  	  	{EndClip++;if (PeakMap[i]==1){PeakBases++;}}
  	  	else
  	  	{break;}
  	}
  	//get the inserted sequence 
  	for (int s = i; s<cigarString.size(); s++)
  	{
  	  	insert = insert+seq[s];
  	}
  	/////////////////////
  	int StartAlign = 0;
  	int j	;
  	for ( j = 0; j< cigarString.size(); j++)
  	{
    		if (cigarString.c_str()[j] != 'H' and cigarString.c_str()[j] != 'S')
    		{StartAlign++;}
    		else
  	  	{break;}
  	}
  	/////////////////////
  	//cout << "pos = " << j << " - " << Positions[j] << endl;
  	pos = Positions[j];
  	Kdepth = (int) qual.c_str()[i]-33; 
  	if (EndClip > 40 && StartAlign > 40)
  	{
    		if ((PeakMap[i-1] or PeakMap[i] or PeakMap[i+1]) and (PeakMap[j-1] or PeakMap[j] or PeakMap[j+1])or PeakBases > 10)
    		{return true;}
  	}


 	 	return false;
 	}  
bool SamRead::StartsWithAlign(int &pos, string &insert)
{
  //cout << "starst with " << endl; 
  ///////////////////
  int EndClip = 0; 
  int i; 
  for ( i = cigarString.size()-1; i>=0; i--)
  {
    if (cigarString.c_str()[i] == 'H' or cigarString.c_str()[i] == 'S')
    {EndClip++;}
    else
    {break;}
  }
  //get the inserted sequence 
  for (int s = i; s<cigarString.size(); s++)
  {
    insert = insert+seq[s]; 
  }
  /////////////////////
  int StartAlign = 0; 
  int j; 
  for ( j = 0; j< cigarString.size(); j++)
  {
    if (cigarString.c_str()[j] != 'H' and cigarString.c_str()[j] != 'S')
    {StartAlign++;}
    else
    {break;}
  } 
  /////////////////////
  //cout << "pos = " << j << " - " << Positions[j] << endl;
  pos = Positions[j]; 
  if (EndClip > 40 && StartAlign > 40)
  {
    
    {return true;}
  }


  return false; 
 }
bool SamRead::EndsWithAlignAtPeak(int &pos, string &insert, int &Kdepth)
{
  //cout << "Running EndsWithAlign" << endl;
  ////////////////////
  int EndAlign = 0;
  int i;
  for ( i = cigarString.size()-1; i>=0; i--)
  {
    if (cigarString.c_str()[i] != 'H' and  cigarString.c_str()[i] != 'S')
    {EndAlign++;}
    else
    {break;}
  }
  pos = Positions[i+1]; 
  //cout << "EndAlign = " << EndAlign << endl;
  ///////////////////
  int StartClip = 0;
  int j;
  int PeakBases = 0; 
  for ( j = 0; j< cigarString.size(); j++)
  {
    if (cigarString.c_str()[j] == 'H' or cigarString.c_str()[j] == 'S')
    {StartClip++; if (PeakMap[j]==1){PeakBases++;}}
    else
    {break;}
  }
  //cout << "StartClip = " << StartClip << endl;
  
  for (int s = 0; s<j; s++)
  {insert = insert+seq[s];}
  Kdepth = (int) qual.c_str()[i]-33;
  if (EndAlign > 40 && StartClip > 40)
  {
    if ((PeakMap[i-1] or PeakMap[i] or PeakMap[i+1]) and (PeakMap[j-1] or PeakMap[j] or PeakMap[j+1])or PeakBases > 10 )
    {return true;}
  }


  return false;
 }
bool SamRead::EndsWithAlign(int &pos, string &insert)
{
  //cout << "Running EndsWithAlign" << endl;
  ////////////////////
  int EndAlign = 0;
  int i;
  for ( i = cigarString.size()-1; i>=0; i--)
  {
    if (cigarString.c_str()[i] != 'H' and  cigarString.c_str()[i] != 'S')
    {EndAlign++;}
    else
    {break;}
  }
   pos = Positions[i+1]; 
  //cout << "EndAlign = " << EndAlign << endl;
  ///////////////////
  int StartClip = 0;
  int j;
  for ( j = 0; j< cigarString.size(); j++)
  {
    if (cigarString.c_str()[j] == 'H' or cigarString.c_str()[j] == 'S')
    {StartClip++; }
    else
    {break;}
  }
  //cout << "StartClip = " << StartClip << endl; 
  for (int s = 0; s<j; s++)
  {insert = insert+seq[s];}

  if (EndAlign > 40 && StartClip > 40)
  {
    
    {return true;}
  }


  return false;
 }

bool CheckTranslocation(SamRead read)
{
	return false; 
}
bool CheckMob(SamRead read)
{
	return false; 
}
bool CheckPolyATail(SamRead read)
{
	return false; 
}
bool CheckLargeInsert(SamRead read)
{
	return false; 
}
string InterpretTargetSize(int size)
{
	if( size == 1)
		return "I"; 
	else if(size == -1)
		return "Y"; 
	else if(size == 0)
		return ""; 
	else if(size == 2)
		return "YY"; 
	else if(size == -2)
		return "DD";
	else if(size >2)
	{
		stringstream ss;
		ss << abs(size) << "Y"; 
		return ss.str();
	}
	else if(size <-2)
	{
		stringstream ss;
		ss << abs(size) << "D";
		return ss.str();
	}
	else
		cout << "ERROR not handeled insert size" << endl; 
	return "ERROR"; 
}
bool checkMobSupAalign(vector <SamRead> R )
{
	//if (reads[i].alignments.size()> 1)
	//{
	//	bool good = false; 
	//	for (int j = 1; j< reads[i].alignments.size(); j++)
	//	{
	//		if (reads[reads[i].alignments[j]].mapQual > 30) 
	//			good = true; 
	//	}
	//	if (good)
	//		return false;
	//	else 
	//		return true; 
	//}
	//else return true; 
	return true; 
}
string GetUnalignedCenter(SamRead A, SamRead B)
{
	//need to test before you call this that the contigs are on the same strand and have the proper clipPattern
	cout << "GetUnalignedCenter" << endl; 
	A.write(); 
	B.write(); 
	bool internal = false;
	bool Afirst = false; 
	bool Bfirst = false; 
	string Return = ""; 
	if (A.seq.size() == B.seq.size())
	{
		int i = 0; 
		for(i = 0; i < A.seq.size(); i++)
		{
			if      ((A.cigarString.c_str()[i] != 'S' && A.cigarString.c_str()[i] != 'H') && (B.cigarString.c_str()[i] == 'S' || B.cigarString.c_str()[i] == 'H'))
			{
				Afirst = true;
				break; 
			}
			else if ((B.cigarString.c_str()[i] != 'S' && B.cigarString.c_str()[i] != 'H') && (A.cigarString.c_str()[i] == 'S' || A.cigarString.c_str()[i] == 'H'))
			{
				Bfirst = true; 
				break; 
			}
		}
		for( i = i; i < A.seq.size(); i++)
		{
			if (Afirst)
			{
				if ((A.cigarString.c_str()[i] != 'S' && A.cigarString.c_str()[i] != 'H') && (B.cigarString.c_str()[i] == 'S' || B.cigarString.c_str()[i] == 'H'))
				{}
				else if ((A.cigarString.c_str()[i] == 'S' || A.cigarString.c_str()[i] != 'H') && (B.cigarString.c_str()[i] == 'S' || B.cigarString.c_str()[i] == 'H'))
				{
					Return += A.seq.c_str()[i]; 
				}
				else
					return Return; 
			}
			else if (Bfirst)
			{
				if ((B.cigarString.c_str()[i] != 'S' && B.cigarString.c_str()[i] != 'H') && (A.cigarString.c_str()[i] == 'S' || A.cigarString.c_str()[i] == 'H'))
				{}
				else if ((B.cigarString.c_str()[i] == 'S' || B.cigarString.c_str()[i] != 'H') && (A.cigarString.c_str()[i] == 'S' || A.cigarString.c_str()[i] == 'H'))
				{
					Return += A.seq.c_str()[i];
				}
				else
					return Return;
			}
			else 
			{
				cout << "WARNING GetUnaligedCenter no one is first" << endl; 
				return "";
			}
		}
	}
	else
	{
		cout << "GetUnalignedCenter WARNING, seq not the same size" << endl; 
	}

	
	return ""; 
}
string InterpretInsertSize(string s)
{
	if (s.size() ==0)
		return "";
	else if (s.size() ==1)
		return "I"; 
	else if (s.size() == 2)
		return "II"; 
	else if (s.size() > 2)
	{
		stringstream a; 
		a << s.size();
		a << "I"; 
		return a.str(); 
	}
	return ""; 
}
bool MobAllA(MobRead R)
{
//	cout << "checking all A" << endl; 
	R.write(); 

	bool allA = true;
	char base = 'Z';
	for (int i = 0; i < R.seq.size(); i++)
	{
//		cout << i << "\t" << R.seq.c_str()[i] << "\t" << R.cigarString.c_str()[i] << "\t" << base; 
		if (R.cigarString.c_str()[i] != 'H' && R.cigarString.c_str()[i] != 'S')
		{
			if (base == 'Z')
				base = R.seq.c_str()[i]; 

			if (base == R.seq.c_str()[i])
			{
//				cout << "yup"; 
			}
			else
			{
//				cout << "nope"; 
				allA = false;
			}
		}
//		cout << endl; 
	}
//	cout << "true = " << true << endl; 
//	cout << "yay done and allA = " << allA << endl; 
	return allA; 


}
float AlignmentAllA(SamRead R)
{
	cout << "checking all A" << endl; 
	R.write(); 

	bool allA = true;
	char base = 'Z';
	int A = 0; 
	int T = 0; 
	float size = 0; 
	for (int i = 0; i < R.seq.size(); i++)
	{
//		cout << i << "\t" << R.seq.c_str()[i] << "\t" << R.cigarString.c_str()[i] << "\t" << base; 
		if (R.cigarString.c_str()[i] != 'H' && R.cigarString.c_str()[i] != 'S')
		{
			size = size +1; 
			if (base == 'Z')
				base = R.seq.c_str()[i]; 

			if (base == R.seq.c_str()[i])
			{
//				cout << "yup"; 
			}
			else
			{
//				cout << "nope"; 
				allA = false;
			}
			if (R.seq.c_str()[i] == 'A')
				A++; 
			else if (R.seq.c_str()[i] == 'T')
				T++; 
		}
//		cout << endl; 
	}
	
	
	cout << "true = " << true << endl; 
	cout << "yay done and allA = " << allA << endl; 
	cout << "A = " << A <<" T = " << T << " Aprop = " << (float) A / size << " Tprop = " << (float) T / size << endl; 
	if (allA)
	{ cout << "returning 1" << endl; return 1; }
	else if ( A > T)
	{cout << "returning " << (float) A / size << endl;return (float) A / size; }
	else
	{cout << "returning " << (float) T / size << endl;return (float) T / size; }
	
}
int MobAligneBases(MobRead M, SamRead R)
{
	cout << "CheckingMobBase" << endl; 
	
	SamRead Read = R; 
	if (GetReadOrientation(R.flag) != GetReadOrientation(M.flag))
	{
		cout << "flipping" << endl; 
		Read.flipRead(); 
	}
	R.write(); 
	M.write(); 
	int MobBase = 0; 
	int MD = 0; 
	int RD = 0; 
	for (int i = 0; i + RD < Read.seq.size() && i+MD < M.seq.size(); i++)
	{
		while (M.seq.c_str()[i+MD] =='-')
		{
			MD++;

		}
 		while (Read.seq.c_str()[i+RD] =='-')
		{
			RD++;

		}
		if (Read.seq.c_str()[i+RD] != M.seq.c_str()[i+MD])
			cout << "out of sync somehow " << Read.seq.c_str()[i+RD] <<"-" << Read.cigarString.c_str()[i+RD] <<  "\t" <<  M.seq.c_str()[i+MD] << "-" << M.cigarString.c_str()[i+MD] << "\t" << i+RD << " - " << i+MD ; 
		else
			cout << "yay"		  << Read.seq.c_str()[i+RD] <<"-" << Read.cigarString.c_str()[i+RD] <<  "\t" <<  M.seq.c_str()[i+MD] << "-" << M.cigarString.c_str()[i+MD] << "\t" << i+RD << " - " << i+MD ;		
		if ( (Read.cigarString.c_str()[i+RD] == 'H' || Read.cigarString.c_str()[i+RD] =='S' ) && (M.cigarString.c_str()[i+MD] !='H' && M.cigarString.c_str()[i+MD] != 'S'))
		{
			MobBase++;
			cout << "\tMOBbase" << endl;
		}
		else
			cout << "\tnone" << endl; 
	}
	Read.write();
	M.write(); 

	cout << "MobBases alibned to this read = " << MobBase << endl; 
	return MobBase; 
}
void FindFirstAndLast(vector<SamRead>& R, int& A, int& B)
{

	cout << "finding flanking aligments for this read " << R.size() << endl; 
	
	long shortest = 30000000000;
	vector<bool> considering; 
	for (int j = 0; j<R.size(); j++)
	{
		if (GetReadOrientation(R[0].flag) != GetReadOrientation(R[j].flag))
		{
			R[j].flipRead();
		}
		if (R[j].seq.size() < shortest)
			shortest = R[j].seq.size(); 
		R[j].write(); 
		if (R[j].sigBreakPoint() > 0)
			considering.push_back(true);
		else
			considering.push_back(false);
		
	}
	

	////find start alignemtn
	for (int i =0; i < shortest; i++)
	{
		for (int j = 0; j<R.size(); j++)
		{
			if (R[j].cigarString.c_str()[i] != 'H' && R[j].cigarString.c_str()[i] != 'S' && considering[j])
			{
				A = j; 
				break; 
			}
		}
		if (A != -1)
			break; 
	}
	cout << "Start alignment is " << A << endl; 
	////flip to find end 
	for (int j = 0; j<R.size(); j++)
	{
		R[j].flipRead();
	
	}

	////find end alignemtn
	
	for (int i =0; i < shortest; i++)
	{
		for (int j = 0; j<R.size(); j++)
		{
			if (R[j].cigarString.c_str()[i] != 'H' && R[j].cigarString.c_str()[i] != 'S' && considering[j])
			{
				B = j;
				break;
			}
		}
		if (B != -1)
			break;
	}
	cout << "end alignemtn is " << B << endl; 
	return; 
	
}
void LastDitch(vector<SamRead>& reads, int i, int A, int B, int& CurrentSVeventID) 
{
	int bp = reads[reads[i].alignments[A]].BreakPoint();
	int sbp = reads[reads[i].alignments[B]].BreakPoint();
	if (bp > 0)
	cout << "Passed SigBreakpoint check LastDitch" << endl; 
	cout << "bp = " << bp << " sbp = " << sbp << endl; 	
	//int start = reads[reads[i].alignments[A]].pos  + bp; 
	CurrentSVeventID++;
	for(int k = 0; k < reads[reads[i].alignments[A]].alignments.size(); k++)
	{reads[reads[reads[i].alignments[A]].alignments[k]].SVeventid = CurrentSVeventID;}
	cout << "here" << endl; 

	string GenotypeField;
	GenotypeField = reads[reads[i].alignments[A]].createStructGenotype(bp);
	stringstream Format; 	
	stringstream alt ;
	cout << "here2 " << endl; 
	Format << "OrphanBND"; 
	Format << "-LC=" << reads[reads[i].alignments[A]].SVCheckParentsForLowCov(bp); 
	string ref = Reff.getSubSequence(reads[reads[i].alignments[A]].chr, reads[reads[i].alignments[A]].pos+bp -1 ,1); 
	reads[reads[i].alignments[A]].BNDid = MaxBND+1; MaxBND++;
	reads[reads[i].alignments[B]].BNDid = MaxBND+1; MaxBND++;
	string SVDES = "";
	if ( reads[reads[i].alignments[A]].clipPattern == "mc")
	{
		cout << "here 3" << endl; 
		ref = Reff.getSubSequence(reads[reads[i].alignments[A]].chr, reads[reads[i].alignments[A]].pos+bp  -1   , 1 );
		string altseq =  Reff.getSubSequence(reads[i].chr, reads[i].pos+bp -1 -1 , 1 );
		if( reads[reads[i].alignments[A]].clipPattern == "mc" && GetReadOrientation(reads[reads[i].alignments[A]].flag) == GetReadOrientation(reads[reads[i].alignments[B]].flag))
		{
			cout << "here 4" << endl; 
			string insertseq = GetUnalignedCenter(reads[reads[i].alignments[A]], reads[reads[i].alignments[B]]);
			alt << altseq << insertseq;
			alt << "[" << reads[reads[i].alignments[B]].chr << ":" << reads[reads[i].alignments[B]].pos+sbp -1 << "[";
			Format <<  "bnd_" << reads[reads[i].alignments[A]].BNDid;
			SVDES = "Translocation";
		}
		else if ( reads[reads[i].alignments[A]].clipPattern == "mc" && GetReadOrientation(reads[reads[i].alignments[A]].flag) !=GetReadOrientation( reads[reads[i].alignments[B]].flag))
		{
			cout << "here 5" << endl; 
			SamRead temp = reads[reads[i].alignments[B]];
			temp.flipRead();
			string insertseq = GetUnalignedCenter(reads[reads[i].alignments[A]], temp);
			alt << altseq << insertseq;
			alt <<  "]" << reads[reads[i].alignments[B]].chr << ":" << reads[reads[i].alignments[B]].pos+sbp -1 << "]";
			Format << "bnd_" << reads[reads[i].alignments[A]].BNDid;
			SVDES = "InvertedTranslocation";
		}
	}
	else if ( reads[reads[i].alignments[A]].clipPattern == "cm")
	{
		cout << "here 3b" << endl; 
		ref = Reff.getSubSequence(reads[reads[i].alignments[A]].chr, reads[reads[i].alignments[A]].pos+bp -1  ,1);
		string altseq = Reff.getSubSequence(reads[reads[i].alignments[A]].chr, reads[reads[i].alignments[A]].pos+bp -1 ,1);
		if ( reads[reads[i].alignments[A]].clipPattern == "cm" && GetReadOrientation(reads[reads[i].alignments[A]].flag) == GetReadOrientation(reads[reads[i].alignments[B]].flag))
		{
			cout << "here 4b" << endl; 
			alt << "]" << reads[reads[i].alignments[B]].chr << ":" << reads[reads[i].alignments[B]].pos + sbp /*check could be +1*/ << "]" << altseq;
			Format << "bnd_" << reads[reads[i].alignments[A]].BNDid ;
			SVDES = "Translocation";
		}
		else if ( reads[reads[i].alignments[A]].clipPattern == "cm" && GetReadOrientation(reads[reads[i].alignments[A]].flag) != GetReadOrientation(reads[reads[i].alignments[B]].flag))
		{
			cout << "here 5b" << endl;
			alt << "[" << reads[reads[i].alignments[B]].chr << ":" << reads[reads[i].alignments[B]].pos+ sbp -1 /*check could be +1*/  << "[" << altseq;
			Format << "bnd_" << reads[reads[i].alignments[A]].BNDid ;
			SVDES = "InvertedTranslocation";
		}
	}
	else
	{

		 	cout << "here 4c" << endl;
                	ref = Reff.getSubSequence(reads[reads[i].alignments[A]].chr, reads[reads[i].alignments[A]].pos+bp -1  ,1);
			string altseq = Reff.getSubSequence(reads[reads[i].alignments[A]].chr, reads[reads[i].alignments[A]].pos+bp -1 ,1);
			string insertseq = GetUnalignedCenter(reads[reads[i].alignments[A]], reads[reads[i].alignments[B]]);
                        alt << altseq << insertseq;
                        alt << "[" << reads[reads[i].alignments[B]].chr << ":" << reads[reads[i].alignments[B]].pos+sbp -1 << "[";
                        Format <<  "bnd_" << reads[reads[i].alignments[A]].BNDid;
                        SVDES = "MessyTranslocations";	 
	}
	cout << "here 6" << endl; 
		
	string FullfilterA = reads[reads[i].alignments[A]].filterSV(); 
	int GMap = 0;
	int minMapQual = 30;
	if (reads[reads[i].alignments[A]].mapQual > minMapQual)
		GMap++;
				
	string InfoFilter = "";
	string Filter = "";
				cout << "here 7"<< endl; 	
	if (reads[reads[i].alignments[A]].SVCheckParentsForLowCov(reads[reads[i].alignments[A]].sigBreakPoint()) >= 1)
	{
		Format << "-Inherited";
		InfoFilter = "Inherited";
		Filter = "LCH"; 
	}
	else if (GMap < 1)
	{
		Format << "-LowMapQual";
		InfoFilter = "LowMapQual";
		Filter = "LMQ"; 
	}
	else if (FullfilterA == "" )
	{
		Format<<"-DeNovo";
		InfoFilter = "Pass";
		Filter = "PASS"; 
	}
	else
	{
		Format<<"-" << FullfilterA; 
		InfoFilter = FullfilterA; 
		Filter = "fail"; 
	}
	//make quality stuff
	int readAmut=0;
	int readApos=0;
	reads[reads[i].alignments[A]].GetQualityHashes(readAmut, readApos, bp);
	float qual = -100; 
	if ((readApos) > 0)
		qual = ((float)readAmut) / ((float)readApos) * 100.0; 
	else
		qual = 0; 
	stringstream info;	
	info << "SVTYPE=BND;MATEID=bnd_" << reads[reads[i].alignments[B]].BNDid << ";";
	string phase="none"; 
	if (reads[reads[i].alignments[A]].phase != "none")
		phase = reads[reads[i].alignments[A]].phase; 
	info << "PH=" << phase << ";"; 
	info << "FEX=" << InfoFilter << ";";
	int SupportingHashes = readAmut;
	int possibleHashes = readApos; 
	info << "FS=" << SupportingHashes << "/" <<possibleHashes << ";"; 
	info << "RN=" << reads[reads[i].alignments[A]].name << ";";
	info << "MQ=" << reads[reads[i].alignments[A]].mapQual << "_and_" << reads[reads[i].alignments[B]].mapQual << ";";
	info << "cigar=" << reads[reads[i].alignments[A]].cigar << "_and_" << reads[reads[i].alignments[B]].cigar << ";";
	info << "SB=" << reads[reads[i].alignments[A]].StrandBias  << ";"; 
	info << "AS=" << reads[reads[i].alignments[A]].AlignmentSegments << "-" << reads[reads[i].alignments[A]].AlignmentSegmentsCigar << "_and_" ;
	stringstream call;
	cout << "in this one" << endl; 
	call << reads[reads[i].alignments[A]].chr << "\t" << reads[reads[i].alignments[A]].pos+bp -1 << "\t" << Format.str() << "\t" << ref <<  "\t" << alt.str() <<  "\t" << qual  << "\t" << Filter << "\t" << info.str() <<  "\t" << "GT:DP:RO:AO\t" << GenotypeField ;
	cout << call.str() << endl ;
	VCFOutFile << call.str();
	VCFOutFile << endl; 
					
				
}
main (int argc, char *argv[])
{

cout << "###########################RUNNING THIS ONE#########################" << endl; 
cout << "Modes	chr	pos	type	reff	alt	MutRef	MutAlt	Par1Ref	Par2Ref" << endl;	
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
  -h [ --help ]  Print help message\n\
  -sam  arg  Path to input SAM file, omit for stdin\n\
  -r    arg  Path to reference file \n\
  -hf   arg  Path to HashFile from RUFUS.build\n\
  -hS   arg  Hash Size\n\
  -o    arg  Output stub\n\
  -m    arg  Maximum varient size: default 1Mb\n\
(Sorry it has to be a num, no 1kb, must be 1000\n\
  -c    arg  Path to sorted.tab file for the parent sample\n\
  -s    arg  Path to sorted.tab file for the subject sample\n\
  -cR   arg  Path to the sorted.tab file fo the parnt sample hashes in the reference\n\
  -sR   arg  Path to the sorted.tab file fo the subject sample hashes in the reference\n\
  -mQ   arg  Minimum map quality to consider varients in\n\
  -mod  arg  Path to the model file from RUFUS.model\n\
  -e    arg  Path to Kmer file to exlude from LowCov check\n\
  -mob  arg  Path to a bam file of the aligned contigs to a mobil element list\n\
  -as   arg  alignemnt segments threshold (default: 10)\n\
";
	
	string MutHashFilePath = "" ;
	string MutHashFilePathReference = "";
	//MaxVarentSize = 1000000;
	string RefFile = ""; 
	string HashListFile = "" ; 	
	string samFile = "stdin"; 
	string outStub= "";
	string ModelFilePath = "";
	string ExcludeFilePath = ""; 
	string MobBam = ""; 
	 SegThreshold = 10; 
	int MinMapQual = 0; 
	for(int i = 1; i< argc; i++)
	{
		cout << i << " = " << argv[i] << endl; 
	}
	cout <<"****************************************************************************************" << endl;
	vector <int> ParentHashFilePaths; 
	vector <int> ParentHashFilePathsReference;
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
		else if (p == "-hs")
		{
			HashSize =  atoi(argv[i+1]);
			i++;
		}
		else if (p == "-m")
		{
			MaxVarentSize =  atoi(argv[i+1]);
			i++;
			cout << "YAAAY added MaxVarSize = " << MaxVarentSize << endl;
		}
		else if (p == "-as")
		{
			SegThreshold = atoi(argv[i+1]); 
			i++;
			cout << "AlignmentSegmentThreshold = " << SegThreshold << endl; 
		}
		else if (p == "-c")
		{
			cout << "Par Hash = " << argv[i+1] << endl;
			ParentHashFilePaths.push_back(i+1);
			i=i+1;
		}
		else if (p == "-cR")
		{
			cout << "Par Ref Hash = " << argv[i+1] << endl;
			ParentHashFilePathsReference.push_back(i+1);
			i=i+1;
		}
		else if (p == "-s")
		{
			cout << "Sub Hash = " << argv[i+1] << endl;
			MutHashFilePath = argv[i+1];
			i+=1;
		}
		else if (p == "-sR")
		{
			cout << "Sub Hash = " << argv[i+1] << endl;
			MutHashFilePathReference = argv[i+1];
			i+=1;
		}
		else if (p == "-mod")
		{
			cout << "model file = " << argv[i+1] << endl;
			ModelFilePath = argv[i+1];
			i+=1;
		}
		else if (p == "-mQ")
		{
			cout << "Min Mapping Qualtiy = " << argv[i+1] << endl;
			MinMapQual = atoi(argv[i+1]);
			i+=1;
		}
		else if(p == "-e")
		{
		  	cout << "Exclue File Path = " << argv[i+1] << endl; 
			ExcludeFilePath = argv[i+1]; 
			i+=1; 
		}
		else if (p == "-mob")
		{
			cout << "Mobil Eelement aligned sam file = " << argv[i+1] << endl; 
			MobBam = argv[i+1]; 
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
	unordered_map <string, MobRead> mobs; 
	if (MobBam != "")
	{	
		cout << "mob aligned contigs provided: " << MobBam << endl; 
		cout << "reading in sam file " << endl; 
		ifstream reader; 
		cout << "1" << endl;
		reader.open(MobBam);
		cout << "2" << endl;
		string line;
		cout << "3" << endl;
		while (getline(reader, line))
		{
			if (line.c_str()[0] == '@')
			{
				cout << " HEADER LINE = " << line << endl;
			}
			else
			{
				MobRead temp;
				
				temp.parse(line);
				if (temp.chr !="*" && MobAllA(temp) == false)
				{
					if (mobs.count(temp.name) > 0)
					{
						if (mobs[temp.name].AS < temp.AS)
						{
							cout << " Writing over mob with better one " << line << endl;
							mobs[temp.name] = temp;
						}
					}
					else
					{
						cout << " Adding Mob Alignment " << line << endl;
						mobs[temp.name] = temp;
					}
				}
			}
		}
	}
	cout << "yaya finished mob " << endl; 
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
	for (int i = 0; i < ParentHashFilePathsReference.size(); i++)
	{
		ifstream reader;
		reader.open (argv[ParentHashFilePathsReference[i]]);
		string line = "";
		unordered_map <unsigned long int, int> hl;
		while (getline(reader, line))
		{
			vector <string> temp = Split(line, ' ');
			unsigned long hash = HashToLong(temp[0]);
			ParentHashes[i][hash] = atoi(temp[1].c_str());
			hash = HashToLong(RevComp(temp[0]));
			ParentHashes[i][hash] = atoi(temp[1].c_str());
		}
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
	
	reader.open (MutHashFilePathReference);
	line = "";
	while (getline(reader, line))
	{

		vector <string> temp = Split(line, ' ');
		unsigned long hash = HashToLong(temp[0]);
		MutantHashes[hash] = atoi(temp[1].c_str());
		hash = HashToLong(RevComp(temp[0]));
		MutantHashes[hash] = atoi(temp[1].c_str());
	}
	reader.close();			
	
	reader.open(ExcludeFilePath); 
	while (getline(reader, line))
	{
		vector <string> temp = Split(line, ' ');
		unsigned long hash = HashToLong(temp[0]);
		//cout << "adding " << temp[0] << " with C=" << temp[1] << endl; 
		ExcludeHashes[hash] =  atoi(temp[1].c_str());
		//hash = HashToLong(RevComp(temp[0]));
		//ExcludeHashes[hash] =  atoi(temp[1].c_str());
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

	ifstream ModelFile; 
	ModelFile.open (ModelFilePath); 
	if (ModelFile.is_open())
	{ cout << "ModelFile is open";}
	else
	{
		cout << "Error no model file given, not worring abou this now" << endl;
		//return -1;
	}

      	ifstream HashList;
      	HashList.open (HashListFile);
      	if ( HashList.is_open())
      	{	 cout << "HashList Open " << HashListFile << endl;}   //cout << "##File Opend\n";
      	else
      	{
      		cout << "Error, HashList could not be opened";
      		return -1;
      	}
	line = "";
	getline(HashList, line);
	cout << "line = " << line << endl; 
	char seperator = '\t';
	vector<string> temp = Split(line, seperator);
	if (temp.size() ==1){
		cout << "separator is not tab" << endl; 
		seperator = ' '; 
		 temp = Split(line, seperator);
	}
	else if (temp.size() == 0)
	{
		cout << "Hash List is empty, aborting" << endl; 
		return -1; 
	}
	else
		cout << "separator is tab" << endl;

	cout << "split = " << temp[0] << " and " << temp[1] << endl;
	HashSize = temp[0].size(); 
	if (temp.size() ==4)
	{ 
		HashSize = temp[3].length(); 
		Hash.insert(pair<string, int>(temp[3], atoi(temp[2].c_str())));
	
		while ( getline(HashList, line))
		{
			vector<string> temp = Split(line, seperator);
			Hash.insert(pair<string, int>(temp[3], atoi(temp[2].c_str())));
			Hash.insert(pair<string, int>(RevComp(temp[3]), atoi(temp[2].c_str()))); 
			//cout << "added pair " << temp[3] << "\t" << temp[2] << endl;
		}
		HashList.close(); 
		cout << "done with HashList" << endl;
	}
	else if (temp.size() ==2)
	{
		HashSize = temp[0].length();
		Hash.insert(pair<string, int>(temp[0], atoi(temp[1].c_str())));
		Hash.insert(pair<string, int>(RevComp(temp[0]), atoi(temp[1].c_str())));
		while ( getline(HashList, line))
		{
			vector<string> temp = Split(line, seperator);
			Hash.insert(pair<string, int>(temp[0], atoi(temp[1].c_str())));
			//cout << "added pair " << temp[3] << "\t" << temp[2] << endl;
		}
		HashList.close();
		cout << "done with HashList" << endl;
	}
	/*else if (temp.size() ==1)
	{
		vector<string> temp = Split(line, ' ');
		HashSize = temp[0].length();
		Hash.insert(pair<string, int>(temp[0], atoi(temp[1].c_str())));
		Hash.insert(pair<string, int>(RevComp(temp[0]), atoi(temp[1].c_str())));
		while ( getline(HashList, line))
		{
			vector<string> temp = Split(line, ' ');
			Hash.insert(pair<string, int>(temp[0], atoi(temp[1].c_str())));
			Hash.insert(pair<string, int>(RevComp(temp[0]), atoi(temp[1].c_str())));
			//cout << "added pair " << temp[3] << "\t" << temp[2] << endl;
		}
		HashList.close();
		cout << "done with HashList" << endl;
	}*/
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
	boom = "Intermediates/" + boom; 
	BEDBigStuff.open(boom+ ".vcf.Big.bed");
	BEDNotHandled.open(boom+ ".vcf.NotHandled.bed");
	Invertions.open(boom+".vcf.invertions.bed");
	Translocations.open(boom+ ".vcf.Translocations");
	Translocationsbed.open(boom+ ".vcf.Translocations.bed");
	Unaligned.open(boom+"vcf.Unaligned");

	//write VCF header
	VCFOutFile << "##fileformat=VCFv4.1" << endl;
	VCFOutFile << "##fileDate=" << time(0) << endl;
	VCFOutFile << "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">" << endl;
	VCFOutFile << "##FORMAT=<ID=AK,Number=1,Type=Integer,Description=\"Alternte Kmer Count\">" << endl;
	VCFOutFile << "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Total Kmer depth across the variant\">" << endl;
	VCFOutFile << "##FORMAT=<ID=RO,Number=1,Type=Integer,Description=\"Mode of reference kmer counts\">" << endl;
	VCFOutFile << "##FORMAT=<ID=AO,Number=1,Type=Integer,Description=\"Mode of alt kmer counts\">" << endl;
	VCFOutFile << "##INFO=<ID=PH,Number=1,Type=String,Description=\"If read backed phasing is possible, the name of the sample that the variant was inherited from\">" << endl;
	VCFOutFile << "##INFO=<ID=FEX,Number=1,Type=String,Description=\"Filters failed and value\">" << endl;
	VCFOutFile << "##INFO=<ID=SB,Number=1,Type=Float,Description=\"Strand Bias of the aassembled contig\">" << endl;
	VCFOutFile << "##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of SV detected\">" << endl;
	VCFOutFile << "##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=\"Length of SV detected\">" << endl; 
	VCFOutFile << "##INFO=<ID=END,Number=1,Type=Integer,Description=\"END of SV detected\">" << endl; 
	VCFOutFile << "##INFO=<ID=AO,Number=1,Type=Integer,Description=\"Alternate allele observations, with partial observations recorded fractionally\">"<<endl;
	VCFOutFile << "##INFO=<ID=HD,Number=.,Type=String,Description=\"Hash counts for each k-mer overlapping the vareint, -1 indicates no info\">"<< endl;
	VCFOutFile << "##INFO=<ID=RN,Number=1,Type=String,Description=\"Name of contig that produced the call\">"<< endl;
	VCFOutFile << "##INFO=<ID=FS,Number=1,Type=String,Description=\"Full score, supporthing kmers \\ possible varient kmers based on sequence\">"<< endl;
	VCFOutFile << "##INFO=<ID=MQ,Number=1,Type=Integer,Description=\"Mapping quality of the contig that created the call\">"<< endl;
	VCFOutFile << "##INFO=<ID=cigar,Number=1,Type=String,Description=\"Cigar string for the contig that created the call\">"<< endl;
	VCFOutFile << "##INFO=<ID=VT,Number=1,Type=String,Description=\"Varient Type\">"<< endl;	
	VCFOutFile << "##INFO=<ID=CVT,Number=1,Type=String,Description=\"Compressed Varient Type\">"<< endl;
	VCFOutFile << "##INFO=<ID=NR,Number=1,Type=Integer,Description=\"Number of total reads in target region\">" << std::endl;
	VCFOutFile << "##INFO=<ID=NH,Number=1,Type=Integer,Description=\"Number of alu heads in target region\">" << std::endl;
	VCFOutFile << "##INFO=<ID=NT,Number=1,Type=Integer,Description=\"Number of polyA tails in target region\">" << std::endl;
	VCFOutFile << "##INFO=<ID=LT,Number=1,Type=Integer,Description=\"Longest polyA tail in target region\">" << std::endl;
	VCFOutFile << "##INFO=<ID=TB,Number=1,Type=Integer,Description=\"Is tail left bound, right bound, or double bound\">" << std::endl;
	VCFOutFile << "##INFO=<ID=AS,Number=1,Type=Integer,Description=\"Number of alignment segments in the contig\">" << std::endl;
	VCFOutFile << "##INFO=<ID=MT,Number=1,Type=String,Description=\"Mobil element sequence inserted\">"<< endl;
	VCFOutFile << "##INFO=<ID=SVID,Number=1,Type=String,Description=\"Uniuqe ID given to an SV event with multiple brekends so it can be quicky identified\">"<< endl;
	VCFOutFile << "##INFO=<ID=SOURCE,Number=1,Type=String,Description=\"Location in the genome where the inserted sequence came from\">"<< endl;
	VCFOutFile << "##INFO=<ID=SVDES,Number=1,Type=String,Description=\"If available RUFUS will interpret the SV type for you\">"<< endl;
	VCFOutFile << "##INFO=<ID=MATEID,Number=1,Type=String,Description=\"If available, the id of the call that is the mate of this one\">"<< endl;
	VCFOutFile << "##FILTER=<ID=PA,Description=\"PoorAlignment\">" << std::endl;	
	VCFOutFile << "##FILTER=<ID=PLC,Description=\"Parents are at low coverage in this region, cannt be sure of genotype\">" << std::endl;
	VCFOutFile << "##FILTER=<ID=LCH,Description=\"Parents have hashes showing variant at low coverage, likely inherited\">" << std::endl;
	VCFOutFile << "##FILTER=<ID=SB,Description=\"Contig fails string bias filter\">" << std::endl; 
	VCFOutFile << "##ALT=<ID=INS:ME:ALU,Description=\"Insertion of ALU element\">" << std::endl;
	VCFOutFile << "##ALT=<ID=INS:ME:L1,Description=\"Insertion of L1 element\">" << std::endl;
	VCFOutFile << "##ALT=<ID=INS:ME:MOB,Description=\"Insertion of ALU or L1element\">" << std::endl; 

	VCFOutFile << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t";
	
	string samplename = outStub.substr(0, outStub.find(".generator")); 
	VCFOutFile << samplename; 
	//VCFOutFile << outStub; 
	for(int i =0; i<ParentHashFilePaths.size(); i++)
	{
		string ParPath = argv[ParentHashFilePaths[i]]; 
		int startpos  = ParPath.find("overlap.asembly.hash.fastq.");
		int endpos = ParPath.find(".generator.Jhash");
		string Par = ParPath.substr(startpos+27, endpos-(startpos+27)); 
		ParNames.push_back(Par); 
		VCFOutFile << "\t" << Par; 
		//VCFOutFile << "\tParent" << i; //argv[ParentHashFilePaths[i]]; 
	}
	VCFOutFile<< endl;

	int lines = 0;

	line="";
 
	unsigned long LongHash;

	cout << "Reading in Sam File"  << endl;
	map <string, int> Names; 
	vector <SamRead> reads; 
	int counter = 0; 
	while (getline(SamFile, line))
	{
		if (line.c_str()[0] == '@')
		{
			cout << " HEADER LINE = " << line << endl; 
			  vector <string> temp = Split(line, '\t');
			cout << temp[0] << endl;
			if (temp[0] == "@SQ")
			{
				cout << temp[1] << endl;
				vector <string> chr = Split(temp[1], ':');
				vector <string> len = Split(temp[2], ':');
	
				cout << "##contig=<ID=" <<  chr[1]<<",length=" << len[1] << ">"<< endl;
				VCFOutFile <<"##contig=<ID=" <<  chr[1]<<",length=" << len[1] << ">" << endl;
			}
		}
		else
		{
			cout << line << endl;
			counter ++; 
			SamRead read; 
			cout << "parse " << endl; 
			read.parse(line);
			//if (read.mapQual > 0)
			if (read.FlagBits[2] != 1) //verify if read is mapped 
			{
				read.parsed = true; 
				read.getRefSeq();
				read.createPeakMap();
				read.checkMob(mobs);
				if (AlignmentAllA(read) > .9)
				{
					read.mapQual = 0; 
					read.AllA = true; 
				}
				int a; 
				string b;
				cout << "Aligned bases = " << read.CheckBasesAligned() << endl; 
				if (read.CheckBasesAligned() > 50 or read.CheckEndsAlign())
				{reads.push_back(read);}
				else 
				{cout << "SKIPPING Alignment" << endl; read.write();}
				if (counter%100 == 0)
					cout << "read " << counter << " entries " << char(13); 
			}
			//else do I want to track unaliged alignments? 
		}
	}
	cout << endl;
	cout << "Read in " << reads.size() << " reads " << endl;
	if (reads.size() == 0)
	{
		cout << "no reads were passed to RUFUS.interpret exiting" << endl;
		return 0; 
	}

	//  VCFOutFile << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t";
	//VCFOutFile << outStub << endl;

	
	cout << "procesing split reads" << endl;
	for (int i = 0; i < reads.size(); i++) 
	{
		cout << "processing read " << reads[i].name << endl;
		if (reads[i].alignments.size() == 0)
		{
			reads[i].alignments.push_back(i);
			int count = 0; 
			for (int j = i+1; j < reads.size(); j++)
			{
				count++; 
				if (strcmp(reads[i].name.c_str(), reads[j].name.c_str()) == 0 && reads[i].pos )// !=reads[j].pos)
				{
					cout << "found mate " << reads[j].name << endl;
					
					reads[i].alignments.push_back(j); 
					reads[j].alignments.push_back(j); 
					reads[j].alignments.push_back(i); 
					reads[j].first = false;
					
				}
				//if (count > 100000)
				//	break; 
			}
		}
		//for (int j = 0; j<reads[i].alignments.size(); j++)
		//	{reads[reads[i].alignments[j]].alignments = reads[i].alignments;} 

	}
	for (int i = 0; i < reads.size(); i++)
	{
		reads[i].LookUpKmers();
		reads[i].CheckPhase(); 
		reads[i].clipPattern = reads[i].ClipPattern();
	}




	cout << "finding multi contig events, " << reads.size() << endl; 
	//find multi contig events insertionsf
	for (int i = 0; i < reads.size()-1; i++)
	{
	 	int pos, pos2, kdep, kdep2;
		string InsStart, InsEnd;
		cout << "########################################STARTING NEW LOOP MOB#############################" << endl; 	
		cout << "i = " << i << endl; 
		reads[i].write();
		//check mobil elements first
		if (reads[i].isSplitRead > 0 && reads[i].MobAligned)
		{
			bool found = false; 
			cout << "Posible multi read mobil element" << endl;
			reads[i].write();
			
			// do I have a breakpoint supporoted by hashes
			int bp = reads[i].sigBreakPoint();
			cout << "HEREbp = " << bp << endl; 
			if (bp > 0)
			{
				cout << "Passed SigBreakpoint check MobLoop" << endl; 
				int start = -2;
                                while (start + i < 0)
                                {start ++;}
				cout << "start = " << start << endl; 
				for ( int j = start; j<= 2 && j+i >= 0 && j+i < reads.size() ; j++)
				{
					cout << "	checking " << j << " = " << reads[i+j].name << " " << reads[i].chr << " == " << reads[i+j].chr  <<" && abs( " <<  reads[i].pos << " - " << reads[i+j].pos << ") < 2000 = " << reads[i].pos-reads[i+j].pos << endl; 
			
					if (j != 0 /*&& reads[i].name != reads[i+j].name*/ &&   reads[i].chr == reads[i+j].chr && abs(reads[i].pos - reads[i+j].pos) < 2000 )	
					{
						vector <SamRead> temp; 
						for (int k = 0; k < reads[i+j].alignments.size(); k++)
						{
							if (reads[i+j].alignments[k] != i+j && reads[reads[i+j].alignments[k]].mapQual > 30)
							{
								temp.push_back(reads[reads[i+j].alignments[k]]);
							}
						}
						int polyAbp =  reads[i+j].isPolyA(temp);
						
						cout << "Poly A pos = " << polyAbp << "read A pos = " << reads[i].pos+bp << " and read B pos = " << reads[i+j].pos + polyAbp  << " with abs " << abs((reads[i].pos+bp) - (reads[i+j].pos+polyAbp))<< endl;
						
						if ( polyAbp > -1 && abs((reads[i].pos+bp) - (reads[i+j].pos+polyAbp)) < 50 )
						{
							if ((reads[i].clipPattern == "cm" &&  reads[i+j].clipPattern == "mc" ) || (reads[i].clipPattern == "mc" &&  reads[i+j].clipPattern == "cm" )  ) 
							{
								cout << "FOUND MOB" << endl;
								if (reads[i].SVeventid ==0)
								{
									CurrentSVeventID++;
								
									int targetsize = 0; 
									if (reads[i].clipPattern == "mc")
									{
										int start = reads[i].pos + bp; 
										int end = reads[i+j].pos + reads[i+j].sigBreakPoint(); 
										targetsize = start - end; 
									}
									else
									{
										int end = reads[i].pos  + bp;
										int start = reads[i+j].pos + reads[i+j].sigBreakPoint();
										targetsize = start - end; 
									}
									//if (targetsize >0)
									{
										
										string GenotypeField;
										if (CheckGenotypes(reads[i].createStructGenotype(reads[i].sigBreakPoint())))
											GenotypeField = reads[i].createStructGenotype(reads[i].sigBreakPoint());
										else 
											GenotypeField = reads[i+j].createStructGenotype(reads[i+j].sigBreakPoint());
										

										stringstream Format ;
										Format << InterpretTargetSize(targetsize); 
										Format << "MOB"; 
										Format << "-"; 
										Format << "LCa-" <<  reads[i].SVCheckParentsForLowCov(reads[i].sigBreakPoint()) << "-LCb-" << reads[i+j].SVCheckParentsForLowCov(reads[i+j].sigBreakPoint()) << "-" ;  
										Format << reads[i].MobAS; 
										string ref = Reff.getSubSequence(reads[i].chr, reads[i].pos+bp -1 ,1); 
										
										stringstream alt ; 
										alt << "<INS:ME:MOB>"; 
		
										string FullfilterA = reads[i].filterSV(); 
										string FullfilterB = "" ; //reads[i+j].filterSV(); 
										int GMap = 0;
										int minMapQual = 30;


										///////// not sure this is a good idea, why am I checking all the alignemtns, im just worried about this one
						//				for(int k = 0; k < reads[i].alignments.size(); k++)
						//				{
						//					if ( reads[reads[i].alignments[k]].mapQual > minMapQual)
						//						GMap++;
						//				}
						//				for(int k = 0; k < reads[i+j].alignments.size(); k++)
						//				{
						//					if ( reads[reads[i+j].alignments[k]].mapQual > minMapQual)
						//						GMap++;
						//				}


										if (reads[i].mapQual > minMapQual)
											GMap++;
										////if (reads[reads[i].alignments[1]].mapQual > minMapQual)
										////	GMap++;
										if ( reads[i+j].mapQual > minMapQual)
											GMap++;
										////if (reads[reads[i+j].alignments[1]].mapQual > minMapQual)
										////	GMap++;
										
										string InfoFilter = "";
										string Filter = "";
										
										if (GMap <= 0)
										{
											Format<< "-LowMapQual";
											InfoFilter = "LowMapQual";
											Filter = "LMQ"; 
										}
										else if (FullfilterA == "" && FullfilterB == "")
										{
											found = true; 
											Format<<"-DeNovo";
											InfoFilter = "Pass";
											Filter = "PASS"; 
                                                                                	for(int k = 0; k < reads[i].alignments.size(); k++)
                                                                                	{reads[reads[i].alignments[k]].SVeventid = CurrentSVeventID;}
                                                                                	for(int k = 0; k < reads[i+j].alignments.size(); k++)
                                                                                	{reads[reads[i+j].alignments[k]].SVeventid = CurrentSVeventID;}	
										}
										else
										{
											Format<<"-" << FullfilterA << "," <<FullfilterB; 
											InfoFilter = FullfilterA; 
											InfoFilter +=FullfilterB;
											Filter = "fail"; 
										}


										//make quality stuff
										int readAmut=0;
										int readApos=0;
										int readBmut=0;
										int readBpos=0;
										reads[i].GetQualityHashes(readAmut, readApos, reads[i].sigBreakPoint());
										reads[i+j].GetQualityHashes(readBmut, readBpos, reads[i+j].sigBreakPoint()); 
										
										float qual = -100; 
										if ((readApos+readBpos) > 0)
											qual = ((float)readAmut+(float)readBmut) / ((float)readApos+(float)readBpos) * 100.0; 
										else
											qual = 0; 
		
										//buildng up info field
										stringstream info; 
										info << "SVTYPE=INS;END=" <<  reads[i].pos+bp -1 << ";";
										info << "MT=" << reads[i].MobContig << ";"; 
										string phase="none"; 
										if (reads[i].phase != "none")
											phase = reads[i].phase; 
										else if (reads[i+j].phase != "none")
											phase = reads[i+j].phase; 
										info << "PH=" << phase << ";"; 
										info << "FEX=" << InfoFilter << ";";
										int SupportingHashes = readAmut+readBmut;
										int possibleHashes = readApos+readBpos; 
										info << "FS=" << SupportingHashes << "/" <<possibleHashes << ";"; 
										 
										info << "RN=" << reads[i].name << "_and_" << reads[i+j].name << ";";
										info << "MQ=" << reads[i].mapQual << "_and_" << reads[i+j].mapQual << ";";
										info << "cigar=" << reads[i].cigar << "_and_" << reads[i+j].cigar << ";";
										info << "SB=" << reads[i].StrandBias << "_and_" << reads[i+j].StrandBias << ";"; 
										info << "AS=" << reads[i].AlignmentSegments << "-" << reads[i].AlignmentSegmentsCigar << "_and_" << reads[i+j].AlignmentSegments << "-" << reads[i+j].AlignmentSegmentsCigar;
										stringstream call;	
										call << reads[i].chr << "\t" << reads[i].pos+bp -1 << "\t" << Format.str() << "\t" << ref <<  "\t" << alt.str() <<  "\t" << qual  << "\t" << Filter << "\t" << info.str() <<  "\t" << "GT:DP:RO:AO\t" << GenotypeField <<  endl;
										cout << call.str();
										VCFOutFile << call.str(); 
									//	break;
									}
									//else
									//	cout << "MOBS cant have deletions, " << targetsize << endl; 
								}
							}
						}
					}
				}
			}
			if (found)
				continue; 
		}
	}
	for (int i = 0; i < reads.size()-1; i++)
	{
		int pos, pos2, kdep, kdep2;
		string InsStart, InsEnd;
		cout << "########################################STARTING NEW LOOP#############################" << endl;
		reads[i].write();
		cout << "checkingbigdel" << endl; 
		reads[i].write(); 
		if (reads[i].alignments.size() == 2 && reads[i].SVeventid == 0)
		{
			cout << "possible big del" << endl; 
			reads[i].write(); 
			reads[reads[i].alignments[1]].write();
			if (reads[i].chr == reads[reads[i].alignments[1]].chr)
			{
				cout << "reads on same chr" << endl; 
				if (reads[i].sigBreakPoint() > 0 || reads[reads[i].alignments[1]].sigBreakPoint() > 0 || BreakpointInUnalignedCenter(reads[i], reads[reads[i].alignments[1]] ))
				{
					cout << "one of them has a sig breakpoint" << endl; 
					if (((reads[reads[i].alignments[1]].pos + reads[reads[i].alignments[1]].BreakPoint()) - (reads[i].pos+reads[i].BreakPoint())) > MaxVarentSize)
					{
						cout << "its big enough " << abs((reads[reads[i].alignments[1]].pos + reads[reads[i].alignments[1]].BreakPoint()) - (reads[i].pos+reads[i].BreakPoint())) << endl; 
						if (GetReadOrientation(reads[i].flag) == GetReadOrientation(reads[reads[i].alignments[1]].flag))
						{
							cout << "same orientation" << endl; 
							if (reads[i].clipPattern == "mc" && reads[reads[i].alignments[1]].clipPattern == "cm" )
							{
								cout << "correct pattern its a del " << endl; 
								//deletion
						
								CurrentSVeventID++;
								//reads[i].SVeventid = CurrentSVeventID;
								//reads[reads[i].alignments[1]].SVeventid = CurrentSVeventID;
									
								string GenotypeField;
								if (CheckGenotypes(reads[i].createStructGenotype(reads[i].sigBreakPoint())))
									GenotypeField = reads[i].createStructGenotype(reads[i].sigBreakPoint());
								else 
									GenotypeField = reads[reads[i].alignments[1]].createStructGenotype(reads[reads[i].alignments[1]].sigBreakPoint());
					
								if (GenotypeField == "")
									cout << "ERRORGeno DUP" << " " << GenotypeField << endl;
								string insertseq = GetUnalignedCenter(reads[i], reads[reads[i].alignments[1]]);	
								int targetsize = ((reads[reads[i].alignments[1]].pos + reads[reads[i].alignments[1]].BreakPoint()) - (reads[i].pos + reads[i].BreakPoint())) *-1;
								stringstream Format ;
								Format << InterpretInsertSize(insertseq);
								Format << InterpretTargetSize(targetsize); ////////////////////////////////////
								
								
								string ref = Reff.getSubSequence(reads[i].chr, reads[i].pos + reads[i].BreakPoint() -1 , 1); 
								
								stringstream alt ; 
								alt << insertseq; 
								alt << "<DEL>"; 

								string FullfilterA = reads[i].filterSV(); 
								int GMap = 0;
								int minMapQual = 40;
								if (reads[i].mapQual > minMapQual)
									GMap++;
								if (reads[reads[i].alignments[1]].mapQual > minMapQual)
									GMap++;
								
								string InfoFilter = "";
								string Filter = "";
								
								if (GMap < 1)
								{
									Format<< "-LowMapQual";
									InfoFilter = "LowMapQual";
									Filter = "LMQ"; 
								}
								else if (FullfilterA == "")
								{
									Format<<"-DeNovo";
									InfoFilter = "Pass";
									Filter = "PASS"; 
									reads[i].SVeventid = CurrentSVeventID;
                                                                	reads[reads[i].alignments[1]].SVeventid = CurrentSVeventID;	
								}
								else
								{
									Format<<"-" << FullfilterA; 
									InfoFilter = FullfilterA; 
									Filter = "fail"; 
								}
								//make quality stuff
								int readAmut=0;
								int readApos=0;
								reads[i].GetQualityHashes(readAmut, readApos, reads[i].BreakPoint());
								
								float qual = -100; 
								if ((readApos) > 0)
									qual = ((float)readAmut) / ((float)readApos) * 100.0; 
								else
									qual = 0; 

								
								//buildng up info field
								stringstream info; 
								info << "SVTYPE=DEL;END=" <<  reads[reads[i].alignments[1]].pos + reads[reads[i].alignments[1]].BreakPoint() << ";";
								info << "SVLEN=" << targetsize*1 << ";"; 
								string phase="none"; 
								if (reads[i].phase != "none")
									phase = reads[i].phase; 
								info << "PH=" << phase << ";"; 
								info << "FEX=" << InfoFilter << ";";
								int SupportingHashes = readAmut;
								int possibleHashes = readApos; 
								info << "FS=" << SupportingHashes << "/" <<possibleHashes << ";"; 
								 
								info << "RN=" << reads[i].name << ";";
								info << "MQ=" << reads[i].mapQual << "_and_" << reads[reads[i].alignments[1]].mapQual << ";";
								info << "cigar=" << reads[i].cigar << "_and_" << reads[reads[i].alignments[1]].cigar << ";";
								info << "SB=" << reads[i].StrandBias << ";"; 
								info << "AS=" << reads[i].AlignmentSegments << "-" << reads[i].AlignmentSegmentsCigar << "_and_" << reads[reads[i].alignments[1]].AlignmentSegments << "-" << reads[reads[i].alignments[1]].AlignmentSegmentsCigar;

								stringstream call;	
								call << reads[i].chr << "\t" << reads[i].pos+ reads[i].BreakPoint()  -1 << "\t" << Format.str() << "\t" << ref <<  "\t" << alt.str() <<  "\t" << qual  << "\t" << Filter << "\t" << info.str() <<  "\t" << "GT:DP:RO:AO\t" << GenotypeField <<  endl;
								cout << call.str();
								VCFOutFile << call.str(); 
								
							
							}
							else if (reads[i].clipPattern == "cm" && reads[reads[i].alignments[1]].clipPattern == "mc")
							{
							//	dup
								cout << "event is dup" << endl; 
								CurrentSVeventID++;
								//reads[i].SVeventid = CurrentSVeventID;
								//reads[reads[i].alignments[1]].SVeventid = CurrentSVeventID;
								cout << "here -1" << endl; 

								string GenotypeField;
								
								int sbpA = reads[i].sigBreakPoint(); 
								int sbpB = reads[reads[i].alignments[1]].sigBreakPoint(); 
								cout << "sbpA = " << sbpA << " sbpB = " << sbpB << endl; 

								if (CheckGenotypes(reads[i].createStructGenotype(reads[i].sigBreakPoint())))
									GenotypeField = reads[i].createStructGenotype(reads[i].sigBreakPoint());
								else 
									GenotypeField = reads[reads[i].alignments[1]].createStructGenotype(reads[reads[i].alignments[1]].sigBreakPoint());
								
								if (GenotypeField == "")
									cout << "ERRORGeno DUP" << " " << GenotypeField << endl; 
								cout << "here" << endl; 
								string insertseq = GetUnalignedCenter(reads[i], reads[reads[i].alignments[1]]);	
								int targetsize = ((reads[reads[i].alignments[1]].pos + reads[reads[i].alignments[1]].BreakPoint()) - (reads[i].pos + reads[i].BreakPoint())) ;
								stringstream Format ;
								Format << InterpretInsertSize(insertseq);
								Format << InterpretTargetSize(targetsize); ////////////////////////////////////
								cout << "here2" << endl; 
								
								string ref = Reff.getSubSequence(reads[i].chr, reads[i].pos + reads[i].BreakPoint() -1 -1 , 1); 
								
								stringstream alt ; 
								alt << insertseq; 
								alt << "<DUP>"; 
								cout << "here3" << endl; 
								string FullfilterA = reads[i].filterSV(); 
								int GMap = 0;
								int minMapQual = 20;
								if (reads[i].mapQual > minMapQual)
									GMap++;
								if (reads[reads[i].alignments[1]].mapQual > minMapQual)
									GMap++;
								cout << "here4" << endl; 
								string InfoFilter = "";
								string Filter = "";
								
								if (GMap < 2)
								{
									Format<< "-LowMapQual";
									InfoFilter = "LowMapQual";
									Filter = "LMQ"; 
								}
								else if (FullfilterA == "")
								{
									Format<<"-DeNovo";
									InfoFilter = "Pass";
									Filter = "PASS"; 
									reads[i].SVeventid = CurrentSVeventID;
                                                                	reads[reads[i].alignments[1]].SVeventid = CurrentSVeventID;	
								}
								else
								{
									Format<<"-"<< FullfilterA; 
									InfoFilter = FullfilterA; 
									Filter = "fail"; 
								}
								cout << "here 5" << endl; 
								//make quality stuff
								int readAmut=0;
								int readApos=0;
								reads[i].GetQualityHashes(readAmut, readApos, reads[i].BreakPoint());
								cout << "here 6" << endl; 
								float qual = -100; 
								if ((readApos) > 0)
									qual = ((float)readAmut) / ((float)readApos) * 100.0; 
								else
									qual = 0; 

								cout << "here 7" << endl; 
								//buildng up info field
								stringstream info; 
								info << "SVTYPE=DUP;END=" <<  reads[reads[i].alignments[1]].pos + reads[reads[i].alignments[1]].BreakPoint() << ";";
								cout << "here 8" << endl; 
								info << "SVLEN=" << targetsize << ";"; 
								string phase="none"; 
								if (reads[i].phase != "none")
									phase = reads[i].phase; 
								info << "PH=" << phase << ";"; 
								info << "FEX=" << InfoFilter << ";";
								int SupportingHashes = readAmut;
								int possibleHashes = readApos; 
								info << "FS=" << SupportingHashes << "/" <<possibleHashes << ";"; 
								 
								info << "RN=" << reads[i].name << ";";
								info << "MQ=" << reads[i].mapQual << "_and_" << reads[reads[i].alignments[1]].mapQual << ";";
								info << "cigar=" << reads[i].cigar << "_and_" << reads[reads[i].alignments[1]].cigar << ";";
								info << "SB=" << reads[i].StrandBias << ";"; 
								info << "AS=" << reads[i].AlignmentSegments << "-" << reads[i].AlignmentSegmentsCigar << "_and_" << reads[reads[i].alignments[1]].AlignmentSegments << "-" << reads[reads[i].alignments[1]].AlignmentSegmentsCigar;

								stringstream call;	
								call << reads[i].chr << "\t" << reads[i].pos+ reads[i].BreakPoint()  -1 << "\t" << Format.str() << "\t" << ref <<  "\t" << alt.str() <<  "\t" << qual  << "\t" << Filter << "\t" << info.str() <<  "\t" << "GT:DP:RO:AO\t" << GenotypeField <<  endl;
								cout << call.str();
								VCFOutFile << call.str(); 
							}
						}
					}
				}
			}
		}
		
		
	}
	cout << "###################done with multi contigs####################" << endl; 
	for (int i = 0; i < reads.size(); i++)
	{
		SamRead read = reads[i];
		cout << "starting work on  " << endl;
		cout << read.name  << endl;
		cout << "alignments = " << read.alignments.size()  << endl;

		if (strcmp(read.chr.c_str(), "*") == 0)
		{
			read.writetofile(Unaligned);
			Unaligned << endl;
			cout << "read unaligned, skipping" << endl; 
		}
		//if  (read.first)
		{
			//if it looks like a simple split read alignment colaps it into a single read

			if (read.first  and read.alignments.size() > 1)
			{
				cout << "picking two best alignments" << endl;
				map <float,int> alignScores;
				for (int j = 0; j < read.alignments.size(); j++){
					float score = (float) reads[read.alignments[j]].AlignScore;
					while (not (alignScores.find(score) == alignScores.end())){
					  score = score * 1.0001;
					}
					alignScores[score] = j;
				}
				vector <int> goodPos;
				std::map<float, int>::reverse_iterator it;
				for ( it = alignScores.rbegin(); it != alignScores.rend(); it++ )
				{
					cout << it->first << " - " << it->second << endl;
					goodPos.push_back(it->second);
				}
				cout << "atempting colaps "<< read.name << endl;
				vector<SamRead> R;

				for(int j =0; j < read.alignments.size(); j++) //read.alignments.size(); j++)
				{
					//these better be sorted by position
					//if (reads[read.alignments[j]].chr == read.chr)
					if (j == goodPos[0] or j == goodPos[1])
					{
						R.push_back(reads[read.alignments[j]]);
						cout << reads[read.alignments[j]].name << endl;
					}
				}
				cout << "made it out o this loop" << endl;
				if (R.size() ==2 && /*R[0].chr == R[1].chr && */ (R[0].mapQual > 0 and R[1].mapQual > 0) && R[0].SVeventid == 0 )
				{
					cout << "straing better way Rsize = " << R.size() << endl;
					read = BetterWay(R);
					cout << "ending better way" << endl;
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
			if (read.mapQual > MinMapQual and read.alignments.size() <=2)
			{
				//reads[i].CheckPhase();
				read.parseMutations(argv, reads);
			}
			else
			{
				cout << "map qual " << read.mapQual << " less than mim map qual of " << MinMapQual << endl; 
			}
		}
	}
	cout << "finding multi contig events, " << reads.size() << endl;
	//Mopping up less likey multi contig optoins
	for (int i = 0; i < reads.size()-1; i++)
		
	{
		// chekcing translocation
		if (reads[i].alignments.size()==2)
		{
			if ((reads[i].sigBreakPoint() > 0 || reads[reads[i].alignments[1]].sigBreakPoint() > 0 || BreakpointInUnalignedCenter(reads[i], reads[reads[i].alignments[1]] )) && reads[i].clipPattern.length()==2 && reads[reads[i].alignments[1]].clipPattern.length()==2)
			{
				cout << "passed First Trans " << endl; 
				reads[i].write();
				//if (reads[i].chr != reads[reads[i].alignments[1]].chr) //check if read I looks like a trans 
				{
					cout << "Posible translocation" << endl;
					reads[reads[i].alignments[1]].write();
					//cout << "done writing possible trans" << endl; 
					int start = -2;
                               	 	while (start + i < 0)
                                        {start ++;}
					for (int j = start; j<3 && i+j < reads.size(); j++)
					{
						cout << "cheking " << j << " with i = " << i << " and max = " << reads.size() << " " << reads[i+j].name << " " << reads[i+j].alignments.size()<< endl; 
						if (reads[i+j].name != reads[i].name && reads[i+j].alignments.size()==2  )
						{

							cout << "	 " << reads[i+j].name << " != " << reads[i].name << "  " << reads[i].clipPattern << " "<< reads[reads[i].alignments[1]].clipPattern<< " && " << reads[i+j].clipPattern << " && " << reads[reads[i+j].alignments[1]].clipPattern << endl; 
							if (reads[i+j].clipPattern.length()==2 && reads[reads[i+j].alignments[1]].clipPattern.length()==2 && /*what the funk and I checking here I thing im saying for the remote alignment, this reads alignemnt needs to be the next one*/ (reads[i+j].alignments[1]-1 == reads[i].alignments[1] || reads[i+j].alignments[1]+1 == reads[i].alignments[1] || reads[i+j].alignments[1]+2 == reads[i].alignments[1] || reads[i+j].alignments[1]-2 == reads[i].alignments[1])) 
							{
								reads[i+j].write();
								reads[reads[i+j].alignments[1]].write(); 
								cout<< "passed alignments check " << reads[i+j].sigBreakPoint() << " " << reads[reads[i+j].alignments[1]].sigBreakPoint() <<endl; 
								if ((reads[i+j].sigBreakPoint() > 0 ||  reads[reads[i+j].alignments[1]].sigBreakPoint() > 0 || BreakpointInUnalignedCenter(reads[i+j], reads[reads[i+j].alignments[1]] ))) // check if read I+J looks like a trans 
								{
									cout << "passig sig break point check " << endl; 
									if ((reads[i+j].chr == reads[i].chr && reads[reads[i+j].alignments[1]].chr ==  reads[reads[i].alignments[1]].chr) ) // chekc if readsa i and reads i+j show the same trans
									{
										cout << "passed chr check" << endl; 

										int breaks = 0; 
										if (reads[i].sigBreakPoint() > 0)
											breaks++; 
										if (reads[reads[i].alignments[1]].sigBreakPoint() > 0)
											breaks++; 
										if ( reads[i+j].sigBreakPoint() > 0)
											breaks++; 
										if (reads[reads[i+j].alignments[1]].sigBreakPoint() > 0)
											breaks++; 
										
										int GMap = 0; 
										int minMapQual = 30; 
	
										if (reads[i].mapQual > minMapQual)
											GMap++;
										if (reads[reads[i].alignments[1]].mapQual > minMapQual)
											GMap++;
										if ( reads[i+j].mapQual > minMapQual)
											GMap++;
										if (reads[reads[i+j].alignments[1]].mapQual > minMapQual)
											GMap++;
	
										if (reads[i].sigBreakPoint() > 0 || reads[i+j].sigBreakPoint() > 0 && breaks >= 3  ) //atleast ons of the reads in this location has to have a sig break point 
										{
											cout << "passed sig breakcheck" << endl; 
											cout << reads[i].chr << " == " << reads[reads[i].alignments[1]].chr << endl;
											if (reads[i].chr != reads[reads[i].alignments[1]].chr)
											{
												//trans chromosomal event
												if (reads[i].SVeventid == reads[i+j].SVeventid)
												{
													if (reads[i].SVeventid ==0)
													{
														CurrentSVeventID++;
														reads[i].BNDid = MaxBND+1; MaxBND++;
														reads[i+j].BNDid = MaxBND+1; MaxBND++;
														reads[reads[i].alignments[1]].BNDid = MaxBND+1; MaxBND++;
														reads[reads[i+j].alignments[1]].BNDid = MaxBND+1; MaxBND++;
													}	
													int targetsize = 0;
													int bp = reads[i].BreakPoint();
													int bpj = reads[i+j].BreakPoint();
													int sbp = reads[reads[i].alignments[1]].BreakPoint(); 
													int sbpj = reads[reads[i+j].alignments[1]].BreakPoint();
													
													if (reads[i].clipPattern == "mc")
													{
														int start = reads[i].pos  + bp;
														int end = reads[i+j].pos + bpj;
														targetsize = start - end;
													}
													else
													{
														int end = reads[i].pos  + bp;
														int start = reads[i+j].pos + bpj;
														targetsize = start - end;
													}
													
													cout << "FOUND TRANSLOCATION" << endl;
													cout << "targetsize = " << targetsize << endl; 	
													stringstream Format; 
													int InsCorrect = targetsize; 
													int DelCorrect = targetsize; 
													if (InsCorrect <0){InsCorrect = 0;}
													if (DelCorrect >0){DelCorrect = 0;}
													cout << "Reff.getSubSequence(" << reads[i].chr << " , " <<reads[i].pos+bp -1 -1 - InsCorrect  << "   , 5 " << ");" << " = " << Reff.getSubSequence( reads[i].chr, reads[i].pos+bp -1 -1 -InsCorrect ,5) << endl; 	
													string ref ;
													string altseq ;
													string insertseq; 
													stringstream alt;
													string SVDES = ""; 
													int offset = bp;
													if ( reads[i].clipPattern == "mc")
													{
														
														offset = bp -1 - InsCorrect; 
														ref = Reff.getSubSequence(reads[i].chr, reads[i].pos+bp -1 -1 -  InsCorrect   , 1 +abs(DelCorrect));
		       										 		altseq = Reff.getSubSequence(reads[i].chr, reads[i].pos+bp -1 -1  - InsCorrect  , 1 +InsCorrect );
														if      ( reads[i].clipPattern == "mc" && GetReadOrientation(reads[i].flag) == GetReadOrientation(reads[reads[i].alignments[1]].flag))
														{
															insertseq = GetUnalignedCenter(reads[i], reads[reads[i].alignments[1]]); 
															alt << altseq <<  insertseq; 
															alt << "[" << reads[reads[i].alignments[1]].chr << ":" << reads[reads[i].alignments[1]].pos+sbp -1 << "["; 
															Format <<InterpretInsertSize(insertseq); 
															Format <<InterpretTargetSize(targetsize) << "_";  
															Format <<  "bnd_" << reads[i].BNDid;
															SVDES = "Translocation"; 
														}
														else if ( reads[i].clipPattern == "mc" && GetReadOrientation(reads[i].flag) !=GetReadOrientation( reads[reads[i].alignments[1]].flag))
														{
															SamRead temp = reads[reads[i].alignments[1]]; 
															temp.flipRead(); 
															insertseq = GetUnalignedCenter(reads[i], temp);
															alt << altseq <<  insertseq;
															alt << "]" << reads[reads[i].alignments[1]].chr << ":" << reads[reads[i].alignments[1]].pos+sbp -1 << "]";
															Format <<InterpretInsertSize(insertseq);
															Format <<InterpretTargetSize(targetsize) << "_";
															Format << "bnd_" << reads[i].BNDid;
															SVDES = "InvertedTranslocation"; 
														}
													}
													else if ( reads[i].clipPattern == "cm")
													{
														offset = bp -1 ; 
														ref = Reff.getSubSequence(reads[i].chr, reads[i].pos+bp -1  ,1);
														altseq = Reff.getSubSequence(reads[i].chr, reads[i].pos+bp -1 ,1);
														if ( reads[i].clipPattern == "cm" && GetReadOrientation(reads[i].flag) == GetReadOrientation(reads[reads[i].alignments[1]].flag))
														{
															alt << "]" << reads[reads[i].alignments[1]].chr << ":" << reads[reads[i].alignments[1]].pos + sbp /*check could be +1*/ << "]" << altseq; 
		       											 		Format << "bnd_" << reads[i].BNDid ;
															SVDES = "Translocation";
														}
														else if ( reads[i].clipPattern == "cm" && GetReadOrientation(reads[i].flag) != GetReadOrientation(reads[reads[i].alignments[1]].flag))
														{										
															alt << "[" << reads[reads[i].alignments[1]].chr << ":" << reads[reads[i].alignments[1]].pos+ sbp -1 /*check could be +1*/  << "[" << altseq;
															Format << "bnd_" << reads[i].BNDid ;
															SVDES = "InvertedTranslocation";
														}	
													}								
													
			
													///Filter is caculate on both reads, if one fails they both fail since its detection relies on both pairs agreeing with eachother. 
													string FullfilterA = reads[i].filterSV();
													string FullfilterB = reads[i+j].filterSV(); 
													string InfoFilter = "";
													string Filter = "";
													if (GMap < 1 || !( (reads[i].mapQual > 0 || reads[i+j].mapQual > 0) && (reads[reads[i].alignments[1]].mapQual > 0 || reads[reads[i+j].alignments[1]].mapQual > 0)))
													{
														Format << "-LowMapQual";
														InfoFilter  = "LowMapQual";
														Filter = "LMQ";
													}
													else if (FullfilterA == "" && FullfilterB == "")
													{
														InfoFilter = "Pass";
														Filter = "PASS";
														Format<<"-DeNovo";
														reads[i].SVeventid = CurrentSVeventID;
                                                                                                                reads[i+j].SVeventid = CurrentSVeventID;
                                                                                                                reads[reads[i].alignments[1]].SVeventid = CurrentSVeventID;
                                                                                                                reads[reads[i+j].alignments[1]].SVeventid = CurrentSVeventID;
													}
													else
													{
														InfoFilter = FullfilterA;
														InfoFilter +=FullfilterB;
														Filter = "fail";
													}
			
													int readAmut=0;
													int readApos=0;
													reads[i].GetQualityHashes(readAmut, readApos, bp);
			
													float qual = -100;
													if ((readApos) > 0)
														qual = ((float)readAmut) / ((float)readApos) * 100.0;
													else
														qual = 0;
			
			
													//buildng up info field
													stringstream info;
										       			info << "SVTYPE=BND;MATEID=bnd_" << reads[reads[i].alignments[1]].BNDid << ";";
													info << "SVID=" << CurrentSVeventID<< ";"; 
													if (SVDES !="")
													{info << "SVDES=" << SVDES << ";";}
													string phase="none";
													if (reads[i].phase != "none")
														phase = reads[i].phase;
													else if (reads[i+j].phase != "none")
														phase = reads[i+j].phase;
													info << "PH=" << phase << ";";
													info << "FEX=" << InfoFilter << ";";
		       									 		int SupportingHashes = readAmut;
		       									 		int possibleHashes = readApos;
		       										 	info << "FS=" << SupportingHashes << "/" <<possibleHashes << ";";
			
													info << "RN=" << reads[i].name << ";";
													info << "MQ=" << reads[i].mapQual << ";";
													info << "cigar=" << reads[i].cigar << ";";
													info << "SB=" << reads[i].StrandBias << ";";
													info << "AS=" << reads[i].AlignmentSegments << "-" << reads[i].AlignmentSegmentsCigar ;
			
			
			
													
													string GenotypeFieldA = reads[i].createStructGenotype(bp);
													stringstream call; 
													call << reads[i].chr << "\t" << reads[i].pos+offset  << "\t" << Format.str()<< "\t" << ref <<  "\t" << alt.str() <<  "\t" << qual  << "\t" << Filter << "\t" << info.str() <<  "\t" << "GT:DP:RO:AO\t" << GenotypeFieldA <<  endl; 
													cout << call.str(); 
													VCFOutFile << call.str(); 
			
			
										//////////////////////////OTHER STRAND//////////////////////////////////////////////////////////////////////////////////////////////
													Format.str(std::string());
													info.str(std::string()); 
													stringstream altj;
													SVDES = ""; 
													ref = ""; 
													altseq = "";
													offset = 0; 
													if ( reads[i+j].clipPattern == "mc")
													{
														offset = bpj -1 -  InsCorrect ; 
														ref = Reff.getSubSequence(reads[i+j].chr, reads[i+j].pos+bpj -1 -1 -  InsCorrect , 1 +abs(DelCorrect));
														altseq = Reff.getSubSequence(reads[i+j].chr, reads[i+j].pos+bpj -1 -1  - InsCorrect  , 1 +InsCorrect );
														if ( reads[i+j].clipPattern == "mc" && GetReadOrientation(reads[i+j].flag) == GetReadOrientation(reads[reads[i+j].alignments[1]].flag))
														{
															altj << altseq << "[" << reads[reads[i+j].alignments[1]].chr << ":" << reads[reads[i+j].alignments[1]].pos+ sbpj  << "[";
															Format << "bnd_" << reads[i+j].BNDid;
															SVDES = "Translocation";
														}
														else if ( reads[i+j].clipPattern == "mc" && GetReadOrientation(reads[i+j].flag) !=GetReadOrientation( reads[reads[i+j].alignments[1]].flag))
														{
															altj << altseq << "]" << reads[reads[i+j].alignments[1]].chr << ":" << reads[reads[i+j].alignments[1]].pos+ sbpj << "]";
															 Format << "bnd_" << reads[i+j].BNDid;
															SVDES = "InvertedTranslocation";
														}
													}
													else if (reads[i+j].clipPattern == "cm")
													{
														offset = bpj ;
														ref = Reff.getSubSequence(reads[i+j].chr, reads[i+j].pos+bpj -1  ,1);
														altseq = Reff.getSubSequence(reads[i+j].chr, reads[i+j].pos+bpj-1 ,1);
														if ( reads[i+j].clipPattern == "cm" && GetReadOrientation(reads[i+j].flag) == GetReadOrientation(reads[reads[i+j].alignments[1]].flag))
														{
															altj << "]" << reads[reads[i+j].alignments[1]].chr << ":" << reads[reads[i+j].alignments[1]].pos+ sbpj << "]" << altseq;
															 Format << "bnd_" << reads[i+j].BNDid;
															SVDES = "Translocation";
														}
														else if	 ( reads[i+j].clipPattern == "cm" && GetReadOrientation(reads[i+j].flag) != GetReadOrientation(reads[reads[i+j].alignments[1]].flag))
														{	
													 	       altj << "[" << reads[reads[i+j].alignments[1]].chr << ":" << reads[reads[i+j].alignments[1]].pos+sbpj << "[" << altseq;
															 Format << "bnd_" << reads[i+j].BNDid;
															SVDES = "InvertedTranslocation";
														}
													}
													
													readAmut=0;
													readApos=0;
													reads[i+j].GetQualityHashes(readAmut, readApos, bpj);
			
		       											 qual = -100;
			       										 if ((readApos) > 0)
		       					  						       qual = ((float)readAmut) / ((float)readApos) * 100.0;
		       											 else
		       												 qual = 0;
		
		
													//buildng up info field
													if (GMap < 1|| !( (reads[i].mapQual > 0 || reads[i+j].mapQual > 0) && (reads[reads[i].alignments[1]].mapQual > 0 || reads[reads[i+j].alignments[1]].mapQual > 0)))
													{
														Format << "-LowMapQual";
														InfoFilter  = "LowMapQual";
														Filter = "LMQ";
													}
													else if (FullfilterA == "" && FullfilterB == "")
													{
														InfoFilter = "Pass";
														Filter = "PASS";
													}
													else
													{
														InfoFilter = FullfilterA;
														InfoFilter +=FullfilterB;
														Filter = "fail";
													}	
													info << "SVTYPE=BND;MATEID=bnd_" << reads[reads[i+j].alignments[1]].BNDid << ";";
													info << "SVID=" << CurrentSVeventID << ";";
													if (SVDES !="")
													{info << "SVDES=" << SVDES << ";";}
		       											phase="none";
													if (reads[i+j].phase != "none")
			      											 phase = reads[i+j].phase;
													else if (reads[i+j].phase != "none")
														phase = reads[i+j].phase;
													info << "PH=" << phase << ";";
													info << "FEX=" << InfoFilter << ";";
													SupportingHashes = readAmut;
													possibleHashes = readApos;
													info << "FS=" << SupportingHashes << "/" <<possibleHashes << ";";
			
													info << "RN=" << reads[i+j].name << ";";
													info << "MQ=" << reads[i+j].mapQual << ";";
													info << "cigar=" << reads[i+j].cigar << ";";
													info << "SB=" << reads[i+j].StrandBias << ";";
													info << "AS=" << reads[i+j].AlignmentSegments << "-" << reads[i+j].AlignmentSegmentsCigar ;
			
			
			
													string GenotypeFieldB = reads[i+j].createStructGenotype(bpj);
												
													call.str(std::string());
													call << reads[i+j].chr << "\t" << reads[i+j].pos+offset  << "\t" << Format.str() << "\t" << ref <<  "\t" << altj.str() <<  "\t" << qual  << "\t" << Filter << "\t" << info.str() <<  "\t" << "GT:DP:RO:AO\t" << GenotypeFieldB <<  endl;
													cout << call.str(); 
													VCFOutFile << call.str(); 
			
													//cout << reads[i].chr << " " << reads[i].pos + readAsig << " to " << reads[reads[i].alignments[1]].chr << " " << reads[reads[i].alignments[1]].sigBreakPoint()+ reads[reads[i].alignments[1]].pos << "\tGenoA " << GenotypeFieldA << "\tGenoB " << GenotypeFieldB << endl; 
													if (j >= 0)
														i = i+j; 
													continue; 
												}
											}
											else //were dealing with an intrachromosomal event, probably a duplication chr6:162,664,265-162,666,099 and 
											{
												cout << "possible interchomosmal event" << endl; 
												int EnterA;
												int ExitA;
												int EnterB;
												int ExitB; 
												int EventPos = -1;
												int TargetSize =0; 
												string insert; 
												string RefSeq; 
												string AltSeq;
												string insertchr ; 
												int insertstart ;
												int insertend ;
												int insertSize;
												bool nope = false; 
												if( reads[i].clipPattern == "mc" && reads[reads[i].alignments[1]].clipPattern == "cm" && reads[i+j].clipPattern == "cm" && reads[reads[i+j].alignments[1]].clipPattern == "mc")
												{
													cout << "yaya typeA" << endl; 
													EnterA = i; 
													ExitB = reads[i].alignments[1]; 
													ExitA = i+j;
													EnterB = reads[i+j].alignments[1];
												}
												else if( reads[i].clipPattern == "cm" && reads[reads[i].alignments[1]].clipPattern == "mc" &&  reads[i+j].clipPattern == "mc" && reads[reads[i+j].alignments[1]].clipPattern == "cm")
												{
													cout << "yay typeB" << endl; 
													ExitA = i; 
													EnterB = reads[i].alignments[1]; 
													EnterA = i+j;
													ExitB = reads[i+j].alignments[1];
												}
												else
												{
													cout << "booo dosnt fit any tyep" << endl; 
													nope = true; 
												}
												if (nope == false )
												{
													if (reads[EnterA].pos + reads[EnterA].BreakPoint() <= reads[ExitA].pos+reads[ExitA].BreakPoint())
													{
														EventPos = reads[EnterA].pos+ reads[EnterA].BreakPoint()-1;
														TargetSize =  (reads[ExitA].pos + reads[ExitA].BreakPoint()) - (reads[EnterA].pos+ reads[EnterA].BreakPoint()); 
													}
													else if (reads[EnterB].pos+ reads[EnterB].BreakPoint() <= reads[ExitB].pos+reads[ExitB].BreakPoint())
													{
														EventPos = reads[EnterB].pos+ reads[EnterB].BreakPoint()-1;
														TargetSize =  (reads[ExitB].pos+reads[ExitB].BreakPoint()) - (reads[EnterB].pos+ reads[EnterB].BreakPoint());
													}
													else
													{TargetSize = -1;}
													if (TargetSize >= 0 && TargetSize < 1000000)
													{
														RefSeq =  Reff.getSubSequence(reads[EnterA].chr, EventPos-1, 1+TargetSize);
														AltSeq = Reff.getSubSequence(reads[EnterA].chr, EventPos-1, 1);
													
													
														if (reads[EnterB].pos+reads[EnterB].BreakPoint() > reads[ExitB].pos + reads[ExitB].BreakPoint())
														{
															insertchr = reads[EnterB].chr; 
															insertstart = reads[ExitB].pos + reads[ExitB].BreakPoint(); 
															insertend  = reads[EnterB].pos + reads[EnterB].BreakPoint();
															insertSize = insertend - insertstart; 
															insert = Reff.getSubSequence(reads[EnterA].chr, insertstart -1 , insertSize );
														}
														else if (reads[EnterA].pos+reads[EnterA].BreakPoint() > reads[ExitA].pos + reads[ExitA].BreakPoint())
														{
															insertchr = reads[EnterA].chr; 
															insertstart = reads[ExitA].pos + reads[ExitA].BreakPoint();
															insertend  = reads[EnterA].pos + reads[EnterA].BreakPoint();
															insertSize = insertend - insertstart;
															insert = Reff.getSubSequence(reads[EnterA].chr, insertstart -1 , insertSize );
														}
														else 
														{
															insertSize = -1; 
														}
														AltSeq += insert; 	
														if (insertSize > 0)
														{
														
	
		
															string FullfilterA = reads[i].filterSV();
															string FullfilterB = reads[i+j].filterSV();
															string InfoFilter = "";
															string Filter = "";	
															stringstream Format ;
															Format << InterpretTargetSize(TargetSize *-1 );
															Format << insert.size();
															Format << "-"<<insertSize; 
															Format << "CopyPaste"; 
															int minMapQual = 30;
															
															
															
															if (GMap < 1|| !( (reads[i].mapQual > 0 || reads[i+j].mapQual > 0) && (reads[reads[i].alignments[1]].mapQual > 0 || reads[reads[i+j].alignments[1]].mapQual > 0)))
															{
																Format << "-LowMapQual";
																InfoFilter  = "LowMapQual";
																Filter = "LMQ";
															}
															else if (FullfilterA == "" && FullfilterB == "")
															{
																InfoFilter = "Pass";
																Filter = "PASS";
															}
															else
															{
																InfoFilter = FullfilterA;
																InfoFilter +=FullfilterB;
																Filter = "fail";
															}
	///////////////////////////////////do I need to add filter stuff to this one? ///////////////////////////////////	
															if (reads[i].SVeventid ==0)
															{
																CurrentSVeventID++;
																for(int k = 0; k < reads[i].alignments.size(); k++)
																{reads[reads[i].alignments[k]].SVeventid = CurrentSVeventID;}
																for(int k = 0; k < reads[i+j].alignments.size(); k++)
																{reads[reads[i+j].alignments[k]].SVeventid = CurrentSVeventID;}	
																
																int readAmut=0;
																int readApos=0;
																int readBmut=0;
																int readBpos=0;
																reads[i].GetQualityHashes(readAmut, readApos, reads[i].sigBreakPoint());
																reads[i+j].GetQualityHashes(readBmut, readBpos, reads[i+j].sigBreakPoint());
								
																float qual = -100;
																if ((readApos+readBpos) > 0)
																	qual = ((float)readAmut+(float)readBmut) / ((float)readApos+(float)readBpos) * 100.0;
																else
																	qual = 0;		
																stringstream info; 
																info << "SVTYPE=COPY:PASTE;" << ";";
																info << "SOURCE=" << insertchr <<":" << insertstart <<"-"<<insertend << ";"; 
																info << "SVID=" << reads[i+j].SVeventid << ";";
																string phase="none";
																if (reads[i+j].phase != "none")
																	 phase = reads[i+j].phase;
																else if (reads[i+j].phase != "none")
																	phase = reads[i+j].phase;
																info << "PH=" << phase << ";";
																info << "FEX=" << InfoFilter << ";";
																int SupportingHashes = readAmut+readBmut;
																int possibleHashes = readApos+readBpos;
																info << "FS=" << SupportingHashes << "/" <<possibleHashes << ";";
			
																info << "RN=" << reads[i].name << "_and_" << reads[i+j].name << ";";
																info << "MQ=" << reads[i].mapQual << "_and_" << reads[i+j].mapQual << ";";
																info << "cigar=" << reads[i].cigar << "_and_" << reads[i+j].cigar << ";";
												       				info << "SB=" << reads[i].StrandBias << "_and_" << reads[i+j].StrandBias << ";";
																info << "AS=" << reads[i].AlignmentSegments << "-" << reads[i].AlignmentSegmentsCigar << "_and_" << reads[i+j].AlignmentSegments << "-" << reads[i+j].AlignmentSegmentsCigar;
			
																string GenotypeField; 
																if (CheckGenotypes(reads[i].createStructGenotype(reads[i].sigBreakPoint())))
																	GenotypeField = reads[i].createStructGenotype(reads[i].sigBreakPoint());
																else if (CheckGenotypes(reads[i+j].createStructGenotype(reads[i+j].sigBreakPoint())))
																	GenotypeField = reads[i+j].createStructGenotype(reads[i+j].sigBreakPoint());
																else if (CheckGenotypes(reads[reads[i].alignments[1]].createStructGenotype(reads[reads[i].alignments[1]].sigBreakPoint())))
																	GenotypeField = reads[reads[i].alignments[1]].createStructGenotype(reads[reads[i].alignments[1]].sigBreakPoint());
																else
																	GenotypeField =  reads[reads[i+j].alignments[1]].createStructGenotype(reads[reads[i+j].alignments[1]].sigBreakPoint());	
																 stringstream call;		
																call << reads[i].chr << "\t" << EventPos  << "\t" << Format.str() << "\t" << RefSeq <<  "\t" << AltSeq <<  "\t" << qual  << "\t" << Filter << "\t" << info.str() <<  "\t" << "GT:DP:RO:AO\t" << GenotypeField <<  endl;
																cout << call.str();
																VCFOutFile << call.str();
																if (j >= 0)
                                                                                                                			i = i+j;
                                                                                                        			continue;	
															}
														}
		
													}
												}


											}
										}
									}
								
							
									else if ((reads[i+j].chr == reads[reads[i+j].alignments[1]].chr  && reads[reads[i+j].alignments[1]].chr ==  reads[i].chr))
									{
										cout << "odd backwords translocation detection, I should do something about this" << endl; 
									}
									cout << "donewith " << j << " with i = " << i << " and max = " << reads.size() << endl;
								}
							}
						}
						cout << "done checking this loop" << endl; 
					}
				}
			}
		}
		cout << "frinsihed trans" << endl; 
		//checking for invertions
		if ( reads[i].alignments.size() == 2 && reads[i].SVeventid ==0)
		{
			if(reads[i].chr == reads[reads[i].alignments[1]].chr  && GetReadOrientation(reads[i].flag) != GetReadOrientation(reads[reads[i].alignments[1]].flag)  && reads[i].sigBreakPoint() > 0)
			{
				cout << "newpossible inversion" << endl;
				reads[i].write(); 
				reads[reads[i].alignments[1]].write(); 
				int start = -2;
                                while (start + i < 0)
                                {start ++;}	
				for ( int j = start ;j<=1 && j+i >= 0 && j+i < reads.size() ; j++)
				{
					if(reads[i+j].alignments.size() > 1 & j != 0 )
					{
						if(reads[i+j].chr == reads[reads[i+j].alignments[1]].chr  && GetReadOrientation(reads[i+j].flag) != GetReadOrientation(reads[reads[i+j].alignments[1]].flag)  && reads[i+j].sigBreakPoint() > 0 &&  (reads[i+j].alignments[1]-1 == reads[i].alignments[1] || reads[i+j].alignments[1]+1 == reads[i].alignments[1]))
						{
							if(reads[i].chr == reads[i+j].chr)
							{
								int positionAa = reads[i].pos + reads[i].sigBreakPoint();
								int positionBa = reads[i+j].pos + reads[i+j].sigBreakPoint();
								int positionAb = reads[reads[i].alignments[1]].pos + reads[reads[i].alignments[1]].sigBreakPoint();
								int positionBb = reads[reads[i+j].alignments[1]].pos + reads[reads[i+j].alignments[1]].sigBreakPoint();
								if (positionAa < positionAb && positionBa < positionBb &&  reads[i].clipPattern != reads[i+j].clipPattern )
								{
									//if (abs(positionAa - positionBa) < HashSize || abs(positionAb -  positionBb) < HashSize || abs(positionAa -  positionBb) < HashSize || abs(positionAb - positionBa) < HashSize )
									

									if (abs(positionAa - positionBa) < HashSize || abs(positionAb -  positionBb) < HashSize )
									{
										CurrentSVeventID++;
										cout << "Found Invertion" << endl; 
										reads[i].write(); 
										reads[i+j].write(); 
										int pos = -1; 
										if(positionAa < positionBa)
											pos = positionAa;
										else
											pos = positionBa;

										int end = 0; 
										if(positionAb > positionBb)
											end = positionAb;
										else
											end = positionBb; 
										
										int startBreak  = 0; 
										if (reads[i].clipPattern == "mc" && reads[i+j].clipPattern == "cm")
											startBreak = positionAa - positionBa;
										else if (reads[i].clipPattern == "cm" && reads[i+j].clipPattern == "mc")
											startBreak = positionBa - positionAa; 
										int endBreak  = 0;
										if (reads[reads[i].alignments[1]].clipPattern == "mc" && reads[reads[i+j].alignments[1]].clipPattern == "cm")
											endBreak = positionAb - positionBb;
										else if (reads[reads[i].alignments[1]].clipPattern == "cm" &&  reads[reads[i+j].alignments[1]].clipPattern == "mc")
											startBreak = positionBb - positionAb;	
										

										int size = end-pos; 

										stringstream call; 




											
										stringstream alt; 
										stringstream ref; 
										stringstream info; 
										string GenotypeField;
										SamRead temp = reads[reads[i].alignments[1]];
										temp.flipRead(); 
										string STARTinsertedseq = GetUnalignedCenter(reads[i], temp); 
										temp = reads[reads[i+j].alignments[1]];
										temp.flipRead();
										string ENDinsertseq = GetUnalignedCenter(reads[i+j], temp); 

										ref << Reff.getSubSequence(reads[i].chr, pos -1 -1, 1); 
										alt << STARTinsertedseq; 
										alt << "<INV>"; 
										alt << ENDinsertseq; 
																			
										int readAmut=0;
										int readApos=0;
										int readBmut=0;
										int readBpos=0;
										reads[i].GetQualityHashes(readAmut, readApos, reads[i].sigBreakPoint());
										reads[i+j].GetQualityHashes(readBmut, readBpos, reads[i+j].sigBreakPoint());

										float qual = -100;
										if ((readApos+readBpos) > 0)
											qual = ((float)readAmut+(float)readBmut) / ((float)readApos+(float)readBpos) * 100.0;
										else
											qual = 0;
										int SupportingHashes = readAmut+readBmut;
										int possibleHashes = readApos+readBpos;

										string phase="none";
										if (reads[i].phase != "none")
											phase = reads[i].phase;
										else if (reads[i+j].phase != "none")
											phase = reads[i+j].phase;	
										

										string FullfilterA = reads[i].filterSV();
										string FullfilterB = reads[i+j].filterSV();

										string InfoFilter = "";
										string Filter = "";
		
										int GMap = 0;
										int minMapQual = 30;
	
										if (reads[i].mapQual > minMapQual)
											GMap++;
										if (reads[reads[i].alignments[1]].mapQual > minMapQual)
											GMap++;
										if ( reads[i+j].mapQual > minMapQual)
											GMap++;
										if (reads[reads[i+j].alignments[1]].mapQual > minMapQual)
											GMap++;
										
										stringstream Format;
										if (startBreak>0)
										{
											Format << abs(startBreak) << "Y";
										}
										else if (startBreak<0)
										{
											Format << abs(startBreak) << "D";
										}
										Format << InterpretInsertSize(STARTinsertedseq); 
										Format <<  size << "V";
										if (endBreak>0)
										{
											Format << abs(endBreak) << "Y";
										}
										else if (endBreak<0)
										{
											Format << abs(endBreak) << "D";
										}	
										Format <<  InterpretInsertSize(ENDinsertseq); 
										
										
										if (GMap < 1)
										{
											Format<< "-LowMapQual";
											InfoFilter = "LowMapQual";
											Filter = "LMQ"; 
										}
										else if (FullfilterA == "" && FullfilterB == "")
										{
											Format << "-DeNovo";
											InfoFilter = "Pass";
											Filter = "PASS";
											for(int k = 0; k < reads[i].alignments.size(); k++)
                                                                                	{reads[reads[i].alignments[k]].SVeventid = CurrentSVeventID;}
                                                                                	for(int k = 0; k < reads[i+j].alignments.size(); k++)
                                                                                	{reads[reads[i+j].alignments[k]].SVeventid = CurrentSVeventID;}	
										}
										else
										{
											Format << "-" << FullfilterA << "," << FullfilterB;
											InfoFilter = FullfilterA;
											InfoFilter +=FullfilterB;
											Filter = "fail";
										}
						
										info << "SVTYPE=INV;END=" << end << ";"; 
										info << "PH=" << phase << ";";
										info << "FEX=" << InfoFilter << ";";
										info << "FS=" << SupportingHashes << "/" <<possibleHashes << ";";

										info << "RN=" << reads[i].name << "_and_" << reads[i+j].name << ";";
										info << "MQ=" << reads[i].mapQual << "_and_" << reads[i+j].mapQual << ";";
										info << "cigar=" << reads[i].cigar << "_and_" << reads[i+j].cigar << ";";
										info << "SB=" << reads[i].StrandBias << "_and_" << reads[i+j].StrandBias << ";";
										info << "AS=" << reads[i].AlignmentSegments << "-" << reads[i].AlignmentSegmentsCigar << "_and_" << reads[i+j].AlignmentSegments << "-" << reads[i+j].AlignmentSegmentsCigar;	

										cout << "read[i] genotype " << reads[i].sigBreakPoint() << " = " << reads[i].createStructGenotype(reads[i].sigBreakPoint()) << endl;
										cout << "read[i+j] genotype " << reads[i+j].sigBreakPoint() << "= " << reads[i+j].createStructGenotype(reads[i+j].sigBreakPoint())<< endl;
										cout << "read[reads[i].alignments[1]] genotype " <<reads[reads[i].alignments[1]].sigBreakPoint() << " = " << reads[reads[i].alignments[1]].createStructGenotype(reads[reads[i].alignments[1]].sigBreakPoint())<< endl;
										cout << "read[reads[i+j].alignments[1]] genotype "<<reads[reads[i+j].alignments[1]].sigBreakPoint() << " = " << reads[reads[i+j].alignments[1]].createStructGenotype(reads[reads[i+j].alignments[1]].sigBreakPoint())<< endl;
									
										if (CheckGenotypes(reads[i].createStructGenotype(reads[i].sigBreakPoint())))
											GenotypeField = reads[i].createStructGenotype(reads[i].sigBreakPoint()); 
										else if (CheckGenotypes(reads[i+j].createStructGenotype(reads[i+j].sigBreakPoint())))
											GenotypeField = reads[i+j].createStructGenotype(reads[i+j].sigBreakPoint());
										else if (CheckGenotypes(reads[reads[i].alignments[1]].createStructGenotype(reads[reads[i].alignments[1]].sigBreakPoint())))
											GenotypeField = reads[reads[i].alignments[1]].createStructGenotype(reads[reads[i].alignments[1]].sigBreakPoint()); 
										else 
											GenotypeField =  reads[reads[i+j].alignments[1]].createStructGenotype(reads[reads[i+j].alignments[1]].sigBreakPoint()); 


										call << reads[i].chr << "\t" << pos -1  << "\t" << Format.str() << "\t" << ref.str() <<  "\t" << alt.str() <<  "\t" << qual  << "\t" << Filter << "\t" << info.str() <<  "\t" << "GT:DP:RO:AO\t" << GenotypeField <<  endl;
										VCFOutFile << call.str(); 
										cout << call.str(); 
									}
									if (j >= 0)
                                                                        	i = i+j;
                                                                        continue; 
								}
								else
								{
									cout << "hopefully already got this one" << endl; 
								}
							}
						}
					}	
				}
			}
		}
		//checking for large insertions
		
		//chekcing for shorter tripple aligned insertions or duplications
		if ( reads[i].alignments.size() ==3 && reads[i].sigBreakPoint() > 0)
		{
			cout << "cheking tripple" << endl; 
			reads[i].write();
			reads[reads[i].alignments[1]].write(); 
			reads[reads[i].alignments[2]].write();
			int start = -1; 
			int mid = -1; 
			int exit = -1;
			//find entering contig
			if (reads[i].clipPattern == "mc")
				start = i; 
			else if (reads[reads[i].alignments[1]].clipPattern == "mc")
				start = reads[i].alignments[1]; 
			else if (reads[reads[i].alignments[2]].clipPattern == "mc")
				start = reads[i].alignments[2];

			//find source contig	
			if (reads[i].clipPattern == "cmc")
				mid = i;
			else if (reads[reads[i].alignments[1]].clipPattern == "cmc")     
			       	mid = reads[i].alignments[1];
			else if (reads[reads[i].alignments[2]].clipPattern == "cmc")
				mid = reads[i].alignments[2];

			 //find exit contig
			if (reads[i].clipPattern == "cm")
				exit = i;
			else if (reads[reads[i].alignments[1]].clipPattern == "cm")     
				exit = reads[i].alignments[1];
			else if (reads[reads[i].alignments[2]].clipPattern == "cm")
				exit = reads[i].alignments[2];
	
			cout << "start = " << start << " mid = " << mid << " exit = " << exit<< endl; 

			if (start > 1 && mid > 1 && exit > 1)
			{
				cout << "yay" << endl; 
				cout << reads[start].chr << " and " << reads[exit].chr << endl; 
				if (reads[start].chr == reads[exit].chr && (reads[exit].sigBreakPoint() > 0 || reads[start].sigBreakPoint() > 0)  && ( reads[exit].mapQual > 0 && reads[start].mapQual > 0 ))
				{
					cout << "yay2" << endl; 
					int TargetSize = ((reads[exit].pos+reads[exit].sigBreakPoint()) - (reads[start].pos+reads[start].sigBreakPoint())) * -1;
					if (reads[start].SVeventid ==0)
					{
						cout << "found tripple" << endl; 
						int pos = reads[start].pos +  reads[start].BreakPoint() -1 ; 
						reads[start].write();
						reads[mid].write();
						reads[exit].write();
						CurrentSVeventID++;
						cout << "here " << endl; 	
						string GenotypeField;
						cout << "start bkreapoing = " << reads[start].BreakPoint() <<" and exit = " << reads[exit].BreakPoint() << endl; 
						cout << "reads start genoptype = " ;
						cout << reads[start].createStructGenotype(reads[start].BreakPoint()) << endl; 
						cout << "reads end genotype = " ; 
						cout << reads[exit].createStructGenotype(reads[mid].BreakPoint()); 
						cout << "boom" << endl; 
						if (CheckGenotypes(reads[start].createStructGenotype(reads[start].BreakPoint())))
							GenotypeField = reads[start].createStructGenotype(reads[start].BreakPoint());
						else if (CheckGenotypes(reads[exit].createStructGenotype(reads[exit].BreakPoint()))) 
							GenotypeField = reads[exit].createStructGenotype(reads[exit].BreakPoint());
						else 
							GenotypeField = reads[exit].createStructGenotype(reads[mid].BreakPoint());
						
						cout << "here 2" << endl; 
						string Format = InterpretTargetSize(TargetSize); 
						Format += "trippleDUP"; 
						cout << "here 3" << endl; 
						string ref = Reff.getSubSequence(reads[start].chr, reads[start].pos+reads[start].BreakPoint()-1-1  ,1); 
						if (TargetSize < 0)
							ref = ref + Reff.getSubSequence(reads[start].chr, reads[start].pos+reads[start].BreakPoint()-1  , TargetSize *-1);

						stringstream alt ; 
						
						alt << Reff.getSubSequence(reads[start].chr, reads[start].pos+reads[start].BreakPoint() -1 -1  ,1);
						if (TargetSize > 0)
							 alt << Reff.getSubSequence(reads[start].chr, reads[start].pos+reads[start].BreakPoint()-1  , TargetSize );
						alt << reads[mid].seq.substr(reads[start].BreakPoint(), reads[exit].BreakPoint() - reads[start].BreakPoint()  ) ;
						//alt << reads[mid].seq.substr(reads[start].BreakPoint(), reads[mid].CountBasesAligned(reads[mid].BreakPoint())+1 +reads[mid].BreakPoint()- reads[start].BreakPoint()  ) ;
						//alt << Reff.getSubSequence(reads[mid].chr, reads[mid].pos+reads[mid].BreakPoint(), reads[mid].CountBasesAligned(reads[mid].BreakPoint())+1) ;
						cout << "here 3" << endl; 
						string FullfilterA = reads[start].filterSV(); 
						string FullfilterB = reads[mid].filterSV(); 
						string FullfilterC = reads[exit].filterSV();
						cout << "here 4" << endl; 
						int GMap = 0;
						int minMapQual = 30;
						if (reads[start].mapQual > minMapQual)
							GMap++;
						if (reads[mid].mapQual > minMapQual)
							GMap++;
						if (reads[exit].mapQual > minMapQual)
							GMap++;
						cout << "here 5"<< endl; 
						string InfoFilter = "";
						string Filter = "";
						
						if (GMap < 1)
						{
							Format+= "-LowMapQual";
							InfoFilter = "LowMapQual";
							Filter = "LMQ"; 
						}
						else if (FullfilterA == "" && FullfilterB == "" && FullfilterC == "")
						{
							Format+="-DeNovo";
							InfoFilter = "Pass";
							Filter = "PASS"; 
							reads[start].SVeventid = CurrentSVeventID;
                                                	reads[mid].SVeventid = CurrentSVeventID;
                                                	reads[exit].SVeventid = CurrentSVeventID;
						}
						else
						{
							Format+="-";
							Format+= FullfilterA;
							Format+= "," ; 
							Format+= FullfilterB; 
							InfoFilter = FullfilterA; 
							InfoFilter +=FullfilterB;
							Filter = "fail"; 
						}
						cout << "here 6" << endl; 
						//make quality stuff
						int readAmut=0;
						int readApos=0;
						int readBmut=0;
						int readBpos=0;
						reads[start].GetQualityHashes(readAmut, readApos, reads[start].BreakPoint());
						reads[exit].GetQualityHashes(readBmut, readBpos, reads[exit].BreakPoint()); 
						
						float qual = -100; 
						if ((readApos+readBpos) > 0)
							qual = ((float)readAmut+(float)readBmut) / ((float)readApos+(float)readBpos) * 100.0; 
						else
							qual = 0; 
		
						cout << "here 7" << endl; 
						//buildng up info field
						stringstream info; 
						info << "SVTYPE=INS;END=" <<  reads[start].pos+reads[start].BreakPoint() -1 << ";";
						info << "SOURCE=" << reads[mid].chr <<":" << reads[mid].pos + reads[mid].BreakPoint() <<"-"<< reads[mid].pos + reads[mid].BreakPoint() + reads[mid].CountBasesAligned(reads[mid].BreakPoint())<< ";";
						string phase="none"; 
						if (reads[start].phase != "none")
							phase = reads[start].phase; 
						else if (reads[exit].phase != "none")
							phase = reads[exit].phase; 
						else if (reads[mid].phase != "none")
							phase = reads[mid].phase;
						info << "PH=" << phase << ";"; 
						info << "FEX=" << InfoFilter << ";";
						int SupportingHashes = readAmut+readBmut;
						int possibleHashes = readApos+readBpos; 
						info << "FS=" << SupportingHashes << "/" <<possibleHashes << ";"; 
						 
						info << "RN=" << reads[start].name << ";";
						info << "MQ=" << reads[start].mapQual << "_and_" << reads[mid].mapQual << "_and_" << reads[exit].mapQual << ";";
						info << "cigar=" << reads[start].cigar << "_and_" << reads[mid].cigar << "_and_" << reads[exit].cigar << ";";
						info << "SB=" << reads[start].StrandBias << ";"; 
						info << "AS=" << reads[start].AlignmentSegments << "-" << reads[start].AlignmentSegmentsCigar << "_and_" << reads[mid].AlignmentSegments << "-" << reads[mid].AlignmentSegmentsCigar << "_and_" << reads[exit].AlignmentSegments << "-" << reads[exit].AlignmentSegmentsCigar;
		
							
						//cout << reads[start].chr << "\t" << reads[start].pos+reads[start].sigBreakPoint() -1 << "\t" << Format << "\t" << ref <<  "\t" << alt.str() <<  "\t" << qual  << "\t" << Filter << "\t" << info.str() <<  "\t" << "GT:DP:RO:AO\t" << GenotypeField <<  endl;
						//VCFOutFile << reads[start].chr << "\t" << pos << "\t" << Format << "\t" << ref <<  "\t" << alt.str() <<  "\t" << qual  << "\t" << Filter << "\t" << info.str() <<  "\t" << "GT:DP:RO:AO\t" << GenotypeField <<  endl;
						stringstream call; 
						call << reads[start].chr << "\t" << pos << "\t" << Format << "\t" << ref <<  "\t" << alt.str() <<  "\t" << qual  << "\t" << Filter << "\t" << info.str() <<  "\t" << "GT:DP:RO:AO\t" << GenotypeField <<  endl;
						VCFOutFile << call.str(); 
						cout << "here 8 " << endl; 
						//break;
					
					}	
				}
			}
		}
		//could add the very large insert tandem stuff here	
		if ( reads[i].alignments.size() >=1 && reads[i].clipPattern == "mc" && reads[i].sigBreakPoint() > 0 && reads[i].SVeventid == 0 )
		{
			cout << "possible large insertion" << endl;
			reads[i].write();
		 	int start = -5;
                        while (start + i < 0)
                        {start ++;}
			for ( int j = start; j <= 5 && j+i >= 0 && j+i < reads.size() ; j++)
			{
				cout << "I = " << i <<  "and  J = " << j << endl; 
				cout << reads[i+j].name << " sig break point = " << reads[i+j].sigBreakPoint() << endl; 
				if (reads[i+j].alignments.size() == 1 && reads[i+j].clipPattern == "cm"  && reads[i+j].sigBreakPoint() > 0 && reads[i].chr == reads[i+j].chr  && reads[i+j].SVeventid ==  0)
				{
					cout << "even more possible large insert" << endl;
					int sbI = reads[i].sigBreakPoint();
					int sbJ = reads[i+j].sigBreakPoint();
					int positionI = reads[i].pos + sbI;
					int positionJ = reads[i+j].pos + sbJ;
					
						//if (abs(positionAa - positionBa) < HashSize || abs(positionAb -  positionBb) < HashSize || abs(positionAa -  positionBb) < HashSize || abs(positionAb - positionBa) < HashSize )
						
						
						if (abs(positionI - positionJ) < 1000000 && reads[i].SVeventid ==0 && reads[i+j].SVeventid == 0 && (reads[i].mapQual > 0 && reads[i+j].mapQual > 0))
						{
							cout << "Found Large Insert" << endl; 
							reads[i].write(); 
							reads[i+j].write(); 
							int pos = -1; 
							if(positionI < positionJ)
								pos = positionI;
							else
								pos = positionJ;

							cout << "pos = " << pos << endl; 
							int end = 0; 
							if(positionI > positionJ)
								end = positionI;
							else
								end = positionJ; 
							cout << "end = " << end << endl; 
							

							int startBreak = positionI - positionJ ; 

							stringstream call; 

							cout << "here 1" << endl; 
							stringstream Format;
							if (startBreak>0)
							{
								Format << abs(startBreak) << "Y";
							}
							else if (startBreak<0)
							{
								Format << abs(startBreak) << "D"; 
							}
							cout << "here 2 " << endl; 
							stringstream alt; 
							stringstream ref; 
							stringstream info; 
							string GenotypeField; 
							//ref << Reff.getSubSequence(reads[i].chr, pos -1 -1, 1); 
							alt << "<INS>"; 
							if (startBreak>0){
								ref << Reff.getSubSequence(reads[i].chr, pos -1 -1  , 1);
								alt << Reff.getSubSequence(reads[i].chr, pos -1 -1 , 1+abs(startBreak)); 
							}
							else if (startBreak<0){
								ref << Reff.getSubSequence(reads[i].chr, pos -1 -1  , 1+abs(startBreak));
								 alt << Reff.getSubSequence(reads[i].chr, pos -1 -1 , 1); 
							}
							cout << "here 3" << endl; 
							//ref << "-" << Reff.getSubSequence(reads[i].chr, pos -1 -1 + abs(startBreak)+1, size - abs(startBreak) - abs(endBreak) ); 
							cout <<" posI = " << sbI << " and pos j = " << sbJ << endl; 
							string leftseq = reads[i].getClippedSequence(sbI, "mc");  
							cout << "here 3.11" << endl; 
							string rightseq =  reads[i+j].getClippedSequence(sbJ, "cm"); 
							cout << "here 3.12" << endl; 
							alt << "-" << reads[i].getClippedSequence(sbI, "mc") << "NNNNNNNNNNNNNNNNNNNN" << reads[i+j].getClippedSequence(sbJ, "cm"); 
							cout << "here 3.1" << endl;  
							Format <<  alt.str().length()<<"+" << "LargeInsert"; 	
							cout << "here 3.2" << endl; 
							int readAmut=0;
							int readApos=0;
							int readBmut=0;
							int readBpos=0;
							reads[i].GetQualityHashes(readAmut, readApos, reads[i].sigBreakPoint());
							cout << "here 3.3" << endl; 
							reads[i+j].GetQualityHashes(readBmut, readBpos, reads[i+j].sigBreakPoint());
							cout << "here 4" << endl; 
							float qual = -100;
							if ((readApos+readBpos) > 0)
								qual = ((float)readAmut+(float)readBmut) / ((float)readApos+(float)readBpos) * 100.0;
							else
								qual = 0;
							int SupportingHashes = readAmut+readBmut;
							int possibleHashes = readApos+readBpos;

							string phase="none";
							if (reads[i].phase != "none")
								phase = reads[i].phase;
							else if (reads[i+j].phase != "none")
								phase = reads[i+j].phase;	
							cout << "here 5" << endl; 
							CurrentSVeventID++;
							for(int k = 0; k < reads[i].alignments.size(); k++)
							{reads[reads[i].alignments[k]].SVeventid = CurrentSVeventID;}
							for(int k = 0; k < reads[i+j].alignments.size(); k++)
							{reads[reads[i+j].alignments[k]].SVeventid = CurrentSVeventID;}	
							string FullfilterA = reads[i].filterSV();
							string FullfilterB = reads[i+j].filterSV();

							string InfoFilter = "";
							string Filter = "";
		
							int GMap = 0;
							int minMapQual = 30;
							cout << "here 6" << endl; 
							if (reads[i].mapQual > minMapQual)
								GMap++;
							if ( reads[i+j].mapQual > minMapQual)
								GMap++;
							if (GMap < 1)
							{
								Format<< "-LowMapQual";
								InfoFilter = "LowMapQual";
								Filter = "LMQ"; 
							}
							else if (FullfilterA == "" && FullfilterB == "")
							{
								Format << "-DeNovo";
								InfoFilter = "Pass";
								Filter = "PASS";
							}
							else
							{
								Format << "-" << FullfilterA << "," << FullfilterB;
								InfoFilter = FullfilterA;
								InfoFilter +=FullfilterB;
								Filter = "fail";
							}
			
							info << "SVTYPE=INS;END=" << end << ";"; 
							info << "PH=" << phase << ";";
							
							info << "FEX=" << InfoFilter << ";";
							info << "FS=" << SupportingHashes << "/" <<possibleHashes << ";";

							info << "RN=" << reads[i].name << "_and_" << reads[i+j].name << ";";
							info << "MQ=" << reads[i].mapQual << "_and_" << reads[i+j].mapQual << ";";
							info << "cigar=" << reads[i].cigar << "_and_" << reads[i+j].cigar << ";";
							info << "SB=" << reads[i].StrandBias << "_and_" << reads[i+j].StrandBias << ";";
							info << "AS=" << reads[i].AlignmentSegments << "-" << reads[i].AlignmentSegmentsCigar << "_and_" << reads[i+j].AlignmentSegments << "-" << reads[i+j].AlignmentSegmentsCigar;	

							cout << "read[i] genotype " << reads[i].sigBreakPoint() << " = " << reads[i].createStructGenotype(reads[i].sigBreakPoint()) << endl;
							cout << "read[i+j] genotype " << reads[i+j].sigBreakPoint() << "= " << reads[i+j].createStructGenotype(reads[i+j].sigBreakPoint())<< endl;
						
							if (CheckGenotypes(reads[i].createStructGenotype(reads[i].sigBreakPoint())))
								GenotypeField = reads[i].createStructGenotype(reads[i].sigBreakPoint()); 
							else
								GenotypeField = reads[i+j].createStructGenotype(reads[i+j].sigBreakPoint());


							call << reads[i].chr << "\t" << pos -1  << "\t" << Format.str() << "\t" << ref.str() <<  "\t" << alt.str() <<  "\t" << qual  << "\t" << Filter << "\t" << info.str() <<  "\t" << "GT:DP:RO:AO\t" << GenotypeField <<  endl;
							VCFOutFile << call.str(); 
							cout << call.str(); 
						}
				}
			}
		}
		if (reads[i].isSplitRead > 0 && reads[i].SVeventid ==0)
		{
			bool found = false; 
			cout << "Posible mobil single read element" << endl;
			reads[i].write();
			for (int j = 1; j < reads[i].alignments.size(); j ++)
			{
				 reads[reads[i].alignments[j]].write(); 
			}
			mobs[reads[i].name].write();	
			// do I have a breakpoint supporoted by hashes
			int bp = reads[i].sigBreakPoint();
			cout << "sig break point is " << bp << endl; 
			if (bp > 0)
			{
				cout << "Passed SigBreakpoint check " << endl; 
				
				//if (checkMobSupAalign(i, reads))
				{
					int maxSupAlign = 0; 
					for( int j = 1; j<reads[i].alignments.size(); j++)
					{
						if (reads[reads[i].alignments[j]].mapQual > maxSupAlign)
							maxSupAlign = reads[reads[i].alignments[j]].mapQual;
					}
					if (reads[i].mapQual >= maxSupAlign)	
					{
						cout << "passed masSUpAlign" << endl; 
						vector <SamRead> temp;
						for( int j = 1; j<reads[i].alignments.size(); j++)
						{
							if (reads[reads[i].alignments[j]].mapQual > 30)
								temp.push_back(reads[reads[i].alignments[j]]);
						}
						cout << "built thing" << endl; 
						int polyAbp = reads[i].isPolyA(temp);
						int myMobBase = MobAligneBases(mobs[reads[i].name], reads[i]); 
						vector <int> secondMobBase;
						int maxSecondMob = 0; 
						for( int j = 1; j<reads[i].alignments.size(); j++)
						{
							secondMobBase.push_back(MobAligneBases(mobs[reads[i].name], reads[reads[i].alignments[j]]));
							if ( secondMobBase[i] > maxSecondMob)
							{maxSecondMob = secondMobBase[i];}
						}
						cout << "checking Distnace" << endl; 
						cout << "size = " << reads[i].alignments.size()<< endl; ; 
						bool checkDistance = true; 
						for( int j = 1; j<reads[i].alignments.size(); j++)
						{
							cout << "checking " << reads[i].chr << " - " << reads[i].pos << " VS " << reads[reads[i].alignments[j]].chr << " - " << reads[reads[i].alignments[j]].pos; 
							if (reads[i].chr == reads[reads[i].alignments[j]].chr && abs(reads[i].pos - reads[reads[i].alignments[j]].pos) < 10000)
								checkDistance = false; 
							cout << " and checkDist = " << checkDistance << endl; 
						}

						if (( polyAbp > -1 ||  (myMobBase > maxSecondMob && myMobBase > 10 ))  && checkDistance)
						{
								cout << "FOUND MOB with no PolyA" << endl;
								if (reads[i].SVeventid ==0 )
								{
									int targetsize = 0; 
										int start = reads[i].pos  + bp; 
										CurrentSVeventID++;
										for(int k = 0; k < reads[i].alignments.size(); k++)
										{reads[reads[i].alignments[k]].SVeventid = CurrentSVeventID;}
										string GenotypeField;
										GenotypeField = reads[i].createStructGenotype(reads[i].sigBreakPoint());
										stringstream Format; 	
										Format << "OrphanBND"; 
										if (polyAbp > -1)
											Format << "-PolyA" << polyAbp; 
										if (myMobBase > 10)
										{
											Format << "-MOB" << myMobBase; 
											for( int j = 0; j<secondMobBase.size(); j++)
                                                                                	{
                                                                                	        Format << "+" << secondMobBase[j] ;
                                                                                	}	
										}
										Format << "-" << reads[i].MobAS; 
										Format << "LC=" << reads[i].SVCheckParentsForLowCov(reads[i].sigBreakPoint()); 
										


										string ref = Reff.getSubSequence(reads[i].chr, reads[i].pos+bp -1 ,1); 
										stringstream alt ; 
										alt << "<INS:ME:MOB>"; 
		
										string FullfilterA = reads[i].filterSV(); 
										int GMap = 0;
										int minMapQual = 30;
										if (reads[i].mapQual > minMapQual)
											GMap++;
										
										
										string InfoFilter = "";
										string Filter = "";
										if (reads[i].SVCheckParentsForLowCov(reads[i].sigBreakPoint()) >= 1)
										{
											Format << "-Inherited";
											InfoFilter = "Inherited";
											Filter = "LCH"; 
										}
										else if (GMap < 1)
										{
											Format << "-LowMapQual";
											InfoFilter = "LowMapQual";
											Filter = "LMQ"; 
										}
										else if (FullfilterA == "" )
										{
											found = true; 
											Format<<"-DeNovo";
											InfoFilter = "Pass";
											Filter = "PASS"; 
										}
										else
										{
											Format<<"-"<<FullfilterA; 
											InfoFilter = FullfilterA; 
											Filter = "fail"; 
										}
										//make quality stuff
										int readAmut=0;
										int readApos=0;
										reads[i].GetQualityHashes(readAmut, readApos, bp);
										
										float qual = -100; 
										if ((readApos) > 0)
											qual = ((float)readAmut) / ((float)readApos) * 100.0; 
										else
											qual = 0; 
		
										
										//buildng up info field
										stringstream info; 
										info << "SVTYPE=INS;END=" <<  reads[i].pos+bp -1 << ";";
										info << "MT=" << reads[i].MobContig << ";"; 
										string phase="none"; 
										if (reads[i].phase != "none")
											phase = reads[i].phase; 
										info << "PH=" << phase << ";"; 
										info << "FEX=" << InfoFilter << ";";
										int SupportingHashes = readAmut;
										int possibleHashes = readApos; 
										info << "FS=" << SupportingHashes << "/" <<possibleHashes << ";"; 
										 
										info << "RN=" << reads[i].name << ";";
										info << "MQ=" << reads[i].mapQual  << ";";
										info << "cigar=" << reads[i].cigar << ";";
										info << "SB=" << reads[i].StrandBias  << ";"; 
										info << "AS=" << reads[i].AlignmentSegments << "-" << reads[i].AlignmentSegmentsCigar << "_and_" ;
	
										cout << "in that one " << endl; 
											
										cout << reads[i].chr << "\t" << reads[i].pos+bp -1 << "\t" << Format.str() << "\t" << ref <<  "\t" << alt.str() <<  "\t" << qual  << "\t" << Filter << "\t" << info.str() <<  "\t" << "GT:DP:RO:AO\t" << GenotypeField <<  endl;
										VCFOutFile << reads[i].chr << "\t" << reads[i].pos+bp -1 << "\t" << Format.str() << "\t" << ref <<  "\t" << alt.str() <<  "\t" << qual  << "\t" << Filter << "\t" << info.str() <<  "\t" << "GT:DP:RO:AO\t" << GenotypeField <<  endl;
										//break;
									
									
										 
								}
							
						}
					}
				}
			}
			if (found)
				continue; 
		}
		if (reads[i].isSplitRead > 0 && reads[i].SVeventid ==0 && reads[i].alignments.size()>1)
		{
			bool found = false; 
			cout << "Last Ditch effort, lests give this something " << endl;
			reads[i].write();
			for (int j = 1; j < reads[i].alignments.size(); j ++)
			{
				 reads[reads[i].alignments[j]].write(); 
			}
			int A= -1; 
			int B= -1;
			vector<SamRead> temp; 
			for( int j = 0; j<reads[i].alignments.size(); j++)
			{
				temp.push_back(reads[reads[i].alignments[j]]);
			}
			FindFirstAndLast(temp, A, B); 
			// do I have a breakpoint supporoted by hashes
			if (A >= 0 && B >= 0 )
			{
				LastDitch( reads,  i,  A,  B, CurrentSVeventID);
				if (A !=B && (reads[reads[i].alignments[A]].sigBreakPoint() > 0 || reads[reads[i].alignments[B]].sigBreakPoint() > 0 || BreakpointInUnalignedCenter(reads[reads[i].alignments[A]], reads[reads[i].alignments[B]]) ))
					LastDitch( reads,  i,  B,  A, CurrentSVeventID);
				
			}
		}
	}

	cout << "Done with Multi contig events" << endl; 


	VCFOutFile.close(); 
	BEDOutFile.close();
	BEDBigStuff.close();
	BEDNotHandled.close();
	Invertions.close(); 
	cout << "finishing RUFUS.Interpret for " << outStub << std::endl;
	return 0; 
}
	
