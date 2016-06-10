//hosdie, how are you 
/*By ANDREW FARRELL 
 * the Marth Lab BC
 * 6.3 = 6.2 + testing sectioning the file 
 * now accepts nfastqd as input
 * 7 = hash based searching method for matches
 */
 
#include <algorithm>
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
#include<sys/stat.h>
#include<sys/types.h>
#include <unordered_map>
#include <map>
#include <sys/time.h>
#include <unistd.h>
using namespace std;
bool FullOut = false;
///Call is BitHashCompare Parent Mutant firstpassfile hashsize

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
string trim (string s)
{
	string newS = "";
	cout << s << endl;
	for (int i = 0; i < s.size(); i++)
	{
		if ((int)s.c_str()[i]>=33)
		{
			newS = newS+s.c_str()[i];
			cout << s.c_str()[i] << endl;
		}
	}
	cout << newS << endl;
return newS;
}

string LongToHash (unsigned long LongHash, int HashSize)
{
	string value = "";
	bitset<64> test (LongHash);
	for (int i = 1; i < HashSize*2; i+=2)
	{
		if (test[i-1] == 0)
		{
			if (test[i] == 0)
			{value = value + "A";}
			else
			{value = value + "C";}
		}
		else
		{
			if (test[i] == 0)
			{value = value + "G";}
			else
			{value = value + "T";}
		}
	}
	return value;
}


string RevComp (string Sequence)
{
	string NewString = "";
	//cout << "Start - " << Sequence << "\n";
	for(int i = Sequence.size()-1; i>=0; i+= -1)
	{
		char C = Sequence.c_str()[i];
	//      cout << C << endl;
		if (C == 'A')
			{NewString += 'T';}
		else if (C == 'C')
		       	{NewString += 'G';}
		else if (C == 'G')
			{NewString += 'C';}
		else if (C == 'T')
			{NewString += 'A';}
		else if (C == 'N')
			{NewString += 'N';}
		else
			{cout << "\nERROR IN RevComp - " << C << "\n";}
	}
	//cout << "end\n";
	return NewString;
}
string RevQual (string Sequence)
{
	string NewString = "";
	for(int i = Sequence.size()-1; i>=0; i+= -1)
	{
		
		unsigned char C = Sequence.c_str()[i];
		NewString += C;
	}
	return NewString;
}

unsigned long HashToLong (string hash, string calledby)
{
	//cout << "booya" << hash << endl;
	bitset<64> HashBits;
	int bitcout = 0;
	for(int i=0; i<hash.size();i++)
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
			cout << "ERROR, invalid character - " << hash.c_str()[i] << ", in hash " << hash << ", called by " << calledby <<  endl;
		}
	}
	//cout <<  HashBits.to_ulong() << "-" << endl;
	return  HashBits.to_ulong();
}
/*	
unsigned long HashToLong (string hash)
{
	//cout << "booya" << hash << endl;
	bitset<64> HashBits;
	//#pragma omp parellel for
	for(int i=0; i<hash.size();i++)
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
			cout << "ERROR, invalid character - " << hash.c_str()[i] << " in " << hash <<  endl;
		}
	}
	//cout <<  HashBits.to_ulong() << "-" << endl;
	return HashBits.to_ulong();
}
*/

int RebuildHashTable(vector<string>& sequenes, int Ai, int SearchHash, unordered_map<unsigned long, vector<int>>& Hashes, int Threads)
{
	cout << "\nDestrying HashTable\n";
	Hashes.clear();
	cout << "HashTable destroyed\n";
	cout << "Rebuilding HashTable - starting at " << Ai << endl;
	int size = sequenes.size();
 	#pragma omp parallel for num_threads(6) shared(Hashes)
	for(int i = Ai; i < size; i++)
	{
		if (i % 10000 > 1 && i % 10000 < Threads  ){cout << "   Hashed " << i << " of " << sequenes.size() << "\r";}
		
		string Sequence =  sequenes[i];
		//cout << Sequence << endl;
		int LoopLimit = Sequence.size()-SearchHash;
		
		for( int j = 0; j < LoopLimit; j++)
		{
			string hash = Sequence.substr(j,SearchHash);
			size_t found = hash.find('N');
  			if (found==std::string::npos)
			{
				unsigned long LongHash = HashToLong(hash, "Rebuild");
				unsigned long  RevHash = HashToLong(RevComp(hash), "Rebuild");
				#pragma omp critical(updateHash)
				{Hashes[LongHash].push_back(i);
				Hashes[RevHash].push_back(i);}
			}
	       	}
	}
	cout << "\nDone Rebulding HashTable\n";
	return 0;
}
int  PrepairSearchList(string A, int Ai, unordered_map<unsigned long, vector<int>>& Hashes,int SearchHash,int ACT, map<int, vector<int>>& array, bool& hitPosLimit, bool& hitIndexLimit, int& NumberPos, int& NumberIndex )
{
	int Alength = A.size();
	map<int,int> Positions;
	int added = 0 ;
	for (int i = 0; i<Alength-SearchHash; i++)
	{
		string hash = A.substr(i,SearchHash);
		size_t found = hash.find('N');
		if (found==std::string::npos)
		{
			unsigned long LongHash = HashToLong(hash, "Align");
			int max = Hashes[LongHash].size();
			for (vector<int>::size_type i = 0; i < max; i++)
			{
				int holder = Hashes[LongHash][i];
				if(holder > Ai+1 )
				{
					if (Positions.count(holder) > 0)
					{
						{Positions[holder]++;added++;}
					}
					else
					{
						{Positions[holder] = 1;added++;}
					}
				}
				if (added > 100000)
				{hitPosLimit = true; NumberPos = added;break;}
			}
		}
	}
									if(FullOut){cout << "done Hasing read" << endl;}
	NumberPos = added;
	map <int,int>::iterator uspos;
	map<double, int> SortedPositions;
	double trick = 0.0;
	for (uspos = Positions.begin(); uspos != Positions.end(); ++uspos)
	{
		trick++;
		SortedPositions[(double)uspos->second + (2.0/trick)] = uspos->first;
	}
									if(FullOut){cout << "found - " << Positions.size() << " possible locations" << endl;}
	map<double, int>::iterator pos;
	vector<int> indexes;
	int sanity = 0;
	for (pos = SortedPositions.end(); pos != SortedPositions.begin(); --pos)
	{
		if (pos->first >= ACT)
		{
			indexes.push_back(pos->second+0);
			sanity++;
			if (sanity > 1000){hitIndexLimit = true; NumberIndex = sanity;break;}
		}
	}
	NumberIndex = sanity;
	#pragma omp critical
	{array[Ai] = indexes;}
	return 1;
									if(FullOut){cout << "       " << indexes.size() << " locations passed filter" << endl;}


}
int Align3 (vector<string>& sequenes, string Ap, string Aq, int Ai, int& overlap, int& BestIndex, float minPercent, bool& PerfectMatch, int MinOverlap, vector<int>& indexes, int Threads, int NumReads)
{
	int QualityOffset = 33;

	bool verbose = false;
	int Alength = Ap.size();
	int bestScore = 0;
	
	#pragma omp parallel for  num_threads(Threads) shared(BestIndex) 	
	for (int booya = 0 ; booya < indexes.size(); booya++)
	{
		string A = Ap;
		int AlengthL = A.size();
		int j = indexes[booya];	
		
		string B, Bq;
		bool localcheck;
	       	B = sequenes[j];

		
		float score = 0 ; 
		int Blength =  B.size();
		int k;
		
		int window = -1;
		int longest = -1;
		bool Asmaller = true;
		if (Blength > AlengthL)
		{
			window = AlengthL;
			longest = Blength;
			Asmaller = false;	
		}
		else 
		{
			Asmaller = true;
			window = Blength;
			longest = AlengthL;
		}	
		int MM = window - (window*minPercent);
		int LbestScore = -1;
		int LBestIndex = -1;
		int Loverlap = 0;

		//checking bits that are within the two reads 
		int Acount = 0;
		int Bcount = 0;
		for (int i = 0; i<= longest - window ;i++)
		{
			score = 0;
			
			
			for (k = 0; k<window; k++) //base compare loop
			{
				if(A.c_str()[k+Acount] == B.c_str()[k+Bcount])
				{if (A.c_str()[k+Acount ] != 'N'){score++;}}
				if ((k-score) > MM)
				{score = -1;break;}
			  
			}	
			if (Asmaller)
			{Acount++;}
			else 
			{Bcount++;}
			
			if(verbose){cout << "  Score = " << score << endl;}
			
			float percent = score/(window);
			if (percent >= minPercent) 
			{
				if (LbestScore < score)
				{
					LbestScore = score;
					LBestIndex = j;
					if (Asmaller)
					{Loverlap =i*-1;}
					else 
					{Loverlap = i;}
					
					if (score == window)
					{		
						PerfectMatch = true;
						break;
					}
				}
			}
		}
		
		//first part of overhang
		//cout << "window = " << window << " Percent = " << minPercent << " MM = " << MM << endl;
		if (PerfectMatch == false)
		{
			for (int i = window-1; i>=MinOverlap; i--)
			{
				if(verbose){cout << "i = " << i << endl;}
				score = 0;
				for(k = 0; k <= i; k++)
				{
					if(verbose){cout << "  k = " << k << " so A = " <<  AlengthL - i + k << " /\\ B = " <<0+k<<  endl;}
					if(verbose){cout << "     A >> " << A.c_str()[ AlengthL - i + k -1] << "=" << B.c_str()[ 0+k] << " << B" << endl;}
				
					if (A.c_str()[ AlengthL - i + k - 1 ] == B.c_str()[ 0+k] )
					{if (B.c_str()[ 0+k] != 'N'){score++;}}
					if ((k-score) > MM)
					{score = -1;break;}
				}
				if(verbose){cout << "  Score = " << score << endl;}
				float percent = score/(k);
				if (percent >= minPercent) 
				{
					if (LbestScore < score)
					{
						LbestScore = score;
						LBestIndex = j;
						Loverlap = i-AlengthL +1;
						
						if (score == i)
						{break;}
					
					}
				}
			}	
			//check other side of overhang
		
			for (int i = window-1; i>=MinOverlap; i--)
			{
				if(verbose){cout << "i = " << i << endl;}
			
				score = 0;
				for(k = 0; k <= i; k++)
				{
					if (B.c_str()[ Blength - i + k - 1 ] == A.c_str()[ 0+k])
					{if (A.c_str()[ 0+k] != 'N'){score++;}}
					 if ((k-score) > MM)
					      	{score = -1;break;}
				}
				if(verbose){cout << "  Score = " << score << endl;}
				
				float percent = score/(k);
				if (percent > minPercent) 
				{
					if (LbestScore < score)
					{
						LbestScore = score;
						LBestIndex = j;
						Loverlap = Blength-i-1;
						
						if (score == i)
						{break;}
					}			
				}
			}
		}
		#pragma omp critical
		{
			if (LbestScore > bestScore)
			{
				bestScore = LbestScore;
				BestIndex = LBestIndex;
				overlap = Loverlap;
			}
		}
	}
	return bestScore;	
}

string ColapsContigs(string A, string B, int k, string Aq, string& Bq, string Ad, string& Bd)
{
	bool verbose = false;
	if (verbose){cout << "Combining; \n" << A << endl << B << endl;} 

	int Asize = A.size();
	int Bsize = B.size();
	
	int Aoffset = 0;
	int Boffset = 0;
	int window;
	string newString = "";
	string newQual = "";
	string newDepth = "";
	if (k>0)
	{Aoffset = k;}else{Boffset = abs(k);}
	if (verbose){cout << "K = " << k << " so Aofset = " << Aoffset << " and Boffset = " << Boffset << endl;}
	for (int i = 0; i<Asize+Bsize;i++)
	{
		char Abase = 'Z';
		char Bbase = 'Z';
		char Aqualc = '!';
		char Bqualc = '!';
		unsigned char Adepc = 0;
		unsigned char Bdepc = 0;

		//Get furrent Bases
		if (((i-Aoffset) >= 0 ) && ((i-Aoffset) < A.size()))
		{
			Abase = A.c_str()[i-Aoffset];
			Aqualc = Aq.c_str()[i-Aoffset];
			Adepc = Ad.c_str()[i-Aoffset];
		}
		else 
		{
			Abase = 'Z';
			Aqualc = '!';
			Adepc = 0;	
		}

		if (i-Boffset >= 0 && i-Boffset < B.size())
		{
			Bbase = B.c_str()[i-Boffset];
			Bqualc = Bq.c_str()[i-Boffset];
			Bdepc = Bd.c_str()[i-Boffset];
		}
		else
		{
			Bbase = 'Z';
			Bqualc = '!';
			Bdepc = 0;
		}
		if (verbose){cout << "I = " << i << " Bi = " << i-Boffset << " Ai = " << i-Aoffset <<  " thus " << Abase << "-" << Bbase<< endl;}
		
		//Compare Current Bases
		if (Abase == Bbase && Abase != 'Z' )
		{
			newString+=Abase;
			if (Aqualc >= Bqualc){newQual+=Aqualc;}else{newQual+=Bqualc;}
			if ((int)Adepc+(int)Bdepc<250){newDepth+=(Adepc+Bdepc);}else{newDepth += (char)250;}
		}
		else if(Abase == 'Z' && Bbase != 'Z')
		{
			newString += Bbase;
			newQual += Bqualc;
			newDepth += Bdepc;
		}
		else if(Abase !='Z' && Bbase == 'Z' )
		{
			newString += Abase;
			newQual += Aqualc;
			newDepth += Adepc;
		}
		else if (Abase !='Z' && Bbase != 'Z')
		{
			/*if (Abase == 'N' && Bbase != 'N')
			{
				newString+=Bbase;
			       	newQual+=Bqual;
		    		newDepth+=Bdep;
				cout << "   A = N" << endl;
			}
			else if (Abase != 'N' && Bbase == 'N')
			{
				newString+=Abase;
				newQual+=Aqual;
				newDepth+=Adep;
				cout << "   B = N" << endl;
			}
			else*/ if ( Adepc > Bdepc)
			{
				newString+=Abase;
				newQual+=Aqualc;
				newDepth+=Adepc;
			}
			else if( Adepc < Bdepc)
			{
                                newString+=Bbase;
                                newQual+=Bqualc;
                                newDepth+=Bdepc;
                        }
			else if ( Aqualc >= Bqualc)
                        {
                                newString+=Abase;
                                newQual+=Aqualc;
                                newDepth+=Adepc;
                        }
			else 
			{
				newString+=Bbase;
				newQual+=Bqualc;
				newDepth+=Bdepc;
			}

		}
		else if (Abase =='Z' && Bbase == 'Z')
		{Bq = newQual; Bd = newDepth;return newString;}
	}
	Bq = newQual;
	Bd = newDepth;
	return newString;

	
}


string TrimNends(string S, string& qual)
{
	bool base = false;
	string NewS = "";
	string NewQ = "";
	for (int i = S.size()-1; i>=0; i--)
	{
		if (base)
		{
			NewS = S.c_str()[i] + NewS;
			NewQ = qual.c_str()[i] +NewQ;
		}
		else if(S.c_str()[i] != 'A' && S.c_str()[i] != 'C' && S.c_str()[i] != 'G' && S.c_str()[i] != 'T')
		{}
		else
		{
			base = true;
			NewS = S.c_str()[i] + NewS;
			NewQ = qual.c_str()[i] +NewQ;
		}
	}
	

	S = NewS;
	qual = NewQ;
	
	base = false;
	NewS = "";
	NewQ = "";
	for (int i = 0; i<S.size(); i++)
	{
		if (base)
		{
			NewS = NewS + S.c_str()[i] ;
			NewQ = NewQ + qual.c_str()[i] ;
		}
		else if(S.c_str()[i] != 'A' && S.c_str()[i] != 'C' && S.c_str()[i] != 'G' && S.c_str()[i] != 'T')
		{}
		else
		{
			base = true;
			NewS = NewS + S.c_str()[i];
			NewQ = NewQ + qual.c_str()[i] ;
		}
	}



	//cout << "Trimmed Laggin N\'s " << endl << S <<endl<< NewS<<endl;
	qual = NewQ;
	return NewS;
}

string TrimLowCoverageEnds(string S,string& quals, string& depth, int cutoff)
{
 	//cout << "S = " << S << endl;
       	bool base = false;
	string NewS = "";
	string NewD = "";
	string NewQ = "";
	//cout << "Starting at left edge" << endl;
	for (int i = S.size()-1; i>=0; i--)
	{
		//cout << "Base " << i << " = " << S.c_str()[i] << " dep = " << (int)depth.c_str()[i] << endl;
		if (base)
		{
		//	cout << "found hihgh base" << endl;
			NewS = S.c_str()[i] + NewS;
			NewD = depth.c_str()[i] + NewD;
			NewQ = quals.c_str()[i] + NewQ;
		}
		else if((int)depth.c_str()[i]>cutoff)
		{
		//	cout << (int)depth.c_str()[i] << " > " << cutoff << "reading rest of bases" << endl;
			base = true;
			NewS = S.c_str()[i] + NewS;
			NewD = depth.c_str()[i] + NewD;
			NewQ = quals.c_str()[i] + NewQ;
		}
		else
		{}//cout << (int)depth.c_str()[i] << " still too low " << endl;}
	}
	if (NewS.size() > 1)
	{
		S = NewS;
		depth = NewD;
		quals = NewQ;
		//cout << "so far trimmed to " << S << endl;
		base = false;
		NewS = "";
		NewD = "";
		NewQ = "";
		for (int i = 0; i<S.size(); i++)
		{
			if (base)
			{
				NewS = NewS + S.c_str()[i] ;
				NewD = NewD + depth.c_str()[i];
				NewQ = NewQ + quals.c_str()[i];
			}
			else if((int)depth.c_str()[i]>cutoff)
			{
				base = true;
				NewS = NewS + S.c_str()[i];
				NewD = NewD + depth.c_str()[i];
				NewQ = NewQ + quals.c_str()[i];
			}
		}
	}
		//cout << "Trimmed Laggin N\'s " << endl << S <<endl<< NewS<<endl;
	depth = NewD;
	quals = NewQ;
	//cout << "returning S = " << NewS << endl <<  "returning D = " << NewD << endl;
	return NewS;
}
string AdjustBases(string sequence, string qual)
{
	int MinQ = 5;
	int QualOffset = 33;
	string NewString = "";
	for(int i=0; i<sequence.size(); i++)
	{
		//cout << qual.c_str()[i]-QualOffset << " " ;
		if (qual.c_str()[i]-QualOffset < MinQ)
		{NewString += 'N';}
		else
		{NewString+=sequence.c_str()[i];}
	}
	if (NewString != sequence)
	//{cout <<"Assinged New String" << endl<<sequence<<endl<<NewString<<endl;}
	return NewString;
}
bool replace(std::string& str, const std::string& from, const std::string& to) {
    size_t start_pos = str.find(from);
    if(start_pos == std::string::npos)
	return false;
    str.replace(start_pos, from.size(), to);
    return true;
}
int main (int argc, char *argv[])
{
	cout << "Testing Trim \n";
	string test = "ANGT";
	string qu   = "4567";
	cout << test << endl << qu << endl;
	test = 	TrimNends(test, qu);
	cout << test << endl << qu << endl << "done" << endl;
	
	cout << "Testing HashToLong \n";
	test = "CACCACCGGCAAGCTGCCCGTGCCCTGCC";
	unsigned long testLong = HashToLong(test, "test");
	string test2 = LongToHash(testLong, test.size());
	cout << "String = " << test << ", hash = " << testLong << endl << "String2= " << test2 << endl;	

//return 0;
	int SearchHash = 30;
	int ACT = 0;
	float MinPercent ;
	int MinOverlap ;
	int MinCoverage ;			
	cout << "you gave "<< argc << " Arguments"  << endl;
	if (argc != 11)
	{cout << "ERROR, wrong numbe of arguemnts\nCall is: FASTQ, MinPercent, MinOverlap, MinCoverage, ReportStub, SearchHashSize, ACT, OutFile LCendTrimEpth Threads"<< endl; return 0;}
		
	
	ifstream fastq;
	fastq.open (argv[1]);
	if ( fastq.is_open())
	{ cout << "Parent File open - " << argv[1] << endl;}	//cout << "##File Opend\n";
	else
	{
		cout << "Error, ParentHashFile could not be opened";
		return 0;
	}
	
	string temp = argv[2];
	MinPercent = atof(temp.c_str());	

	temp = argv[3];
	MinOverlap = atoi(temp.c_str());

	temp = argv[4];
	MinCoverage = atoi(temp.c_str());
	
	temp = argv[6];
	SearchHash = atoi(temp.c_str());
	
	temp = argv[7];
	ACT = atoi(temp.c_str());
	
	temp = argv[9];
	int TrimLCcuttoff = atoi(temp.c_str());
	temp = argv[10];
	int Threads = atoi(temp.c_str());
	int Buffer = 100 * Threads;
	ofstream report;
	std::stringstream ss;
	string FirstPassFile = argv[1];
	ss << argv[8] << ".fastq";
	FirstPassFile = ss.str();
	report.open (FirstPassFile.c_str());
	if (report.is_open())
	{}
	else
	{
		cout << "ERROR, Mut-Output file could not be opened - " << FirstPassFile << endl;
		return 0;
	}
	
	ofstream Depreport;
	FirstPassFile += "d";
	Depreport.open (FirstPassFile.c_str());
	if (report.is_open())
	{}
	else
	{
		cout << "ERROR, Mut-Output file could not be opened - " << FirstPassFile << endl;
		return 0;
	}

	ofstream good;
	FirstPassFile = ss.str(); 
	FirstPassFile += "good.fastq";
	good.open (FirstPassFile.c_str());
	if (good.is_open())
	{}
	else
	{
		cout << "ERROR, Mut-Output file could not be opened - " << FirstPassFile << endl;
		return 0;
	}
	ofstream bad;
	FirstPassFile = ss.str();
	FirstPassFile += "bad.fastq";
	bad.open (FirstPassFile.c_str());
	if (bad.is_open())
	{}
	else
	{
		cout << "ERROR, Mut-Output file could not be opened - " << FirstPassFile << endl;
		return 0;
	}
	
	
	string line;
	
	vector<string> sequenes;// = new vector<string>;
	vector<string> qual;//= new vector<string>;
	vector<string> depth;// = new vector<string>;
		
	std::unordered_map<unsigned long, vector<int>> Hashes;
	

	int lines = -1;
	int goodlines = 0; 
	int dup = 0; 
	string L1;
	string L2;
	string L3;
	string L4;
	string L5;
	string L6;
	int Rejects = 0;
	
	//read in fastq
	//int readcount = 0;
	string Fastqd = argv[1];
	size_t found = Fastqd.find(".fastqd");
	if (found !=string::npos)	
	{
		
		int counter =0;
		cout << "ATTENTION - Fastq+depth input detected, reading in FASTQD file \n";
		while (getline(fastq, L1))
		{
			counter++;
			if (counter % 100 == 1)
			{
				cout << "Read in " << counter << " lines and rejected " << Rejects << " reads\r";
			}

			getline(fastq, L2);
			getline(fastq, L3);
			getline(fastq, L4);
			getline(fastq, L5);
			getline(fastq, L6);
			string depths = "";
			int ReadSize = L2.size();
				
			bool Multiple = false;
		       	vector<string> temp = Split(L6, ' ');
			 		
			for (vector<string>::size_type i = 0; i < temp.size(); i++)
			{unsigned char C = atoi(temp[i].c_str());depths += C;if((int)C>1){Multiple = true;}}
			
			Multiple = false;	
			if (Multiple == true){L2 = TrimLowCoverageEnds(L2, L4, depths,TrimLCcuttoff);}
			
			if (L2.size() > SearchHash+1)
			{
				lines++;	
				sequenes.push_back(L2);
				qual.push_back(L4);
				depth.push_back(depths);
				ReadSize = L2.size();
			}
			else
			{Rejects++;}
		}
		

	}
	else
	{
		vector <string> DupCheck; 
		cout << "Reading in raw fastq \n"; 
		int counter = 0;	
		while (getline(fastq, L1))
		{
			counter++;
			if (counter % 100 == 1)
			{
				cout << "Read in " << counter << " lines, rejected " << Rejects << " reads with " << dup << " duplicates\r";
			}
			//if (lines == 10000)
			//{break;}
			
			getline(fastq, L2);
			getline(fastq, L3);
			getline(fastq, L4);
			int Ns = 0; 
			for (int i =0; i<L2.size(); i++)
			{
				if (L2.c_str()[i] == 'N')
					{Ns++;}
			}
			
			lines++;
			int ReadSize = L2.size();
				
			string depths = "";	
			ReadSize = L2.size()-1;
			bool found = false;
			
			bool RunDupCheck = true; 
			if (RunDupCheck)
			{
				#pragma omp parallel for  num_threads(Threads) shared(DupCheck, L2, found)
				for(int i = 0; i < DupCheck.size(); i++) 
				{
					if (L2.size() == DupCheck[i].size())
					{
						bool AllBasesMatch = true;
						for (int k = 0; k < L2.size(); k++)
						{
							if (L2.c_str()[k] == 'N' or DupCheck[i].c_str()[k] == 'N'){}
							else if (L2.c_str()[k] == DupCheck[i].c_str()[k]){}
							else
							{
								AllBasesMatch = false;
								break;
							}
						}

						if (AllBasesMatch)
						{
							#pragma omp critical
							{found = true;}
						}
					}
				}
			
				if ((double)Ns/(double) L2.size() < 0.20)
				{DupCheck.push_back(L2);}	
			}
			if (found ==  false)
			{
				L2 = AdjustBases(L2, L4);
				L2 = TrimNends(L2, L4);
				if ((double)L2.size()/(double)ReadSize > .6)
				{
					goodlines++;
					sequenes.push_back(L2);
					qual.push_back(L4);
					unsigned char C = 1;
					for (int i=0; i<=ReadSize;i++)
					{depths += C;}
					depth.push_back(depths);
					good << L1 << endl << L2 << endl << L3 << endl << L4 << endl;
				} 
				else
	       			{Rejects++;}	
			}
			else
			{
				dup++;
				bad << L1 << endl << L2 << endl << L3 << endl << L4 << endl;
			}					
		}
		DupCheck.clear();
	}	
	good.close();
	bad.close();    
	cout << "done reading "<<endl;
	int NumReads = sequenes.size();
	cout << "\nDone reading in \n     Read in a total of " << lines << " and rejected "<<Rejects  << " with " << dup << " duplicate reads detected for a total of " << goodlines << "good reads" << endl; 
	
	RebuildHashTable(sequenes, 0, SearchHash, Hashes, Threads);
	
	clock_t St,Et;
	int FoundMatch = 0;
	struct timeval start, end;	
	gettimeofday(&start, NULL);
	int LinesSinceLastBuild = 1;

	int NumberHitPosLimit = 0;
	int NumberHitSanityLimit = 0;
	double AverageFPos = 0.0;
	double AverageFSanity = 0.0;
	double AverageRPos = 0.0;
	double AverageRSanity = 0.0;
	
	for(std::vector<string>::size_type b = 0; b < sequenes.size(); b+=Buffer)
	{
		LinesSinceLastBuild+= Buffer;
		if ( LinesSinceLastBuild > 1000000 )
		{
			RebuildHashTable(sequenes, b, SearchHash, Hashes, Threads);
			LinesSinceLastBuild = 0;
		}
		vector<string>ToAddHashes;
		vector<int>ToAddPos;
		map<int, vector<int>>Forwards;
		map<int, vector<int>>Revs;
		int max = b+Buffer;
	
		if (max > sequenes.size())
		{max = sequenes.size();}
								if (FullOut){cout << "Bulding list to align" << endl;}	
		#pragma omp parallel for  num_threads(Threads) shared(Hashes, Forwards) 
		for (int i = b; i<max; i++)
		{
			string A = sequenes[i];
			bool posLimit = false;
			bool sanityLimit = false;
			int NumPos = 0;
			int NumSanity = 0;
			PrepairSearchList(A, i, Hashes, SearchHash, ACT, Forwards, posLimit, sanityLimit, NumPos, NumSanity);
			if (posLimit){NumberHitPosLimit++;}
			if(sanityLimit){NumberHitSanityLimit++;}
			AverageFPos = ((AverageFPos *(double) b)+(double)NumPos)/((double)b+1.0);
			AverageFSanity = ((AverageFSanity * (double)b)+(double)NumSanity)/((double)b+1.0);
		}
		#pragma omp parallel for  num_threads(Threads) shared( Hashes, Revs) 
		for (int i = b; i<max; i++)
		{
			string A = RevComp(sequenes[i]);
			bool posLimit = false;
			bool sanityLimit = false;
			int NumPos = 0;
			int NumSanity = 0;
			PrepairSearchList(A, i, Hashes, SearchHash, ACT, Revs, posLimit, sanityLimit, NumPos, NumSanity);
			if (posLimit){NumberHitPosLimit++;}
			if(sanityLimit){NumberHitSanityLimit++;};
			AverageRPos = ((AverageRPos * (double)b)+(double)NumPos)/((double)b+1.0);
			AverageRSanity = ((AverageRSanity * (double)b)+(double)NumSanity)/((double)b+1.0);
		}
								if (FullOut){cout << "Done Bulding List"<< endl;}
		for (int i = b; i<max; i++)
		{
			string A, Aqual, Adep;
			A = sequenes[i];
			Aqual = qual[i];
			Adep = depth[i];
			int k;
			bool PerfectMatch = false;
			int bestIndex = -1;
			
			float Dt;
								if(FullOut){Et = clock();
								Dt = ((double) (Et - St)) / CLOCKS_PER_SEC;
								cout << "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<**************************************************************>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"<< endl;}


			if (i % 100 == 5)
			{
				gettimeofday(&end, NULL);
				float Dt = end.tv_sec - start.tv_sec;
				cout << "done " << i << " of " << NumReads << ", \% = " << ((double)i/(double)NumReads)*100.00 <<", TT= " << Dt << " , S/R = " << Dt/i << ", \% match= " << ((double)FoundMatch/(double)i)*100.00  << " %Plimit = " << NumberHitPosLimit/2 << ", " << ((double)NumberHitPosLimit/2.0)/(double)i << " AvP= " << (AverageFPos +AverageRPos)/2.0 << " % ILmit = " << NumberHitSanityLimit << ", " << ((double)NumberHitSanityLimit/2.0)/(double)i << " AvI= " << (AverageRSanity+AverageFSanity)/2.0  << "\r";
			}	

	
			int booya = Align3(sequenes, A, Aqual , i,k, bestIndex,  MinPercent , PerfectMatch, MinOverlap, Forwards[i], Threads, NumReads);
								if(FullOut){cout << "best forward score is " << booya << " k is " << k  << " index = " << bestIndex << endl; }
		
			if (!(PerfectMatch))
			{
				string revA = RevComp(A);
				string revAqual = RevQual(Aqual);
				string revAdep = RevQual(Adep);
				int revk = -1;
				int revbestIndex = -1;
								if(FullOut){cout << "Checking Reverse\n";}
				int revbooya =Align3(sequenes, revA, revAqual, i,revk, revbestIndex, MinPercent, PerfectMatch, MinOverlap,Revs[i], Threads, NumReads);
								if(FullOut){cout << "best reverse score is " << revbooya << " k is " << revk  << " index = " << revbestIndex  << endl;} 
			
				if (revbooya>booya)
				{
					A = revA;
					Aqual = revAqual;
					Adep = revAdep;
					k = revk;
					booya = revbooya;
					bestIndex = revbestIndex;
				}
			}
			else
								{if(FullOut){cout<<"Perfect Match Found, Skipping Reverse Search"<<endl;}}
			
			if (booya < MinOverlap )
			{
								if(FullOut){cout << "No good match found, skipping"<< endl;}
			}
			else	
			{
				string B, Bqual, Bdep;
				B = sequenes[bestIndex];
				Bqual = qual[bestIndex];
				Bdep = depth[bestIndex];
				FoundMatch++;
				
												if(FullOut){if(k>0){cout << "found match at "<< k  << endl; for( int z = 0; z<k; z++)cout << "+";cout << A << endl << B  << endl;for (int z = 0; z < Bdep.size();z++){int bam = Bdep.c_str()[z];cout << bam;}cout << endl;}else {cout << "found match at "<< k  << endl;cout << A << endl;for( int z = 0; z<abs(k); z++)cout << "-";cout << B  << endl;for( int z = 0; z<abs(k); z++)cout << "-";for (int z = 0; z < Bdep.size();z++){int bam = Bdep.c_str()[z];cout << bam;}cout << endl;}}
				
				string combined = ColapsContigs(A, B, k, Aqual, Bqual, Adep, Bdep);
				
							if (Bqual.size() != combined.size()){cout << "ERRRORRR ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^" << endl;}	
				qual[bestIndex] = Bqual;
				depth[bestIndex] = Bdep;
				sequenes[bestIndex] = combined;
				sequenes[i] = "moved";
			
				//cout << "staring lookup" << endl;
			
				#pragma omp parallel for num_threads(Threads) shared(Hashes)
				for(  int j = 0; j < A.size()-SearchHash; j++)
				{
					string hash = A.substr(j,SearchHash);
					size_t found = hash.find('N');
					if (found==std::string::npos)
					{
						bool found = false;
						int k = 0; 
                                                for ( k = 0; k < Hashes[HashToLong(hash, "update")].size(); k++)
                                                {
                                                	if (Hashes[HashToLong(hash, "update")][k] == bestIndex)
							{
								found = true; 
								break; 		
							}
						}
						//cout << "took " << k << "lookups" << endl; 
						if (found = false)
						{
							#pragma omp critical
							{	
								Hashes[HashToLong(hash, "update")].push_back(bestIndex);
							//	Hashes[HashToLong(RevComp(hash), "upsate")].push_back(bestIndex);
							}
						}
						found = false;
                                                for ( k = 0; k < Hashes[HashToLong(RevComp(hash), "update")].size(); k++)
                                                {
                                                        if (Hashes[HashToLong(RevComp(hash), "update")][k] == bestIndex)        
                                                        {
                                                                found = true;
                                                                break;
                                                        }
                                                }
						//cout << "took " << k << "lookups" << endl;
						if (found = false)
						{
                                                	#pragma omp critical
                                                	{       
                                                	//      Hashes[HashToLong(hash, "update")].push_back(bestIndex);
                                                		Hashes[HashToLong(RevComp(hash), "upsate")].push_back(bestIndex);
                                                	}
						}
					}
				}							if(FullOut){cout << combined << endl;cout << Bqual << endl;for (int z = 0; z < Bdep.size();z++){int bam = Bdep.c_str()[z];cout << bam;}cout << endl;}			
			}
		}
	}
	
	cout << "\nRESULTS\n";
	int count = 0;
	for (int i = 0; i < sequenes.size(); i++)
	{
		
		//cout << ">>>>" << sequenes[i] << endl;
		if (sequenes[i] !="moved" && sequenes[i].size() >= 95)
		{
			string rDep = depth[i];
			int maxDep = -1;
			for (int z = 0; z < rDep.size();z++)
			{
				unsigned char bam = rDep.c_str()[z];
				if ((int)bam > maxDep){maxDep = (int)bam;} 
			}
			if (maxDep >=MinCoverage)
			{
				//if (!(sequenes[i].size() == qual[i].size() && sequenes[i].size() == depth[i].size()))
				//{cout <<"ERRROR lengths arnt equal \n   " << "@NODE_" << i << "_L=" <<sequenes[i].size() << "_D=" << maxDep <<  "\n   "<< sequenes[i] << "\n   " << qual[i]<< endl;}
				count ++;
				report << "@NODE_" << i << "_L=" <<sequenes[i].size() << "_D=" << maxDep <<  endl;
				report << sequenes[i] << endl;
				report << "+" << endl;
				report << qual[i] << endl;

				Depreport << "@NODE_" << i <<  "_L=" << sequenes[i].size() << "_D=" << maxDep << endl;
				Depreport << sequenes[i] << endl;
				Depreport << "+" << endl;
				Depreport << qual[i] << endl;
				Depreport << '+' << endl;
				unsigned char C = depth[i].c_str()[0];
				int booya = C;
				Depreport << booya;
				for (int w = 1; w < depth[i].size(); w++)
				{
			
					C = depth[i].c_str()[w];
					booya = C;
					Depreport << ' ' << booya;
				}
				Depreport << endl;
			}
			//else{cout << "Rejected" << endl;}
		}
	}
	cout << "Wrote " << count << " sequences" << endl;
	report.close();
	Depreport.close();
	
}
	
