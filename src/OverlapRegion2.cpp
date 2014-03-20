/*By ANDREW FARRELL 
 * the Marth Lab BC
 * 6.3 = 6.2 + testing sectioning the file 
 * now accepts fastqd as input
 * 7 = hash based searching method for matches
 */

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
#include <sys/time.h>
#include <unordered_map>
#include <map>
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
	return HashBits.to_ulong();
}
string LongToHash (unsigned long LongHash)
{
	bitset<64> HashBits;
//	char test [33];
//	itoa (4,test,2);
//	cout << test << endl;
//	cint >> test;
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
string RevQual (string Sequence)
{
	string NewString = "";
	for(int i = Sequence.length()-1; i>=0; i+= -1)
	{
		unsigned char C = Sequence.c_str()[i];
		if (C != '\0')
		{NewString += C;}
	}
	return NewString;
}


int Align3 (vector<string>& sequenes, vector<string>& quals, string Ap, string Aqp, int Ai, int& overlap, int& index, float minPercentpassed, bool& PerfectMatch, int MinOverlapPassed, int Threads)
{
	//cout << "started\n";
	int QualityOffset = 64;
	int MinQual = 20; 

	bool verbose = false;
	//string B;
	//string Bq;
	int bestScore = 0;
	//float score = 0;
	//int k, j;
	int NumReads = sequenes.size();
	//int count = 0;
	int start = Ai+1;
	int end = sequenes.size();
	
	#pragma omp parallel for  shared(index, overlap, bestScore) num_threads(Threads)  
	for(int j = start; j < end; j++)
	{
		int MinOverlap = MinOverlapPassed;
		float minPercent = minPercentpassed;
		int LocalBestScore = 0;
		int LocalIndex = -1;
		int LocalOverlap = 0; 
		string A = Ap;
		int Alen = A.length();
		string Aq = Aqp;
		string B;
        	string Bq;
		int Blength = -1;
		int Alength = Alen;
		int k;
		//count ++;
		//cout << "comapring to " << j << "\r" ;
		//#pragma omp critical(read)
		//{
			B = sequenes[j];
			Bq = quals[j];
                	Blength = B.length();
		//}
		//cout  << B << endl << Bq<< endl;
		
		int window = -1;
		int longest = -1;
		bool Asmaller = true;
		if (Blength > Alength)
		{
			window = Alength;
			longest = Blength;
			Asmaller = false;	
		}
		else 
		{
			Asmaller = true;
			window = Blength;
			longest = Alength;
		}	
		int MM = window - (window*minPercent);
		
		//checking bits that are within the two reads 
		int Acount = 0;
		int Bcount = 0;
		//#pragma omp parallel for
		for (int i = 0; i<= longest - window ;i++)
		{
			float score = 0;
			
			
			for (k = 0; k<window; k++) //base compare loop
			{
				if(A.c_str()[k+Acount ] == B.c_str()[k+Bcount])
				{if ( B.c_str()[k+Bcount] != 'N'){score++;}}
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
				if (LocalBestScore < score)
				{
					LocalBestScore = score;
                                        LocalIndex = j;
                                        if (Asmaller)
                                        {LocalOverlap =i*-1;}
                                        else
                                        {LocalOverlap = i;}
					/*#pragma omp critical(updateCounts)
            				{
						bestScore = score;
						index = j;
						if (Asmaller)
						{overlap =i*-1;}
						else 
						{overlap = i;}
					}*/
					
				}
				if (score == window)
				{
					//cout << "Found exct match, exiting" << endl;
					PerfectMatch = true;
					break;
					//return bestScore;
				}
			}
		}
		
		//first part of overhang
		//cout << "window = " << window << " Percent = " << minPercent << " MM = " << MM << endl;
		//#pragma omp parallel for
		if (PerfectMatch == false)
		{
			for (int i = window-1; i>=MinOverlap; i--)
			{
				if(verbose){cout << "i = " << i << endl;}
					
				float score = 0;
				for(k = 0; k <= i; k++)
				{
					if(verbose){cout << "  k = " << k << " so A = " <<  Alength - i + k << " /\\ B = " <<0+k<<  endl;}
					if(verbose){cout << "     A >> " << A.c_str()[ Alength - i + k -1] << "=" << B.c_str()[ 0+k] << " << B" << endl;}
					if (A.c_str()[ Alength - i + k - 1 ] == B.c_str()[ 0+k] )
					{if ( B.c_str()[ 0+k] != 'N'){score++;}}
					if ((k-score) > MM)
					{score = -1;break;}
				}
				if(verbose){cout << "  Score = " << score << endl;}
				
				float percent = score/(k);
				//cout << score << "/" << k << " = " << percent << endl;
				if (percent >= minPercent) 
				{
					//				cout << score << "/" << k << " = " << percent << endl;
					if (LocalBestScore < score)
					{
							LocalBestScore = score;
	                                                LocalIndex = j;
	                                                LocalOverlap = i-Alength +1;
						/*#pragma omp critical(updateCounts)
	            				{
							bestScore = score;
							index = j;
							overlap = i-Alength +1;
						}*/
						if (score == i)
						{break;}
						
					}
				}
			}
			//check other side of overhang
			//#pragma omp parallel for
			for (int i = window-1; i>=MinOverlap; i--)
			{
				if(verbose){cout << "i = " << i << endl;}
				
				float score = 0;
				for(k = 0; k <= i; k++)
				{
					if (B.c_str()[ Blength - i + k - 1 ] == A.c_str()[ 0+k] )
					{if (A.c_str()[ 0+k] != 'N'){score++;}}
					 if ((k-score) > MM)
	                               	{score = -1;break;}
				}
				if(verbose){cout << "  Score = " << score << endl;}
					
				float percent = score/(k);
				if (percent >= minPercent) 
				{
		//			cout << score << "/" << k << " = " << percent << endl;
					if (LocalBestScore < score)
					{
							LocalBestScore = score;
	                                                LocalIndex = j;
	                                                LocalOverlap = Blength-i-1;
	
						/*#pragma omp critical(updateCounts)
	            				{
							bestScore = score;
							index = j;
							overlap = Blength-i-1;
						}*/
						if (score == i)
	                                       	{break;}
					}
				}
			}
		}
		#pragma omp critical(updateCounts)
                {
			if (bestScore < LocalBestScore)
			{
                 		bestScore = LocalBestScore;
                       	 	index = LocalIndex;
                        	overlap = LocalOverlap;
			}
                }
	}
	return bestScore;
}

string ColapsContigs(string A, string B, int k, string Aq, string& Bq, string& Ad, string& Bd)
{
	bool verbose = false;
	if (verbose){cout << "Conbining; \n" << A << endl << B << endl;} 

	int Asize = A.length();
	int Bsize = B.length();
	
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
		char Aqual = '!';
		char Bqual = '!';
		unsigned char Adep = 0;
		unsigned char Bdep = 0;

		//Get furrent Bases
		if (((i-Aoffset) >= 0 ) && ((i-Aoffset) < A.length()))
		{
			Abase = A.c_str()[i-Aoffset];
			Aqual = Aq.c_str()[i-Aoffset];
			Adep = Ad.c_str()[i-Aoffset];
		}
		else 
		{
			Abase = 'Z';
			Aqual = '!';
			Adep = 0;	
		}

		if (i-Boffset >= 0 && i-Boffset < B.length())
		{
			Bbase = B.c_str()[i-Boffset];
			Bqual = Bq.c_str()[i-Boffset];
			Bdep = Bd.c_str()[i-Boffset];
		}
		else
		{
			Bbase = 'Z';
			Bqual = '!';
			Bdep = 0;
		}
		if (verbose){cout << "I = " << i << " Bi = " << i-Boffset << " Ai = " << i-Aoffset <<  " thus " << Abase << "-" << Bbase<< endl;}
		
		//Compare Current Bases
		if (Abase == Bbase && Abase != 'Z')
		{
			newString+=Abase;
			if (Aqual >= Bqual){newQual+=Aqual;}else{newQual+=Bqual;}
			if ((int)Adep+(int)Bdep<250){newDepth+=(Adep+Bdep);}else{newDepth += 250;}
		}
		else if(Abase == 'Z' && Bbase != 'Z')
		{
			newString += Bbase;
			newQual += Bqual;
			newDepth += Bdep;
		}
		else if(Abase !='Z' && Bbase == 'Z')
		{
			newString += Abase;
			newQual += Aqual;
			newDepth += Adep;
		}
		else if (Abase !='Z' && Bbase != 'Z')
                {
                        if (Abase == 'N' && Bbase != 'N')
                        {
                                newString+=Bbase;
                                newQual+=Bqual;
                                newDepth+=Bdep;
                        }
                        else if (Abase != 'N' && Bbase == 'N')
                        {
                                newString+=Abase;
                                newQual+=Aqual;
                                newDepth+=Adep;
                        }
                        else if ( Aqual >= Bqual)
                        {
                                newString+=Abase;
                                newQual+=Aqual;
                                newDepth+=Adep;
                        }
                        else
                        {
                                newString+=Bbase;
                                newQual+=Bqual;
                                newDepth+=Bdep;
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
                //      cout << "found hihgh base" << endl;
                        NewS = S.c_str()[i] + NewS;
                        NewD = depth.c_str()[i] + NewD;
                        NewQ = quals.c_str()[i] + NewQ;
                }
                else if((int)depth.c_str()[i]>cutoff)
                {
                //      cout << (int)depth.c_str()[i] << " > " << cutoff << "reading rest of bases" << endl;
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

/*string TrimLowCoverageEnds(string S,string& quals, string& depth, int cutoff)
{
 	cout << "S = " << S << endl;
       	bool base = false;
        string NewS = "";
        string NewD = "";
	string NewQ = "";
	cout << "Starting at left edge" << endl;
        for (int i = S.size()-1; i>=0; i--)
        {
		cout << "Base " << i << " = " << S.c_str()[i] << " dep = " << (int)depth.c_str()[i] << endl;
                if (base)
                {
			cout << "found hihgh base" << endl;
                        NewS = S.c_str()[i] + NewS;
                        NewD = depth.c_str()[i] + NewD;
			NewQ = quals.c_str()[i] + NewQ;
                }
                else if((int)depth.c_str()[i]>cutoff)
                {
                	cout << (int)depth.c_str()[i] << " > " << cutoff << "reading rest of bases" << endl;
		        base = true;
                        NewS = S.c_str()[i] + NewS;
                        NewD = depth.c_str()[i] + NewD;
			NewQ = quals.c_str()[i] + NewQ;
                }
		else
		{cout << (int)depth.c_str()[i] << " still too low " << endl;}
        }
	S = NewS;
	depth = NewD;
	quals = NewQ;
	cout << "so far trimmed to " << S << endl;
	base = false;
        NewS = "";
        NewD = "";
	NewQ = "";
        for (int i = 0; i<=S.size()-1; i++)
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
        cout << "Trimmed Laggin N\'s " << endl << S <<endl<< NewS<<endl;
	depth = NewD;
	quals = NewQ;
        cout << "returning S = " << NewS << endl <<  "returning D = " << NewD << endl;
	
	return NewS;
}*/
string AdjustBases(string sequence, string qual)
{
	int MinQ = 10;
	int QualOffset = 32;
	string NewString = "";
	for(int i=0; i<sequence.length(); i++)
	{
//		cout << qual.c_str()[i]-QualOffset << " " ;
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
    str.replace(start_pos, from.length(), to);
    return true;
}
bool validateFASTQD(string& L1, string& L2, string& L3, string& L4, string& L5, string& L6)
{
	if (L1.c_str()[0] != '@')
	{
		cout << "error header problems - " << L1 << endl;
		return false;
	}
	if (L2.size() != L4.size())
	{
		cout << "error sequence and qual  problems - \n" << L2 << endl << L4 << endl;
		return false;
	}
	vector<string> temp = Split(L6, ' ');
	if (temp.size() != L2.size())
	{
		cout << "error counts  problems - " << L2.size() << " != " << temp.size()  << endl;
		return false;
	}
	return true;
}
int main (int argc, char *argv[])
{
	cout << "USING THIS ONE" << endl;		
	float MinPercent ;
	int MinOverlap ;
	int MinCoverage ;			
	cout << "you gave "<< argc << " Arguments"  << endl;
	if (argc != 9)
	{cout << "ERROR, wrong numbe of arguemnts\nCall is: FASTQ, MinPercent, MinOverlap, MinCoverage, ReportStub, NodeStub LCcutoff Threads"<< endl << "     You Gave\n     File = " << argv[1] << "\n     MinPercent = " 
		<< argv[2] << "\n     MinOVerlap = " << argv[3] << "\n     MinCoverage = " << argv[4] << "\n     FileStub = " << argv[5] << "\n     NodeStub = " << argv[6] <<  endl; return 0;}
	
	cout << "     You Gave\n     File = " << argv[1] << "\n     MinPercent = "
                << argv[2] << "\n     MinOVerlap = " << argv[3] << "\n     MinCoverage = " << argv[4] << "\n     FileStub = " << argv[5] << "\n     NodeStub = " << argv[6] << "\n     Threads = " << argv[7] <<  endl;	
	
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
	
	temp = argv[7];
	int LCcutoff = atoi(temp.c_str());
	
	temp = argv[8];
	int Threads = atoi(temp.c_str());
	
	ofstream report;
	string FirstPassFile = argv[1]; 

	std::stringstream ss;
	//	ss << FirstPassFile << "_mp" << MinPercent<< "mo" << MinOverlap << "mc" << MinCoverage << "_"<< argv[5] << ".fastq";
	ss << argv[5] << ".fastq";
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
	
	
	string line;
	
	std::vector<string> sequenes;
	std::vector<string> qual;
	std::vector<string> depth;
	
	

	int lines = -1;

	string L1;
	string L2;
	string L3;
	string L4;
	string L5;
	string L6;
	unsigned long LongHash;
	int Rejects = 0;
	
	//read in fastq
	string Fastqd = argv[1];
	size_t found = Fastqd.find(".fastqd");
        if (found !=string::npos)	
	{
		
		int counter = -1;
                cout << "ATTENTION - Fastq+depth input detected, reading in FASTQD file \n";
                while (getline(fastq, L1))
                {
                        //lines++;
                        getline(fastq, L2);
                        getline(fastq, L3);
                        getline(fastq, L4);
                        getline(fastq, L5);
                        getline(fastq, L6);
                        string depths = "";
                        int ReadSize = L2.size();
			if (validateFASTQD(L1, L2, L3, L4, L5, L6))
			{
                        	bool Multiple = false;
                       	 	std::vector<string> temp = Split(L6, ' ');
                        	for (std::vector<string>::size_type i = 0; i < temp.size(); i++)
                        	{
					int e = atoi(temp[i].c_str());
					unsigned char C = e;
					depths += C;
					if((int)C>1){Multiple = true;}
				}
				
                        	if (Multiple == true){L2 = TrimLowCoverageEnds(L2, L4, depths, LCcutoff);}

                                if (L2.size() > 60)
                                {
                                        lines++;
                                        sequenes.push_back(L2);
                                        qual.push_back(L4);
                                        depth.push_back(depths);
                                        ReadSize = L2.size();
                                }
                	}	
			else
			{
				cout << "ERROR in FASTQD file \n   " << L1 << "\n   " << L2 << "\n   "<< L3 << "\n   "<< L4 << "\n   "<< L5 << "\n   "<<L6<< endl;
				return 1;
			}	
		}
	}
	else
	{
		cout << "Reading in raw fastq \n"; 
		int counter = 0;	
		while (getline(fastq, L1))
        	{
			counter++;
			if (counter % 100 == 1)
                	{
                        	cout << "Read in " << counter << " lines and rejected " << Rejects << " reads\r";
                	}
			//if (lines == 100)
			//{break;}
			
			getline(fastq, L2);
			getline(fastq, L3);
			getline(fastq, L4);
			int ReadSize = L2.size();
			//size_t found = L2.find("N");
                        //if (found ==string::npos)
                        //{
				
			//L2 = AdjustBases(L2, L4);
			L2 = TrimNends(L2, L4);
			 if ((double)L2.size()/(double)ReadSize > .6)
                        {
                                ReadSize = L2.size();
                                lines++;
                                sequenes.push_back(L2);
                                qual.push_back(L4);
                                string depths = "";
                                unsigned char C = 1;
                                for (int i=0; i<L2.length();i++)
                                {depths += C;}
                                depth.push_back(depths);


                        }
                        else
                        {Rejects++;}	

		}
	}               
	cout << endl;
	int NumReads = sequenes.size();
	cout << "\nDone reading in \n     Read in a total of " << NumReads+Rejects << " and rejected "<<Rejects << endl; 
	//for(std::vector<string>::size_type i = 0; i < sequenes.size(); i++)
	//	cout << sequenes[i]<< endl;
	//cout << "done" << endl;
	clock_t St,Et;
	float Dt;
	
	struct timeval start, end;
	 gettimeofday(&start, NULL);
	int FoundMatch = 0;	
	St = clock();
	//int breakcount = 0;
	for(std::vector<string>::size_type i = 0; i < sequenes.size(); i++)
	{
                        //
		//breakcount++;
		//if (breakcount==20){return 0;}
		string A = sequenes[i];
		string Aqual = qual[i];
		string Adep = depth[i];
		
									if(FullOut){cout << "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<**************************************************************>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"<< endl;}
		Et = clock();
	 	//if ( (int)i%10 < 1 )
                //{
			gettimeofday(&end, NULL);
                        float Dt = end.tv_sec - start.tv_sec;
	       		cout << "aligning " << i << " of " << NumReads << ", \% done = " << ((double)i/(double)NumReads)*100.00 <<", TotalTime= " << Dt << " , second per read = " << Dt/i<< ", \% finding match = " << ((double)FoundMatch/(double)i)*100.00  << "\r";
                //}
									if(FullOut){Dt = ((double) (Et - St)) / CLOCKS_PER_SEC;cout << "aligning " << i << " of " << NumReads << "\% done = " << ((double)i/(double)NumReads)*100.00 <<", TotalTime= " << Dt << " , second per read = " << Dt/i << ", \% finding match = " << ((double)FoundMatch/(double)i)*100.00  << endl <<A << endl;
									for (int z = 0; z < Adep.length();z++){int bam = Adep.c_str()[z];cout << bam;}
									cout << endl;}
		int k = -1;
		int bestIndex = -1;
		bool PerfectMatch = false;
		int booya = Align3(sequenes, qual, A, Aqual , i,k, bestIndex, MinPercent , PerfectMatch, MinOverlap, Threads);
		if(FullOut){cout << "best forward score is " << booya << " k is " << k  << endl; }
	
		if (!(PerfectMatch))
		{
			string revA = RevComp(A);
			string revAqual = RevQual(Aqual);
			string revAdep = RevQual(Adep);
			int revk = -1;
			int revbestIndex = -1;
			int revbooya =Align3(sequenes, qual, revA, revAqual, i,revk, revbestIndex, MinPercent, PerfectMatch, MinOverlap, Threads);
			if(FullOut){cout << "best reverse score is " << revbooya << " k is " << revk  << endl;} 
		
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
		{if(FullOut){cout<<"Perfect Match Found, Skipping Referse Search"<<endl;}}
		if (booya < MinOverlap )
		{if(FullOut){cout << "No good match found, skipping"<< endl;}}
		
		else
		{
			//cout << "here " << bestIndex<< endl;

			FoundMatch++;
			string B = sequenes[bestIndex];
                 	string Bqual = qual[bestIndex];
			//cout << sequenes.size() << " - " << qual.size() << " - " << depth.size() << endl;
			string Bdep = depth[bestIndex];
			//cout << "here 2" << endl;
			if(k>0)
			{
				if(FullOut){cout << "found match at "<< k  << endl; 
				for( int z = 0; z<k; z++)
				cout << "+";
				cout << A << endl << B  << endl;
				 for (int z = 0; z < Bdep.length();z++){int bam = Bdep.c_str()[z];cout << bam;}cout << endl;}
			}
			else {
				if(FullOut){cout << "found match at "<< k  << endl;
				cout << A << endl;
				for( int z = 0; z<abs(k); z++)
					cout << "-";
				cout << B  << endl;
				for( int z = 0; z<abs(k); z++)
                                        cout << "-";
				for (int z = 0; z < Bdep.length();z++){int bam = Bdep.c_str()[z];cout << bam;}cout << endl;}
			}
			
			if (i == bestIndex )
			{ cout << "ERROR ____________________ SAME READS " << endl;}
			if (A.size() != Adep.size() && B.size() != Bdep.size()){cout <<" ERRPR somethis the wrong size\n  A= " << A.size() << " Ad = " << Adep.size() << " B= " << B.size() << " Bd = " << Bdep.size()<<endl;}
			string combined = ColapsContigs(A, B, k, Aqual, Bqual, Adep, Bdep) ;
			 if (combined.size() != Bdep.size()){cout << " ERRPR combined is the wrong size\n  C= " << combined.size() << " Bd = " << Bdep.size()<<endl;}
			sequenes[bestIndex] = combined;
			qual[bestIndex] = Bqual;
			depth[bestIndex] = Bdep;
			sequenes[i] = "moved";
			if(FullOut){cout << combined << endl;
			cout << Bqual << endl;
			 for (int z = 0; z < Bdep.length();z++){int bam = Bdep.c_str()[z];cout << bam;}
			cout << endl;
			}
			//cout << "Done colapsing\n";
		}
		//cout << "done with this read \n";
	}
	
	cout << "\n\nRESULTS\n";
	int count = 0;
	for (int i = 0; i < sequenes.size(); i++)
	{
		
		if(FullOut){cout << ">>>>" << sequenes[i] << endl;}
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
				if (sequenes[i].size() !=  qual[i].size() &&  qual[i].size() != depth[i].size()){cout << "ERROR, read " << "@NODE-" << argv[6] << "_" << i << "_L=" <<sequenes[i].size() << "_D=" << maxDep << " Has the wrong size, Seq = " << sequenes[i].size() << " Qual = " << qual[i].size() << " Dep = " <<  depth[i].size() << endl;}
				count ++;
				report << "@NODE-" << argv[6] << "_" << i << "_L=" <<sequenes[i].size() << "_D=" << maxDep <<  endl;
				report << sequenes[i] << endl;
				report << "+" << endl;
				report << qual[i] << endl;

				Depreport << "@NODE_" << argv[6] << "_" << i <<  "_L=" << sequenes[i].size() << "_D=" << maxDep << endl;
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
					Depreport << " " << booya;
				}
				Depreport << endl;
			}
		}
	}
	cout << "\nWrote " << count << " sequences" << endl;
	report.close();
	
}
	
