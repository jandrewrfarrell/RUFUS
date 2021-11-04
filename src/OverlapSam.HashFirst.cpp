/*By ANDREW FARRELL
 * OverlapSam.cpp
 * -------------------------------------------------- 
 * Assembles k-mers containing variation into contigs
 * that represent the variant sequence
 * --------------------------------------------------
 */

#include <bitset>
#include <fstream>
#include <iostream>
#include <map>
#include <math.h>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <sys/stat.h>
#include <sys/types.h>
#include <sys/time.h>
#include <time.h>
#include <unordered_map>
#include <vector>


#include "Util.h"

using namespace std;
unordered_map<unsigned long, int> Mutations;
int HashSize = -1; 
bool FullOut = false;
unordered_map<string, bool> DupCheck; 
int Align3(vector<string>& sequenes, vector<string>& quals, vector<string>& strand, string phase, string Ap, string Aqp, int Ai, int& overlap, int& index, float minPercentpassed, bool& PerfectMatch, int MinOverlapPassed, int Threads) 
{
	int QualityOffset = 33; //=64; 
	int MinQual = 20;
	bool verbose = false;
	int bestScore = 0;
	int NumReads = sequenes.size();
	int start = Ai + 1;
	int end = start + 10; 
	if (end > sequenes.size()) 
	{
		end = sequenes.size();
	}

	#pragma omp parallel for shared(Ap, Aqp, index, overlap, bestScore) num_threads(Threads)
	for (int j = start; j < end; j++) 
	{
		bool run = false; 

		size_t found = strand[j].find(".");
		
		if (phase == "wh" && found == string::npos)
		{run = true;}
		if (phase == "nh" && found != string::npos)
		{run = true;}
		if (phase == "all")
		{run = true;}
		

			if (run)
			{
			int MinOverlap = MinOverlapPassed;
			float minPercent = minPercentpassed;
			int LocalBestScore = 0;
			int LocalIndex = -1;
			int LocalOverlap = 0;
			string A;
			int Alen;
			string Aq;
			//#pragma omp critical
			{
				A = Ap;
				Alen = A.length();
				Aq = Aqp;
			}
			string B;
			string Bq;
			int Blength = -1;
			int Alength = Alen;
			int k;
	
			//#pragma omp critical
			{
				B = sequenes[j];
				Bq = quals[j];
			}
			Blength = B.length();
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
				window = Blength;
				longest = Alength;
				Asmaller = true;
			}
	
			int MM = window - (window * minPercent);
			int Acount = 0;
			int Bcount = 0;
	
			for (int i = 0; i <= longest - window; i++) 
			{
				float score = 0;
				//first check where the reads completely overlap
				for (k = 0; k < window; k++) 
				{
					if (A.c_str()[k + Acount] == B.c_str()[k + Bcount]) 
					{
						if (B.c_str()[k + Bcount] != 'N' && (int)Aq.c_str()[k + Acount] > 5 && (int)Bq.c_str()[k + Bcount] > 5) 
						{
							score++;
						}
					}
					if ((k - score) > MM) 
					{
						score = -1;
						break;
					}
				}
	
				if (Asmaller) 
				{
					Acount++;
				} else {
					Bcount++;
				}
	
				if (verbose) 
				{
					cout << "	Score = " << score << endl;
				}
	
				float percent = score / (window);
				if (percent >= minPercent) 
				{
					if (LocalBestScore < score) 
					{
						LocalBestScore = score;
						LocalIndex = j;
						if (Asmaller) 
						{
							LocalOverlap = i * -1;
						} else {
							LocalOverlap = i;
						}
					}
					if (score == window) {
						PerfectMatch = true;
						break;
					}
				}
			}
	
			if (PerfectMatch == false) 
			{
				for (int i = window - 1; i >= MinOverlap; i--) 
				{
					if (verbose) {cout << "i = " << i << endl;}
	
					float score = 0;
					for (k = 0; k <= i; k++) 
					{
						if (verbose) {cout << "	k = " << k << " so A = " << Alength - i + k<< " /\\ B = " << 0 + k << endl;}
						if (verbose) {cout << "		 A >> " << A.c_str()[Alength - i + k - 1] << "="<< B.c_str()[0 + k] << " << B" << endl;}
						
						if (A.c_str()[Alength - i + k - 1] == B.c_str()[0 + k]) 
						{
							if (B.c_str()[0 + k] != 'N' && (int)Aq.c_str()[Alength - i + k - 1] > 5 &&(int)Bq.c_str()[0 + k] > 5) 
							{
								score++;
							}
						}
						if ((k - score) > MM) {
							score = -1;
							break;
						}
					}
					if (verbose) {cout << "	Score = " << score << endl;}
	
					float percent = score / (k);
					if (percent >= minPercent) 
					{
						if (LocalBestScore < score) 
						{
							LocalBestScore = score;
							LocalIndex = j;
							LocalOverlap = i - Alength + 1;
							if (score == i) 
							{
								break;
							}
						}
					}
				}
	
				for (int i = window - 1; i >= MinOverlap; i--) 
				{
					if (verbose) {	cout << "i = " << i << endl;}
					float score = 0;
					for (k = 0; k <= i; k++) {
						if (B.c_str()[Blength - i + k - 1] == A.c_str()[0 + k]) 
						{
							if (A.c_str()[0 + k] != 'N' && (int)Aq.c_str()[0 + k] > 5 && (int)Bq.c_str()[Blength - i + k - 1] > 5) 
							{
								score++;
							}
						}
						if ((k - score) > MM) {
							score = -1;
							break;
						}
					}
	
					if (verbose) {cout << "	Score = " << score << endl;}
	
					float percent = score / (k);
					if (percent >= minPercent) 
					{
						if (LocalBestScore < score) 
						{
							LocalBestScore = score;
							LocalIndex = j;
							LocalOverlap = Blength - i - 1;
							if (score == i) 
							{
								break;
							}
						}
					}
				}
			}
			#pragma omp critical(updateCounts)
			{
				if (bestScore < LocalBestScore) {
					bestScore = LocalBestScore;
					index = LocalIndex;
					overlap = LocalOverlap;
				}
			}
		}
	}
	return bestScore;
}

string ColapsContigs(string A, string B, int k, string Aq, string& Bq,
										 string& Ad, string& Bd, string& As, string& Bs) {
	bool verbose = false;
	if (verbose) {
		cout << "Combinding; \n" << A << endl << B << endl;
	}

	int Asize = A.length();
	int Bsize = B.length();

	int Aoffset = 0;
	int Boffset = 0;
	int window;
	string newString = "";
	string newQual = "";
	string newDepth = "";

	if (k > 0) {
		Aoffset = k;
	} else {
		Boffset = abs(k);
	}

	if (verbose) {
		cout << "K = " << k << " so Aofset = " << Aoffset
				 << " and Boffset = " << Boffset << endl;
	}

	for (int i = 0; i < Asize + Bsize; i++) {
		char Abase = 'Z';
		char Bbase = 'Z';
		char Aqual = '!';
		char Bqual = '!';
		unsigned char Adep = 0;
		unsigned char Bdep = 0;

		if (((i - Aoffset) >= 0) && ((i - Aoffset) < A.length())) {
			Abase = A.c_str()[i - Aoffset];
			Aqual = Aq.c_str()[i - Aoffset];
			Adep = Ad.c_str()[i - Aoffset];
		} else {
			Abase = 'Z';
			Aqual = '!';
			Adep = 0;
		}

		if (i - Boffset >= 0 && i - Boffset < B.length()) {
			Bbase = B.c_str()[i - Boffset];
			Bqual = Bq.c_str()[i - Boffset];
			Bdep = Bd.c_str()[i - Boffset];
		} else {
			Bbase = 'Z';
			Bqual = '!';
			Bdep = 0;
		}

		if (verbose) {
			cout << "I = " << i << " Bi = " << i - Boffset << " Ai = " << i - Aoffset
					 << " thus " << Abase << "-" << Bbase << endl;
		}

		if (Abase == Bbase && Abase != 'Z') {
			newString += Abase;

			if (Aqual >= Bqual) {
				newQual += Aqual;
			} else {
				newQual += Bqual;
			}

			if ((int)Adep + (int)Bdep < 250) {
				newDepth += (Adep + Bdep);
			} else {
				newDepth += 250;
			}

		} else if (Abase == 'Z' && Bbase != 'Z') {
			newString += Bbase;
			newQual += Bqual;
			newDepth += Bdep;
		} else if (Abase != 'Z' && Bbase == 'Z') {
			newString += Abase;
			newQual += Aqual;
			newDepth += Adep;
		} else if (Abase != 'Z' && Bbase != 'Z') {

			if (Abase == 'N' && Bbase != 'N') {
				newString += Bbase;
				newQual += Bqual;
				newDepth += Bdep;
			} else if (Abase != 'N' && Bbase == 'N') {
				newString += Abase;
				newQual += Aqual;
				newDepth += Adep;
			} else if (Aqual >= Bqual) {
				newString += Abase;
				newQual += Aqual;
				newDepth += Adep;
			} else {
				newString += Bbase;
				newQual += Bqual;
				newDepth += Bdep;
			}

		} else if (Abase == 'Z' && Bbase == 'Z') {
			Bq = newQual;
			Bd = newDepth;
			break;
		} 
	}
	Bq = newQual;
	Bd = newDepth;
	Bs += As;
	return newString;
}

string TrimNends(string S, string& qual) {
	bool base = false;
	string NewS = "";
	string NewQ = "";
	for (int i = S.size() - 1; i >= 0; i--) {
		if (base) {
			NewS = S.c_str()[i] + NewS;
			NewQ = qual.c_str()[i] + NewQ;
		} else if (S.c_str()[i] != 'A' && S.c_str()[i] != 'C' &&
							 S.c_str()[i] != 'G' && S.c_str()[i] != 'T') {
		} else {
			base = true;
			NewS = S.c_str()[i] + NewS;
			NewQ = qual.c_str()[i] + NewQ;
		}
	}
	qual = NewQ;
	return NewS;
}

string TrimLowCoverageEnds(string S, string& quals, string& depth, int cutoff) {
	bool base = false;
	string NewS = "";
	string NewD = "";
	string NewQ = "";

	for (int i = S.size() - 1; i >= 0; i--) {
		if (base) {
			NewS = S.c_str()[i] + NewS;
			NewD = depth.c_str()[i] + NewD;
			NewQ = quals.c_str()[i] + NewQ;
		} else if ((int)depth.c_str()[i] > cutoff) {
			base = true;
			NewS = S.c_str()[i] + NewS;
			NewD = depth.c_str()[i] + NewD;
			NewQ = quals.c_str()[i] + NewQ;
		} 
	}

	if (NewS.size() > 1) {
		S = NewS;
		depth = NewD;
		quals = NewQ;
		base = false;
		NewS = "";
		NewD = "";
		NewQ = "";

		for (int i = 0; i < S.size(); i++) {
			if (base) {
				NewS = NewS + S.c_str()[i];
				NewD = NewD + depth.c_str()[i];
				NewQ = NewQ + quals.c_str()[i];
			} else if ((int)depth.c_str()[i] > cutoff) {
				base = true;
				NewS = NewS + S.c_str()[i];
				NewD = NewD + depth.c_str()[i];
				NewQ = NewQ + quals.c_str()[i];
			}
		}
	}
	depth = NewD;
	quals = NewQ;
	return NewS;
}

string AdjustBases(string sequence, string qual) {
	int MinQ = 10;
	int QualOffset = 32;
	string NewString = "";
	for (int i = 0; i < sequence.length(); i++) {
		if (qual.c_str()[i] - QualOffset < MinQ) {
			NewString += 'N';
		} else {
			NewString += sequence.c_str()[i];
		}
	}

	if (NewString != sequence) {
		return NewString;
	}
}

bool replace(std::string& str, const std::string& from, const std::string& to) {
	size_t start_pos = str.find(from);
	if (start_pos == std::string::npos) {
		return false;
	}
	str.replace(start_pos, from.length(), to);
	return true;
}

bool validateFASTQD(string& L1, string& L2, string& L3, string& L4, string& L5,
										string& L6) {
	if (L1.c_str()[0] != '@') {
		cout << "error header problems - " << L1 << endl;
		return false;
	}
	if (L2.size() != L4.size()) {
		cout << "error sequence and qual	problems - \n" << L2 << endl
				 << L4 << endl;
		return false;
	}
	vector<string> temp = Util::Split(L6, ' ');
	if (temp.size() != L2.size()) {
		cout << "error counts	problems - " << L2.size() << " != " << temp.size()
				 << endl;
		return false;
	}
	return true;
}

bool IsBitSet(int num, int bit) { return 1 == ((num >> bit) & 1); }

int GetReadOrientation(int flag) {
	bool is_set = IsBitSet(flag, 4);
	cout << "flag is " << flag << endl;
	cout << "orientation is " << is_set << endl;
	return is_set;
}

string FlipStrands(string strand) {
	string NewStrand = "";
	for (int i = 0; i < strand.size(); i++) {
		if (strand.c_str()[i] == '+')
			NewStrand += "-";
		else if (strand.c_str()[i] == '-')
			NewStrand += "+";
		else if(strand.c_str()[i] == '.')
			 NewStrand += ".";
	}
	return NewStrand;
}

void compresStrand(string S, int& F, int& R) {
        for (int i = 0; i < S.size(); i++) {
                if (S.c_str()[i] == '+')
                        F++;
                else if (S.c_str()[i] == '-')
                        R++;
        }
        return;
}
int CountHashes(string seq)
{
	int count=0; 
	for (int i = 0; i < seq.size()-HashSize; i++)
	{
		string hash = seq.substr(i, HashSize); 
		size_t found = hash.find("N");
		if (found == string::npos) 
		{
			if(Mutations.count(Util::HashToLong( hash)) > 0)
			{count++;}
		}
	}
	return count; 
}
int NumLowQbases(string qual, int min)
{
	int count = 0; 
	for (int i =0; i < qual.size(); i++)
	{
		if( int(qual.c_str()[i]) - 33 < min)
			count++;
	}
	return count; 
}
int main(int argc, char* argv[]) {
	float MinPercent;
	int MinOverlap;
	int MinCoverage;
	cout << "you gave " << argc << " Arguments" << endl;
	if (argc != 10) {
		cout << "ERROR, wrong numbe of arguemnts\nCall is: SAM, MinPercent, "
						"MinOverlap, MinCoverage, ReportStub, NodeStub LCcutoff Threads"
				 << endl
				 << "		You Gave\n		 File = " << argv[1]
				 << "\n		MinPercent = " << argv[2]
				 << "\n	MinOVerlap = " << argv[3]
				 << "\n	MinCoverage = " << argv[4]
				 << "\n	FileStub = " << argv[5] 
				 << "\n	NodeStub = " << argv[6]
				 << "\n	MinCov = " << argv[7]
				 << "\n	HashPath = " << argv[8]
				 << "\n	Threads = " << argv[9]
				 << endl;
		return 0;
	}
			cout << "YAY,right numbe of arguemnts\nCall is: SAM, MinPercent, "
						"MinOverlap, MinCoverage, ReportStub, NodeStub LCcutoff Threads"
				 << endl
				 << "	   You Gave\n	       File = " << argv[1]
				 << "\n	 MinPercent = " << argv[2]
				 << "\n MinOVerlap = " << argv[3]
				 << "\n MinCoverage = " << argv[4]
				 << "\n FileStub = " << argv[5] 
				 << "\n NodeStub = " << argv[6] 
				 << "\n MinCovTrimCov = " << argv[7] 
				 << "\n HashPath = " << argv[8]
				 << "\n Threads = " << argv[9]
				 << endl;	
	
	ifstream fastq;
	fastq.open(argv[1]);
	if (fastq.is_open()) {
		cout << "File open - " << argv[1] << endl;
	}	else {
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
	string HashPath = argv[8]; 
	temp = argv[9];
	int Threads = atoi(temp.c_str());
	ofstream report;
	string FirstPassFile = argv[1];
	std::stringstream ss;
	ss << argv[5] << ".fastq";
	FirstPassFile = ss.str();
	report.open(FirstPassFile.c_str());

	if (report.is_open()) {
	} else {
		cout << "ERROR, Mut-Output file could not be opened - " << FirstPassFile << endl;
		return 0;
	}

	ofstream Depreport;
	FirstPassFile += "d";
	Depreport.open(FirstPassFile.c_str());

	if (report.is_open()) {
	} else {
		cout << "ERROR, Mut-Output file could not be opened - " << FirstPassFile << endl;
		return 0;
	}

	string line;
	std::vector<string> sequenes;
	std::vector<string> qual;
	std::vector<string> depth;
	std::vector<string> strand;
	std::vector<string> Unsequenes;
	std::vector<string> Unqual;
	std::vector<string> Undepth;
	std::vector<string> Unstrand;
	int lines = -1;
	string L1;
	string L2;
	string L3;
	string L4;
	string L5;
	string L6;
	unsigned long LongHash;
	int Rejects = 0;
	string Fastqd = argv[1];
	cout << "Reading in SAM \n";
	int counter = 0;
	int unalignedCounter = 0;
	int goodreads = 0;
	int lowMapQual = 0; 
	int other = 0; 
	cout << "reading in hash list" << endl; 
	ifstream MutHashFile;
	MutHashFile.open(HashPath);
	if (MutHashFile.is_open()) {
		cout << "Parent File open - " << argv[1] << endl;
	}	// cout << "##File Opend\n";
	else {
		cout << "Error, ParentHashFile could not be opened";
		return 0;
	}


	while (getline(MutHashFile, L1)) {
		vector<string> temp;
		temp = Util::Split(L1, '\t');

		if (temp.size() == 2) {
			unsigned long b = Util::HashToLong(temp[0]);
			unsigned long revb = Util::HashToLong(Util::RevComp(temp[0]));
			Mutations.insert(pair<unsigned long, int>(b, 0));
			Mutations.insert(pair<unsigned long, int>(revb, 0));
			HashSize = temp[0].size();
		} else if (temp.size() == 4) {
			unsigned long b = Util::HashToLong(temp[3]);
			unsigned long revb = Util::HashToLong(Util::RevComp(temp[3]));
			Mutations.insert(pair<unsigned long, int>(b, 0));
			Mutations.insert(pair<unsigned long, int>(revb, 0));
			HashSize = temp[3].size();
		}
		if (temp.size() == 1) {
			temp = Util::Split(L1, ' ');
			unsigned long b = Util::HashToLong(temp[0]);
			unsigned long revb = Util::HashToLong(Util::RevComp(temp[0]));
			Mutations.insert(pair<unsigned long, int>(b, 0));
			Mutations.insert(pair<unsigned long, int>(revb, 0));
			HashSize = temp[0].size();
		}
	}

	if (HashSize == -1)
	{
		cout << "ERROR Hash Size could not be determined by the HashFile" << endl;
		return -1; 
	}
	cout << "HashSize = " << HashSize << endl; 

	while (getline(fastq, L1)) {
		counter++;
		//cout << L1 << endl; 
		if (counter % 100 == 1) {
			cout << "Read in " << counter << " reads, with " << goodreads << " aligned reads, " << unalignedCounter<< " unaligned reads, " << lowMapQual << "low map qual and " << other <<" other with rejected " << Rejects
					 << " reads\r";
		}

		vector<string> temp = Util::Split(L1, '\t');
		int ReadSize = temp[10].size();
		bool b[16];
		int v = atoi(temp[1].c_str()); 

		//TODO: understand this syntax
		for (int j = 0; j < 16; ++j) 
		{
			b[j] = 0 != (v & (1 << j));
		}
		//if (DupCheck.count(temp[9]) > 0)
		//{
		//	cout << "skipping exact match sequence " << L1 << endl; 
		//}
		//else 
		{
			int lowq = NumLowQbases(temp[10], 20); 
			DupCheck[temp[9]] == true; 
			if (b[8] or b[11] or b[10] or temp[9].length() < 100 or lowq > 50) 
			{
				//cout << "rejected" << endl; 
				//cout << L1 << endl; 
				Rejects++;
			} 
			else if ( b[2] )//or atoi(temp[4].c_str())<5) 
			{
				if (b[2])
					unalignedCounter++;
				else if (atoi(temp[4].c_str())<5)
					lowMapQual++;
				else
					other++; 
				string L4 = temp[10];
				string L2 = temp[9]; ////////sequence//////////
				L2 = TrimNends(L2, L4);
				int hashes = CountHashes(L2);
				//cout << "poorly mapped read " << L1 << endl; 
				//cout << "with Hash = " << hashes << endl; 
				if ((double)L2.size() / (double)ReadSize > .6) 
				{
					ReadSize = L2.size();
					lines++;
					Unsequenes.push_back(L2);
					Unqual.push_back(L4);
					if (hashes > 0)
					{
						if (b[0]== 0)
						{
							Unstrand.push_back(".");
						}
						else if (b[4] == 0) 
						{
							Unstrand.push_back("+");
						} 
						else if (b[4] == 1) 
						{
							Unstrand.push_back("-");
						}
					}
					else
						Unstrand.push_back("."); 

					string depths = "";
					unsigned char C = 1;
	
					for (int i = 0; i < L2.length(); i++) 
					{
						depths += C;
					}
	
					Undepth.push_back(depths);
	
				} 
				else 
				{
					Rejects++;
				}
	
			} 
			else 
			{
				goodreads++; 
				//cout << "good alignment" << endl; 
				string L4 = temp[10];
				string L2 = temp[9];
				L2 = TrimNends(L2, L4);
				int hashes = CountHashes(L2);
				
				//cout << "Hash = " << hashes << endl;
				if ((double)L2.size() / (double)ReadSize > .6) 
				{
					ReadSize = L2.size();
					lines++;
					sequenes.push_back(L2);
					qual.push_back(L4);
					string depths = "";
	
					if (hashes > 0)
					{
						if (b[0]== 0)
                                                {
                                                        strand.push_back(".");
                                                }
                                                else if (b[4] == 0) 
						{
							strand.push_back("+");
						}
						else if (b[4] == 1) 
						{
							strand.push_back("-");
						}
					}
					else
						strand.push_back(".");
	
					unsigned char C = 1;
	
					for (int i = 0; i < L2.length(); i++) 
					{
						depths += C;
					}
	
					depth.push_back(depths);
	
				} 
				else 
				{
					Rejects++;
				}
			}
		}
	}

	cout << endl;
	int NumReads = sequenes.size();
	cout << "\nDone reading in \n		 Read in a total of " << NumReads + Rejects
			 << " and rejected " << Rejects << endl;
	clock_t St, Et;
	float Dt;

	struct timeval start, end;
	gettimeofday(&start, NULL);
	int FoundMatch = 0;
	St = clock();

	for (std::vector<string>::size_type i = 0; i < sequenes.size(); i++) 
	{

		string A = sequenes[i];
		string Aqual = qual[i];
		string Adep = depth[i];
		string Astr = strand[i];
		size_t found = Astr.find(".");
                if (found == string::npos)
		{
			if (FullOut) {
				cout << "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<**************************************************************>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>" << endl;
			}
	
			Et = clock();
	
			if ((int)i % 10 < 1){
				gettimeofday(&end, NULL);
				float Dt = end.tv_sec - start.tv_sec;
				cout << "aligning " << i << " of " << NumReads << ", \% done = " << ((double)i / (double)NumReads) * 100.00<< ", TotalTime= " << Dt << " , second per read = " << Dt / i<< ", \% finding match = "<< ((double)FoundMatch / (double)i) * 100.00 << "\r";
			}
	
			if (FullOut) {
				Dt = ((double)(Et - St)) / CLOCKS_PER_SEC;
				cout << "aligning " << i << " of " << NumReads
						 << "\% done = " << ((double)i / (double)NumReads) * 100.00
						 << ", TotalTime= " << Dt << " , second per read = " << Dt / i
						 << ", \% finding match = "
						 << ((double)FoundMatch / (double)i) * 100.00 << endl
						 << A << endl;
	
				for (int z = 0; z < Adep.length(); z++) {
					int bam = Adep.c_str()[z];
					cout << bam;
				}
				cout << endl;
			}
	
			int k = -1;
			int bestIndex = -1;
			bool PerfectMatch = false;
			int booya = Align3(sequenes, qual, strand, "wh", A, Aqual, i, k, bestIndex, MinPercent, PerfectMatch, MinOverlap, Threads);
			if (FullOut) {
				cout << "best forward score is " << booya << " k is " << k << endl;
			}
	
			if (!(PerfectMatch)) {
				string revA = Util::RevComp(A);
				string revAqual = Util::RevQual(Aqual);
				string revAdep = Util::RevQual(Adep);
				string revAstr = FlipStrands(Astr);
				int revk = -1;
				int revbestIndex = -1;
				int revbooya =	Align3(sequenes, qual, strand, "wh",  revA, revAqual, i, revk, revbestIndex,
	
					 MinPercent, PerfectMatch, MinOverlap, Threads);
				if (FullOut) {
					cout << "best reverse score is " << revbooya << " k is " << revk
							 << endl;
				}
	
				if (revbooya > booya) {
					A = revA;
					Aqual = revAqual;
					Adep = revAdep;
					Astr = revAstr;
					k = revk;
					booya = revbooya;
					bestIndex = revbestIndex;
				}
	
			} else {
				if (FullOut) {
					cout << "Perfect Match Found, Skipping Referse Search" << endl;
				}
			}
	
			if (booya < MinOverlap) {
				if (FullOut) {
					cout << "No good match found, skipping" << endl;
				}
			} else {
				FoundMatch++;
				string B = sequenes[bestIndex];
				string Bqual = qual[bestIndex];
				string Bdep = depth[bestIndex];
				string Bstr = strand[bestIndex];
	
				if (k > 0) {
	
					if (FullOut) {
						cout << "found match at " << k << endl;
	
						for (int z = 0; z < k; z++) {
				cout << "+";
			}
	
						cout << A << endl << B << endl;
	
						for (int z = 0; z < Bdep.length(); z++) {
							int bam = Bdep.c_str()[z];
							cout << bam;
						}
	
						cout << endl;
					}
				} else {
					if (FullOut) {
						cout << "found match at " << k << endl;
						cout << A << endl;
	
						for (int z = 0; z < abs(k); z++) {
				cout << "-";
			}
	
						cout << B << endl;
	
						for (int z = 0; z < abs(k); z++) { 
				cout << "-";
			}
	
						for (int z = 0; z < Bdep.length(); z++) {
							int bam = Bdep.c_str()[z];
							cout << bam;
						}
	
						cout << endl;
					}
				}
	
				if (i == bestIndex) {
					cout << "ERROR ____________________ SAME READS " << endl;
				}
				if (A.size() != Adep.size() && B.size() != Bdep.size()) {
					cout << " ERRPR somethis the wrong size\n	A= " << A.size()
							 << " Ad = " << Adep.size() << " B= " << B.size()
							 << " Bd = " << Bdep.size() << endl;
				}
	
				string combined =
		ColapsContigs(A, B, k, Aqual, Bqual, Adep, Bdep, Astr, Bstr);
	
				if (combined.size() != Bdep.size()) {
					cout << " ERRPR combined is the wrong size\n	C= " << combined.size()
							 << " Bd = " << Bdep.size() << endl;
				}
	
				sequenes[bestIndex] = combined;
				qual[bestIndex] = Bqual;
				depth[bestIndex] = Bdep;
				strand[bestIndex] = Bstr;
				sequenes[i] = "moved";
	
				if (FullOut) {
					cout << combined << endl;
	
					for (int z = 0; z < Bdep.length(); z++) {
						int bam = Bdep.c_str()[z];
						cout << bam;
					}
				}
			}
		}
	}

	cout << "\n\nRESULTS\n";
	int count = 0;
	for (int i = 0; i < sequenes.size(); i++) {

		if (FullOut) {
			cout << ">>>>" << sequenes[i] << endl;
		}

		if (sequenes[i] != "moved" && sequenes[i].size() >= 95) {
			string rDep = depth[i];
			int maxDep = -1;

			for (int z = 0; z < rDep.size(); z++) {
				unsigned char bam = rDep.c_str()[z];
				if ((int)bam > maxDep) {
					maxDep = (int)bam;
				}
			}

			if (maxDep >= MinCoverage && maxDep >= 2) {

				if (sequenes[i].size() != qual[i].size() && qual[i].size() != depth[i].size()) {
						cout << "ERROR, read "
							 << "@NODE_" << argv[6] << "_" << i << "_L=" << sequenes[i].size()
							 << "_D=" << maxDep
							 << " Has the wrong size, Seq = " << sequenes[i].size()
							 << " Qual = " << qual[i].size() << " Dep = " << depth[i].size()
							 << endl;
				}

				count++;
				int F = 0;
                                int R = 0;
                                compresStrand(strand[i], F, R);
				report << "@NODE_" << argv[6] << "_" << i << "_L=" << sequenes[i].size() << "_D=" << maxDep  << ":" << F << ":" << R << ":" << endl;
				report << sequenes[i] << endl;
				report << "+" << endl;
				report << qual[i] << endl;

				Depreport << "@NODE_" << argv[6] << "_" << i<< "_L=" << sequenes[i].size() << "_D=" << maxDep  << ":" << F << ":" << R << ":" << endl;
				Depreport << sequenes[i] << endl;
				Depreport << "+" << endl;
				Depreport << qual[i] << endl;
				Depreport << strand[i] << endl;
				unsigned char C = depth[i].c_str()[0];
				int booya = C;
				Depreport << booya;

				for (int w = 1; w < depth[i].size(); w++) {
					C = depth[i].c_str()[w];
					booya = C;
					Depreport << " " << booya;
				}

				Depreport << endl;
			}
		}
	}
	if (MinCoverage <= 1  )
	{
		for (int i = 0; i < Unsequenes.size(); i++) {
	
			if (FullOut) {
				cout << ">>>>" << Unsequenes[i] << endl;
			}
	
			if (Unsequenes[i] != "moved" && Unsequenes[i].size() >= 95) {
				int maxDep = -1;
	
				if (Unsequenes[i].size() != Unqual[i].size() &&
						Unqual[i].size() != Undepth[i].size()) {
					cout << "ERROR, read "
							 << "@NODE_" << argv[6] << "_" << i << "_L=" << Unsequenes[i].size()
							 << "_D=" << maxDep
							 << " Has the wrong size, Seq = " << Unsequenes[i].size()
							 << " Qual = " << Unqual[i].size() << " Dep = " << Undepth[i].size()
							 << endl;
				}
	
				count++;
				report << "@NODE_" << argv[6] << "_" << i << "_L=" << Unsequenes[i].size()
							 << "_D" << maxDep << endl;
				report << Unsequenes[i] << endl;
				report << "+" << endl;
				report << Unqual[i] << endl;
	
				Depreport << "@NODE_" << argv[6] << "_" << i
									<< "_L=" << Unsequenes[i].size() << "_D" << maxDep << endl;
				Depreport << Unsequenes[i] << endl;
				Depreport << "+" << endl;
				Depreport << Unqual[i] << endl;
				Depreport << Unstrand[i] << endl;
				unsigned char C = Undepth[i].c_str()[0];
				int booya = C;
				Depreport << booya;
	
				for (int w = 1; w < Undepth[i].size(); w++) {
					C = Undepth[i].c_str()[w];
					booya = C;
					Depreport << " " << booya;
				}
	
				Depreport << endl;
			}
		}
	}
	else
		cout << "min coverage = " << MinCoverage << " skipping Unaligned sequences" << endl; 
	cout << "\nWrote " << count << " sequences" << endl;
	report.close();
}
