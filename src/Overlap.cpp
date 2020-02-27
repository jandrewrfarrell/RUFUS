/*By ANDREW FARRELL
 * Overlap.cpp
 * --------------------------------------------------
 * Assembles k-mers containing variation into contigs
 * that represent the variant sequence
 * --------------------------------------------------
 */

#include <algorithm>
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
#include <unistd.h>
#include <unordered_map>
#include <vector>
#include <omp.h> 

#include "Util.h"

using namespace std;

bool FullOut = false;

int RebuildHashTable(vector<string>& sequenes, int Ai, int SearchHash, unordered_map<unsigned long, vector<int>>& Hashes, int Threads, unordered_map<unsigned long, int>& Hashesize) 
{
	cout << "\nDestrying HashTable\n";
	Hashes.clear();
	cout << "HashTable destroyed\n";
	cout << "Rebuilding HashTable - starting at " << Ai << endl;
	int size = sequenes.size();

	#pragma omp parallel for num_threads(12) shared(Hashes)
	for (int i = Ai; i < size; i++) {

		if (i % 10000 > 1 && i % 10000 < Threads) {
			//pragma omp critical (sequenes)
			{cout << "	 Hashed " << i << " of " << sequenes.size() << "\r";}
		
		}
		string Sequence; 
		//pragma omp critical (sequenes)
		{Sequence = sequenes[i];}
		int LoopLimit = Sequence.size() - SearchHash;

		for (int j = 0; j < LoopLimit; j++) {
			string hash = Sequence.substr(j, SearchHash);
			size_t found = hash.find('N');
			
			if (found == std::string::npos) {
				unsigned long LongHash = Util::HashToLong(hash);
				unsigned long RevHash = Util::HashToLong(Util::RevComp(hash));
				#pragma omp critical(updateHash)
				{
					Hashes[LongHash].push_back(i);
					Hashes[RevHash].push_back(i);
				} //end pragma??
			} //end if
		} 
	}
	Hashesize.clear(); 
	for (auto it = Hashes.begin(); it != Hashes.end(); it++)
	{
		Hashesize[it->first] = it->second.size(); 
	}
	cout << "\nDone Rebulding HashTable\n";
	return 0;
}

int PrepairSearchList(string A, int Ai,	unordered_map<unsigned long, vector<int>>& Hashes,int SearchHash, int ACT, map<int, vector<int>>& array,bool& hitPosLimit, bool& hitIndexLimit, int& NumberPos,	int& NumberIndex , unordered_map<unsigned long, int>& Hashesize) 
{
	int Alength = A.size();
	map<int, int> Positions;
	int added = 0;

	for (int i = 0; i < Alength - SearchHash; i++) {
		string hash = A.substr(i, SearchHash);
		size_t found = hash.find('N');

		if (found == std::string::npos) {
			unsigned long LongHash = Util::HashToLong(hash);
			int max=0; 
			
			#pragma omp atomic 
				max += Hashesize[LongHash];
			for (vector<int>::size_type i = 0; i < max; i++) 
			{
				int holder = 0; 
				#pragma omp atomic 
				holder += Hashes[LongHash][i];
				if (holder > Ai + 1) {

					if (Positions.count(holder) > 0) {
			Positions[holder]++;
			added++;
		} else {
			Positions[holder] = 1;
			added++;
		}
				}

				if (added > 100000) {
					hitPosLimit = true;
					NumberPos = added;
					break;
				}
			}
		}
	}

	if (FullOut) {
		cout << "done Hashing read" << endl;
	}

	NumberPos = added;
	map<int, int>::iterator uspos;
	map<double, int> SortedPositions;
	double trick = 0.0;

	for (uspos = Positions.begin(); uspos != Positions.end(); ++uspos) {
		trick++;
		SortedPositions[(double)uspos->second + (2.0 / trick)] = uspos->first;
	}

	if (FullOut) {
		cout << "found - " << Positions.size() << " possible locations" << endl;
	}

	map<double, int>::iterator pos;
	vector<int> indexes;
	int sanity = 0;

	for (pos = SortedPositions.end(); pos != SortedPositions.begin(); --pos) {

		if (pos->first >= ACT) {
			indexes.push_back(pos->second + 0);
			sanity++;

			if (sanity > 1000) {
				hitIndexLimit = true;
				NumberIndex = sanity;
				break;
			}
		}
	}

	NumberIndex = sanity;
	#pragma omp critical (array)
	{ array[Ai] = indexes; }

	if (FullOut) {
		cout << "			 " << indexes.size() << " locations passed filter" << endl;
	}
	return 1;
}

int Align3(vector<string>& sequenes, string Ap, string Aq, int Ai, int& overlap,int& BestIndex, float minPercent, bool& PerfectMatch, int MinOverlap,vector<int>& indexes, int Threads, int NumReads) 
{
	int QualityOffset = 33;
	bool verbose = false;
	int Alength = Ap.size();
	int bestScore = 0;

	#pragma omp parallel for num_threads(Threads) shared(BestIndex)
	for (int booya = 0; booya < indexes.size(); booya++) 
	{
		string A; 
		int AlengthL; 
		int j; 
		//pragma omp critical (A)
		{
			A = Ap;
			AlengthL = A.size();
			j = indexes[booya];
		}
		string B;
		bool localcheck;
		//pragma omp critical (sequenes)
		{B = sequenes[j];}
		float score = 0;
		int Blength = B.size();
		int k;
		int window = -1;
		int longest = -1;
		bool Asmaller = true;

		if (Blength > AlengthL) 
		{
			window = AlengthL;
			longest = Blength;
			Asmaller = false;
		} else 
		{
			Asmaller = true;
			window = Blength;
			longest = AlengthL;
		}

		int MM = window - (window * minPercent);
		int LbestScore = -1;
		int LBestIndex = -1;
		int Loverlap = 0;
		int Acount = 0;
		int Bcount = 0;

		for (int i = 0; i <= longest - window; i++) 
		{
			score = 0;

			for (k = 0; k < window; k++)	// base compare loop
			{
				if (A.c_str()[k + Acount] == B.c_str()[k + Bcount]) 
				{
					if (A.c_str()[k + Acount] != 'N') 
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

			if (Asmaller) {
				Acount++;
			} else {
				Bcount++;
			}

			if (verbose) {cout << "	Score = " << score << endl;}

			float percent = score / (window);
			if (percent >= minPercent) 
			{
				if (LbestScore < score) 
				{
						LbestScore = score;
						LBestIndex = j;
					
					if (Asmaller) {
						Loverlap = i * -1;
					} else {
						Loverlap = i;
					}

					if (score == window) 
					{
						PerfectMatch = true;
						break;
					}
				}
			}
		}

		if (PerfectMatch == false) 
		{
			for (int i = window - 1; i >= MinOverlap; i--) 
			{
				if (verbose) {cout << "i = " << i << endl;}
				score = 0;
				for (k = 0; k <= i; k++) 
				{
					if (verbose) {cout << "	k = " << k << " so A = " << AlengthL - i + k<< " /\\ B = " << 0 + k << endl;}
					if (verbose) {cout << "		 A >> " << A.c_str()[AlengthL - i + k - 1] << "="<< B.c_str()[0 + k] << " << B" << endl;}

					if (A.c_str()[AlengthL - i + k - 1] == B.c_str()[0 + k]) 
					{
						if (B.c_str()[0 + k] != 'N') 
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
					if (LbestScore < score) 
					{
							LbestScore = score;
							LBestIndex = j;
							Loverlap = i - AlengthL + 1;
						if (score == i) 
						{
							break;
						}
					}
				}
			}

			for (int i = window - 1; i >= MinOverlap; i--) 
			{
				if (verbose) {cout << "i = " << i << endl;}

				score = 0;
				for (k = 0; k <= i; k++) 
				{
					if (B.c_str()[Blength - i + k - 1] == A.c_str()[0 + k]) 
					{
						if (A.c_str()[0 + k] != 'N') 
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

				if (verbose) {cout << "	Score = " << score << endl;}

				float percent = score / (k);

				if (percent > minPercent) 
				{
					if (LbestScore < score) 
					{
							LbestScore = score;
							LBestIndex = j;
							Loverlap = Blength - i - 1;
						if (score == i) {
							break;
						}
					}
				}
			}
		}
		#pragma omp critical (best)
		{
			if (LbestScore > bestScore) {
				bestScore = LbestScore;
				BestIndex = LBestIndex;
				overlap = Loverlap;
			}
		}
	}
	return bestScore;
}

string ColapsContigs(string A, string B, int k, string Aq, string& Bq,string Ad, string& Bd, string As, string& Bs) {
	bool verbose = false;
	if (verbose) {cout << "Combining; \n" << A << endl << B << endl;}

	int Asize = A.size();
	int Bsize = B.size();
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

	if (verbose) {cout << "K = " << k << " so Aofset = " << Aoffset<< " and Boffset = " << Boffset << endl;}

	for (int i = 0; i < Asize + Bsize; i++) 
	{
		char Abase = 'Z';
		char Bbase = 'Z';
		char Aqualc = '!';
		char Bqualc = '!';
		unsigned char Adepc = 0;
		unsigned char Bdepc = 0;

		if (((i - Aoffset) >= 0) && ((i - Aoffset) < A.size())) 
		{
			Abase = A.c_str()[i - Aoffset];
			Aqualc = Aq.c_str()[i - Aoffset];
			Adepc = Ad.c_str()[i - Aoffset];
		} else {
			Abase = 'Z';
			Aqualc = '!';
			Adepc = 0;
		}

		if (i - Boffset >= 0 && i - Boffset < B.size()) {
			Bbase = B.c_str()[i - Boffset];
			Bqualc = Bq.c_str()[i - Boffset];
			Bdepc = Bd.c_str()[i - Boffset];
		} else {
			Bbase = 'Z';
			Bqualc = '!';
			Bdepc = 0;
		}

		if (verbose) {cout << "I = " << i << " Bi = " << i - Boffset << " Ai = " << i - Aoffset<< " thus " << Abase << "-" << Bbase << endl;}

		if (Abase == Bbase && Abase != 'Z') {

			newString += Abase;
			if (Aqualc >= Bqualc) {
				newQual += Aqualc;
			} else {
				newQual += Bqualc;
			}
			if ((int)Adepc + (int)Bdepc < 250) {
				newDepth += (Adepc + Bdepc);
			} else {
				newDepth += (char)250;
			}
		} else if (Abase == 'Z' && Bbase != 'Z') {
			newString += Bbase;
			newQual += Bqualc;
			newDepth += Bdepc;
		} else if (Abase != 'Z' && Bbase == 'Z') {
			newString += Abase;
			newQual += Aqualc;
			newDepth += Adepc;

		} else if (Abase != 'Z' && Bbase != 'Z') {
			 if (Adepc > Bdepc) {
				newString += Abase;
				newQual += Aqualc;
				newDepth += Adepc;
			} else if (Adepc < Bdepc) {
				newString += Bbase;
				newQual += Bqualc;
				newDepth += Bdepc;
			} else if (Aqualc >= Bqualc) {
				newString += Abase;
				newQual += Aqualc;
				newDepth += Adepc;
			} else {
				newString += Bbase;
				newQual += Bqualc;
				newDepth += Bdepc;
			}

		} else if (Abase == 'Z' && Bbase == 'Z') {
			Bq = newQual;
			Bd = newDepth;
			break;
		}	
	}
	Bs += As;
	Bq = newQual;
	Bd = newDepth;
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

	S = NewS;
	qual = NewQ;
	base = false;
	NewS = "";
	NewQ = "";

	for (int i = 0; i < S.size(); i++) {

		if (base) {
			NewS = NewS + S.c_str()[i];
			NewQ = NewQ + qual.c_str()[i];
		} else if (S.c_str()[i] != 'A' && S.c_str()[i] != 'C' &&
							 S.c_str()[i] != 'G' && S.c_str()[i] != 'T') {
		} else {
			base = true;
			NewS = NewS + S.c_str()[i];
			NewQ = NewQ + qual.c_str()[i];
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
 

	string NewS2 = "";
	string NewD2 = "";
	string NewQ2 = "";
	if (NewS.size() > 1) {
		base = false;
		for (int i = 0; i < NewS.size(); i++) {
			if (base) {
				NewS2 = NewS2 + NewS.c_str()[i];
				NewD2 = NewD2 + NewD.c_str()[i];
				NewQ2 = NewQ2 + NewQ.c_str()[i];
			} else if ((int)NewD.c_str()[i] > cutoff) {
				base = true;
				NewS2 = NewS2 + NewS.c_str()[i];
				NewD2 = NewD2 + NewD.c_str()[i];
				NewQ2 = NewQ2 + NewQ.c_str()[i];
			}
		}
	}

	depth = NewD2;
	quals = NewQ2;
	return NewS2;
}


string AdjustBases(string sequence, string qual) {
	int MinQ = 5;
	int QualOffset = 33;
	string NewString = "";

	for (int i = 0; i < sequence.size(); i++) {
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
	str.replace(start_pos, from.size(), to);
	return true;
}

string FlipStrands(string strand) {
	string NewStrand = "";

	for (int i = 0; i < strand.size(); i++) {
		if (strand.c_str()[i] == '+') {
			NewStrand += "-";
		} else if (strand.c_str()[i] == '-') {
			NewStrand += "+";
		}
	}
	return NewStrand;
}
void compresStrand(string S, int& F, int& R) {
	for (int i = 0; i < S.size(); i++) {
		if (S.c_str()[i] == '+')
			F++;
		else
			R++;
	}
	return;
}

int main(int argc, char* argv[]) {
	int SearchHash = 30;
	int ACT = 0;
	float MinPercent;
	int MinOverlap;
	int MinCoverage;
	cout << "you gave " << argc << " Arguments" << endl;

	if (argc != 11) {
		cout << "ERROR, wrong numbe of arguemnts\nCall is: FASTQ, MinPercent, "
						"MinOverlap, MinCoverage, ReportStub, SearchHashSize, ACT, OutFile "
						"LCendTrimEpth Threads"
				 << endl;
		return 0;
	}

	ifstream fastq;
	fastq.open(argv[1]);

	if (fastq.is_open()) {
		cout << "Parent File open - " << argv[1] << endl;
	}	// cout << "##File Opend\n";
	else {
		cout << "Error, ParentHashFile could not be opened";
		return 0;
	}

	//TODO: Factor out temps and cast to string
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
	report.open(FirstPassFile.c_str());

	if (report.is_open()) {
	} else {
		cout << "ERROR, Mut-Output file could not be opened - " << FirstPassFile
				 << endl;
		return 0;
	}

	ofstream Depreport;
	FirstPassFile += "d";
	Depreport.open(FirstPassFile.c_str());
	if (report.is_open()) {
	} else {
		cout << "ERROR, Mut-Output file could not be opened - " << FirstPassFile
				 << endl;
		return 0;
	}

	ofstream good;
	FirstPassFile = ss.str();
	FirstPassFile += "good.fastq";
	good.open(FirstPassFile.c_str());

	if (good.is_open()) {
	} else {
		cout << "ERROR, Mut-Output file could not be opened - " << FirstPassFile
				 << endl;
		return 0;
	}

	ofstream bad;
	FirstPassFile = ss.str();
	FirstPassFile += "bad.fastq";
	bad.open(FirstPassFile.c_str());

	if (bad.is_open()) {
	} else {
		cout << "ERROR, Mut-Output file could not be opened - " << FirstPassFile
				 << endl;
		return 0;
	}

	string line;
	vector<string> sequenes;	// = new vector<string>;
	vector<string> qual;			//= new vector<string>;
	vector<string> depth;		 // = new vector<string>;
	vector<string> strand;
	std::unordered_map<unsigned long, vector<int>> Hashes;
	std::unordered_map<unsigned long, int> Hashesize; 
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
	string Fastqd = argv[1];
	size_t found = Fastqd.find(".fastqd");

	if (found != string::npos) {
		int counter = 0;
		cout << "ATTENTION - Fastq+depth input detected, reading in FASTQD file \n";

		while (getline(fastq, L1)) {
			counter++;
			if (counter % 100 == 1) {
				cout << "Read in " << counter << " lines and rejected " << Rejects
						 << " reads\r";
			}

			getline(fastq, L2);
			getline(fastq, L3);
			getline(fastq, L4);
			getline(fastq, L5);
			getline(fastq, L6);
			string depths = "";
			int ReadSize = L2.size();
			bool Multiple = false;
			vector<string> temp = Util::Split(L6, ' ');

			for (vector<string>::size_type i = 0; i < temp.size(); i++) {
				unsigned char C = atoi(temp[i].c_str());
				depths += C;
				if ((int)C > 1) {
					Multiple = true;
				}
			}

			//Multiple = false;
			if (Multiple == true) {
				L2 = TrimLowCoverageEnds(L2, L4, depths, TrimLCcuttoff);
			}

			if (L2.size() > SearchHash + 1) {
				lines++;
				sequenes.push_back(L2);
				qual.push_back(L4);
				depth.push_back(depths);
				strand.push_back(L5);
				ReadSize = L2.size();
			} else {
				Rejects++;
			}
		}
	} else {
		vector<string> DupCheck;
		cout << "Reading in raw fastq \n";
		int counter = 0;

		while (getline(fastq, L1)) {
			counter++;
			if (counter % 100 == 1) {
				cout << "Read in " << counter << " lines, rejected " << Rejects
						 << " reads with " << dup << " duplicates\r";
			}
			getline(fastq, L2);
			getline(fastq, L3);
			getline(fastq, L4);
			int Ns = 0;

			for (int i = 0; i < L2.size(); i++) {
				if (L2.c_str()[i] == 'N') {
					Ns++;
				}
			}

			lines++;
			int ReadSize = L2.size();
			string depths = "";
			ReadSize = L2.size() - 1;
			bool found = false;
			bool RunDupCheck = true;

			if (RunDupCheck) {

#pragma omp parallel for num_threads(Threads) shared(DupCheck, L2, found)
				for (int i = 0; i < DupCheck.size(); i++) {
					if (L2.size() == DupCheck[i].size()) {
						bool AllBasesMatch = true;

						for (int k = 0; k < L2.size(); k++) {
							if (L2.c_str()[k] == 'N' or DupCheck[i].c_str()[k] == 'N') {
							} else if (L2.c_str()[k] == DupCheck[i].c_str()[k]) {
							} else {
								AllBasesMatch = false;
								break;
							}
						}

						if (AllBasesMatch) {
							#pragma omp critical (found)
							{ found = true; }
						}
					}
				}

				if ((double)Ns / (double)L2.size() < 0.20) {
					DupCheck.push_back(L2);
				}
			}

			if (found == false) {
				L2 = AdjustBases(L2, L4);
				L2 = TrimNends(L2, L4);

				if ((double)L2.size() / (double)ReadSize > .6) {
					goodlines++;
					sequenes.push_back(L2);
					qual.push_back(L4);
					strand.push_back("+");
					unsigned char C = 1;

					for (int i = 0; i <= ReadSize; i++) {
						depths += C;
					}

					depth.push_back(depths);
					good << L1 << endl << L2 << endl << L3 << endl << L4 << endl;
				} else {
					Rejects++;
				}
			} else {
				dup++;
				bad << L1 << endl << L2 << endl << L3 << endl << L4 << endl;
			}
		}
		DupCheck.clear();
	}

	good.close();
	bad.close();
	cout << "done reading " << endl;
	int NumReads = sequenes.size();
	cout << "\nDone reading in \n		 Read in a total of " << lines
			 << " and rejected " << Rejects << " with " << dup
			 << " duplicate reads detected for a total of " << goodlines
			 << "good reads" << endl;

	RebuildHashTable(sequenes, 0, SearchHash, Hashes, Threads, Hashesize);
	clock_t St, Et;
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

	for (std::vector<string>::size_type b = 0; b < sequenes.size(); b += Buffer) {
		LinesSinceLastBuild += Buffer;

		if (LinesSinceLastBuild > 1000000) {
			RebuildHashTable(sequenes, b, SearchHash, Hashes, Threads, Hashesize);
			LinesSinceLastBuild = 0;
		}

		vector<string> ToAddHashes;
		vector<int> ToAddPos;
		map<int, vector<int>> Forwards;
		map<int, vector<int>> Revs;
		int max = b + Buffer;

		if (max > sequenes.size()) {
			max = sequenes.size();
		}

		if (FullOut) {
			cout << "Bulding list to align" << endl;
		}

		
#pragma omp parallel for num_threads(Threads) shared(Hashes, Forwards)
		for (int i = b; i < max; i++) 
		{
			string A = sequenes[i];
			bool posLimit = false;
			bool sanityLimit = false;
			int NumPos = 0;
			int NumSanity = 0;
			PrepairSearchList(A, i, Hashes, SearchHash, ACT, Forwards, posLimit,sanityLimit, NumPos, NumSanity, Hashesize);
			if (posLimit) {
				NumberHitPosLimit++;
			}
			if (sanityLimit) {
				NumberHitSanityLimit++;
			}
			AverageFPos =((AverageFPos * (double)b) + (double)NumPos) / ((double)b + 1.0);
			AverageFSanity = ((AverageFSanity * (double)b) + (double)NumSanity) /((double)b + 1.0);
		}

#pragma omp parallel for num_threads(Threads) shared(Hashes, Revs)
		for (int i = b; i < max; i++) 
		{
			string A = Util::RevComp(sequenes[i]);
			bool posLimit = false;
			bool sanityLimit = false;
			int NumPos = 0;
			int NumSanity = 0;
			PrepairSearchList(A, i, Hashes, SearchHash, ACT, Revs, posLimit,sanityLimit, NumPos, NumSanity, Hashesize);
			if (posLimit) {
				NumberHitPosLimit++;
			}
			if (sanityLimit) {
				NumberHitSanityLimit++;
			};
			AverageRPos =((AverageRPos * (double)b) + (double)NumPos) / ((double)b + 1.0);
			AverageRSanity = ((AverageRSanity * (double)b) + (double)NumSanity) /((double)b + 1.0);
		}

		if (FullOut) {
			cout << "Done Bulding List" << endl;
		}

		for (int i = b; i < max; i++) {
			string A, Aqual, Adep, Astr;
			A = sequenes[i];
			Aqual = qual[i];
			Adep = depth[i];
			Astr = strand[i];
			int k;
			bool PerfectMatch = false;
			int bestIndex = -1;
			float Dt;

			if (FullOut) {
				Et = clock();
				Dt = ((double)(Et - St)) / CLOCKS_PER_SEC;
				cout << "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<**************************************************************>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>" << endl;
			}

			if (i % 100 == 5) {
				gettimeofday(&end, NULL);
				float Dt = end.tv_sec - start.tv_sec;
				cout << "done " << i << " of " << NumReads
						 << ", \% = " << ((double)i / (double)NumReads) * 100.00
						 << ", TT= " << Dt << " , S/R = " << Dt / i
						 << ", \% match= " << ((double)FoundMatch / (double)i) * 100.00
						 << " %Plimit = " << NumberHitPosLimit / 2 << ", "
						 << ((double)NumberHitPosLimit / 2.0) / (double)i
						 << " AvP= " << (AverageFPos + AverageRPos) / 2.0
						 << " % ILmit = " << NumberHitSanityLimit << ", "
						 << ((double)NumberHitSanityLimit / 2.0) / (double)i
						 << " AvI= " << (AverageRSanity + AverageFSanity) / 2.0 << "\r";
			}

			int booya =Align3(sequenes, A, Aqual, i, k, bestIndex, MinPercent, PerfectMatch, MinOverlap, Forwards[i], Threads, NumReads);

			if (FullOut) {
				cout << "best forward score is " << booya << " k is " << k
						 << " index = " << bestIndex << endl;
			}

			if (!(PerfectMatch)) {
				string revA = Util::RevComp(A);
				string revAqual = Util::RevQual(Aqual);
				string revAdep = Util::RevQual(Adep);
				string revAstr = FlipStrands(Astr);
				int revk = -1;
				int revbestIndex = -1;

				if (FullOut) {
					cout << "Checking Reverse\n";
				}

				int revbooya =	Align3(sequenes, revA, revAqual, i, revk, revbestIndex, MinPercent, PerfectMatch, MinOverlap, Revs[i], Threads, NumReads);
				if (FullOut) {
					cout << "best reverse score is " << revbooya << " k is " << revk
							 << " index = " << revbestIndex << endl;
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
					cout << "Perfect Match Found, Skipping Reverse Search" << endl;
				}
			}

			if (booya < MinOverlap) {
				if (FullOut) {
					cout << "No good match found, skipping" << endl;
				}
			} else {
				string B, Bqual, Bdep, Bstr;
				B = sequenes[bestIndex];
				Bqual = qual[bestIndex];
				Bdep = depth[bestIndex];
				Bstr = strand[bestIndex];
				FoundMatch++;

				if (FullOut) {
					if (k > 0) {
						cout << "found match at " << k << endl;

						for (int z = 0; z < k; z++) {
				cout << "+";
			}

						cout << A << endl << B << endl;

						for (int z = 0; z < Bdep.size(); z++) {
							int bam = Bdep.c_str()[z];
							cout << bam;
						}

						cout << endl;

					} else {
						cout << "found match at " << k << endl;
						cout << A << endl;

						for (int z = 0; z < abs(k); z++) {
				cout << "-";
			}
						cout << B << endl;
						for (int z = 0; z < abs(k); z++) {
				cout << "-"; 
			}
						for (int z = 0; z < Bdep.size(); z++) {
							int bam = Bdep.c_str()[z];
							cout << bam;
						}
			
						cout << endl;
					}
				}

				string combined =
		ColapsContigs(A, B, k, Aqual, Bqual, Adep, Bdep, Astr, Bstr);

				if (Bqual.size() != combined.size()) {
					cout << "ERRRORRR "
									"^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^"
									"^^^^"
							 << endl;
				}

				qual[bestIndex] = Bqual;
				depth[bestIndex] = Bdep;
				sequenes[bestIndex] = combined;
				strand[bestIndex] = Bstr;
				sequenes[i] = "moved";

#pragma omp parallel for num_threads(Threads) shared(Hashes)
				for (int j = 0; j < A.size() - SearchHash; j++) {
					string hash = A.substr(j, SearchHash);
					size_t found = hash.find('N');

					if (found == std::string::npos) {
						bool found = false;
						int k = 0;

						for (k = 0; k < Hashes[Util::HashToLong(hash)].size(); k++) {
							if (Hashes[Util::HashToLong(hash)][k] == bestIndex) {
								found = true;
								break;
							}
						}

						if (found = false) {
							#pragma omp critical (Hashes)
							{
								Hashes[Util::HashToLong(hash)].push_back(bestIndex);
							}
						}
						found = false;

						for (k = 0;
								 k < Hashes[Util::HashToLong(Util::RevComp(hash))].size();
								 k++) {
							if (Hashes[Util::HashToLong(Util::RevComp(hash))][k] ==	bestIndex) {
								found = true;
								break;
							}
						}

						if (found = false) {
							#pragma omp critical (Hashes)
							{
								Hashes[Util::HashToLong(Util::RevComp(hash))].push_back(bestIndex);
							}
						}
					}
				}

				if (FullOut) {
					cout << combined << endl;
					cout << Bqual << endl;

					for (int z = 0; z < Bdep.size(); z++) {
						int bam = Bdep.c_str()[z];
						cout << bam;
					}

					cout << endl;
				}
			}
		}
	}

	cout << "\nRESULTS\n";
	int count = 0;

	for (int i = 0; i < sequenes.size(); i++) {

		if (sequenes[i] != "moved" && sequenes[i].size() >= 95) {
			string rDep = depth[i];
			int maxDep = -1;

			for (int z = 0; z < rDep.size(); z++) {
				unsigned char bam = rDep.c_str()[z];
				if ((int)bam > maxDep) {
					maxDep = (int)bam;
				}
			}

			if (maxDep >= MinCoverage) {
				count++;
			 int F = 0;
				int R = 0;
				compresStrand(strand[i], F, R); 
	//report << "@NODE_" << i << "_L=" << sequenes[i].size()			 << "_D=" << maxDep << endl;
	report << "@NODE_" << argv[6] << "_" << i << "_L" << sequenes[i].size()<< "_D" << maxDep << ":" << F << ":" << R << ":" << endl;
				report << sequenes[i] << endl;
				report << "+" << endl;
				report << qual[i] << endl;

				Depreport << "@NODE_" << argv[6] << "_" << i << "_L" << sequenes[i].size()<< "_D" << maxDep << ":" << F << ":" << R << ":" << endl;
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
					Depreport << ' ' << booya;
				}

				Depreport << endl;
			}
		}
	}
	cout << "Wrote " << count << " sequences" << endl;
	report.close();
	Depreport.close();
}
