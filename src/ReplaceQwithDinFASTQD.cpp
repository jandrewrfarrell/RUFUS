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
	int MinQ = 10;
	int QualOffset = 32;
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
	int SearchHash = 30;
	int ACT = 0;
	float MinPercent ;
	int MinOverlap ;
	int MinCoverage ;			
	//cout << "you gave "<< argc << " Arguments"  << endl;
	if (argc != 2)
	{cout << "ERROR, wrong numbe of arguemnts\nCall is: FASTQD " << endl; return 0;}
		
	
	ifstream fastq;
	fastq.open (argv[1]);
	if ( fastq.is_open())
	{}	//cout << "##File Opend\n";
	else
	{
		cout << "Error, ParentHashFile could not be opened";
		return 0;
	}
	
	int lines = -1;

	string L1;
	string L2;
	string L3;
	string L4;
	string L5;
	string L6;
	int Rejects = 0;
	
	//read in fastq
	//int readcount = 0;
		int counter =0;
		//cout << "ATTENTION - Fastq+depth input detected, reading in FASTQD file \n";
		while (getline(fastq, L1))
                {
			counter++;
			if (counter % 100 == 1)
                        {
                        //        cout << "Read in " << counter << " lines and rejected " << Rejects << " reads\r";
                        }

                        getline(fastq, L2);
                        getline(fastq, L3);
                        getline(fastq, L4);
			getline(fastq, L5);
			getline(fastq, L6);
			string depths = "";
			string adjustedDepths = "";
			int ReadSize = L2.size();
				
			bool Multiple = false;
                       	vector<string> temp = Split(L6, ' ');
			 		
			for (vector<string>::size_type i = 0; i < temp.size(); i++)
                        {unsigned char C = atoi(temp[i].c_str());depths += C;}
			
			for (int i =0; i< depths.size(); i++)
			{
				int dep = (unsigned char )depths.c_str()[i]; 
				if ((dep + 33) > 126)
				{ adjustedDepths += (unsigned char) 126; }
				else
				{adjustedDepths += char ((unsigned char )depths.c_str()[i] + 33 );}
			}
			
			
			//if (Multiple == true){L2 = TrimLowCoverageEnds(L2, L4, depths,TrimLCcuttoff);}
			
			//if (L2.size() > SearchHash+1)
			//{
				lines++;	
			cout << L1 << endl;
                        cout << L2 << endl;
			cout << L3 << endl;
			cout << adjustedDepths << endl;
			cout << L5 << endl;
			cout << L6 << endl;
			//}
                }

}
	
