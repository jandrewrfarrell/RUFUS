
/*This vertion imcorporates both copy number and mutation detection in 1
 *  * it needs to be run in two sptes, first the build, then the filter
 *   * it is split up to allow distribution to a cluster */

#include <unistd.h>
#include <ios>

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <string.h>
#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <map>

using namespace std;


const vector<string> Split(const string& line, const char delim) {
    vector<string> tokens;
    stringstream lineStream(line);
    string token;
    while ( getline(lineStream, token, delim) )
        tokens.push_back(token);
    return tokens;
}

int main (int argc, char *argv[])
{
        ifstream SamIn;         
	SamIn.open ("/dev/stdin");

	ofstream ChrOut;
	ChrOut.open (argv[1]);
	if (ChrOut.is_open())
	{}
	else
	{	
		cout << "ERROR, Output file could not be opened -" << argv[3] << endl;
		return 0;
	}

	
	string L1;
	string current = "notachr";
	//string chr = "";
	

	while (getline(SamIn, L1))
	{
		//string name = "";
		int i = 0; 
		
		const char* L1_array = L1.c_str(); 		
		int start =i;
		const char* tmpName = L1_array + i ; 	
		while (L1_array[i] != '\t')
		{
			//name+=L1_array[i]; 
			++i;
		}
		int lenName = i-start; 
		++i;
		//string flag = "";
		while (L1_array[i] != '\t')
                {
                        ++i;
                }
		++i;
		//string chr = "";
		start = i; 
		const char* tmpChr = L1_array + i ;
		while (L1_array[i] != '\t')
                {
                	//chr+=L1_array[i];
		        ++i;
                }
		int lenChr = i-start; 
		++i;
		//string pos = "";
		while (L1_array[i] != '\t')
                {
                        ++i;
                }
                ++i;
		//string something = "";
                while (L1_array[i] != '\t')
                {
                        ++i;
                }
		++i;
		//string cigar = "";
                while (L1_array[i] != '\t')
                {
                        ++i;
                }
		++i;
		//string something2 = "";
                while (L1_array[i] != '\t')
                {
                        ++i;
                }
		++i;
		//string something3=  ""; 
                while (L1_array[i] != '\t')
                {
                        ++i;
                }
		++i;
                //string something4=  "";
                while (L1_array[i] != '\t')
                {
                        ++i;
                }
		++i;
		//string seq = "";
		start = i; 
		const char* tmpSeq=L1_array + i ;
                while (L1_array[i] != '\t')
                {
		//	seq += L1_array[i];
                        ++i;
                }
		int lenSeq = i-start; 
		++i;
                //string qual = "";
		start = i;
		const char* tmpQual = L1_array + i ;
                while (L1_array[i] != '\t')
                {
                 //       qual += L1_array[i];
                        ++i;
                }
		int lenQual = i - start; 
		

		
		if (strncmp( tmpChr,  current.c_str(), lenChr) != 0 or lenChr != current.size() )
		{
			ChrOut << current  << endl;
			
			current = string(tmpChr, lenChr); 
		}
		//cout << "@" << name << endl << seq << endl << "+" << endl << qual << endl;
		cout.write("@", 1);
		cout.write(tmpName, lenName);
		cout << endl; 
		cout.write(tmpSeq, lenSeq);
		cout << endl << "+" << endl; 
		cout.write(tmpQual, lenQual); 
		cout << endl;  
	}
       ChrOut << current << endl;
	ChrOut << "booya" << endl;	
	SamIn.close();
return 0; 
}
	
