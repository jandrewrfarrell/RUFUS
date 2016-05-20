
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

	
	string line;
	map  <string, unsigned long> Chrs;


	string L1;
	while (getline(SamIn, L1))
	{
		vector <string> temp;
		temp = Split(L1, '\t');
		Chrs[temp[2]]++;
		cout << temp[0] << endl << temp[9] << endl << "+" << endl << temp[10] << endl;
	}
	SamIn.close();

	typedef std::map<std::string, unsigned long>::iterator it_type;
	for(it_type iterator = Chrs.begin(); iterator != Chrs.end(); iterator++) 
	{
	    ChrOut << iterator->first << "\t" << iterator->second << endl;
	}
	ChrOut.close(); 	
	
}
	
