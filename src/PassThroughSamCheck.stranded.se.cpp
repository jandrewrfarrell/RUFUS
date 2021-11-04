

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
      cout << "ERROR, Output file could not be opened -" << argv[1] << endl;
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
      //field 1
      while (L1_array[i] != '\t')
	{
	  //name+=L1_array[i]; 
	  ++i;
	}
      int lenName = i-start; 
      ++i;
      //string flag = "";
      //field 2
      int flagStart = i; 
      while (L1_array[i] != '\t')
	{
	  ++i;
	}
      int flagEnd = i;  
      ++i;
      //string chr = "";
      start = i; 
      const char* tmpChr = L1_array + i ;
      //field 3
      while (L1_array[i] != '\t')
	{
	  //chr+=L1_array[i];
	  ++i;
	}
      int lenChr = i-start; 
      ++i;
      //string pos = "";
      //field 4
      while (L1_array[i] != '\t')
	{
	  ++i;
	}
      ++i;
      //string something = "";
      //field 5
      while (L1_array[i] != '\t')
	{
	  ++i;
	}
      ++i;
      //string cigar = "";
      //field 6
      while (L1_array[i] != '\t')
	{
	  ++i;
	}
      ++i;
      //string something2 = "";
      //field 7
      while (L1_array[i] != '\t')
	{
	  ++i;
	}
      ++i;
      //string something3=  ""; 
      //field 8
      while (L1_array[i] != '\t')
	{
	  ++i;
	}
      ++i;
      //string something4=  "";
      //field 9
      while (L1_array[i] != '\t')
	{
	  ++i;
	}
      ++i;
      //string seq = "";
      start = i; 
      int seqStart = i; 
      const char* tmpSeq=L1_array + i ;
      //field 10
      while (L1_array[i] != '\t')
	{
	  //seq += L1_array[i];
	  ++i;
	}
      int seqEnd = i; 
      int lenSeq = i-start; 
      ++i;
      //string qual = "";
      start = i;
      int qualStart = i; 
      const char* tmpQual = L1_array + i ;
      //field 11
      while (L1_array[i] != '\t')
	{
	  //       qual += L1_array[i];
	  ++i;
	}
      int qualEnd = i; 
      int lenQual = i - start; 
      

      
      if (strncmp( tmpChr,  current.c_str(), lenChr) != 0 or lenChr != current.size() )
	{
	  ChrOut << current  << endl;
	  
	  current = string(tmpChr, lenChr); 
	}
      //cout << "@" << name << endl << seq << endl << "+" << endl << qual << endl;
      
      string flagstring = "";

      for (int i = flagStart; i< flagEnd; i++)
	{
	  flagstring+=L1_array[i];
	}
      
      //int v = atoi(flagstring.c_str());  // flag to dissect
      //int strand = 0 != (v & (1 << 4));
      if ( 0 != (atoi(flagstring.c_str()) & (1 << 4))){
	cout.write("@", 1);
	cout.write(tmpName, lenName);
	cout << endl;
	for (int j =seqEnd-1; j>=seqStart; j--){
	  switch (L1_array[j]){
	  case 'A' : cout << 'T'; break;
	  case 'C' : cout << 'G'; break;
	  case 'G' : cout << 'C'; break; 
	  case 'T' : cout << 'A'; break; 
	  case 'N' : cout << 'N'; break;
	  }
	}
	cout << endl; 
	cout << "+" << endl;
	for (int j =qualEnd-1; j>=qualStart; j--){
	  cout << L1_array[j]; 
	}
	cout << endl; 
      }
      else
	{
	  cout.write("@", 1);
	  cout.write(tmpName, lenName);
	  cout << endl; 
	  cout.write(tmpSeq, lenSeq);
	  cout << endl << "+" << endl; 
	  cout.write(tmpQual, lenQual); 
	  cout << endl;  
	}
    }
  ChrOut << current << endl;
  SamIn.close();
  return 0; 
}
