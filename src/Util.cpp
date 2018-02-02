#include <unistd.h>
#include <ios>
#include <bitset>
#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <math.h>
#include <sys/resource.h>
#include <stdio.h>
#include <stdlib.h>


#include "Util.h"

using namespace std;

bool Util::fncomp(char lhs, char rhs) 
{
  return lhs<rhs;
}

 const vector<string> Util::Split(const string& line, const char delim)
{
  vector<string> tokens;
  stringstream lineStream(line);
  string token;
  while(getline(lineStream, token, delim)){
    tokens.push_back(token);
  }
  return tokens;
}

string Util::trim(string s)
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

unsigned long Util::HashToLong(string hash)
{
  bitset<64> HashBits;
  //#pragma omp parellel for                                                                                                                                                                    
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
  //cout <<  HashBits.to_ulong() << "-" << endl;                                                                                                                                                
  return HashBits.to_ulong();
}

unsigned long HashToLongTB(string hash, string calledby)
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

string Util::LongToHash(unsigned long LongHash, int HashSize)
{
  string value = "";
  //cout << LongHash << endl;                                                                                                                                                                   
  bitset<64> test (LongHash);
  for (int i = 1; i < HashSize*2; i+=2)
    {
      if (test[i-1] == 0)
	{
	  if (test[i] == 0)
	    {
	      value += "A";
	    }
	  else
	    {
	      value += "C";
	    }
	}
      else
	{
	  if (test[i] == 0)
	    {
	      value += "G";
	    }
	  else
	    {
	      value += "T";
	    }
	}
    }
  return value;
}

/*string Util::RevComp(string Sequence)
{
  string NewString = "";
  for(int i = Sequence.length(); i > 0; --i)
    {
      switch(char(Sequence.c_str()[i]))
	{
	case 'A':
	  NewString += 'T';
	case 'C':
	  NewString += 'G';
	case 'G':
	  NewString += 'C';
	case 'T':
	  NewString += 'A';
	case 'N':
	  NewString += 'N';
	default: 
	  cout << "ERROR IN RevComp - \n" << Sequence.c_str()[i] << endl;
	}
    }
  return NewString;
  }*/

  string Util::RevComp (string Sequence)
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

string Util::RevQual(string Sequence)
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

void Util::process_mem_usage(double& vm_usage, double& resident_set, double& MAXvm, double& MAXrss)
{
  using std::ios_base;
  using std::ifstream;
  using std::string;

  vm_usage     = 0.0;
  resident_set = 0.0;

  // 'file' stat seems to give the most reliable results                                                                                                               

  ifstream stat_stream("/proc/self/stat",ios_base::in);

  // dummy vars for leading entries in stat that we don't care about                                                                                                   
  //                                                                                                                                                                   
  string pid, comm, state, ppid, pgrp, session, tty_nr;
  string tpgid, flags, minflt, cminflt, majflt, cmajflt;
  string utime, stime, cutime, cstime, priority, nice;
  string O, itrealvalue, starttime;

  // the two fields we want                                                                                                                                            

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

