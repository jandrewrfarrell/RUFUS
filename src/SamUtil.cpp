#include <cmath>
#include <ctime>
#include <errno.h>
#include <fcntl.h>
#include <fstream>
#include <ios>
#include <iostream>
#include <map>
#include <math.h>
#include <sstream>
#include <stack>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <sys/mman.h>
#include <sys/resource.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <time.h>
#include <unordered_map>
#include <unistd.h>
#include <vector>

#include "SamUtil.h"

using namespace std;

int SamUtil::GetReadOrientation(int flag) {
  bool b[16];

  for (int j = 0;  j < 16;  ++j){
    b [j] =  0 != (flag & (1 << j));
  }

  cout << "bits are " << b[0] << " " << b[1] << " " << b[2] << " " << b[3]
       << " " << b[4] << " " << b[5] << " " << b[6] << " " << b[7] << " "
       << b[8] << " " << b[9] << " " << b[10] << " " << b[11] << " " << b[12]
       << endl;
  cout << "flag is " << flag << endl;
  cout << "orientation is " << b[4] << endl;
  return b[4];
}

int SamUTil::getHash(string seq, int , int HashSize) {
  int bases = 0;
  string NewHash = "";
  if (j <seq.size()) {
    if (seq.c_str()[j] == 'A' or seq.c_str()[j] == 'C' or
        seq.c_str()[j] == 'G' or seq.c_str()[j] == 'T') {
      while (bases<HashSize and j<seq.size()) {
        if (seq.c_str()[j] == 'A' or seq.c_str()[j] == 'C' or
            seq.c_str()[j] == 'G' or seq.c_str()[j] == 'T') {
          NewHash+=seq.c_str()[j];
          bases++;
        }
        j++;
      }
    }
  }
  return NewHash;
}

string SamUtil::compressVar(string line, int start, string& StructCall)
{
  cout << "compressing var" << endl;
  char current = line.c_str()[0];
  int currentCount = 1;
  string CV = "";
  for (int i = 1; i< line.size(); i++)
    {
      cout << current << endl;
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
                  cout << "YAAAY STRUCT" << endl;
                  StructCall = "SVTYPE=DUP;END=";
                  StructCall += convertEND.str();
                  StructCall += ";SVLEN=";
                  StructCall += convert.str();
                  StructCall += ";";
                  cout << StructCall << endl;
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
              cout << "ERROR in compress " << current << " " << currentCount << endl;
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
          cout << "YAAAY STRUCT" << endl;
          StructCall = "SVTYPE=DUP:TANDEM;END=";
          StructCall += convertEND.str();
          StructCall += ";SVLEN=";
          StructCall += convert.str();
          StructCall += ";";
          cout << StructCall << endl;
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
      cout << "ERROR in compress " << current << " " << currentCount << endl;
    }
  return CV;
}
