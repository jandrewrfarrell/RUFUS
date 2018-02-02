/*By ANDREW FARRELL
 * ReplaceQwithDinFASTQ.cpp
 * TODO: describe funciton of file
 */

#include <algorithm>
#include <bitset>
#include <fstream>
#include <iostream>
#include <math.h>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <sys/stat.h>
#include <sys/types.h>
#include <time.h>
#include <vector>

#include "Util.h"

using namespace std;

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

int main(int argc, char* argv[]) {
  int SearchHash = 30;
  int ACT = 0;
  float MinPercent;
  int MinOverlap;
  int MinCoverage;
  // cout << "you gave "<< argc << " Arguments"  << endl;
  if (argc != 2) {
    cout << "ERROR, wrong numbe of arguemnts\nCall is: FASTQD " << endl;
    return 0;
  }

  ifstream fastq;
  fastq.open(argv[1]);

  if (fastq.is_open()) {
  } else {
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

  while (getline(fastq, L1)) {
    getline(fastq, L2);
    getline(fastq, L3);
    getline(fastq, L4);
    getline(fastq, L5);
    getline(fastq, L6);
    string depths = "";
    string adjustedDepths = "";
    int ReadSize = L2.size();
    bool Multiple = false;
    vector<string> temp = Util::Split(L6, ' ');

    for (vector<string>::size_type i = 0; i < temp.size(); i++) {
      unsigned char C = atoi(temp[i].c_str());
      depths += C;
    }

    for (int i = 0; i < depths.size(); i++) {
      int dep = (unsigned char)depths.c_str()[i];
      if ((dep + 33) > 126) {
        adjustedDepths += (unsigned char)126;
      } else {
        adjustedDepths += char((unsigned char)depths.c_str()[i] + 33);
      }
    }

    lines++;
    cout << L1 << endl;
    cout << L2 << endl;
    cout << L3 << endl;
    cout << adjustedDepths << endl;
    cout << L5 << endl;
    cout << L6 << endl;
  }
}
