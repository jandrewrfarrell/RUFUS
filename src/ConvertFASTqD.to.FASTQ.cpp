/*By ANDREW FARRELL
 * ConvertFASTqD.to.FASTQ.cpp
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

int main(int argc, char *argv[]) {
  int SearchHash = 30;
  int ACT = 0;
  float MinPercent;
  int MinOverlap;
  int MinCoverage;

  if (argc != 2) {
    cout << "ERROR, wrong numbe of arguemnts\nCall is: FASTQD " << endl;
    return 0;
  }

  ifstream fastq;
  fastq.open(argv[1]);
  if (fastq.is_open()) {
  }  // cout << "##File Opend\n";
  else {
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
    cout << L1 << endl;
    cout << L2 << endl;
    cout << L3 << endl;
    cout << L4 << endl;
  }
}	
