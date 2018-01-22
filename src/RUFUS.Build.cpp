/* By ANDREW FARRELL
 * RUFUS.1kg.filter.cpp
 * ---------------------------------------------------
 * Builds intersections between k-mer frequency tables
 * to identify k-mers that correspond to a mutation or
 * copy number variation
 * ---------------------------------------------------
 */

#include <bitset>
#include <fstream>
#include <ios>
#include <iostream>
#include <math.h>
#include <omp.h>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <stack>
#include <sys/resource.h>
#include <sys/time.h>
#include <time.h>
#include <unistd.h>
#include <vector>

#include "Util.h"

using namespace std;

ifstream ParentHash[100];

int main(int argc, char *argv[]) {
  double z = 0;
  double one = 1;
  cout << "0/1 = " << z / 1 << endl;
  double MaxCount = 65534;
  double vm, rss, MAXvm, MAXrss;
  MAXvm = 0;
  MAXrss = 0;
  Util::process_mem_usage(vm, rss, MAXvm, MAXrss);
  string helptext;
  helptext =
      \
"\
RUFUS.interpret: converts RUFUS aligned contigs into a VCF \n\
By Andrew Farrell\n\
   The Marth Lab\n\
\n\
options:\
  -h [ --help ]          Print help message\n\
  -c  arg                Control or parent file pair, Hash file then dist file separeted by space\n\
  -s  arg                mutant or subject file paie, Hash file then dist file separeted by space \n\
  -o  arg                Output stub\n\
  -hs arg Hash size used in kmer counting\n\
  -mS arg minimum coverage in the subject\n\
  -mC arg  maximum coverage in the control samples to consider 0\n\
  -max arg Maximum coverage to consider\n\
  -t  arg  Number of threds to use\n\
  -d  arg                Deliminator for parent file, default = tab\n\
";

  vector<string> ParentHashFilePaths;
  string MutHashFilePath = "";
  string OutFileStub = "";
  int HashSize = -1;
  double MinCovSubject = -1.0;  // probability of a mutation;
  double MaxCovControl = -1.0;  // probablitliy of a copy nubmer event;
  double Pvalue = -1.0;  // Pvlaue threshold
  int threads = 1;
  long MaxCoverage = 100000000;
  char delim = '\t';

  for (int i = 1; i < argc; i++) {
    string p = argv[i];
    cout << i << " = " << argv[i] << endl;
    if (p == "-h") {
      // print help
      cout << helptext << endl;
      return 0;
    } else if (p == "-c") {
      cout << "Par Hash = " << argv[i + 1] << endl;
      ParentHashFilePaths.push_back(argv[i + 1]);
      i = i + 1;
    } else if (p == "-s") {
      cout << "Sub Hash = " << argv[i + 1] << endl;
      cout << "Sub Dist = " << argv[i + 2] << endl;
      MutHashFilePath = argv[i + 1];
      i += 1;
    } else if (p == "-o") {
      OutFileStub = argv[i + 1];
      i++;
    } else if (p == "-hs") {
      HashSize = atoi(argv[i + 1]);
      cout << "HS = " << HashSize << endl;
      i++;
    } else if (p == "-mS") {
      MinCovSubject = atof(argv[i + 1]);
      cout << "MinCovSubject = " << MinCovSubject << endl;
      i++;
    } else if (p == "-mC") {
      MaxCovControl = atof(argv[i + 1]);
      cout << "MaxCovControl = " << MaxCovControl << endl;
      i++;
    } else if (p == "-max") {
      MaxCoverage = atoi(argv[i + 1]);
      i++;
    } else if (p == "-t") {
      threads = atoi(argv[i + 1]);
      i++;
    } else if (p == "-d") {
      delim = argv[i + 1][0];
      i++;
    } else {
      cout << "ERROR: unkown command line paramater -" << argv[i] << "-"
           << endl;
      return 0;
    }
  }

  if (ParentHashFilePaths.size() == 0) {
    cout << "Error in control file inputs, either empty or number of file not "
            "the same for both"
         << endl;
    return 0;
  }
  if (MutHashFilePath == "") {
    cout << "Error subject file required, curnt = " << MutHashFilePath << endl;
    return 0;
  }
  if (OutFileStub == "") {
    cout << "not out file provided, using Subject Hash path as prefix" << endl;
    OutFileStub = MutHashFilePath;
  }
  if (HashSize == -1) {
    cout << "Error Hash Size required" << endl;
    return 0;
  }
  if (MinCovSubject < 0) {
    cout << "Error Minimum coverage in the subject required and must be "
            "greather than 0"
         << endl;
    return 0;
  }
  if (MaxCovControl < -1) {
    cout << "Error Maximum Coverage in the control must be set and must be "
            "greaterthan or equal to 0"
         << endl;
    return 0;
  }

  string LastLine = "";

  for (int c = 0; c < HashSize; c++) {
    LastLine += "T";
  }

  LastLine += "0";

  for (int pi = 0; pi < ParentHashFilePaths.size(); pi++) {
    ParentHash[pi].open(ParentHashFilePaths[pi].c_str());
    if (ParentHash[pi].is_open()) {
      cout << "Parent Hash File open - " << ParentHashFilePaths[pi] << endl;
    }  // cout << "##File Opend\n";
    else {
      cout << "Error, ParentHashFile could not be opened" << endl
           << ParentHashFilePaths[pi] << endl;
      return 0;
    }
  }

  ifstream MutHashFile;
  MutHashFile.open(MutHashFilePath.c_str());

  if (MutHashFile.is_open()) {
  } else {
    cout << "Error, MutHashFile could not be opened";
    return 0;
  }

  ofstream MutMutHashTable;
  string MutHashTablePath = OutFileStub;
  MutMutHashTable.open(MutHashTablePath.c_str());

  if (MutMutHashTable.is_open()) {
  } else {
    cout << "Error MutHashTableFilter file couldnt be opened - "
         << MutHashTablePath << endl;
    return 0;
  }

  ofstream Trace;
  string TracePath = OutFileStub;
  TracePath += ".Buld.Log";
  Trace.open(TracePath.c_str());

  if (Trace.is_open()) {
  } else {
    cout << "Error MutHashTableFilter file couldnt be opened - " << TracePath
         << endl;
    return 0;
  }

  string line;
  long MutCount = 0;
  long lines = 0;
  string ML1, PL1;
  clock_t St, Et;
  struct timeval start, end;
  gettimeofday(&start, NULL);
  bool notdone = true;
  int r = 0;
  cout << "starting " << endl;
  string MutLine;
  long int HetCount = 0;
  vector<string> ParLines;

  for (int pi = 0; pi < ParentHashFilePaths.size(); pi++) {
    string line;
    getline(ParentHash[pi], line);
    ParLines.push_back(line);
  }

  while (getline(MutHashFile, MutLine)) {
    lines++;
    vector<string> Mline = Util::Split(MutLine, '\t');
    int Acount = atoi(Mline[1].c_str());

    if (Acount >= MinCovSubject and Acount <= MaxCoverage) {

      for (int pi = 0; pi < ParentHashFilePaths.size(); pi++) {
        vector<string> Pline = Util::Split(ParLines[pi], delim);

        while (Pline[0] < Mline[0]) {
          string last = ParLines[pi];
          getline(ParentHash[pi], ParLines[pi]);

          if (ParLines[pi].size() == 0) {
            cout << "Parent " << ParentHashFilePaths[pi] << " died at " << last
                 << endl;
            ParLines[pi] = LastLine;
          }

          Pline = Util::Split(ParLines[pi], delim);

        }  // will need to add here check to see if at eof
      }

      stringstream Bcounts;
      int TotalParentDepth = 0;

      for (int pi = 0; pi < ParentHashFilePaths.size(); pi++) {
        vector<string> Pline = Util::Split(ParLines[pi], delim);

        if (Pline[0] == Mline[0]) {
          int Bcount = atoi(Pline[1].c_str());
          TotalParentDepth += Bcount;
          Bcounts << Bcount << "-";
          stringstream ss;
          ss << Pline[0] << '\t' << Bcount << '\t' << Acount << '\t';
        }
      }

      if (Acount >= MinCovSubject and TotalParentDepth <= MaxCovControl) {
#pragma omp critical(h)
        {
          MutMutHashTable << Util::HashToLong(Mline[0].c_str()) << "\t"
                          << TotalParentDepth << "\t" << Acount << "\t"
                          << Mline[0] << endl;
          HetCount++;
        }
      }
    }

    if ((lines) % 10000 == 0) {
      gettimeofday(&end, NULL);
      double Dt = end.tv_sec - start.tv_sec;
      Util::process_mem_usage(vm, rss, MAXvm, MAXrss);
      cout << "Read in " << lines << " lines "
           << "; VM: " << vm << "; RSS: " << rss << " muts = " << HetCount
           << " , lines per second = " << 1.0 / (Dt / (lines)) << "\r";
    }
  }

  cout << "\nDone reading Parent File\n";
  int who = RUSAGE_SELF;
  struct rusage usage;
  int b = getrusage(RUSAGE_SELF, &usage);
  cout << "I am using " << usage.ru_maxrss << endl;
  MutMutHashTable.close();
  Trace.close();
  cout << "\nreally done\n";
}

