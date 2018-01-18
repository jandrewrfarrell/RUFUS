/*By ANDREW FARRELL
 * RUFUS.1kg.filter.cpp
 * TODO: describe funciton of file
 */

#include <bitset>
#include <errno.h>
#include <fcntl.h>
#include <fstream>
#include <ios>
#include <iostream>
#include <unistd.h>
#include <math.h>
#include <omp.h>
#include <sstream>
#include <stack>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <sys/mman.h>
#include <sys/resource.h>
#include <sys/stat.h>
#include <sys/time.h>
#include <sys/types.h>
#include <time.h>
#include <unordered_map>
#include <vector>

#include "Util.h"

using namespace std;

int main(int argc, char *argv[]) {
  cout << "Call is Parent_hash Mutant_Hash FilterLiestOutFile hashsize "
          "MaxParDepth MinMutDepth MaxDepth Threads"
       << endl;

  double vm, rss, MAXvm, MAXrss;
  MAXvm = 0;
  MAXrss = 0;
  Util::process_mem_usage(vm, rss, MAXvm, MAXrss);
  string temp = argv[4];
  int HashSize = atoi(temp.c_str());
  temp = argv[5];
  int MaxParDepth = atoi(temp.c_str());
  temp = argv[6];
  int MinMutDepth = atoi(temp.c_str());
  temp = argv[7];
  int TooHigh = atoi(temp.c_str());
  temp = argv[8];
  int Threads = atoi(temp.c_str());
  int BufferSize = 10000;

  cout << "Paramaters are:\n  Parent_Hash = " << argv[1]
       << "\n  Mutant_Hash = " << argv[2]
       << "\n  FilterLiestOutFile = " << argv[3] << "\n  HashSize = " << argv[4]
       << "\n  Max Parent Depth = " << argv[5]
       << "\n  Min Mut Depth = " << argv[6] << "\n  Maximum count = " << argv[7]
       << "\n  Threads = " << argv[8] << endl;

  ifstream MutHashFile;
  MutHashFile.open(argv[2]);

  if (MutHashFile.is_open()) {
  } else {
    cout << "Error, MutHashFile could not be opened";
    return 0;
  }

  ofstream HashTable;
  string HashTablePath = argv[3];
  HashTable.open(HashTablePath.c_str());

  if (HashTable.is_open()) {
  } else {
    cout << "Error HashTableFilter file couldnt be opened - " << HashTablePath
         << endl;
    return 0;
  }

  string line;
  unordered_map<unsigned long, int> good;
  unsigned long int lines = 0;
  string L1;
  string L2;
  string L3;
  string L4;
  string ML1, PL1;
  clock_t St, Et;
  struct timeval start, end;
  gettimeofday(&start, NULL);
  bool notdone = true;
  int r = 0;
  cout << "Reading in Mutant Hash Table" << endl;

  while (notdone == true) {
    lines++;
    vector<string> Mstack;

    while (Mstack.size() < BufferSize && getline(MutHashFile, ML1)) {
      Mstack.push_back(ML1);
    }

    for (int bam = 0; bam < Mstack.size(); bam++) {
      string B1 = Mstack[bam];
      vector<string> stuff;
      stuff = Util::Split(B1, '\t');
      int Bcount = atoi(stuff[1].c_str());
      unsigned long MutantLongHash = Util::HashToLong(stuff[0]);
      good.insert(pair<unsigned long, int>(MutantLongHash, Bcount)); 
    }

    if (Mstack.size() == 0) {
      notdone = false;
    }

    if ((lines * BufferSize) % 1000 == 0) {
      gettimeofday(&end, NULL);
      float Dt = end.tv_sec - start.tv_sec;
      Util::process_mem_usage(vm, rss, MAXvm, MAXrss);
      cout << "Read in " << BufferSize * lines << " lines "
           << "; VM: " << vm << "; RSS: " << rss << " good = " << good.size()
           << ", TotalTime= " << Dt
           << " , second per line = " << Dt / (lines * BufferSize) << "\r";
    }
  }
  cout << "Done reading in mut table\nStarting filter\n";
  lines = 0;

#pragma omp parallel for num_threads(Threads)
  for (int p = 0; p < Threads; p++) {
    struct timeval start, end;
    gettimeofday(&start, NULL);
    long int fd;
    struct stat sb;
    char *fName = argv[1];
    char *data;
    char *lineMap;
    fd = open(fName, O_RDONLY);
    fstat(fd, &sb);
    long int pageSize;
    pageSize = sysconf(_SC_PAGE_SIZE);
    long int off = 0;
    long int count = 0;
    char *fileptr = NULL;
    string foo = "";
    unordered_map<char, bool> chars;
    chars.insert(pair<char, bool>('A', true));
    chars.insert(pair<char, bool>('C', true));
    chars.insert(pair<char, bool>('G', true));
    chars.insert(pair<char, bool>('T', true));
    long int ChunkSize = sb.st_size / pageSize / Threads;
    off = ChunkSize * p * pageSize;

    cout << "total File size = " << sb.st_size << endl;
    cout << "Chunc Size = " << ChunkSize << endl;
    cout << "thread " << p + 1 << endl;
    cout << "I go from " << off << "    to    "
         << sb.st_size / pageSize / Threads * (p + 1) << endl;

    long pages = 0;

    while (off <= ChunkSize * pageSize * (p + 1)) {
      pages++;
      fileptr = (char *)mmap64(NULL, pageSize, PROT_READ,
                               MAP_PRIVATE | MAP_POPULATE, fd, off);
      data = fileptr;
      for (int i = 0; i < pageSize; i++) {
        if (chars.count(data[i]) > 0) {
          foo += data[i];
        } else {
          if (foo.length() == HashSize) {
            unsigned long ParHash = Util::HashToLong(foo);
            int Count;
            Count = good.count(ParHash);
            if (Count > 0) {
              gettimeofday(&end, NULL);
              float Dt = end.tv_sec - start.tv_sec;

#pragma omp critical(boom)
              {
                cout << "   FoundHash -" << foo << "- in thread " << p << " "
                     << foo << " - " << ParHash << " - % done ="
                     << ((float)pages / (float)ChunkSize) * 100.0
                     << ", Pages/sec = " << ((float)pages * (float)Threads) / Dt
                     << endl;
              }
	      good[ParHash] = -1;
            }
          }
          foo = "";
        }
      }
      
      munmap(fileptr, pageSize);
      off += pageSize;
    }
  }

  cout << "\nDone reading Parent File\n";
  int who = RUSAGE_SELF;
  struct rusage usage;
  int b = getrusage(RUSAGE_SELF, &usage);
  cout << "I am using " << usage.ru_maxrss << endl;
  unordered_map<unsigned long, int>::iterator pos;

  for (pos = good.begin(); pos != good.end(); ++pos) {
    if (pos->second > MinMutDepth and pos->second < TooHigh) {
      HashTable << pos->first << "\t"
                << "0"
                << "\t" << pos->second << "\t"
                << Util::LongToHash(pos->first, HashSize) << endl;
    }
  }

  HashTable.close();
  cout << "\nreally done\n";
}
