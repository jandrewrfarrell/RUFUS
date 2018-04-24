/*By ANDREW FARRELL
 * RUFUS.interpret.cpp
 * TODO: describe funciton of file
 */

#include <bitset>
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

#include "externals/fastahack/Fasta.h"
#include "SamRead.h"
#include "Util.h"

#define NUMINTS  (1000)
#define FILESIZE (NUMINTS * sizeof(int))

using namespace std;

vector <unordered_map <unsigned long int, int> > ParentHashes;
unordered_map  <unsigned long int, int> MutantHashes; 
FastaReference Reff;
int HashSize = 25; 
int totalDeleted; 
int totalAdded;
int MaxVarentSize = 1000; 
ofstream VCFOutFile;
ofstream BEDOutFile;
ofstream BEDBigStuff;
ofstream BEDNotHandled;
ofstream Invertions; 
ofstream Translocations; 
ofstream Translocationsbed; 
ofstream Unaligned; 
map <string, int> Hash; 

int checkPage( char *data, string hash, long int pageSize, string line) {
  bool firstNew = false;
  for (int i = 0; i <  pageSize; i++) {
    if (data[i] == '\n') {
      if (firstNew != true) {
	firstNew = true;
      } else {
	vector<string> stuff;
	stuff = Util::Split(line, '\t');
	string PageHash = stuff[0];
	if (hash == PageHash) {
	  return atoi(stuff[1].c_str());
	}
	line = "";
      }
    }
    else if (firstNew == true) {
      line += data[i];
    }
  }
  return 0;
}

void ProcessPage(char *data, string& PageFirstHash, string& PageLastHash, 
		 long int pageSize) {
  string line = "";
  bool firstNew = false;
  
  for (int i = 0; i < pageSize; i++) {
    if (data[i] == '\n') {
      if (firstNew == true) {
	break;
      }
      else {
	firstNew = true;
      }
    }
    else if (firstNew == true) {
      line += data[i];
    }
  }
  
  vector<string> stuff;
  stuff = Util::Split(line, '\t');
  PageFirstHash = stuff[0];
  firstNew = false;
  line = "";
  for (int i = pageSize-1; i > 0; i+=-1) {
    if ( data[i] == '\n') {
      if (firstNew == true) {
	break;
      } else {
	firstNew = true;
      }
    }
    else if(firstNew == true) {
      line = data[i] + line;
    }
  }
  stuff = Util::Split(line, '\t');
  PageLastHash = stuff[0];
}

int search(long int& fd, string hash, char* fileptr) {
  char *data;
  struct stat sb;
  fstat(fd, &sb);
  long int pageSize;
  pageSize = sysconf(_SC_PAGE_SIZE);
  long int NumPages = sb.st_size/pageSize;
  long int off = 0;
  long int firstPos;
  long int lastPos;
  firstPos = 0;
  fileptr = (char*)mmap64(NULL, pageSize, PROT_READ, MAP_PRIVATE | 
			  MAP_POPULATE, fd, firstPos*pageSize);
  data =  fileptr;
  string FirstPageFirstHash;
  string FirstPageLastHash;
  ProcessPage(data, FirstPageFirstHash, FirstPageLastHash, pageSize);
  
  if (munmap(fileptr, pageSize) == -1) {
    perror("Error un-mmapping the file");
  }
  
  if (hash >= FirstPageFirstHash and hash <= FirstPageLastHash) {
    fileptr = (char*)mmap64(NULL, pageSize, PROT_READ, MAP_PRIVATE | 
			    MAP_POPULATE, fd, firstPos*pageSize);
    int val =  checkPage(data, hash, pageSize, "");
    if (munmap(fileptr, pageSize) == -1) {
      perror("Error un-mmapping the file");
    }
    return val; 
  }

  if (hash < FirstPageFirstHash) {
      cout << "HASH NOT IN FILE " << hash << endl;
      return 0; 
  }
  
  fileptr = (char*)mmap64(NULL, pageSize, PROT_READ, MAP_PRIVATE |
			  MAP_POPULATE, fd, pageSize*(NumPages-1));
  data = fileptr;
  lastPos = NumPages-1;
  string LastPageFirstHash;
  string LastPageLastHash;
  ProcessPage(data, LastPageFirstHash, LastPageLastHash, pageSize);
  
  if (munmap(fileptr, pageSize) == -1) {
    perror("Error un-mmapping the file");
  }
  
  if (hash >= LastPageFirstHash and hash <= LastPageLastHash) {
    fileptr = (char*)mmap64(NULL, pageSize, PROT_READ, MAP_PRIVATE | 
			    MAP_POPULATE, fd, pageSize*(NumPages-1));
    int val =  checkPage(data, hash, pageSize, "");
    if (munmap(fileptr, pageSize) == -1) {
      perror("Error un-mmapping the file");
    }
    return val; 
  }

  if (hash > LastPageLastHash) {
    cout << "HASH NOT IN FILE " << hash << endl;
    return 0;
  }

  int counter = 0;

  while (true) {
    counter++;
    long int currentPage = lastPos - ((lastPos-firstPos)/2);

    if (currentPage == lastPos or currentPage == firstPos or 
	lastPos - firstPos < 3) {
      string extra = "";
      fileptr = (char*)mmap64(NULL, pageSize*5, PROT_READ, MAP_PRIVATE 
			      | MAP_POPULATE, fd, pageSize * (firstPos-1));
      data =  fileptr;
      int val = checkPage(data, hash, pageSize*5, extra);

      if (munmap(fileptr, pageSize*5) == -1) {
	perror("Error un-mmapping the file");
      }

      return val; 
    }
    
    fileptr = (char*)mmap64(NULL, pageSize, PROT_READ, MAP_PRIVATE | 
			    MAP_POPULATE, fd, pageSize*currentPage);
    data =  fileptr;
    string CurrentPageFirstHash;
    string CurrentPageLastHash;
    ProcessPage(data, CurrentPageFirstHash, CurrentPageLastHash, pageSize);
    
    if (munmap(fileptr, pageSize) == -1) {
      perror("Error un-mmapping the file");
    }
    
    if (hash >= CurrentPageFirstHash and hash <= CurrentPageLastHash) {
      fileptr = (char*)mmap64(NULL, pageSize, PROT_READ, MAP_PRIVATE | 
			      MAP_POPULATE, fd, pageSize*currentPage);
      int val = checkPage(data, hash, pageSize, "");
      if (munmap(fileptr, pageSize) == -1) {
	perror("Error un-mmapping the file");
      }
      return val;
    } else {
      if (hash < CurrentPageFirstHash) {
	lastPos = currentPage;
	LastPageFirstHash = CurrentPageFirstHash;
	LastPageLastHash = CurrentPageLastHash;
      } else if (hash > CurrentPageLastHash) {
	firstPos = currentPage;
	FirstPageFirstHash = CurrentPageFirstHash;
	FirstPageLastHash = CurrentPageLastHash;
      }
    }
  }
  close(fd); 
}

//return value of 0 represents forward strand
//return value of 1 represents reverse strand
int GetReadOrientation(int flag) {
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

string getHash(string seq, int j, int HashSize) {
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

string compressVar(string line, int start, string& StructCall) {
  cout << "compressing var" << endl;
  char current = line.c_str()[0];
  int currentCount = 1;
  string CV = "";
  for (int i = 1; i < line.size(); i++) {
    cout << current << endl;
    if (line.c_str()[i] == current) {
      currentCount++;
    } else {
      if (currentCount > 2) {
        ostringstream convert;
        convert << currentCount;
        CV += convert.str();
        CV += current;

        ostringstream convertEND;
        int end = currentCount + start;
        convertEND << end;

        if (current == 'Y') {
          cout << "YAAAY STRUCT" << endl;
          StructCall = "SVTYPE=DUP;END=";
          StructCall += convertEND.str();
          StructCall += ";SVLEN=";
          StructCall += convert.str();
          StructCall += ";";
          cout << StructCall << endl;
        }
      } else if (currentCount == 2) {
        CV += current;
        CV += current;
      } else if (currentCount == 1) {
        CV += current;
      } else {
        cout << "ERROR in compress " << current << " " << currentCount << endl;
      }

      current = line.c_str()[i];
      currentCount = 1;
    }
  }
  if (currentCount > 2) {
    ostringstream convert;
    convert << currentCount;
    CV += convert.str();
    CV += current;

    ostringstream convertEND;
    int end = currentCount + start;
    convertEND << end;

    if (current == 'Y') {
      cout << "YAAAY STRUCT" << endl;
      StructCall = "SVTYPE=DUP:TANDEM;END=";
      StructCall += convertEND.str();
      StructCall += ";SVLEN=";
      StructCall += convert.str();
      StructCall += ";";
      cout << StructCall << endl;
    }
  } else if (currentCount == 2) {
    CV += current;
    CV += current;
  } else if (currentCount == 1) {
    CV += current;
  } else {
    cout << "ERROR in compress " << current << " " << currentCount << endl;
  }
  return CV;
}

int findBreak(SamRead& read) {
  char Afirst = read.cigarString.c_str()[0];

  cout << "Afirst = " << Afirst << endl;
  cout << "starting A check " << endl;
  if (Afirst == 'H' or Afirst == 'S') {
    cout << "forward" << endl;
    for (int i = 0; i < read.seq.size(); i++) {
      if (read.cigarString.c_str()[i] == 'H' or
          read.cigarString.c_str()[i] == 'S') {
        // keep going
        cout << i << " == " << read.cigarString.c_str()[i] << ' '
             << read.seq.c_str()[i] << endl;
      } else {
        cout << i << " == " << read.cigarString.c_str()[i] << ' '
             << read.seq.c_str()[i] << endl;
        cout << "fond break at " << i << endl;
        return i;
      }
    }
  } else {
    cout << "reverse" << endl;
    for (int i = read.seq.size() - 1; i >= 0; i += -1) {
      if (read.cigarString.c_str()[i] == 'H' or
          read.cigarString.c_str()[i] == 'S') {
        cout << i << " == " << read.cigarString.c_str()[i] << ' '
             << read.seq.c_str()[i] << endl;
        // keep going
      } else {
        cout << i << " == " << read.cigarString.c_str()[i] << ' '
             << read.seq.c_str()[i] << endl;
        cout << "fond break at " << i << endl;
        return i;
      }
    }
  }
}

SamRead BetterWay(vector<SamRead> reads) {
  cout << "In BetterWay" << endl;
  int A = 0;
  int B = 1;
  // cout << "B = " << B << endl;
  cout << "Working on " << reads[A].name
       << ", size = " << reads[B].pos - reads[A].pos << endl;
   vector<vector<int>> AlignmentPos;
  vector<vector<string>> AlignmentChr;
  int ALastGoodRef = -1;
  string ALastGoodChr = "nope";
  int BLastGoodRef = -1;
  string BLastGoodChr = "nope";

  vector<string> NewSeqs;
  vector<string> NewQuals;
  vector<string> NewRefs;
  vector<string> NewCigars;

  for (int i = 0; i < reads.size(); i++) {
    NewSeqs.push_back("");
    NewQuals.push_back("");
    NewRefs.push_back("");
    NewCigars.push_back("");
  }

  for (int i = 0; i < reads.size(); i++) {
    cout << "Pre Lining up reads " << i << endl;
    reads[i].write();
  }
  int Acount = 0;
  int Bcount = 0;
  if (B == 1 and
      GetReadOrientation(reads[A].flag) != GetReadOrientation(reads[B].flag)) {
    cout << "FLIPPING reads not on the same strand";
    reads[B].flipRead();
  }
  while (Acount < reads[A].cigarString.size() and
         Bcount < reads[B].cigarString.size()) {
    vector<int> currentPos;
    vector<string> currentChr;
    // need to get all the reads lined up with the same number of bases, taking
    // account of I's and D's that change length//
    if (reads[A].cigarString.c_str()[Acount] == 'D' and
        reads[B].cigarString.c_str()[Acount] != 'D') {
      currentPos.push_back(reads[A].Positions[Acount]);
      currentPos.push_back(-1);

      currentChr.push_back(reads[A].ChrPositions[Acount]);
      currentChr.push_back("nope");

      ALastGoodRef = reads[A].Positions[Acount];
      ALastGoodChr = reads[A].ChrPositions[Acount];

      NewSeqs[A] += reads[A].seq.c_str()[Acount];
      NewQuals[A] += reads[A].qual.c_str()[Acount];
      NewRefs[A] += reads[A].RefSeq.c_str()[Acount];
      NewCigars[A] += reads[A].cigarString.c_str()[Acount];

      Acount++;

      NewSeqs[B] += '-';
      NewQuals[B] += '!';
      NewRefs[B] += "-";
      NewCigars[B] += 'R';
    } else if (reads[A].cigarString.c_str()[Acount] != 'D' and
               reads[B].cigarString.c_str()[Acount] == 'D') {
      currentPos.push_back(-1);
      currentPos.push_back(reads[B].Positions[Bcount]);

      currentChr.push_back("nope");
      currentChr.push_back(reads[B].ChrPositions[Bcount]);

      BLastGoodRef = reads[B].Positions[Bcount];
      BLastGoodChr = reads[B].ChrPositions[Bcount];

      NewSeqs[B] += reads[B].seq.c_str()[Bcount];
      NewQuals[B] += reads[B].qual.c_str()[Bcount];
      NewRefs[B] += reads[B].RefSeq.c_str()[Bcount];
      NewCigars[B] += reads[B].cigarString.c_str()[Bcount];

      Bcount++;

      NewSeqs[A] += '-';
      NewQuals[A] += '!';
      NewRefs[A] += '-';
      NewCigars[A] += 'R';
    } else {
      if (reads[A].cigarString.c_str()[Acount] == 'H' or
          reads[A].cigarString.c_str()[Acount] == 'S') {
        currentPos.push_back(-1);
        currentChr.push_back("nope");

        NewSeqs[A] += reads[A].seq.c_str()[Acount];
        NewQuals[A] += reads[A].qual.c_str()[Acount];
        NewRefs[A] += reads[A].RefSeq.c_str()[Acount];
        NewCigars[A] += reads[A].cigarString.c_str()[Acount];

        Acount++;

      } else if (reads[A].cigarString.c_str()[Acount] == 'M' or
                 reads[A].cigarString.c_str()[Acount] == 'X' or
                 reads[A].cigarString.c_str()[Acount] == 'D') {
        currentPos.push_back(reads[A].Positions[Acount]);
        currentChr.push_back(reads[A].ChrPositions[Acount]);
        ALastGoodRef = reads[A].Positions[Acount];
        ALastGoodChr = reads[A].ChrPositions[Acount];

        NewSeqs[A] += reads[A].seq.c_str()[Acount];
        NewQuals[A] += reads[A].qual.c_str()[Acount];
        NewRefs[A] += reads[A].RefSeq.c_str()[Acount];
        NewCigars[A] += reads[A].cigarString.c_str()[Acount];

        Acount++;
      } else if (reads[A].cigarString.c_str()[Acount] == 'I') {
        currentPos.push_back(ALastGoodRef);
        currentChr.push_back(ALastGoodChr);

        NewSeqs[A] += reads[A].seq.c_str()[Acount];
        NewQuals[A] += reads[A].qual.c_str()[Acount];
        NewRefs[A] += reads[A].RefSeq.c_str()[Acount];
        NewCigars[A] += reads[A].cigarString.c_str()[Acount];

        Acount++;
      } else
        cout << "WTF, cigar = " << reads[A].cigarString.c_str()[Acount];

      if (reads[B].cigarString.c_str()[Bcount] == 'H' or
          reads[B].cigarString.c_str()[Bcount] == 'S') {
        currentPos.push_back(-1);
        currentChr.push_back("nope");

        NewSeqs[B] += reads[B].seq.c_str()[Bcount];
        NewQuals[B] += reads[B].qual.c_str()[Bcount];
        NewRefs[B] += reads[B].RefSeq.c_str()[Bcount];
        NewCigars[B] += reads[B].cigarString.c_str()[Bcount];

        Bcount++;
      } else if (reads[B].cigarString.c_str()[Bcount] == 'M' or
                 reads[B].cigarString.c_str()[Bcount] == 'X' or
                 reads[B].cigarString.c_str()[Bcount] == 'D') {
        currentPos.push_back(reads[B].Positions[Bcount]);
        currentChr.push_back(reads[B].ChrPositions[Bcount]);
        BLastGoodRef = reads[B].Positions[Bcount];
        BLastGoodChr = reads[B].ChrPositions[Bcount];

        NewSeqs[B] += reads[B].seq.c_str()[Bcount];
        NewQuals[B] += reads[B].qual.c_str()[Bcount];
        NewRefs[B] += reads[B].RefSeq.c_str()[Bcount];
        NewCigars[B] += reads[B].cigarString.c_str()[Bcount];

        Bcount++;
      } else if (reads[B].cigarString.c_str()[Bcount] == 'I') {
        currentPos.push_back(BLastGoodRef);
        currentChr.push_back(BLastGoodChr);

        NewSeqs[B] += reads[B].seq.c_str()[Bcount];
        NewQuals[B] += reads[B].qual.c_str()[Bcount];
        NewRefs[B] += reads[B].RefSeq.c_str()[Bcount];
        NewCigars[B] += reads[B].cigarString.c_str()[Bcount];

        Bcount++;
      } else
        cout << "WTF, cigar = " << reads[B].cigarString.c_str()[Bcount];
    }
    AlignmentPos.push_back(currentPos);
    AlignmentChr.push_back(currentChr);
  }

  for (int i = 0; i < reads.size(); i++) {
    reads[i].Positions.clear();
    reads[i].ChrPositions.clear();
    reads[i].seq = NewSeqs[i];
    reads[i].qual = NewQuals[i];
    reads[i].RefSeq = NewRefs[i];
    reads[i].cigarString = NewCigars[i];
    for (int j = 0; j < AlignmentPos.size(); j++) {
      reads[i].Positions.push_back(AlignmentPos[j][i]);
    }
    for (int j = 0; j < AlignmentPos.size(); j++) {
      reads[i].ChrPositions.push_back(AlignmentChr[j][i]);
    }
  }

  bool deletion = true;
  string NewCigar = "";
  string NewSeq = "";
  string NewQual = "";
  string NewRef = "";
  vector<int> NewPos;
  vector<string> NewChr;

  char LastAlignedQ = ' ';
  int LastAlignedPos = -1;
  string LastAlignedChr = "nope";

  // set LastAlignedPos to the first base with an aligned base
  bool notfound = true;
  int base = 0;
  while (notfound) {
    for (int i = 0; i < reads.size(); i++) {
      if (reads[i].Positions[base] > -1) {
        LastAlignedPos = reads[i].Positions[base];
        LastAlignedChr = reads[i].ChrPositions[base];
        notfound = false;
        break;
      }
    }
    if (notfound) base++;
  }
  for (int i = 0; i < base; i++) {
    NewCigar += reads[A].cigarString.c_str()[i];
    NewSeq += reads[A].seq.c_str()[i];
    NewQual += reads[A].qual.c_str()[i];
    NewRef += reads[A].RefSeq.c_str()[i];
    NewPos.push_back(reads[A].Positions[i]);
    NewChr.push_back(reads[A].ChrPositions[i]);
  }
  // corect qualites so everyone has the same ones, H will produce no quality
  cout << "checking quals" << endl;
  string bestQual = reads[0].qual;
  for (int i = 0; i < reads.size(); i++) {
    bool h = false;
    cout << "read " << reads[i].name << endl;
    for (int j = 0; j < reads[i].RefSeq.size(); j++) {
      if (reads[i].RefSeq.c_str()[j] == 'H') {
        h = true;
      }
    }
    if (h) {
      cout << "contains H" << reads[i].RefSeq << endl;
    } else {
      cout << "does not contain H " << reads[i].RefSeq << endl;
      bestQual = reads[i].qual;
    }
  }
  for (int i = 0; i < reads.size(); i++) {
    reads[i].qual = bestQual;
  }

  for (int i = 0; i < reads.size(); i++) {
    reads[i].createPeakMap();
  }

  cout << "Post adjustment" << endl;
  reads[A].write();
  reads[B].write();

  cout << "************INTO**********" << endl;
  if (B == 1 and GetReadOrientation(reads[A].flag) ==
      GetReadOrientation(reads[B].flag))  // if they are on the
     {
      for (int i = base; i < reads[A].seq.size(); i++) {
	if (reads[A].Positions[i] > -1)  // if this base is aligned in A
	  {
	    if (reads[A].Positions[i] - LastAlignedPos > 1)  // indicates a deletion
	      {
		if (LastAlignedQ == '!' or reads[A].qual.c_str()[i] == '!') {
		  cout << "Unexpected behavior in RUFUS.interpret.cpp A" << endl;
		  return reads[A];
		}
		if (reads[A].chr == reads[B].chr and
		    abs(reads[A].Positions[i] - LastAlignedPos) <
		    MaxVarentSize)  // must be on the same chromosome
		  {
		    cout << "wel this dosnt make any sense"
			 << endl;  // reads are in order in the bam so A should always
		    // be downstream of B, theus the deletion shoould be
		    // detected in B
		    BEDBigStuff << reads[A].chr << "\t" << LastAlignedPos << "\t"
				<< reads[A].Positions[i] << "\t"
				<< "Deletion" << endl;
		    for (int j = LastAlignedPos; j < reads[A].Positions[i] - 1; j++) {
		      NewCigar += 'D';
		      NewSeq += '-';
		      NewQual += LastAlignedQ;
              NewRef +=
		toupper(Reff.getSubSequence(reads[A].chr, i, 1).c_str()[0]);
              NewPos.push_back(j);
              NewChr.push_back(reads[A].ChrPositions[i]);
		    }
		  } else {
		  // if( reads[A].ChrPositions[i] == LastAlignedChr and
		  // abs(reads[A].Positions[i] -LastAlignedPos ) >= MaxVarentSize )
		  if (reads[A].chr == reads[B].chr and
		      abs(reads[A].Positions[i] - LastAlignedPos) >= MaxVarentSize) {
		    int Abreak = findBreak(reads[A]);
		    int Bbreak = findBreak(reads[B]);
		    cout << " if( " << reads[A].PeakMap[Abreak] << " == 1 or "
			 << reads[A].PeakMap[Abreak - 1] << " == 1 or "
			 << reads[B].PeakMap[Bbreak] << " == 1 or "
			 << reads[B].PeakMap[Bbreak - 1] << " == 1)";
		    if (/*1==1 or*/ (reads[A].PeakMap[Abreak] == 1 or
				     reads[A].PeakMap[Abreak - 1] == 1) and
			(reads[B].PeakMap[Bbreak] == 1 or
			 reads[B].PeakMap[Bbreak - 1] == 1) and
			Abreak > 0 and Bbreak > 0) {
		      cout << "INVERSION written to file" << endl;
		      Translocations << "Too Big, Same strand and chr "
				     << abs(reads[A].Positions[i] - LastAlignedPos)
				     << endl;
		      reads[A].writetofile(Translocations);
		      reads[B].writetofile(Translocations);
		      Translocations << endl << endl;
                Translocationsbed
		  << reads[A].chr << "\t" << reads[A].Positions[Abreak] - 200
		  << "\t" << reads[A].Positions[Abreak] + 200 << endl
		  << reads[B].chr << "\t" << reads[B].Positions[Bbreak] - 200
		  << "\t" << reads[B].Positions[Bbreak] + 200 << endl;
		    } else {
		      cout << "INVERSION skipped" << endl;
		    }
		  } else {
		    if ((reads[A].chr == "hs37d5" and reads[B].chr != "hs37d5") or
			(reads[A].chr != "hs37d5" and reads[B].chr == "hs37d5")) {
		      int Abreak = findBreak(reads[A]);
		      int Bbreak = findBreak(reads[B]);
		      cout << " if( " << reads[A].PeakMap[Abreak] << " == 1 or "
			   << reads[A].PeakMap[Abreak - 1] << " == 1 or "
			   << reads[B].PeakMap[Bbreak] << " == 1 or "
			   << reads[B].PeakMap[Bbreak - 1] << " == 1)";
		      if (/*1==1 or*/ (reads[A].PeakMap[Abreak] == 1 or
				       reads[A].PeakMap[Abreak - 1] == 1) and
			  (reads[B].PeakMap[Bbreak] == 1 or
			   reads[B].PeakMap[Bbreak - 1] == 1) and
			  Abreak > 0 and Bbreak > 0) {
			cout << "INVERSION written to file" << endl;
			Translocations << "Possible mob event "
				       << abs(reads[A].Positions[i] - LastAlignedPos)
				       << endl;
			reads[A].writetofile(Translocations);
			reads[B].writetofile(Translocations);
			Translocations << endl << endl;
			Translocationsbed << reads[A].chr << "\t"
					  << reads[A].Positions[Abreak] - 200 << "\t"
					  << reads[A].Positions[Abreak] + 200 << endl
					  << reads[B].chr << "\t"
					  << reads[B].Positions[Bbreak] - 200 << "\t"
					  << reads[B].Positions[Bbreak] + 200 << endl;
		      } else {
			cout << "INVERSION skipped" << endl;
		      }
		    } else {
		      int Abreak = findBreak(reads[A]);
		      int Bbreak = findBreak(reads[B]);
		      cout << " if( " << reads[A].PeakMap[Abreak] << " == 1 or "
			   << reads[A].PeakMap[Abreak - 1] << " == 1 or "
			   << reads[B].PeakMap[Bbreak] << " == 1 or "
			   << reads[B].PeakMap[Bbreak - 1] << " == 1)";
		      if (/*1==1 or*/ (reads[A].PeakMap[Abreak] == 1 or
				       reads[A].PeakMap[Abreak - 1] == 1) and
			  (reads[B].PeakMap[Bbreak] == 1 or
			   reads[B].PeakMap[Bbreak - 1] == 1) and
			  Abreak > 0 and Bbreak > 0) {
			cout << "INVERSION written to file" << endl;
			Translocations << "Translocataion, same strand "
				       << abs(reads[A].Positions[i] - LastAlignedPos)
				       << endl;
			reads[A].writetofile(Translocations);
			reads[B].writetofile(Translocations);
			Translocations << endl << endl;
			Translocationsbed << reads[A].chr << "\t"
					  << reads[A].Positions[Abreak] - 200 << "\t"
					  << reads[A].Positions[Abreak] + 200 << endl
					  << reads[B].chr << "\t"
					  << reads[B].Positions[Bbreak] - 200 << "\t"
					  << reads[B].Positions[Bbreak] + 200 << endl;
		      } else {
			cout << "INVERSION skipped" << endl;
		      }
		    }
		  }
		}
	      } else if (reads[A].Positions[i] - LastAlignedPos < 0 and
			 abs(reads[A].Positions[i] - LastAlignedPos) <
			 MaxVarentSize)  // indicates a possible insertion or
	      // tandem duplication
	      {
		if (LastAlignedQ == '!' or reads[A].qual.c_str()[i] == '!') {
		  cout << "Unexpected behavior in RUFUS.interpret.cpp B" << endl;
		  return reads[A];
		}
		cout << "this could be one A, last = " << LastAlignedPos
		     << " Current = " << reads[A].Positions[i]
		     << " chr = " << LastAlignedChr << " and "
		     << reads[A].ChrPositions[i] << endl;
		if (reads[A].chr == reads[B].chr) {
		  cout << "This is an insertion A at base " << i << endl;
		  BEDBigStuff << reads[A].chr << "\t" << reads[A].Positions[i] << "\t"
			      << LastAlignedPos << "\tTandemDup" << endl;
		  cout << "tadem dup" << endl;
		  int j = 0;
		  for (j = i; j < reads[A].seq.size() and
			 reads[A].Positions[j] < LastAlignedPos;
		       j++) {
		    NewCigar += 'Y';  //'I';
		    NewSeq += reads[A].seq.c_str()[j];
		    NewQual += reads[A].qual.c_str()[j];
		    NewRef += '-';
		    NewPos.push_back(reads[A].Positions[i]);
		    NewChr.push_back(reads[A].ChrPositions[i]);
		  }
		  i = j;
		  // find the last base that was aligned, htere can be novel insertion
		  // stuff so you cant jsut take the last base
		  int k;
		  for (k = reads[A].Positions.size() - 1; k >= 0; k += -1) {
		    if (reads[A].Positions[k] + 1 > 1) break;
		  }
		  for (j = reads[A].Positions[k] + 1; j < LastAlignedPos; j++) {
		    NewCigar += 'Y';  //'I';
              NewSeq +=
		toupper(Reff.getSubSequence(reads[A].chr, j, 1).c_str()[0]);
              NewQual += '!';
              NewRef += '-';
              NewPos.push_back(j);
              NewChr.push_back(reads[A].chr);
		  }
		} else {
		  if ((reads[A].chr == "hs37d5" and reads[B].chr != "hs37d5") or
		      (reads[A].chr != "hs37d5" and reads[B].chr == "hs37d5")) {
		    int Abreak = findBreak(reads[A]);
		    int Bbreak = findBreak(reads[B]);
		    cout << " if( " << reads[A].PeakMap[Abreak] << " == 1 or "
			 << reads[A].PeakMap[Abreak - 1] << " == 1 or "
			 << reads[B].PeakMap[Bbreak] << " == 1 or "
			 << reads[B].PeakMap[Bbreak - 1] << " == 1)";
		    if (/*1==1 or*/ (reads[A].PeakMap[Abreak] == 1 or
				     reads[A].PeakMap[Abreak - 1] == 1) and
			(reads[B].PeakMap[Bbreak] == 1 or
			 reads[B].PeakMap[Bbreak - 1] == 1) and
			Abreak > 0 and Bbreak > 0) {
		      cout << "INVERSION written to file" << endl;
		      Translocations << "Possible mob event "
				     << abs(reads[A].Positions[i] - LastAlignedPos)
				     << endl;
		      reads[A].writetofile(Translocations);
		      reads[B].writetofile(Translocations);
		      Translocations << endl << endl;
                Translocationsbed
		  << reads[A].chr << "\t" << reads[A].Positions[Abreak] - 200
		  << "\t" << reads[A].Positions[Abreak] + 200 << endl
		  << reads[B].chr << "\t" << reads[B].Positions[Bbreak] - 200
		  << "\t" << reads[B].Positions[Bbreak] + 200 << endl;
		    } else {
		      cout << "INVERSION skipped" << endl;
		    }
		  }
		  int Abreak = findBreak(reads[A]);
		  int Bbreak = findBreak(reads[B]);
		  cout << " if( " << reads[A].PeakMap[Abreak] << " == 1 or "
		       << reads[A].PeakMap[Abreak - 1] << " == 1 or "
		       << reads[B].PeakMap[Bbreak] << " == 1 or "
		       << reads[B].PeakMap[Bbreak - 1] << " == 1)";
		  if (/*1==1 or*/ (reads[A].PeakMap[Abreak] == 1 or
				   reads[A].PeakMap[Abreak - 1] == 1) and
		      (reads[B].PeakMap[Bbreak] == 1 or
		       reads[B].PeakMap[Bbreak - 1] == 1) and
		      Abreak > 0 and Bbreak > 0) {
		    cout << "INVERSION written to file" << endl;
		    Translocations << "Translocation, same strand "
				   << abs(reads[A].Positions[i] - LastAlignedPos)
				   << endl;
		    reads[A].writetofile(Translocations);
		    reads[B].writetofile(Translocations);
		    Translocations << endl << endl;
              Translocationsbed
		<< reads[A].chr << "\t" << reads[A].Positions[Abreak] - 200
		<< "\t" << reads[A].Positions[Abreak] + 200 << endl
		<< reads[B].chr << "\t" << reads[B].Positions[Bbreak] - 200
		<< "\t" << reads[B].Positions[Bbreak] + 200 << endl;
		  } else {
		    cout << "INVERSION skipped" << endl;
		  }
		}

	      } else if (reads[A].Positions[i] - LastAlignedPos < 0 and
			 abs(reads[A].Positions[i] - LastAlignedPos) >=
			 MaxVarentSize) {
	      if (LastAlignedQ == '!' or reads[A].qual.c_str()[i] == '!') {
		cout << "Unexpected behavior in RUFUS.interpret.cpp C" << endl;

		return reads[A];
	      }
	      int Abreak = findBreak(reads[A]);
	      int Bbreak = findBreak(reads[B]);
	      cout << " if( " << reads[A].PeakMap[Abreak] << " == 1 or "
		   << reads[A].PeakMap[Abreak - 1] << " == 1 or "
		   << reads[B].PeakMap[Bbreak] << " == 1 or "
		   << reads[B].PeakMap[Bbreak - 1] << " == 1)";
	      if (/*1==1 or*/ (reads[A].PeakMap[Abreak] == 1 or
			       reads[A].PeakMap[Abreak - 1] == 1) and
		  (reads[B].PeakMap[Bbreak] == 1 or
		   reads[B].PeakMap[Bbreak - 1] == 1) and
		  Abreak > 0 and Bbreak > 0) {
		cout << "INVERSION written to file" << endl;
		if (reads[A].chr == reads[B].chr)
		  Translocations << "TOO BIG 3 "
				 << abs(reads[A].Positions[i] - LastAlignedPos)
				 << endl;
		else
		  Translocations << "Translocation 3 "
				 << abs(reads[A].Positions[i] - LastAlignedPos)
				 << endl;
		reads[A].writetofile(Translocations);
		reads[B].writetofile(Translocations);
		Translocations << endl << endl;
            Translocationsbed
	      << reads[A].chr << "\t" << reads[A].Positions[Abreak] - 200
	      << "\t" << reads[A].Positions[Abreak] + 200 << endl
	      << reads[B].chr << "\t" << reads[B].Positions[Bbreak] - 200
	      << "\t" << reads[B].Positions[Bbreak] + 200 << endl;
	      } else {
		cout << "INVERSION skipped" << endl;
	      }
	    }

	    if (i < reads[A].cigarString.size()) {
	      NewCigar += reads[A].cigarString.c_str()[i];
	      NewSeq += reads[A].seq.c_str()[i];
	      NewQual += reads[A].qual.c_str()[i];
	      NewRef += reads[A].RefSeq.c_str()[i];
	      NewPos.push_back(reads[A].Positions[i]);
	      NewChr.push_back(reads[A].ChrPositions[i]);
	      LastAlignedQ = reads[A].qual.c_str()[i];
	      LastAlignedPos = reads[A].Positions[i];
	      LastAlignedChr = reads[A].ChrPositions[i];
	    }
	  } else if (reads[B].Positions[i] > -1) {
	  if (reads[B].Positions[i] - LastAlignedPos > 1) {
	    if (LastAlignedQ == '!' or reads[B].qual.c_str()[i] == '!') {
	      cout << "Unexpected behavior in RUFUS.interpret.cpp D" << endl;

	      return reads[A];
	    }
	    if (reads[B].chr == reads[A].chr and
		abs(reads[B].Positions[i] - LastAlignedPos) < MaxVarentSize) {
	      cout << "striahgtup deletion, size = "
		   << abs(reads[B].Positions[i] - LastAlignedPos) << " at base "
		   << i << " from Position " << LastAlignedPos << " to "
		   << reads[B].Positions[i] << endl;
	      BEDBigStuff << reads[B].chr << "\t" << LastAlignedPos << "\t"
			  << reads[B].Positions[i] << "\t"
			  << "Deletion" << endl;
	      for (int j = LastAlignedPos; j < reads[B].Positions[i] - 1; j++) {
		// cout << j << " - " << j - LastAlignedPos<< endl;
		NewCigar += 'D';
		NewSeq += '-';
		NewQual += LastAlignedQ;
              NewRef +=
		toupper(Reff.getSubSequence(reads[B].chr, j, 1).c_str()[0]);
              NewPos.push_back(j);
              NewChr.push_back(reads[B].ChrPositions[i]);
	      }
	    } else {
	      if (reads[A].chr == reads[B].chr and
		  abs(reads[B].Positions[i] - LastAlignedPos) >= MaxVarentSize) {
		int Abreak = findBreak(reads[A]);
		int Bbreak = findBreak(reads[B]);
		cout << " if( " << reads[A].PeakMap[Abreak] << " == 1 or "
		     << reads[A].PeakMap[Abreak - 1] << " == 1 or "
		     << reads[B].PeakMap[Bbreak] << " == 1 or "
		     << reads[B].PeakMap[Bbreak - 1] << " == 1)";
		if (/*1==1 or*/ (reads[A].PeakMap[Abreak] == 1 or
				 reads[A].PeakMap[Abreak - 1] == 1) and
		    (reads[B].PeakMap[Bbreak] == 1 or
		     reads[B].PeakMap[Bbreak - 1] == 1) and
		    Abreak > 0 and Bbreak > 0) {
		  cout << "INVERSION written to file" << endl;
		  Translocations << "TOO BIG 2 "
				 << abs(reads[A].Positions[i] - LastAlignedPos)
				 << endl;
		  reads[A].writetofile(Translocations);
		  reads[B].writetofile(Translocations);
		  Translocations << endl << endl;
                Translocationsbed
		  << reads[A].chr << "\t" << reads[A].Positions[Abreak] - 200
		  << "\t" << reads[A].Positions[Abreak] + 200 << endl
		  << reads[B].chr << "\t" << reads[B].Positions[Bbreak] - 200
		  << "\t" << reads[B].Positions[Bbreak] + 200 << endl;
		} else {
		  cout << "INVERSION skipped" << endl;
		}
	      } else {
		if ((reads[A].chr == "hs37d5" and reads[B].chr != "hs37d5") or
		    (reads[A].chr != "hs37d5" and reads[B].chr == "hs37d5")) {
		  int Abreak = findBreak(reads[A]);
		  int Bbreak = findBreak(reads[B]);
		  cout << " if( " << reads[A].PeakMap[Abreak] << " == 1 or "
		       << reads[A].PeakMap[Abreak - 1] << " == 1 or "
		       << reads[B].PeakMap[Bbreak] << " == 1 or "
		       << reads[B].PeakMap[Bbreak - 1] << " == 1)";
		  if (/*1==1 or*/ (reads[A].PeakMap[Abreak] == 1 or
				   reads[A].PeakMap[Abreak - 1] == 1) and
		      (reads[B].PeakMap[Bbreak] == 1 or
		       reads[B].PeakMap[Bbreak - 1] == 1) and
		      Abreak > 0 and Bbreak > 0) {
		    cout << "INVERSION written to file" << endl;
		    Translocations << "Possible mob event "
				   << abs(reads[A].Positions[i] - LastAlignedPos)
				   << endl;
		    reads[A].writetofile(Translocations);
		    reads[B].writetofile(Translocations);
		    Translocations << endl << endl;
		    Translocationsbed << reads[A].chr << "\t"
				      << reads[A].Positions[Abreak] - 200 << "\t"
				      << reads[A].Positions[Abreak] + 200 << endl
				      << reads[B].chr << "\t"
				      << reads[B].Positions[Bbreak] - 200 << "\t"
				      << reads[B].Positions[Bbreak] + 200 << endl;
		  } else {
		    cout << "INVERSION skipped" << endl;
		  }
		} else {
		  int Abreak = findBreak(reads[A]);
		  int Bbreak = findBreak(reads[B]);
		  cout << " if( " << reads[A].PeakMap[Abreak] << " == 1 or "
		       << reads[A].PeakMap[Abreak - 1] << " == 1 or "
		       << reads[B].PeakMap[Bbreak] << " == 1 or "
		       << reads[B].PeakMap[Bbreak - 1] << " == 1)";
		  if (/*1==1 or*/ (reads[A].PeakMap[Abreak] == 1 or
				   reads[A].PeakMap[Abreak - 1] == 1) and
		      (reads[B].PeakMap[Bbreak] == 1 or
		       reads[B].PeakMap[Bbreak - 1] == 1) and
		      Abreak > 0 and Bbreak > 0) {
		    cout << "INVERSION written to file" << endl;
		    Translocations << "Translocation 2 "
				   << abs(reads[A].Positions[i] - LastAlignedPos)
				   << endl;
		    reads[A].writetofile(Translocations);
		    reads[B].writetofile(Translocations);
		    Translocations << endl << endl;
		    Translocationsbed << reads[A].chr << "\t"
				      << reads[A].Positions[Abreak] - 200 << "\t"
				      << reads[A].Positions[Abreak] + 200 << endl
				      << reads[B].chr << "\t"
				      << reads[B].Positions[Bbreak] - 200 << "\t"
				      << reads[B].Positions[Bbreak] + 200 << endl;
		  } else {
		    cout << "INVERSION skipped" << endl;
		  }
		}
	      }
	    }

	  } else if (reads[B].Positions[i] - LastAlignedPos < 0 and
		     abs(reads[B].Positions[i] - LastAlignedPos) <
		     MaxVarentSize)  // indicates a possible insertion or
	    // tandem duplication
	    {
	      if (LastAlignedQ == '!' or reads[B].qual.c_str()[i] == '!') {
		cout << "Unexpected behavior in RUFUS.interpret.cpp E" << endl;

		return reads[A];
	      }
	      if (reads[B].chr == reads[A].chr) {
		cout << "This is an insertion B at base " << i << endl;

		BEDBigStuff << reads[B].chr << "\t" << reads[B].Positions[i] << "\t"
			    << LastAlignedPos << "\tTandemDup" << endl;
		int j = 0;
		for (j = i; j < reads[B].seq.size() and
		       reads[B].Positions[j] <= LastAlignedPos;
		     j++) {
		  NewCigar += 'Y';
		  NewSeq += reads[B].seq.c_str()[j];
		  NewQual += reads[B].qual.c_str()[j];
		  NewRef += '-';
		  NewPos.push_back(reads[B].Positions[i]);
		  NewChr.push_back(reads[B].ChrPositions[i]);
		}
		i = j;
		// need to find last base that was alinged, insetion can mess this
		// up so you can just take the last base pos
		int k;
		for (k = reads[B].Positions.size() - 1; k >= 0; k += -1) {
		  if (reads[B].Positions[k] + 1 > 1) break;
		}
		// cout << "j= " << reads[B].Positions[k]+1 << " < " <<
		// LastAlignedPos << " - " << endl;
		for (j = reads[B].Positions[k] + 1; j < LastAlignedPos; j++) {
		  NewCigar += 'Y';
              NewSeq +=
		toupper(Reff.getSubSequence(reads[B].chr, j, 1).c_str()[0]);
              NewQual += '!';
              NewRef += '-';
              NewPos.push_back(j);
              NewChr.push_back(reads[B].chr);
		}

	      } else {
		cout << "we got a translocation" << endl;
		if ((reads[A].chr == "hs37d5" and reads[B].chr != "hs37d5") or
		    (reads[A].chr != "hs37d5" and reads[B].chr == "hs37d5")) {
		  cout << "mobil element " << endl;
		}
	      }

	    } else if (reads[B].Positions[i] - LastAlignedPos < 0 and
		       abs(reads[B].Positions[i] - LastAlignedPos) >=
                       MaxVarentSize) {
	    if (LastAlignedQ == '!' or reads[B].qual.c_str()[i] == '!') {
	      cout << "Unexpected behavior in RUFUS.interpret.cpp F" << endl;
	      return reads[A];
	    }
	    int Abreak = findBreak(reads[A]);
	    int Bbreak = findBreak(reads[B]);

	    cout << " if( " << reads[A].PeakMap[Abreak] << " == 1 or "
		 << reads[A].PeakMap[Abreak - 1] << " == 1 or "
		 << reads[B].PeakMap[Bbreak] << " == 1 or "
		 << reads[B].PeakMap[Bbreak - 1] << " == 1)";
	    if (/*1==1 or*/ (reads[A].PeakMap[Abreak] == 1 or
			     reads[A].PeakMap[Abreak - 1] == 1) and
		(reads[B].PeakMap[Bbreak] == 1 or
		 reads[B].PeakMap[Bbreak - 1] == 1) and
		Abreak > 0 and Bbreak > 0) {
	      cout << "INVERSION written to file" << endl;
	      Translocations << "TOO BIG 1 "
			     << abs(reads[A].Positions[i] - LastAlignedPos)
			     << endl;
	      reads[A].writetofile(Translocations);
	      reads[B].writetofile(Translocations);
	      Translocations << endl << endl;
            Translocationsbed
	      << reads[A].chr << "\t" << reads[A].Positions[Abreak] - 200
	      << "\t" << reads[A].Positions[Abreak] + 200 << endl
	      << reads[B].chr << "\t" << reads[B].Positions[Bbreak] - 200
	      << "\t" << reads[B].Positions[Bbreak] + 200 << endl;
	    } else {
	      cout << "INVERSION skipped" << endl;
	    }
	  }

	  if (i < reads[B].cigarString.size()) {
	    NewCigar += 'M';  // reads[B].cigarString.c_str()[i];
	    NewSeq += reads[B].seq.c_str()[i];
	    NewQual += reads[B].qual.c_str()[i];
	    NewRef += reads[B].RefSeq.c_str()[i];
	    NewPos.push_back(reads[B].Positions[i]);
	    NewChr.push_back(reads[B].ChrPositions[i]);

	    LastAlignedQ = reads[B].qual.c_str()[i];
	    LastAlignedPos = reads[B].Positions[i];
	    LastAlignedChr = reads[B].ChrPositions[i];
	    // cout << "B " << reads[B].seq.c_str()[i] << " " <<
	    //reads[B].cigarString.c_str()[i] << " " << reads[B].RefSeq.c_str()[i]
	    //<<" " << reads[B].ChrPositions[i] << " " << reads[B].Positions[i] <<
	    //NewSeq << endl;
	  }
	} else {
	  if (reads[A].cigarString.c_str()[i] == 'S') {
	    NewCigar += reads[A].cigarString.c_str()[i];
	    NewSeq += reads[A].seq.c_str()[i];
	    NewQual += reads[A].qual.c_str()[i];
	    NewRef += reads[A].RefSeq.c_str()[i];
	    NewPos.push_back(reads[A].Positions[i]);
	    NewChr.push_back(reads[A].ChrPositions[i]);

	  } else if (reads[B].cigarString.c_str()[i] == 'S') {
	    NewCigar += reads[B].cigarString.c_str()[i];
	    NewSeq += reads[B].seq.c_str()[i];
	    NewQual += reads[B].qual.c_str()[i];
	    NewRef += reads[B].RefSeq.c_str()[i];
	    NewPos.push_back(reads[B].Positions[i]);
	    NewChr.push_back(reads[B].ChrPositions[i]);

	  } else
	    cout << "no base works :(" << endl;
	}
      }

      // fix inernal S bases
      int First = -1;
      int Last = -1;
      string NewNewCigar = "";
      for (int i = 0; i < NewCigar.size(); i++) {
	if (NewCigar.c_str()[i] != 'S' and NewCigar.c_str()[i] != 'H') {
	  First = i;
	  break;
	}
      }
      for (int i = NewCigar.size() - 1; i >= 0; i--) {
	if (NewCigar.c_str()[i] != 'S' and NewCigar.c_str()[i] != 'H') {
	  Last = i;
	  break;
	}
      }

      for (int i = 0; i < NewCigar.size(); i++) {
	if (i > First and i < Last) {
	  if (NewCigar.c_str()[i] == 'S' or NewCigar.c_str()[i] == 'H')
	    NewNewCigar += 'I';
	  else
	    NewNewCigar += NewCigar.c_str()[i];
	} else
	  NewNewCigar += NewCigar.c_str()[i];
      }
      NewCigar = NewNewCigar;

      int UnalignedCount = 0;
      for (int i = 0; i < NewCigar.size(); i++) {
	if (NewCigar.c_str()[i] == 'H' || NewCigar.c_str()[i] == 'S')
	  UnalignedCount++;
      }
      cout << "Unaligned Baises = " << UnalignedCount
	   << " %= " << (double)UnalignedCount / (double)NewCigar.size() << endl;
      if (UnalignedCount < 10) {
	NewCigar = NewNewCigar;
	reads[A].first = true;
	int temp = reads[A].alignments[0];
	reads[A].alignments.clear();
	reads[A].alignments.push_back(temp);
	reads[A].cigarString = NewCigar;
	reads[A].seq = NewSeq;
	reads[A].qual = NewQual;
	reads[A].RefSeq = NewRef;
	reads[A].Positions.clear();
	reads[A].ChrPositions.clear();
	reads[A].Positions = NewPos;
	reads[A].ChrPositions = NewChr;
	reads[A].combined = true;
	reads[B].combined = true;
      } else {
	cout << "Skipping" << endl;
      }
      cout << "Done with " << reads[A].name << endl;
     } else {
    
    if (reads[A].chr == reads[B].chr) {
      // need to check each base and find the breakpoint
      cout << "found possible Inversion " << endl;
      
      int Abreak = findBreak(reads[A]);
      int Bbreak = findBreak(reads[B]);
      
      cout << " if( " << reads[A].PeakMap[Abreak] << " == 1 or "
	   << reads[A].PeakMap[Abreak - 1] << " == 1 or "
	   << reads[B].PeakMap[Bbreak] << " == 1 or "
	   << reads[B].PeakMap[Bbreak - 1] << " == 1)";
      if (/*1==1 or*/ (reads[A].PeakMap[Abreak] == 1 or
		       reads[A].PeakMap[Abreak - 1] == 1) and
	  (reads[B].PeakMap[Bbreak] == 1 or
	   reads[B].PeakMap[Bbreak - 1] == 1) and
	  Abreak > 0 and Bbreak > 0) {
	cout << "INVERSION written to file" << endl;
	Translocations << "INVERSION" << endl;
	reads[A].writetofile(Translocations);
	    reads[B].writetofile(Translocations);
	    Translocations << endl << endl;
	    Translocationsbed << reads[A].chr << "\t"
			      << reads[A].Positions[Abreak] - 200 << "\t"
			      << reads[A].Positions[Abreak] + 200 << endl
			      << reads[B].chr << "\t"
			      << reads[B].Positions[Bbreak] - 200 << "\t"
			      << reads[B].Positions[Bbreak] + 200 << endl;
      } else {
	cout << "INVERSION skipped" << endl;
      }
      string Acig = "";
      string Bcig = "";
      for (int i = 0; i < reads[A].seq.size(); i++) {
	char Ab = reads[A].cigarString.c_str()[i];
	char Bb = reads[B].cigarString.c_str()[i];
	if ((reads[A].cigarString.c_str()[i] == 'M' or
	     reads[A].cigarString.c_str()[i] == 'X') and
	    (reads[B].cigarString.c_str()[i] == 'S' or
	     reads[B].cigarString.c_str()[i] == 'H'))
	  Bb = 'U';
	if ((reads[B].cigarString.c_str()[i] == 'M' or
	     reads[B].cigarString.c_str()[i] == 'X') and
	    (reads[A].cigarString.c_str()[i] == 'S' or
	     reads[A].cigarString.c_str()[i] == 'H'))
	  Ab = 'U';
	Acig += Ab;
	Bcig += Bb;
      }
      reads[A].cigarString = Acig;
      reads[B].cigarString = Bcig;
      cout << "invertion adjust string";
      reads[A].write();
      reads[B].write();
    } else if ((reads[A].chr == "hs37d5" and reads[B].chr != "hs37d5") or
	       (reads[A].chr != "hs37d5" and reads[B].chr == "hs37d5")) {
      int Abreak = findBreak(reads[A]);
      int Bbreak = findBreak(reads[B]);
      
      cout << " if( " << reads[A].PeakMap[Abreak] << " == 1 or "
	   << reads[A].PeakMap[Abreak - 1] << " == 1 or "
	   << reads[B].PeakMap[Bbreak] << " == 1 or "
	   << reads[B].PeakMap[Bbreak - 1] << " == 1)";
      if (/*1==1 or*/ (reads[A].PeakMap[Abreak] == 1 or
		       reads[A].PeakMap[Abreak - 1] == 1) and
	    (reads[B].PeakMap[Bbreak] == 1 or
	     reads[B].PeakMap[Bbreak - 1] == 1) and
	  Abreak > 0 and Bbreak > 0) {
	Translocations << "mobil elemnt inverted" << endl;
	reads[A].writetofile(Translocations);
	reads[B].writetofile(Translocations);
	Translocations << endl << endl;
	Translocationsbed << reads[A].chr << "\t"
			  << reads[A].Positions[Abreak] - 200 << "\t"
			  << reads[A].Positions[Abreak] + 200 << endl
			  << reads[B].chr << "\t"
			  << reads[B].Positions[Bbreak] - 200 << "\t"
			  << reads[B].Positions[Bbreak] + 200 << endl;
      }
    } else {
      int Abreak = findBreak(reads[A]);
      int Bbreak = findBreak(reads[B]);
      
      cout << " if( " << reads[A].PeakMap[Abreak] << " == 1 or "
	   << reads[A].PeakMap[Abreak - 1] << " == 1 or "
	   << reads[B].PeakMap[Bbreak] << " == 1 or "
	   << reads[B].PeakMap[Bbreak - 1] << " == 1)";
      if (/*1==1 or*/ (reads[A].PeakMap[Abreak] == 1 or
		       reads[A].PeakMap[Abreak - 1] == 1) and
	  (reads[B].PeakMap[Bbreak] == 1 or
	   reads[B].PeakMap[Bbreak - 1] == 1) and
	  Abreak > 0 and Bbreak > 0) {
	Translocations << "we got a translocation and invertion" << endl;
	reads[A].writetofile(Translocations);
	reads[B].writetofile(Translocations);
	Translocations << endl << endl;
	Translocationsbed << reads[A].chr << "\t"
			  << reads[A].Positions[Abreak] - 200 << "\t"
			  << reads[A].Positions[Abreak] + 200 << endl
			  << reads[B].chr << "\t"
			  << reads[B].Positions[Bbreak] - 200 << "\t"
			  << reads[B].Positions[Bbreak] + 200 << endl;
      }
    }
    
    cout << "*********SKIPPING************\ndifference strands" << endl;
    BEDNotHandled << "Different strands" << endl;
    BEDNotHandled << reads[A].chr << "\t" << reads[A].pos << "\t"
		  << reads[A].pos + reads[A].seq.size() << "\t" << reads[A].name
		  << "\t" << reads[A].cigar << endl;
    reads[A].writetofile(BEDNotHandled);
    BEDNotHandled << reads[B].chr << "\t" << reads[B].pos << "\t"
		  << reads[B].pos + reads[B].seq.size() << "\t" << reads[B].name
		  << "\t" << reads[B].cigar << endl;
    reads[B].writetofile(BEDNotHandled);
    BEDNotHandled << endl << endl;
    Invertions << reads[A].chr << "\t" << reads[A].pos << "\t" << reads[B].pos
		 << "\t" << reads[B].pos - reads[A].pos << endl;
  }
  
  reads[A].LookUpKmers();
  cout << "ReAdjustedKmers" << endl;
  reads[A].writeVertical();
  return reads[A];
  cout << "Out Of BetterWay " << endl;
}


int main (int argc, char *argv[]) {
  cout << "ModeschrpostypereffaltMutRefMutAltPar1RefPar2Ref" << endl;
   string helptext;
  helptext = \
"\
RUFUS.interpret: converts RUFUS aligned contigs into a VCF \n\
By Andrew Farrell\n\
   The Marth Lab\n\
\n\
options:\
  -h [ --help ]  Print help message\n\
  -sam  argPath to input SAM file, omit for stdin\n\
  -r    argPath to reference file \n\
  -hf   argPath to HashFile from RUFUS.build\n\
  -hS   argHash Size\n\
  -o    argOutput stub\n\
  -m    argMaximum varient size: default 1Mb\n\
(Sorry it has to be a num, no 1kb, must be 1000\n\
  -c    argPath to sorted.tab file for the parent sample\n\
  -s    arg Path to sorted.tab file for the subject sample\n\
  -cR   argPath to the sorted.tab file fo the parnt sample hashes in the reference\n\
  -sR   argPath to the sorted.tab file fo the subject sample hashes in the reference\n\
  -mQ   argMinimum map quality to consider varients in\n\
  -mod argPath to the model file from RUFUS.model\n\
";
 
  string MutHashFilePath = "" ;
  string MutHashFilePathReference = "";
  MaxVarentSize = 1000000;
  string RefFile = ""; 
  string HashListFile = "" ; 
  string samFile = "stdin"; 
  string outStub= "";
  string ModelFilePath = "";
  int MinMapQual = 0; 
  for(int i = 1; i < argc; i++) {
      cout << i << " = " << argv[i] << endl; 
    }
  cout <<"****************************************************************************************" << endl;
  vector <int> ParentHashFilePaths; 
  vector <int> ParentHashFilePathsReference;
  for(int i = 1; i< argc; i++) {
    string p = argv[i];
    cout << i << " = " << argv[i]<< endl;
    if (p == "-h") {
      cout << helptext << endl;
      return 0; 
    } else if (p == "-r") {
      RefFile = argv[i+1];
      i=i+1;
      cout << "Added RefFile = " << RefFile << endl;
    } else if (p == "-sam") {
      samFile =  argv[i+1];
      i++;
    } else if (p == "-o") {
      outStub =  argv[i+1];
      i++;
    } else if (p == "-hf") {
      HashListFile =  argv[i+1];
      i++;
    } else if (p == "-hs") {
      HashSize =  atoi(argv[i+1]);
      i++;
    } else if (p == "-m") {
      MaxVarentSize =  atoi(argv[i+1]);
      i++;
      cout << "Added MaxVarSize = " << MaxVarentSize << endl;
    } else if (p == "-c") {
      cout << "Par Hash = " << argv[i+1] << endl;
      ParentHashFilePaths.push_back(i+1);
      i=i+1;
    } else if (p == "-cR") {
      cout << "Par Ref Hash = " << argv[i+1] << endl;
      ParentHashFilePathsReference.push_back(i+1);
      i=i+1;
    } else if (p == "-s") {
      cout << "Sub Hash = " << argv[i+1] << endl;
      MutHashFilePath = argv[i+1];
      i+=1;
    } else if (p == "-sR") {
      cout << "Sub Hash = " << argv[i+1] << endl;
      MutHashFilePathReference = argv[i+1];
      i+=1;
    } else if (p == "-mod") {
      cout << "model file = " << argv[i+1] << endl;
      ModelFilePath = argv[i+1];
      i+=1;
    } else if (p == "-mQ") {
      cout << "Min Mapping Qualtiy = " << argv[i+1] << endl;
      MinMapQual = atoi(argv[i+1]);
      i+=1;
    } else {
      cout << "ERROR: unkown command line paramater -" <<  argv[i] << "-"<< endl;
      return 0; 
    }
  }
  //check values 
  if (RefFile == "") {
    cout << "ERROR Reference required" << endl;
    return 0; 
  }
  if (HashListFile == "") {
    cout << "Error HashList required" << endl;
    return 0; 
  }
  if (outStub == "") {
    if (samFile != "stdin") {
      outStub = samFile; 
    } else {
      cout << "ERROR out file stub required " << endl;
      return -1; 
    }
  }
  
  for (int i = 0; i < ParentHashFilePaths.size(); i++) {
    ifstream reader; 
    reader.open (argv[ParentHashFilePaths[i]]); 
    string line = "";
    unordered_map <unsigned long int, int> hl;  
    while (getline(reader, line)) {
      vector <string> temp = Util::Split(line, ' '); 
      unsigned long hash = Util::HashToLong(temp[0]); 
      hl[hash] = atoi(temp[1].c_str());
      hash = Util::HashToLong(Util::RevComp(temp[0])); 
      hl[hash] = atoi(temp[1].c_str()); 
    }
    ParentHashes.push_back(hl); 
    reader.close(); 
  }

  for (int i = 0; i < ParentHashFilePathsReference.size(); i++) {
    ifstream reader;
    reader.open (argv[ParentHashFilePathsReference[i]]);
    string line = "";
    unordered_map <unsigned long int, int> hl;
    while (getline(reader, line)) {
      vector <string> temp = Util::Split(line, ' ');
      unsigned long hash = Util::HashToLong(temp[0]);
      ParentHashes[i][hash] = atoi(temp[1].c_str());
      hash = Util::HashToLong(Util::RevComp(temp[0]));
      ParentHashes[i][hash] = atoi(temp[1].c_str());
    }
    reader.close();
  }
  
  for(int i =0; i < ParentHashes.size(); i++) {
    cout << "sample " << i << endl;
  }
  
  ifstream reader;
  reader.open (MutHashFilePath);
  string line = "";
  
  while (getline(reader, line)) {
    vector <string> temp = Util::Split(line, ' ');
    unsigned long hash = Util::HashToLong(temp[0]);
    MutantHashes[hash] = atoi(temp[1].c_str());
    hash = Util::HashToLong(Util::RevComp(temp[0]));
    MutantHashes[hash] = atoi(temp[1].c_str());
  }

  reader.close(); 
  reader.open (MutHashFilePathReference);
  line = "";

  while (getline(reader, line)) {
    vector <string> temp = Util::Split(line, ' ');
    unsigned long hash = Util::HashToLong(temp[0]);
    MutantHashes[hash] = atoi(temp[1].c_str());
    hash = Util::HashToLong(Util::RevComp(temp[0]));
    MutantHashes[hash] = atoi(temp[1].c_str());
  }
  
  reader.close();
  double vm, rss, MAXvm, MAXrss;
  MAXvm = 0;
  MAXrss = 0;
  Util::process_mem_usage(vm, rss, MAXvm, MAXrss);
  cout << "VM: " << vm << "; RSS: " << rss << endl;
  int BufferSize = 1000;
  Reff.open(RefFile);
  ifstream ModelFile; 
  ModelFile.open (ModelFilePath); 

  if (ModelFile.is_open()) {
    cout << "ModelFile is open";
  } else {
    cout << "Error no model file given" << endl;
    return -1;
  }
  
  ifstream HashList;
  HashList.open (HashListFile);

  if ( HashList.is_open()) {
    cout << "HashList Open " << HashListFile << endl;
  } else {
    cout << "Error, HashList could not be opened";
    return -1;
  }
  
  line = "";
  getline(HashList, line);
  cout << "line = " << line << endl; 
  char seperator = '\t';
  vector<string> temp = Util::Split(line, seperator);
  if (temp.size() ==1) {
    cout << "separator is not tab" << endl; 
    seperator = ' '; 
    temp = Util::Split(line, seperator);
  } else {
    cout << "separator is tab" << endl;
  }
  
  cout << "split = " << temp[0] << " and " << temp[1] << endl;
  HashSize = temp[0].size(); 
  if (temp.size() == 4) { 
    HashSize = temp[3].length(); 
    Hash.insert(pair<string, int>(temp[3], atoi(temp[2].c_str())));
    
    while (getline(HashList, line)) {
      vector<string> temp = Util::Split(line, seperator);
      Hash.insert(pair<string, int>(temp[3], atoi(temp[2].c_str())));
      Hash.insert(pair<string, int>(Util::RevComp(temp[3]), atoi(temp[2].c_str()))); 
    }
    HashList.close(); 
    cout << "done with HashList" << endl;
  }
  else if (temp.size() == 2) {
    HashSize = temp[0].length();
    Hash.insert(pair<string, int>(temp[0], atoi(temp[1].c_str())));
    Hash.insert(pair<string, int>(Util::RevComp(temp[0]), atoi(temp[1].c_str())));
    while ( getline(HashList, line)) {
      vector<string> temp = Util::Split(line, seperator);
      Hash.insert(pair<string, int>(temp[0], atoi(temp[1].c_str())));
    }
    HashList.close();
    cout << "done with HashList" << endl;
  }
 
  ifstream SamFile;
  if (samFile == "stdin") {
    cout << "Sam File is STDIN" << endl;
    SamFile.open ("/dev/stdin");
  } else {
    cout << "Sam File is " << samFile << endl;
    SamFile.open (samFile);
  }
  if (SamFile.is_open()) {
    cout << "Sam File Opend\n";
  } else {
    cout << "Error, SamFile could not be opened";
    return 0;
  }
  
  string boom = outStub;
  VCFOutFile.open(boom+ ".vcf"); 
  BEDOutFile.open(boom+ ".vcf.bed"); 
  BEDBigStuff.open(boom+ ".vcf.Big.bed");
  BEDNotHandled.open(boom+ ".vcf.NotHandled.bed");
  Invertions.open(boom+".vcf.invertions.bed");
  Translocations.open(boom+ ".vcf.Translocations");
  Translocationsbed.open(boom+ ".vcf.Translocations.bed");
  Unaligned.open(boom+"vcf.Unaligned");
  
  //write VCF header
  VCFOutFile << "##fileformat=VCFv4.1" << endl;
  VCFOutFile << "##fileDate=" << time(0) << endl;
  VCFOutFile << "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">" << endl;
  VCFOutFile << "##FORMAT=<ID=AK,Number=1,Type=Integer,Description=\"Alternte Kmer Count\">" << endl;
  VCFOutFile << "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Total Kmer depth across the variant\">" << endl;
  VCFOutFile << "##FORMAT=<ID=RO,Number=1,Type=Integer,Description=\"Mode of reference kmer counts\">" << endl;
  VCFOutFile << "##FORMAT=<ID=AO,Number=1,Type=Integer,Description=\"Mode of alt kmer counts\">" << endl;
  VCFOutFile << "##FORMAT=<ID=LP,Number=1,Type=Integer,Description=\"Number of lowcoverage parent bases\">" << endl;
  VCFOutFile << "##FORMAT=<ID=PC,Number=1,Type=Integer,Description=\"Mode of parents coverage\">" << endl;
  VCFOutFile << "##FORMAT=<ID=SB,Number=1,Type=Float,Description=\"StrandBias\">" << endl;
  VCFOutFile << "##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of SV detected\">" << endl;
  VCFOutFile << "##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=\"Length of SV detected\">" << endl; 
  VCFOutFile << "##INFO=<ID=END,Number=1,Type=Integer,Description=\"END of SV detected\">" << endl; 
  VCFOutFile << "##INFO=<ID=AO,Number=1,Type=Integer,Description=\"Alternate allele observations, with partial observations recorded fractionally\">"<<endl;
  VCFOutFile << "##INFO=<ID=HD,Number=.,Type=String,Description=\"Hash counts for each k-mer overlapping the vareint, -1 indicates no info\">"<< endl;
  VCFOutFile << "##INFO=<ID=RN,Number=1,Type=String,Description=\"Name of contig that produced the call\">"<< endl;
  VCFOutFile << "##INFO=<ID=MQ,Number=1,Type=Integer,Description=\"Mapping quality of the contig that created the call\">"<< endl;
  VCFOutFile << "##INFO=<ID=cigar,Number=1,Type=String,Description=\"Cigar string for the contig that created the call\">"<< endl;
  VCFOutFile << "##INFO=<ID=VT,Number=1,Type=String,Description=\"Varient Type\">"<< endl;
  VCFOutFile << "##INFO=<ID=CVT,Number=1,Type=String,Description=\"Compressed Varient Type\">"<< endl;
  VCFOutFile << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t";
  VCFOutFile << outStub << endl;
  
  int lines = 0;
  line="";
  unsigned long LongHash;
  cout << "Reading in Sam File"  << endl;
  map <string, int> Names; 
  vector <SamRead> reads; 
  int counter = 0; 
  
  while (getline(SamFile, line)) {
    if (line.c_str()[0] == '@') {
      cout << " HEADER LINE = " << line << endl; 
      vector <string> temp = Util::Split(line, '\t');
      cout << temp[0] << endl;
      if (temp[0] == "@SQ") {
	cout << temp[1] << endl;
	vector <string> chr = Util::Split(temp[1], ':');
	vector <string> len = Util::Split(temp[2], ':');
	
	cout << "##contig=<ID=" <<  chr[1]<<",length=" << len[1] << ">"<< endl;
	VCFOutFile <<"##contig=<ID=" <<  chr[1]<<",length=" << len[1] << ">" << endl;
      }
    } else {
      cout << line << endl;
      counter ++; 
     SamRead read; 
     cout << "parse " << endl; 
     read.parse(line);

     if (read.FlagBits[2] != 1)  {
       cout << "RefSeq" << endl;
       read.getRefSeq();
       cout << "peak" << endl;
       read.createPeakMap();
       reads.push_back(read);
       if (counter%100 == 0) {
	 cout << "read " << counter << " entries " << char(13); 
       }
     }
    }
  }
  cout << endl;
  cout << "Read in " << reads.size() << " reads " << endl;
  cout << "procesing split reads" << endl;

  for (int i = 0; i < reads.size(); i++) {
    cout << "processing read " << reads[i].name << endl;
    if (reads[i].alignments.size() == 0) {
      reads[i].alignments.push_back(i);
      int count = 0; 
      for (int j = i+1; j < reads.size(); j++) {
	count++; 
	if (strcmp(reads[i].name.c_str(), reads[j].name.c_str()) == 0) {
	  cout << "found mate " << reads[j].name << endl;
	  reads[i].alignments.push_back(j); 
	  reads[j].first = false;
	}
	if (count > 100000) {
	  break; 
	}
      }
    }
    
    for (int j = 0; j<reads[i].alignments.size(); j++) {
      reads[reads[i].alignments[j]].alignments = reads[i].alignments;
    } 
  }
  
  for (int i = 0; i < reads.size(); i++) {
    SamRead read = reads[i]; 
    cout << read.name  << endl;
    if (strcmp(read.chr.c_str(), "*") == 0) {
      read.writetofile(Unaligned);
      Unaligned << endl;
    }
    //if  (read.first)
     {
       //if it looks like a simple split read alignment colaps it into a single read
       if (read.first and read.alignments.size() == 2 ) {
	 cout << "atempting colaps" << endl;
	 cout << read.name << endl;
	 vector<SamRead> R; 
	 for(int j =0; j< read.alignments.size(); j++) {
	   //if (reads[read.alignments[j]].chr == read.chr)
	   {
	     R.push_back(reads[read.alignments[j]]);
	     cout << reads[read.alignments[j]].name << endl;
	   }
	 }
	 if (R.size() == 2 & /*R[0].chr == R[1].chr & */ R[0].mapQual > 0 & R[1].mapQual > 0) {
	   read = BetterWay(R); 
	 }
       }
       else if(read.first and read.alignments.size() > 2) {
	 BEDNotHandled << "too many alignments" << endl;
	 BEDNotHandled << read.chr << "\t" << read.pos << "\t" << read.pos+read.seq.size() << "\t" << read.name << "\t" << read.cigar << endl;
	 cout << "too many alignments" << endl;
	 cout << read.chr << "\t" << read.pos << "\t" << read.pos+read.seq.size() << "\t" << read.name << "\t" << read.cigar << endl;
	 for (int j = 0; j< read.alignments.size(); j++) {
	   SamRead mate = reads[read.alignments[j]]; 
	   cout << j << "\t" << mate.chr << "\t" << mate.pos << "\t" << mate.pos+mate.seq.size() << "\t" << mate.name << "\t" << mate.cigar << endl;
	   BEDNotHandled << mate.chr << "\t" << mate.pos << "\t" << mate.pos+mate.seq.size() << "\t" << mate.name << "\t" << mate.cigar << endl;
	   mate.writetofile(BEDNotHandled);
	 }
	 BEDNotHandled << endl << endl;
       }
       
       if (read.mapQual > MinMapQual and read.alignments.size() <= 2) {
	 read.parseMutations(argv); 
       }
     }
  }
  VCFOutFile.close(); 
  BEDOutFile.close();
  BEDBigStuff.close();
  BEDNotHandled.close();
  Invertions.close(); 
  cout << "\nDone with RUFUS.interpret\n";
  return 0; 
}
