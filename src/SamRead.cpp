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
#include "SamUtil.h"
#include "Util.h"

using namespace std;

extern vector <unordered_map <unsigned long int, int> > ParentHashes;
extern unordered_map  <unsigned long int, int> MutantHashes;
extern FastaReference Reff;
extern int HashSize;
extern int MaxVarentSize;
extern ofstream VCFOutFile;
extern  ofstream BEDOutFile;
extern ofstream BEDBigStuff;
extern ofstream BEDNotHandled;
extern map <string, int> Hash;

int SamRead::CheckParentCov(int &mode) {
  bool good = true;
  int lowC = 0;
  vector<int> cov;

  for (int pi = 0; pi < RefRefCounts.size(); pi++){
    for (int i = 0; i < RefRefCounts[pi].size(); i++){
      if (RefKmers[i] != ""){
        int ParRef = 0;
        int ParAlt = 0;
        if (RefAltCounts[pi][i] > 0) {
          ParAlt = RefAltCounts[pi][i];
	}
        if (RefRefCounts[pi][i] > 0) {
          ParRef = RefRefCounts[pi][i];
	}
        cov.push_back(ParRef+ParAlt);
        if (ParRef+ParAlt > 0 && ParRef+ParAlt < 10) {
          lowC++;
	}
      }
    }
  }
  
  if(cov.size()>1){
    sort (cov.begin(), cov.end());
    mode = cov[cov.size()/2];
  } else {
    mode = -1;
  }

  return lowC;
}

void SamRead::flipRead() {
  cout <<"FLIPPING reads not on the same strand";
  write();
  string FlipSeq = ""  ;
  string FlipQual = ""   ;
  string FlipRefSeq = "";
  string FlipCigarString = "";
  string FlipStrand = "";
  vector<bool> FlipPeakMap ;
  vector<int> FlipPos;
  vector<string> FlipChrPos;

  for (int i = seq.size() -1;  i >=0; i--) {
      FlipQual += qual.c_str()[i];
      FlipCigarString += cigarString.c_str()[i];
      FlipStrand += '-';
      FlipPos.push_back(Positions[i]);
      FlipChrPos.push_back(ChrPositions[i]);
      FlipPeakMap.push_back(PeakMap[i]);
    }
  
  FlipSeq = Util::RevComp(seq);
  FlipRefSeq = Util::RevComp(RefSeq);
  seq = FlipSeq;
  qual = FlipQual;
  RefSeq = FlipRefSeq;
  cigarString = FlipCigarString;
  strand = FlipStrand;
  Positions = FlipPos;
  ChrPositions = FlipChrPos;
  PeakMap=FlipPeakMap;
  write();
}

void SamRead::write() {
  cout << name << endl;
  cout << "   flag = " << flag << endl;
  cout << "   mapQual = " << mapQual << endl;
  cout << "   Strand = " << SamUtil::GetReadOrientation(flag) << endl;
  cout << "   Alignments = " << alignments.size() << endl;
  cout << "   chr " << chr << " - " << pos << " qual = " << mapQual <<  endl;
  cout << "   cigar  = " << cigar << endl;
  cout << "   Seq    = " << seq << endl;
  cout << "   Qual   = " << qual << endl;
  cout << "   Cigar  = " << cigarString << endl;
  cout << "   RefSeq = " << RefSeq << endl;
  cout << "   strand = " << strand << endl;
  cout << "  PeakMap = ";

  for (int i =0; i < PeakMap.size(); i++) {
    cout << PeakMap[i]; 
  }
  cout << endl;
  cout << "   RefPositions: ";

  for (int i =0; i < Positions.size(); i++) {
    cout << Positions[i] << " \t";
  }

  cout << endl;
  cout << "   RefChromoso: ";

  for (int i =0; i < ChrPositions.size(); i++) {
    cout << ChrPositions[i] << " \t";
  }
  cout << endl;
}

void SamRead::writetofile(ofstream &out) {
  out << name << endl;
  out << "   flag = " << flag << endl;
  out << "   mapQual = " << mapQual << endl;
  out << "   Strand = " << SamUtil::GetReadOrientation(flag) << endl;
  out << "   Alignments = " << alignments.size() << endl;
  out << "   chr " << chr << " - " << pos << " qual = " << mapQual <<  endl;
  out << "   cigar  = " << cigar << endl;
  out << "   Seq    = " << seq << endl;
  out << "   Qual   = " << qual << endl;
  out << "   Cigar  = " << cigarString << endl;
  out << "   RefSeq = " << RefSeq << endl;
  out << "   PeakMap= ";

  for (int i =0; i < PeakMap.size(); i++) {
    out <<  PeakMap[i] ;
  }
  out << endl;
  out << "   PMSize = " << PeakMap.size() << endl;
}

void SamRead::writeVertical() {
  cout << "ParentHashes size = " <<  ParentHashes.size() << "RefAltCounts size " 
       << RefAltCounts.size() << endl;
  cout << name << endl;
  cout << "   flag = " << flag << endl;
  cout << "   chr " << chr << " - " << pos << " qual = " << mapQual <<  endl;
  cout << "   cigar  = " << cigar << endl;
  cout << "   Alignments = " << alignments.size() << endl;

  for (int i =0; i < seq.size(); i++) {
    cout << seq.c_str()[i] << "\t" << qual.c_str()[i] << "\t" << cigarString.c_str()[i]
	 << "\t" << RefSeq.c_str()[i] << "\t" << ChrPositions[i] << "\t" << Positions[i] 
	 << "\t" << MutAltCounts[i] << "\t" << MutRefCounts[i] << "\t" 
	 << MutHashListCounts[i] << "\t" ;
    cout << "\tParents";

    for (int pi=0; pi < RefAltCounts.size(); pi++){
      cout << "\t" << RefAltCounts[pi][i] << "\t" << RefRefCounts[pi][i];
    }

    cout << "\t" << RefKmers[i] << "\t" << AltKmers[i];
    cout<< endl;
  }
}

void SamRead::createPeakMap() {
  vector<bool> tempPeakMap;
  int i =0;
  int max = -1;
  int maxSpot = -1;
  int last = 0;

  for (int i =0; i< qual.size(); i++) {
    if (qual[i] <='!') {
      tempPeakMap.push_back(0);
    } else {
      int j = i;
      max = qual[j];
      
      while (j < qual.size() and qual[j] > '!') {
	if (max < qual[j]) {
	  max = qual[j];
	}
	j++;
      }
      
      cout << endl;
      j = j-1;
      
      for (int k = i; k < qual.size() and k <= j ; k++) {
	if (qual[k]==max and cigarString[k] != 'H') {
	  tempPeakMap.push_back(1);
	} else {
	  tempPeakMap.push_back(0);
	}
      }

      i = j;
    }
  }

  for (int i =0; i< qual.size(); i++) {
    if (seq[i] == '-'){
      tempPeakMap[i] == tempPeakMap[i-1];
    }
  }

  PeakMap.clear();
  PeakMap = tempPeakMap;
}

void SamRead::parseMutations( char *argv[]) {
  cout << "In Parsing Mutations " << endl;
  createPeakMap();
  write(); 
  string StructCall = ""; 
  /////////////////Building up varHash and hash lists ///////////// 
  cout << "Building up varHash" << endl; 
  vector <string> hashes;
  vector <string> hashesRef;
  vector <bool> varHash; 

  for (int i = 0; i <  seq.size() - HashSize; i++) {
    string newHash = "";
    string newHashRef = "";
    newHash += seq.c_str()[i];
    newHashRef += RefSeq.c_str()[i]; 
    int count = 0;
    if ((cigarString.c_str()[i] != 'D' and cigarString.c_str()[i] != 'R'and
	 cigarString.c_str()[i] != 'H')) {

	  for (int j = 1; j<seq.size() - i and count < HashSize-1 ; j++) {
	    if (cigarString.c_str()[i+j] != 'D' and 
		cigarString.c_str()[i+j] != 'R' and 
		cigarString.c_str()[i+j] != 'H') {
	      newHash += seq.c_str()[i+j]; 
	      newHashRef += RefSeq.c_str()[i+j];
	      count++;
	    }
	  }
    }
    hashes.push_back(newHash);
    hashesRef.push_back(newHashRef);

    if (Hash.count(newHash) > 0 or Hash.count(Util::RevComp(newHash)) > 0) {
      varHash.push_back(true); 
    } else {
      varHash.push_back(false); 
    }
  }
  ///////////////////building up parent hash counts //////////////////
  cout << "Bulding Par hash counts" << endl; 
  vector <vector<int>> parentCounts;
  vector <vector<int>> parentCountsReference;
  
  for (int pi = 0; pi<ParentHashes.size(); pi++) {
    vector <int> counts;
    vector <int> countsRef;
    for (int i = 0; i< hashes.size(); i++) {
	  string hash =  hashes[i];
	  string hashRef = hashesRef[i]; 
	  bool checkHash = true; 

	  for (int j = 0; j < HashSize; j++) {
	    if (!(hash[j] == 'A' or hash[j] == 'C' or 
		  hash[j] == 'G' or hash[j] == 'T')) {
		  checkHash = false; 
		  break; 
	    }
	  }
	  
	  if (checkHash) {
	    unsigned long int LongHash = Util::HashToLong(hash);

	    if (ParentHashes[pi].count(LongHash) >0) {
	      counts.push_back(ParentHashes[pi][LongHash]);
	    } else {
	      counts.push_back(0);
	    }
	    
	    unsigned long int LongHashRef = Util::HashToLong(hashRef);
	    if (ParentHashes[pi].count(LongHashRef) >0) {
	      countsRef.push_back(ParentHashes[pi][LongHashRef]);
	    } else {
	      countsRef.push_back(0);
	    }
	  } else {
	    counts.push_back(-1);
	    countsRef.push_back(-1);
	  } 
    }
    parentCounts.push_back(counts);
    parentCountsReference.push_back(countsRef);
  }
  /////////////////////bulid Mut counts/////////////////////////////
  cout << "bulding mut counts" << endl; 
  vector <int> mutCounts;
  vector <int> mutCountsRef;
  cout << hashes.size() << endl;
  cout << hashesRef.size() << endl; 

  for(int i = 0; i < hashes.size(); i++) {
    cout << i<< endl;
    string hash =  hashes[i];
    cout << "   hash = " << hash << endl; 
    string hashRef = hashesRef[i];
    cout << "RefHash = " << hashRef << endl; 
    bool checkHash = true;

    for (int j = 0; j < HashSize; j++) {
      if (!(hash[j] == 'A' or hash[j] == 'C' or hash[j] == 'G' or hash[j] == 'T')) {
	checkHash = false;
	break;
      }
    }
    
    cout << "check hash = " << checkHash << endl; 
    if (checkHash) {
      unsigned long int LongHash = Util::HashToLong(hash);
      
      if (MutantHashes.count(LongHash) > 0) {
	mutCounts.push_back(MutantHashes[LongHash]);
      } else {
	mutCounts.push_back(0);
      }
      
      unsigned long int LongHashRef = Util::HashToLong(hashRef);
      if (MutantHashes.count(LongHashRef) >0) {
	mutCountsRef.push_back(MutantHashes[LongHashRef]);
      } else {
	mutCountsRef.push_back(0);
      }
    } else {
      mutCounts.push_back(-1);
      mutCountsRef.push_back(-1);
    }
  }
  ////////////////////write out vertical table/////////////////////////
  cout << "writing hashes out vert" << endl;
  
  for(int i = 0; i < hashes.size(); i++) {
    cout << i << "\t" << hashes[i] << "\t" << varHash[i] << "\t" 
	 << PeakMap[i] << "\t" << (int) qual.c_str()[i]-33;
    cout << "\t" << "MutVar-" << mutCounts[i]; 
    
    for (int j = 0; j < parentCounts.size(); j++) {
      cout << "\t" <<  parentCounts[j][i];
    }
    cout << "\t" << "MutRef-" << mutCountsRef[i];
    for (int j = 0; j < parentCountsReference.size(); j++) {
      cout << "\t" <<  parentCountsReference[j][i];
    }

    cout << endl;
  }

  string reff = "";
  string alt = "";
  string varType = "";
  for(int i = 0; i < cigarString.size(); i++) {
    reff = ""; 
    alt = "";
    varType = ""; 
    if ((cigarString.c_str()[i] == 'X' or cigarString.c_str()[i] == 'I' or
	 cigarString.c_str()[i] == 'D' or cigarString.c_str()[i] == 'Y') and
	RefSeq.c_str()[i] != 'N'){
      cout << "found a " << cigarString.c_str()[i] << endl;
      int size = -1; 
      int startPos = i; 
      bool  AnyBasesOver0  = false; 
      string Denovo = "inherited"; 
      
      if (qual.c_str()[i] > '!') {
	AnyBasesOver0 = true;
      }
      if (PeakMap[i] == 1) {
	Denovo = "DeNovo";
      }
      
      for(int j = 0; j< cigarString.size() - i; j++) {
	if(cigarString.c_str()[i+j] == 'X' or cigarString.c_str()[i+j] == 'D' or
	   cigarString.c_str()[i+j] == 'I' or cigarString.c_str()[i+j] == 'Y') {
	  size = j; 
	  
	  if (qual.c_str()[i+j] > '!') {
	    AnyBasesOver0 = true;
	  }
	  if (PeakMap[i+j] == 1) {
	    Denovo = "DeNovo";
	  }
	  
	} else {
	  break;
	} 
      }

      cout << "size =" << size<< endl;
      /////////////////////Parent Low Coverage Check ////////////////// 
      bool LowCov = false; 
      int lowCount = 0; 
      int low = i - HashSize - 5; 
      if (low < 0) {
	low = 0;
      }
      cout << "checking bases " << low << " to " << i+size+5 << endl;
      for(int j = low; j < i+size+5 and j < hashes.size(); j++) {
	bool AllPar0 = true; 

	for (int k = 0; k < parentCounts.size(); k++) {
	  if (parentCounts[k][j] != 0) {
	    AllPar0 = false;
	  }
	}

	if (AllPar0) {
	  varHash[j] = true; 
	}
	
	if (hashes[j].size() == HashSize and varHash[j] == false) {
	  for (int k = 0; k < parentCounts.size(); k++) {
	    if (parentCounts[k][j] <= 5 and parentCounts[k][j] > 0) {
	      LowCov = true;
	      lowCount++; 
	    }
	  }
	}
      }
      
      //enabling this will only report varites covered by hashes 
      if (AnyBasesOver0) {
	if (cigarString.c_str()[i] == 'I' or cigarString.c_str()[i] == 'D' or
	    cigarString.c_str()[i] == 'Y') {
	  for (int k = 1; i-k >= 0; k++) {
	    if (ChrPositions[i-k] == "nope") {
	    } else {
	      reff+=RefSeq.c_str()[i-k]; 
	      alt+=seq.c_str()[i-k];
	      startPos = i-k;
	      break;
	    }
	  }  
	}
	/////////build the alleles and var type/////////
	for(int j = 0; j<= size; j++) {

	  if (RefSeq.c_str()[i+j] == 'A' or RefSeq.c_str()[i+j] == 'C' or
	      RefSeq.c_str()[i+j] == 'G' or RefSeq.c_str()[i+j] == 'T') {
	    reff+=RefSeq.c_str()[i+j]; 
	  }
	  
	  if (seq.c_str()[i+j] == 'A' or seq.c_str()[i+j] == 'C' or 
	      seq.c_str()[i+j] == 'G' or seq.c_str()[i+j] == 'T') {
	    alt+=seq.c_str()[i+j]; 
	  }

	  varType += cigarString.c_str()[i+j]; 
	}

	//***********Build up hash depth nehborhod********************
	int lower = i-HashSize; 
	if (lower < 0) {
	  lower =0;
	}
	int upper = i+alt.length()+reff.length();//-1;
	cout << i<<" + "<< alt.length() << " + " << reff.length() << endl;
	if (upper > MutRefCounts.size()) {
	  cout << "this is going to break " << upper << " > " << MutRefCounts.size() << endl; 
	  upper = MutRefCounts.size(); 
	}
	//////////////chekcing allele frequencies ///////////
	vector <int> HashCounts; 
	vector <int> HashCountsOG; 
	vector <int> varMutRefCounts;
	vector <int> varMutAltCounts; 
	vector<vector <int>> varParRefCounts;
	vector<vector <int>> varParAltCounts;  
	vector <int> temp; 

	for(int pi = 0; pi<ParentHashes.size(); pi++){
	  varParRefCounts.push_back(temp);
	  varParAltCounts.push_back(temp); 
	}

	string nonspecific = "nonspecific";
	vector <float> freqs;  
	for (int j = lower; j<upper; j++) {
	  if(MutRefCounts[j] > 0) {
	    varMutRefCounts.push_back(MutRefCounts[j]);
	  }
	  if (MutAltCounts[j]>0) {
	    varMutAltCounts.push_back(MutAltCounts[j]);
	  }
	  if (MutRefCounts[j]>0 and MutAltCounts[j]>0) {
	    cout << "Ref = " << MutRefCounts[j] << " Mut = " << MutAltCounts[j] 
		 << " sum = " << MutRefCounts[j] +  MutAltCounts[j] ;
	    
	    if (MutRefCounts[j] +  MutAltCounts[j] <= 59 and MutRefCounts[j] + 
		MutAltCounts[j] >= 18){
	      nonspecific = "DeNovo";
	      cout << "  yay DeNovo"; 
	      freqs.push_back(( MutRefCounts[j] +  MutAltCounts[j])/MutAltCounts[j]); 
	    } else if (MutRefCounts[j] +  MutAltCounts[j] < 18) {
	      cout <<"   boo too low"; 
	    }
	    
	    cout << endl;
	  } //end of for loop

	  for (int pi=0; pi < varParRefCounts.size(); pi++) {
	    if (RefRefCounts[pi][j] > 0) {
	      varParRefCounts[pi].push_back(RefRefCounts[pi][j]); 
	    }
	    if (RefAltCounts[pi][j] > 0) {
	      varParAltCounts[pi].push_back(RefAltCounts[pi][j]);
	    }
	  }
	  
	  if (Hash.count(AltKmers[j]) > 0 and Hash[AltKmers[j]] > 0) {
	    HashCountsOG.push_back(Hash[AltKmers[j]]);
	  } else if (Hash.count(Util::RevComp(AltKmers[j])) > 0 and 
		     Hash[Util::RevComp(AltKmers[j])]) {
	    HashCountsOG.push_back(Hash[Util::RevComp(AltKmers[j])]);
	  }
	  if (Hash[AltKmers[j]] > 0) {
	    HashCounts.push_back(Hash[AltKmers[j]]);
	  } else {
	    HashCounts.push_back(-1); 
	  }
	} //end of for loop
	
	float freq = 0; 
	if (freqs.size() > 0) {
	  for (int i =0; i < freqs.size(); i++) {
	    freq+=freqs[i]; 
	  }
	  freq = freq/freqs.size(); 
	}
	
	cout << "<><><><><>MutRef<><><><><><>" << endl ;
	for (int s =0; s<varMutRefCounts.size(); s++) {
	  cout << varMutRefCounts[s] << " " ;
	}

	cout << endl << "<><><><><>MutAlt<><><><><><>" << endl;
	for (int s = 0; s < varMutAltCounts.size(); s++) {
	  cout << varMutAltCounts[s] << " " ;
	}

	cout << endl;
	sort (varMutRefCounts.begin(), varMutRefCounts.end());
	sort (varMutAltCounts.begin(), varMutAltCounts.end());

	cout << "<><><><><>MutRefSorted<><><><><><>" << endl;
	for (int s =0; s<varMutRefCounts.size(); s++){
	  cout << varMutRefCounts[s] << " " ; 
	}

	cout << endl << "<><><><><>MutAltSorted<><><><><><>" << endl; 
	for (int s = 0; s < varMutAltCounts.size(); s++) {
	  cout << varMutAltCounts[s] << " " ;
	}
	cout << endl; 

	for(int pi = 0; pi<varParRefCounts.size(); pi++){
	  sort (varParRefCounts[pi].begin(), varParRefCounts[pi].end());
	  cout << "<><><><><>Ref" << pi << "<><><><><<><>" << endl;
	  for (int s =0; s<varParRefCounts[pi].size(); s++){
	    cout << varParRefCounts[pi][s] << " " ;
	  }
	  cout << endl;
	}
	
	int MutRefMode; 
	if (varMutRefCounts.size() > 1) {
	  MutRefMode = varMutRefCounts[(varMutRefCounts.size())/2];
	} else {
	  MutRefMode = -1;
	}
	
	int MutAltMode; 
	if (varMutAltCounts.size() > 2) {
	  MutAltMode= varMutAltCounts[(varMutAltCounts.size()-2)/2];
	} else {
	  MutAltMode=-1;
	}
	
	vector <int> ParModes; 
	for(int pi = 0; pi<varParRefCounts.size(); pi++){
	  if (varParRefCounts[pi].size() > 1) {
	    ParModes.push_back(varParRefCounts[pi][((
						     varParRefCounts[pi].size())/2)]); 
	  } else {
	    ParModes.push_back(-1);
	  }
	}
	//***********check that the alleles contain only valid basis**************
	bool good = true; 
	for (int j = 0; j<reff.size(); j++) {
	  if (reff.c_str()[j] != 'A' or reff.c_str()[j] != 'C' or
	      reff.c_str()[j] != 'G' or reff.c_str()[j] != 'T') {
	    good = false;
	  }
	}

	for (int j = 0; j < alt.size(); j++) {
	  if (alt.c_str()[j] != 'A' or alt.c_str()[j] != 'C' or alt.c_str()[j] != 'G' 
	      or alt.c_str()[j] != 'T'){good = false;
	  }
	}
	
	if (good = false) {
	  cout <<"ERROR in SNP detect"<< endl;
	  cout <<endl<< chr << "\t" << pos+i << "\t" << reff << "\t" << alt << endl;
	  write();
	}
		
	//*****************Prototype Genotyper*********************
	cout << "AF = " 
	     << (double) MutAltMode / ((double) MutRefMode + (double) MutAltMode) 
	     << endl;
	string Genotype; 
	if ((double) MutAltMode / ((double) MutRefMode + (double) MutAltMode)  >.9) {
	  Genotype = "1/1";
	} else if ((double) MutAltMode / ((double) MutRefMode + (double) MutAltMode)  <.1) {
	  Genotype = "0/0";
	} else {
	  Genotype = "0/1";
	}
	
	string CompressedVarType = SamUtil::compressVar(varType, Positions[startPos], StructCall); 
	cout <<  chr << "\t" << pos+i << "\t" << CompressedVarType /*"."*/ 
	     << "\t" << reff << "\t" << alt << "\t" << HashCountsOG.size() 
	     << "\t" << varType << "\t" << "." << "\t" << "." << "\t" << "." 
	     << endl;
	
	////////////////check that parents have enough coverage////////////////////
	for(int pi = 0; pi<ParentHashes.size(); pi++) {
	  if (varParRefCounts[pi].size() > 1) {
	    if (varParRefCounts[pi][(varParRefCounts[pi].size())/2] < 10) {
	      Denovo = "ParLowCov";
	    }
	    cout << "Par" << pi << " = " << varParRefCounts[pi][(varParRefCounts[pi].size())/2] << endl;
	  }
	}
	
	if (LowCov) {
	  cout << "LOW COVERAGE" << endl;
	  Denovo = "LowCov";
	  stringstream ss;
	  ss << lowCount;
	  Denovo += ss.str(); 
	} else {
	  cout << "GOOD COVERAGE" << endl;
	}
	
	cout << "startpos = " << startPos << " chrsize = " << ChrPositions.size() << endl; 
	cout << ChrPositions[startPos] << "\t" << endl;
	cout << Positions[startPos] << "\t"  << endl;
	cout << CompressedVarType <<"-"  << endl;
	cout << Denovo /*"."*/  << "\t"  << endl;
	cout << reff << "\t"  << endl;
	cout << alt << "\t"  << endl;
	cout << HashCountsOG.size() << "\t"  << endl;
	cout << "." << "\t"  << endl;
	cout << StructCall  << endl;
	cout <<"RN=" << name  << endl;
	cout << ";MQ=" << mapQual  << endl;
	cout << ";cigar=" << cigar  << endl;
	cout << ";" << "CVT=" << CompressedVarType << ";HD=" << endl;

	if (StrandBias >= 0){
	  if (StrandBias >0.9 or StrandBias < 0.1) {
		  Denovo = "StrandBias"; 
	  }
	}
	
	////////////////////////Writing var out to file/////////////////////////
	cout << ChrPositions[startPos] << "\t" <<Positions[startPos] << "\t" 
	     << CompressedVarType <<"-" <<Denovo /*"."*/  << "\t" << reff 
	     << "\t" << alt << "\t" << HashCountsOG.size() << "\t" << "." 
	     << "\t" << StructCall <<"RN=" << name << ";MQ=" << mapQual 
	     << ";cigar=" << cigar << ";" << "CVT=" << CompressedVarType << ";HD="; 
	      
	      VCFOutFile << ChrPositions[startPos] << "\t" <<Positions[startPos] 
			 << "\t" << CompressedVarType <<"-" <<Denovo /*"."*/  
			 << "\t" << reff << "\t" << alt << "\t" << HashCountsOG.size() 
			 << "\t" << "." << "\t" << StructCall <<"RN=" << name << ";MQ="
			 << mapQual << ";cigar=" << cigar << ";" << "CVT=" 
			 << CompressedVarType << ";HD="; 
	      
	      for (int j = 0; j < HashCounts.size(); j++) {
		cout << HashCounts[j] << "_";
		VCFOutFile << HashCounts[j] << "_"; 
	      }

	      if (HashCountsOG.size() > 1) {
		std::sort (HashCountsOG.begin(), HashCountsOG.end());
		cout       << ";AO=" << HashCountsOG[HashCountsOG.size()/2];
		VCFOutFile << ";AO=" << HashCountsOG[HashCountsOG.size()/2];
	      } else {
		cout        << ";AO=" << "-1";
		VCFOutFile << ";AO=" << "-1";
	      }
	      
	      cout << ";VT=" <<  varType << "\t" ;
	      VCFOutFile <<  ";VT=" <<  varType << "\t" ;
	      cout << "GT:DP:RO:AO:LP:PC:SB" << "\t" << Genotype << ":" 
		   << MutRefMode + MutAltMode << ":" << MutRefMode << ":" 
		   << MutAltMode ;
	      VCFOutFile << "GT:DP:RO:AO:LP:PC:SB" << "\t" << Genotype 
			 << ":" << MutRefMode + MutAltMode << ":" 
			 << MutRefMode << ":" << MutAltMode ; 

	      int ParentMode; 
	      int lowC = CheckParentCov(ParentMode); 
	      cout << ":" << lowC << ":" << ParentMode << ":" << StrandBias;; 
	      VCFOutFile << ":" << lowC << ":" << ParentMode << ":" << StrandBias;;
	      cout << endl; 
	      VCFOutFile << endl;
	      BEDOutFile << chr << "\t" << pos+i << "\t" <<  pos+i+size << "\t" 
			 << chr << ":" << pos+i << ":" 
			 << (int)(reff.length() - alt.length()) << ":" 
			 << HashCountsOG.size() << endl;
	      i+=size;
	      
	      cout << "\nModes\t" << ChrPositions[startPos] << "\t" 
		   << Positions[startPos] << "\t" << CompressedVarType 
		   <<"-" <<Denovo /*"."*/  << "\t" << reff << "\t" 
		   << alt <<"\t" << MutRefMode <<"\t" << MutAltMode;

	      for(int pi = 0; pi < ParentHashes.size(); pi++) {
		cout << "\t" << ParModes[pi];
	      }

	      cout << endl;
      } 
    }
  }
  cout << "Out parsingMutations" << endl;
}

void SamRead::processCigar() {
  string num = "";
  for (int i = 0; i < cigar.length(); i++) {
    if (cigar.c_str()[i] >= 48 and cigar.c_str()[i] <= 57) {
        num = num + cigar.c_str()[i];
    } else {
      int number = atoi(num.c_str());
      for(int j = 0; j < number; j++) {
	cigarString += cigar.c_str()[i];
      }
      num = "";
    }
  }
}

void SamRead::FixTandemRef() {
  cout << "FOUND TANDEM" << endl;
  write();
  writeVertical();
  string lastChr = "nope";
  int lastPos = -1;
  string NewRef = "";
  for (int i = 0; i < seq.size(); i++) {
    cout << i << " " << seq[i] << " " <<  cigarString[i] << " " << RefSeq[i] << endl;
    if (cigarString[i]=='Y') {
      NewRef+= toupper(Reff.getSubSequence(lastChr, lastPos+1, 1).c_str()[0]);
      lastPos++;
    } else {
      NewRef+=RefSeq[i];
      lastChr = ChrPositions[i];
      lastPos = Positions[i];
    }
  }
  RefSeq = NewRef;
  cout << "done tandtem fix" << endl;
  write();
}

void SamRead::getRefSeq() {
  originalSeq = seq;
  originalQual = qual;
  RefSeq = "";
  string NewSeq = "";
  string NewQual = "";
  string NewCigar = "";
  string NewStrand = "";
  vector<int> NewPositions;
  vector<string> NewChromosome;
  int InsOffset = 0;
  
  if (Reff.sequenceNameStartingWith(chr) == "") {
    cout << "ERROR chr " << chr << " not found\n";
    return;
  }
  
  // correct star position of the read to account for Hard and soft clipped
  // bases as we are counting those now
  for (int i = 0; i < cigarString.size(); i++) {
    if (cigarString.c_str()[i] == 'H') {
    } else {
      pos = pos - i;
      break;
    }
  }
  for (int i = 0; i < cigarString.size(); i++) {
    if (cigarString.c_str()[i] == 'S') {
    } else {
      pos = pos - i;
      break;
    }
  }

  int Roffset = 0;
  int Coffset = 0;
  for (int i = 0; i < cigarString.length(); i++) {
    if (cigarString.c_str()[i] == 'M') {
      RefSeq += toupper(
			Reff.getSubSequence(chr, i + pos - 1 + Roffset, 1).c_str()[0]);
      NewSeq += seq.c_str()[i - Coffset];
      NewQual += qual.c_str()[i - Coffset];
      NewPositions.push_back(pos + i - InsOffset);
      NewChromosome.push_back(chr);
      if (toupper(Reff.getSubSequence(chr, i + pos - 1 + Roffset, 1)
		  .c_str()[0]) == seq.c_str()[i - Coffset]) {
        NewCigar += 'M';
      } else {
	NewCigar += 'X';
      }

    } else if (cigarString.c_str()[i] == 'I') {
      RefSeq += '-';
      Roffset += -1;
      NewSeq += seq.c_str()[i - Coffset];
      NewQual += qual.c_str()[i - Coffset];
      NewCigar += 'I';
      InsOffset++;
      NewPositions.push_back(pos + i - InsOffset);
      NewChromosome.push_back(chr);
    } else if (cigarString.c_str()[i] == 'D') {
      NewSeq += '-';
      NewQual += ' ';
      Coffset++;
      RefSeq += toupper(
			Reff.getSubSequence(chr, i + pos - 1 + Roffset, 1).c_str()[0]);
      NewCigar += 'D';
      NewPositions.push_back(pos + i - InsOffset);
      NewChromosome.push_back(chr);
    }

    else if (cigarString.c_str()[i] == 'H') {
      RefSeq += 'H';
      NewSeq += 'H';
      NewQual += ' ';
      Coffset++;
      NewCigar += 'H';
      NewPositions.push_back(-1);
      NewChromosome.push_back("nope");
    } else if (cigarString.c_str()[i] == 'S') {
      RefSeq += '-';
      NewSeq += seq.c_str()[i - Coffset];
      NewQual += qual.c_str()[i - Coffset];
      NewCigar += 'S';
      NewPositions.push_back(pos + i -
                             InsOffset);
      NewChromosome.push_back(chr); 
    } else {
      cout << "Unexpected behavior in SamRead::getRefSeq()" << endl;
    }
  }
  seq = NewSeq;
  cigarString = NewCigar;
  qual.clear();
  char lastQ = ' ';

  for (int i = 0; i < NewQual.size(); i++) {
    if (NewQual.c_str()[i] == ' ') {
      if (NewCigar.c_str()[i] == 'D')
        qual += lastQ;  //'!';
      else
        qual += '!';
    } else {
      qual += NewQual.c_str()[i];
      lastQ = NewQual.c_str()[i];
    }
  }
  // build up the string that will indidcate if the read has to be flipped
  for (int i = 0; i < qual.size(); i++) {
    strand += "+";
  }
  Positions = NewPositions;
  ChrPositions = NewChromosome;
  //************************Lookup Kmer counts ***************************//
  LookUpKmers();
  vector<long> blank;
  cout << "After getRefSeq";
  writeVertical();
}

void SamRead::LookUpKmers() {
  cout << "SeqSize = " << seq.size() << " RefSize = " << RefSeq.size() << endl;
  vector<long> blank;
  RefAltCounts.clear();
  RefRefCounts.clear();
  MutHashListCounts.clear();
  MutAltCounts.clear();
  MutRefCounts.clear();
  RefKmers.clear();
  AltKmers.clear();

  for (int pi = 0; pi < ParentHashes.size(); pi++) {
    RefAltCounts.push_back(blank);
    RefRefCounts.push_back(blank);
  }

  for (int j = 0; j < seq.size(); j++) {
    string hash = SamUtil::getHash(seq, j, HashSize);
    string Refhash = SamUtil::getHash(RefSeq, j, HashSize);
    RefKmers.push_back(Refhash);
    AltKmers.push_back(hash);
    if (hash != "") {
      if (hash == Refhash) {
        MutAltCounts.push_back(0);
        for (int pi = 0; pi < ParentHashes.size(); pi++) {
          RefAltCounts[pi].push_back(0);
        }
      } else {
        if (MutantHashes.count(Util::HashToLong(hash)) > 0) {
          MutAltCounts.push_back(MutantHashes[Util::HashToLong(hash)]);
        } else {
          MutAltCounts.push_back(-1);
	}
	for (int pi = 0; pi < ParentHashes.size(); pi++) {
          if (ParentHashes[pi].count(Util::HashToLong(hash)) > 0) {
            RefAltCounts[pi].push_back(ParentHashes[pi][Util::HashToLong(hash)]);
          } else {
            RefAltCounts[pi].push_back(-1);
	  }
        }
      }

      if (Hash.count(hash) > 0) {
        MutHashListCounts.push_back(Hash[hash]);
      } else {
        MutHashListCounts.push_back(-1);
      }

    } else {
      MutAltCounts.push_back(-3);
      MutHashListCounts.push_back(-3);
      for (int pi = 0; pi < ParentHashes.size(); pi++) {
        RefAltCounts[pi].push_back(-3);
      }
    }

    if (Refhash != "") {
      if (MutantHashes.count(Util::HashToLong(Refhash)) > 0) {
        MutRefCounts.push_back(MutantHashes[Util::HashToLong(Refhash)]);
      } else { 
        MutRefCounts.push_back(-1);
      }

      for (int pi = 0; pi < ParentHashes.size(); pi++) {
        if (ParentHashes[pi].count(Util::HashToLong(Refhash)) > 0) {
          RefRefCounts[pi].push_back(ParentHashes[pi][Util::HashToLong(Refhash)]);
        } else {
          RefRefCounts[pi].push_back(-1);
	}
      }

    } else {
      MutRefCounts.push_back(-3);

      for (int pi = 0; pi < ParentHashes.size(); pi++) {
        RefRefCounts[pi].push_back(-3);
      }
    }
  }
}

void SamRead::parse(string read) {
  cout << "parsing " << read << endl;
  vector<string> temp = Util::Split(read, '\t');
  cout << "boom" << endl;
  name = temp[0];
  cout << "name " << name;
  flag = atoi(temp[1].c_str());
  chr = temp[2];
  pos = atoi(temp[3].c_str());
  mapQual = atoi(temp[4].c_str());
  cigar = temp[5];
  seq = temp[9];

  if (temp[10] == "*") {
    cout << "correcting missing qualtiy" << endl;
    string newQual = "";
    for (int i = 0; i < seq.size(); i++) {
      newQual += '5';
    }
    temp[10] = newQual;
  }

  qual = temp[10];
  alignments.clear();
  first = true;
  combined = false;
  string NewSeq = "";

  for (int i = 0; i < seq.size(); i++) {
    NewSeq += toupper(seq.c_str()[i]);
  }

  seq = NewSeq;
  processCigar();
  cout << "working on strand bias" << endl;
  cout << temp[0] << endl;
  vector<string> temp2 = Util::Split(name, ':');

  if (temp2.size() >= 2) {
    cout << "break it down " << endl;
    strands = temp2[1];
    cout << "strands = " << strands << endl;
    forward = 0;
    reverse = 0;
    forward = atoi(temp2[1].c_str());
    reverse = atoi(temp2[2].c_str());
    cout << "forward = " << forward << endl;
    cout << "reverse = " << reverse << endl;
    StrandBias = ((float)forward) / ((float)forward + (float)reverse);
    cout << "strand bias = " << StrandBias << endl;
  } else {
    cout << "no strand data " << temp2.size() << endl;
    strands = "";
    StrandBias = -1;
    forward = -1;
    reverse = -1;
  }
  cout << "getting flag bits " << endl;

  for (int j = 0; j < 16; ++j) {
    FlagBits[j] = 0 != (flag & (1 << j));
  }

  cout << "Read Pared = " << FlagBits[0] << endl;
  cout << "read mapped in proper pair = " << FlagBits[1] << endl;
  cout << "read unmapped = " << FlagBits[2] << endl;
  cout << "mate unmapped = " << FlagBits[3] << endl;
  cout << "read reverse strand = " << FlagBits[4] << endl;
  cout << "mate referse strand = " << FlagBits[5] << endl;
  cout << "first in pair =" << FlagBits[6] << endl;
  cout << "second in pair =" << FlagBits[7] << endl;
  cout << "not primary alignment =" << FlagBits[8] << endl;
  cout << "read fails platform or vendor quality checks =" << FlagBits[9]
       << endl;
  cout << "read is PCR or optical duplicate =" << FlagBits[10] << endl;
  cout << "supplementary alignment =" << FlagBits[11] << endl;
}
