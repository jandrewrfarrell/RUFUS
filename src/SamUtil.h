#ifndef __SOURCE_SAMUTIL_H__
#define __SOURCE_SAMUTIL_H__

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

using namespace std;

class SamUtil{
 public:
  static int GetReadOrientation(int);
  static string getHash(string, int, int);
  static string compressVar(string, int, string&);
  
 private:

};

#endif //__SOURCE_SAMUTIL_H__
