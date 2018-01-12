#ifndef __SOURCE_UTIL_H__
#define __SOURCE_UTIL_H__

#include <vector>
#include <string>

struct classcomp{
  bool operator() (const char& lhs, const char& rhs) const
  {return lhs<rhs;}
};

class Util{
 public: 
  static bool fncomp(char, char);
  static const std::vector<std::string> Split(const std::string&, const char);
  static std::string trim(std::string);
  static unsigned long HashToLong(std::string);
  static unsigned long HashToLongTB(std::string, std::string);
  static std::string LongToHash(unsigned long, int);
  static std::string RevComp(std::string);
  static std::string RevQual(std::string);
  static void process_mem_usage(double&, double&, double&, double&);
};

#endif // __SOURCE_UTIL_H__
