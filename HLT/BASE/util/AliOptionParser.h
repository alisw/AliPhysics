#ifndef ALIOPTIONPARSER_H
#define ALIOPTIONPARSER_H

//blame: Mikolaj Krzewicki, mkrzewic@cern.ch
//simple option parser class

#include <vector>
#include <string>

#include "TString.h"

class AliOptionParser {
public:
  typedef std::vector<std::pair<std::string, std::string> > aliStringVec;
  AliOptionParser() {}
  virtual ~AliOptionParser() {}
  //implement this to process one option at a time
  virtual int ProcessOption(TString /*option*/, TString /*value*/) {return 0;}

  //call this to parse the args
  int ProcessOptionString(TString arguments);
  int ProcessOptionString(int argc, char** argv) { return ProcessOptionString(GetFullArgString(argc,argv)); }

  //convert argc/argv into a TString of options
  static TString GetFullArgString(int argc, char** argv);
  static aliStringVec* TokenizeOptionString(const TString str);
};

#endif
