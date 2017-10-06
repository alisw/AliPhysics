//blame: Mikolaj Krzewicki, mkrzewic@cern.ch
//simple option parser class

#include "AliOptionParser.h"
#include "TString.h"
#include "TObjString.h"
#include "TPRegexp.h"
#include "TObjArray.h"

//_______________________________________________________________________________________
TString AliOptionParser::GetFullArgString(int argc, char** argv)
{
  TString argString;
  TString argument="";
  if (argc>0) {
    for (int i=1; i<argc; i++) {
      argument=argv[i];
      if (argument.IsNull()) continue;
      if (!argString.IsNull()) argString+=" ";
      argString+=argument;
    }
  }
  return argString;
}

//______________________________________________________________________________
int AliOptionParser::ProcessOptionString(TString arguments)
{
  //process passed options, return number of processed valid options
  aliStringVec* options = TokenizeOptionString(arguments);
  int nOptions=0;
  for (aliStringVec::iterator i=options->begin(); i!=options->end(); ++i)
  {
    //printf("  %s : %s\n", i->first.data(), i->second.data());
    if (ProcessOption(i->first,i->second)<0)
    {
      nOptions=-1;
      break;
    }
    nOptions++;
  }
  delete options; //tidy up

  return nOptions;
}

//______________________________________________________________________________
AliOptionParser::aliStringVec* AliOptionParser::TokenizeOptionString(const TString strIn)
{
  //options have the form:
  // -option value
  // -option=value
  // -option
  // --option value
  // --option=value
  // --option
  // option=value
  // option value
  // (value can also be a string like 'some string')
  //
  // options can be separated by ' ' arbitrarily combined, e.g:
  //"-option option1=value1 --option2 value2, -option4=\'some string\'"

  //optionRE by construction contains a pure option name as 3rd submatch (without --,-, =)
  //valueRE does NOT match options
  TPRegexp optionRE("(?:(-{1,2})|((?='?[^=]+=?)))"
                    "((?(2)(?:(?(?=')'(?:[^'\\\\]++|\\.)*+'|[^ =]+))(?==?))"
                    "(?(1)[^ =]+(?=[= $])))");
  TPRegexp valueRE("(?(?!(-{1,2}|[^ =]+=))"
                   "(?(?=')'(?:[^'\\\\]++|\\.)*+'"
                   "|[^ =]+))");

  aliStringVec* options = new aliStringVec;

  //first split in lines (by newline) and ignore comments
  TObjArray* lines = strIn.Tokenize("\n\r");
  TIter nextLine(lines);
  while (TObjString* objString = (TObjString*)nextLine())
  {
  TString line = objString->String();
  if (line.BeginsWith("#")) continue;
  if (line.BeginsWith("//")) continue;
  TArrayI pos;
  const TString mods="";
  Int_t start = 0;
  while (1) {
    Int_t prevStart=start;
    TString optionStr="";
    TString valueStr="";

    //check if we have a new option in this field
    Int_t nOption=optionRE.Match(line,mods,start,10,&pos);
    if (nOption>0)
    {
      optionStr = line(pos[6],pos[7]-pos[6]);
      optionStr=optionStr.Strip(TString::kBoth,'\n');
      optionStr=optionStr.Strip(TString::kBoth,'\'');
      optionStr=optionStr.Strip(TString::kLeading,'-');
      start=pos[1]; //update the current character to the end of match
    }

    //check if the next field is a value
    Int_t nValue=valueRE.Match(line,mods,start,10,&pos);
    if (nValue>0)
    {
      valueStr = line(pos[0],pos[1]-pos[0]);
      valueStr=valueStr.Strip(TString::kBoth,'\n');
      valueStr=valueStr.Strip(TString::kBoth,'\'');
      start=pos[1]; //update the current character to the end of match
    }

    //skip empty entries
    if (nOption>0 || nValue>0)
    {
      options->push_back(std::make_pair(optionStr.Data(),valueStr.Data()));
    }

    if (start>=line.Length()-1 || start==prevStart ) break;
  }

  }//while(nextLine())
  lines->Delete();
  delete lines;

  return options;
}

