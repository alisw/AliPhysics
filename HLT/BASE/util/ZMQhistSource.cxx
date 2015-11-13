#include "zmq.h"
#include <iostream>
#include "AliHLTDataTypes.h"
#include "AliHLTComponent.h"
#include "AliHLTMessage.h"
#include "TClass.h"
#include "TMap.h"
#include "TPRegexp.h"
#include "TObjString.h"
#include "TList.h"
#include "TMessage.h"
#include "TH1F.h"
#include "TF1.h"
#include <cmath>
#include <time.h>
#include <string>
#include <map>
#include <unistd.h>
#include "AliZMQhelpers.h"

void* fZMQout = NULL;
void* fZMQcontext = NULL;
TString fZMQconfigOUT = "PUSH>tcp://localhost:60210";
int fZMQsocketModeOUT = -1;

float fSleep = 1e6; //in microseconds
int fCount = 0;
int fNentries = 1;

TH1F* fHistogram;
TString fHistName = "histogram";
TString fHistDistribution = "exp(-0.5*((x-0.)/0.1)**2)";
float fHistRangeLow = -0.5;
float fHistRangeHigh = 0.5;
int fHistNBins = 100;
    
TString GetFullArgString(int argc, char** argv);
int ProcessOptionString(TString arguments);
stringMap* TokenizeOptionString(const TString str);
int ProcessOption(TString option, TString value);

//_______________________________________________________________________________________
int main(int argc, char** argv)
{
  int mainReturnCode=0;
  int rc = 0;

  //process args
  if (ProcessOptionString(GetFullArgString(argc,argv)) != 0) return 1;
  
  //ZMQ init
  fZMQcontext = zmq_ctx_new();
  fZMQsocketModeOUT = alizmq_socket_init(fZMQout, fZMQcontext, fZMQconfigOUT.Data(), -1);
  if (fZMQsocketModeOUT < 0) return 1;

  TF1 formula("histDistribution",fHistDistribution);
  fHistogram = new TH1F(fHistName, fHistName, fHistNBins, fHistRangeLow, fHistRangeHigh);

  //main loop
  int iterations=0;
  while(fCount==0 || iterations++<fCount)
  {
    fHistogram->Reset();
    fHistogram->FillRandom("histDistribution",fNentries);
    
    AliHLTDataTopic topic = kAliHLTDataTypeTObject;
    
    if (fZMQsocketModeOUT==ZMQ_REP)
    {
      //receive the request - could be multipart
      aliZMQmsgStr request;
      rc = alizmq_msg_recv(&request, fZMQout, 0);
    }
    
    rc = alizmq_msg_send(topic, fHistogram, fZMQout, 0, 0);
    if (rc<0) printf("unable to send\n");

    unsigned int microseconds;
    microseconds = fSleep;
    if (fZMQsocketModeOUT!=ZMQ_REP) usleep(microseconds);
  }//main loop

  zmq_close(fZMQout);
  zmq_ctx_destroy(fZMQcontext);
  return mainReturnCode;
}

//______________________________________________________________________________
int ProcessOption(TString option, TString value)
{
  //process option
  //to be implemented by the user
  
  //if (option.EqualTo("ZMQpollIn"))
  //{
  //  fZMQpollIn = (value.EqualTo("0"))?kFALSE:kTRUE;
  //}
 
  if (option.EqualTo("name")) 
  {
    fHistName = value;
  }
  else if (option.EqualTo("out"))
  {
    fZMQconfigOUT = value;
  }
  else if (option.EqualTo("sleep"))
  {
    fSleep = round(value.Atof()*1e6);
  }
  else if (option.EqualTo("distribution"))
  {
    fHistDistribution = value;
  }
  else if (option.EqualTo("range"))
  {
    TString lowString = value(0,value.Index(','));
    TString highString = value(value.Index(',')+1,999);
    fHistRangeLow = lowString.Atof();
    fHistRangeHigh = highString.Atof();
  }
  else if (option.EqualTo("nbins"))
  {
    fHistNBins = value.Atoi();
  }
  else if (option.EqualTo("count"))
  {
    fCount = value.Atoi();
  }
  else if (option.EqualTo("entries"))
  {
    fNentries = value.Atoi();
  }
  return 0; 
}

////////////////////////////////////////////////////////////////////////////////
//_______________________________________________________________________________________
TString GetFullArgString(int argc, char** argv)
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
int ProcessOptionString(TString arguments)
{
  //process passed options
  int rc = 0;
  stringMap* options = TokenizeOptionString(arguments);
  for (stringMap::iterator i=options->begin(); i!=options->end(); ++i)
  {
    Printf("  %s : %s", i->first.data(), i->second.data());
    rc = ProcessOption(i->first,i->second);
    if (rc != 0) break;
  }
  delete options; //tidy up

  return rc; 
}

//______________________________________________________________________________
stringMap* TokenizeOptionString(const TString str)
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

  stringMap* options = new stringMap;

  TArrayI pos;
  const TString mods="";
  Int_t start = 0;
  while (1) {
    Int_t prevStart=start;
    TString optionStr="";
    TString valueStr="";
    
    //check if we have a new option in this field
    Int_t nOption=optionRE.Match(str,mods,start,10,&pos);
    if (nOption>0)
    {
      optionStr = str(pos[6],pos[7]-pos[6]);
      optionStr=optionStr.Strip(TString::kBoth,'\'');
      start=pos[1]; //update the current character to the end of match
    }

    //check if the next field is a value
    Int_t nValue=valueRE.Match(str,mods,start,10,&pos);
    if (nValue>0)
    {
      valueStr = str(pos[0],pos[1]-pos[0]);
      valueStr=valueStr.Strip(TString::kBoth,'\'');
      start=pos[1]; //update the current character to the end of match
    }
    
    //skip empty entries
    if (nOption>0 || nValue>0)
    {
      (*options)[optionStr.Data()] = valueStr.Data();
    }
    
    if (start>=str.Length()-1 || start==prevStart ) break;
  }

  //for (stringMap::iterator i=options->begin(); i!=options->end(); ++i)
  //{
  //  printf("%s : %s\n", i->first.data(), i->second.data());
  //}
  return options;
}

