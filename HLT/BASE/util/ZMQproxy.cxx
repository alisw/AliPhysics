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
#include "TRint.h"
#include "TApplication.h"
#include <time.h>
#include <string>
#include <map>
#include "AliZMQhelpers.h"

//this is meant to become a class, hence the structure with global vars etc.
//Also the code is rather flat - it is a bit of a playground to test ideas.
//TODO structure this at some point, e.g. introduce a SIMPLE unified way of handling
//zmq payloads, maybe a AliZMQmessage class which would by default be multipart and provide
//easy access to payloads based on topic or so (a la HLT GetFirstInputObject() etc...)

typedef std::map<std::string,std::string> stringMap;

//methods
TString GetFullArgString(int argc, char** argv);
int ProcessOptionString(TString arguments);
stringMap* TokenizeOptionString(const TString str);
int ProcessOption(TString option, TString value);
int InitZMQ();
void* work(void* param);
int Run();

//configuration vars
Bool_t  fVerbose = kFALSE;
TString fZMQconfigIN   = "PULL@tcp://*:60201";
TString fZMQconfigOUT  = "PUB@tcp://*:60211";
TString fZMQconfigMON  = "";

Bool_t  fSendOnMerge = kTRUE;
Bool_t  fResetOnSend = kFALSE;

//ZMQ stuff
void* fZMQcontext = NULL;    //ze zmq context

void* fZMQmon = NULL;        //the request-reply socket, here we request the merged data
void* fZMQout = NULL;        //the monitoring socket, here we publish a copy of the data
void* fZMQin  = NULL;        //the in socket - entry point for the data to be merged.

void* work(void* /*param*/)
{
  return NULL;
}

//_______________________________________________________________________________________
int Run()
{
  int rc = 0;

  rc = zmq_proxy(fZMQin, fZMQout, fZMQmon);

  return rc;
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
 
  if (option.EqualTo("ZMQconfigIN") || option.EqualTo("in"))
  {
    fZMQconfigIN = value;
  }
  else if (option.EqualTo("ZMQconfigOUT") || option.EqualTo("out"))
  {
    fZMQconfigOUT = value;
  }
  else if (option.EqualTo("ZMQconfigMON") || option.EqualTo("mon"))
  {
    fZMQconfigMON = value;
  }
  return 1; 
}

//_______________________________________________________________________________________
int InitZMQ()
{
  //init or reinit stuff
  int rc = 0;
  rc += alizmq_socket_init(fZMQin,  fZMQcontext, fZMQconfigIN.Data(), 0);
  rc += alizmq_socket_init(fZMQout, fZMQcontext, fZMQconfigOUT.Data(), 0);
  rc += alizmq_socket_init(fZMQmon, fZMQcontext, fZMQconfigMON.Data(), 0);
  return rc;
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
  stringMap* options = TokenizeOptionString(arguments);
  for (stringMap::iterator i=options->begin(); i!=options->end(); ++i)
  {
    //Printf("  %s : %s", i->first.data(), i->second.data());
    ProcessOption(i->first,i->second);
  }
  delete options; //tidy up

  return 1; 
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
  // options can be separated by ' ' or ',' arbitrarily combined, e.g:
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

//_______________________________________________________________________________________
int main(int argc, char** argv)
{
  int mainReturnCode=0;

  //process args
  TString argString = GetFullArgString(argc,argv);
  ProcessOptionString(argString);

  //globally enable schema evolution for serializing ROOT objects
  TMessage::EnableSchemaEvolutionForAll(kTRUE);
  //the context
  fZMQcontext = zmq_ctx_new();

  //init stuff
  if (InitZMQ()<0) {
    Printf("failed init");
    return 1;
  }

  Run();

  //destroy ZMQ sockets
  zmq_close(fZMQmon);
  zmq_close(fZMQin);
  //zmq_close(fZMQout);
  zmq_ctx_destroy(fZMQcontext);
  return mainReturnCode;
}

