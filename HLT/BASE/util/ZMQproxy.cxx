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

//methods
int ProcessOptionString(TString arguments);
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

int fZMQtimeout = -1;
int fSleep = 0;
std::string reqTopic;
std::string reqBody;

const char* fUSAGE =
  "ZMQproxy: a simple monitored ZMQ proxy\n"
  "options:\n"
  " -in : socket in\n"
  " -out : socket out\n"
  " -mon : monitor socket\n"
  " -sleep : sleep between polls\n"
  " -timeout : timeout for a poll\n"
  " -requestTopic : request topic\n"
  " -requestBody : request body\n"
  ;

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

//______________________________________________________________________________
int ProcessOptionString(TString arguments)
{
  //process passed options
  aliStringVec* options = AliOptionParser::TokenizeOptionString(arguments);
  int nOptions = 0;
  for (aliStringVec::iterator i=options->begin(); i!=options->end(); ++i)
  {
    const TString& option = i->first;
    const TString& value = i->second;
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
    else
    {
      nOptions=-1;
      break;
    }
    nOptions++;
  }
  delete options; //tidy up

  return nOptions; 
}

//_______________________________________________________________________________________
int main(int argc, char** argv)
{
  int mainReturnCode=0;

  //process args
  TString argString = AliOptionParser::GetFullArgString(argc,argv);
  if (ProcessOptionString(argString)<=0)
  {
    printf("%s",fUSAGE);
    return 1;
  }

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

