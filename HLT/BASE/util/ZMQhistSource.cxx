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
    
int ProcessOptionString(TString arguments);

//_______________________________________________________________________________________
int main(int argc, char** argv)
{
  int mainReturnCode=0;
  int rc = 0;

  //process args
  if (ProcessOptionString(AliOptionParser::GetFullArgString(argc,argv))<0) return 1;
  
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
int ProcessOptionString(TString arguments)
{
  //process passed options
  int nOptions = 0;
  stringMap* options = AliOptionParser::TokenizeOptionString(arguments);
  for (stringMap::iterator i=options->begin(); i!=options->end(); ++i)
  {
    Printf("  %s : %s", i->first.data(), i->second.data());
    const TString& option = i->first;
    const TString& value = i->second;
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

