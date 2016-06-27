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
#include "AliHLTMessage.h"
#include "TH1F.h"
#include "TF1.h"
#include <cmath>
#include <time.h>
#include <string>
#include <map>
#include <unistd.h>
#include "AliHLTZMQhelpers.h"
#include <sstream>
#include <vector>
#include "TRandom.h"
#include "TTimeStamp.h"
#include "TSystem.h"
#include "TObjArray.h"
#include "AliHLTObjArray.h"
#include "AliOptionParser.h"
#include "AliHLTexampleMergeableContainer.h"

using namespace AliZMQhelpers;

void* fZMQout = NULL;
void* fZMQcontext = NULL;
TString fZMQconfigOUT = "PUSH>tcp://localhost:60210";
int fZMQsocketModeOUT = -1;

float fSleep = 1e6; //in microseconds
int fCount = 0;
int fNentries = 1;
int fRunNumber = -1;

vector<TH1F*> fHistograms;
int fNHistos = 1;
TString fHistName = "histogram";
TString fHistDistribution = "exp(-0.5*((x-0.)/0.1)**2)";
float fHistRangeLow = -0.5;
float fHistRangeHigh = 0.5;
int fHistNBins = 100;
aliZMQrootStreamerInfo* fSchema = NULL;
bool fVerbose = false;
int fCompression = 0;
TList* fCollection = NULL;
AliHLTObjArray* fAnalContainer = NULL;
TObjArray* fAnalComponentContainer = NULL;
AliHLTexampleMergeableContainer* fCustomContainer = NULL;
TH1F* fTemplateHist = NULL;
Int_t fNumberOfTemplateEntries = 5000000;

const char* fUSAGE =
    "ZMQhstSource: send a randomly filled ROOT histogram\n"
    "options: \n"
    " -out : data out\n"
    " -name : name of the histogram\n"
    " -sleep : [ms] how long to sleep between sending a new one\n"
    " -distribution : the pdf of the distribution\n"
    " -range : the range of the histogram, comma separated, e.g. -12.,12.\n"
    " -nbins : how many bins\n"
    " -count : how many messages to send before quitting (0 is never quit)\n"
    " -entries : how many entries in the histogram before sending\n"
    " -histos : how many histograms per message\n"
    " -schema : include the streamer infos in the message\n"
    " -run : run number\n"
    " -collection : wrap all histograms in a TObjArray\n"
    " -analysisContainer : wrap the collection in an AliHLTObjArray inside TObjArray like online\n"
    " -customContainer : wrap the collection in an AliHLTexampleMergeableContainer\n"
    //" -compression : compression level (0|1)\n"
    ;

int ProcessOptionString(TString arguments);

//_______________________________________________________________________________________
int main(int argc, char** argv)
{
  int mainReturnCode=0;
  int rc = 0;

  //process args
  if (ProcessOptionString(AliOptionParser::GetFullArgString(argc,argv))<=0)
  {
    printf("%s", fUSAGE);
    return 1;
  }

  //ZMQ init
  fZMQcontext = zmq_ctx_new();
  fZMQsocketModeOUT = alizmq_socket_init(fZMQout, fZMQcontext, fZMQconfigOUT.Data(), -1);
  if (fZMQsocketModeOUT < 0) return 1;

  TF1 formula("histDistribution",fHistDistribution);
  for (int i = 0; i < fNHistos; i++)
  {
    stringstream ss;
    ss << fHistName.Data();
    if (i>0) ss << i;
    TH1F* hist = new TH1F(ss.str().c_str(), ss.str().c_str(), fHistNBins, fHistRangeLow, fHistRangeHigh);
    hist->SetXTitle("x title");
    fHistograms.push_back(hist);
  }

  printf("initializing the distribution (%s) \n", formula.GetExpFormula("P").Data());
  TH1F* templateHist = new TH1F("templateHist", "templateHist", fHistNBins, fHistRangeLow, fHistRangeHigh);
  templateHist->FillRandom("histDistribution", fNumberOfTemplateEntries);
  printf("...done\n");

  if (fCollection)
  {
    for (int i = 0; i < fNHistos; i++)
    {
      fCollection->Add(fHistograms[i]);
    }
  }

  if (fCollection && fAnalContainer)
  {
    fAnalContainer->Add(fCollection);
  }

  if (fCustomContainer)
  {
    for (int i = 0; i < fNHistos; i++)
    {
      if (fVerbose) printf("adding histogram to custom container\n");
      fCustomContainer->Add(fHistograms[i]);
    }
  }

  if (fSchema) {
    if (fVerbose) printf("enabling schema for AliHLTMessage\n");
  }

  TTimeStamp time;
  gRandom->SetSeed(time.GetNanoSec()+gSystem->GetPid());

  //main loop
  int iterations=0;
  while(fCount==0 || iterations++<fCount)
  {
    for (int i = 0; i < fNHistos; i++)
    {
      fHistograms[i]->Reset();
      fHistograms[i]->FillRandom(templateHist, fNentries);
    }

    AliHLTDataTopic topic = kAliHLTDataTypeTObject;

    if (fZMQsocketModeOUT==ZMQ_REP)
    {
      //receive the request - could be multipart
      aliZMQmsgStr request;
      rc = alizmq_msg_recv(&request, fZMQout, 0);
    }

    if (fRunNumber>=0)
    {
      TString runInfo = "run=";
      runInfo+=fRunNumber;
      rc=alizmq_msg_send("INFO",runInfo.Data(), fZMQout, ZMQ_SNDMORE);
    }

    aliZMQmsg message;
    if (fCollection && !fAnalComponentContainer && !fCustomContainer)
    {
      if (fVerbose) printf("adding collection\n");
      rc = alizmq_msg_add(&message, &topic, fCollection, fCompression, fSchema);
      if (rc < 0)
        printf("unable to send\n");
    }
    else if (fCollection && fAnalComponentContainer)
    {
      if (fVerbose) printf("adding analysis container\n");
      rc = alizmq_msg_add(&message, &topic, fAnalComponentContainer, fCompression, fSchema);
      if (rc < 0)
        printf("unable to send\n");
    }
    else if (fCustomContainer)
    {
      if (fVerbose) printf("adding custom container\n");
      rc = alizmq_msg_add(&message, &topic, fCustomContainer, fCompression, fSchema);
      if (rc < 0)
        printf("unable to send\n");
    }
    else
    {
      for (int i = 0; i < fNHistos; i++)
      {
        if (fVerbose) printf("adding histogram directly\n");
        rc = alizmq_msg_add(&message, &topic, fHistograms[i], fCompression, fSchema);
        if (rc < 0)
          printf("unable to send\n");
      }
    }

    if (fSchema) alizmq_msg_prepend_streamer_infos(&message,fSchema);

    alizmq_msg_send(&message, fZMQout, 0);
    alizmq_msg_close(&message);
    if (fVerbose) printf("sent!\n");

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
  aliStringVec* options = AliOptionParser::TokenizeOptionString(arguments);
  for (aliStringVec::iterator i=options->begin(); i!=options->end(); ++i)
  {
    //Printf("  %s : %s", i->first.data(), i->second.data());
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
      fSleep = round(value.Atof()*1e3);
    }
    else if (option.EqualTo("distribution"))
    {
      fHistDistribution = value;
    }
    else if (option.EqualTo("run"))
    {
      fRunNumber = value.Atoi();
    }
    else if (option.EqualTo("collection"))
    {
      if (!fCollection) fCollection = new TList();
      fCollection->SetOwner(kTRUE);
    }
    else if (option.EqualTo("analysisContainer"))
    {
      fAnalContainer = new AliHLTObjArray(1);
      fAnalComponentContainer = new TObjArray(1);
      fAnalComponentContainer->Add(fAnalContainer);
      if (!fCollection) fCollection = new TList();
      fCollection->SetOwner(kTRUE);
    }
    else if (option.EqualTo("customContainer"))
    {
      fCustomContainer = new AliHLTexampleMergeableContainer("test");

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
    else if (option.EqualTo("histos")) {
      fNHistos = value.Atoi();
    }
    else if (option.EqualTo("schema")) {
      fSchema = new aliZMQrootStreamerInfo;
    }
    else if (option.EqualTo("Verbose")) {
      fVerbose=true;
    }
    else if (option.EqualTo("compression")) {
      fCompression=value.Atoi();
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

