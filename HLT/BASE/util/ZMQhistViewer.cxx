#include "zmq.h"
#include <algorithm>
#include <iostream>
#include "AliHLTDataTypes.h"
#include "AliHLTComponent.h"
#include "AliHLTMessage.h"
#include "TClass.h"
#include "TCanvas.h"
#include "TMap.h"
#include "TPRegexp.h"
#include "TObjString.h"
#include "TDirectory.h"
#include "TList.h"
#include "AliZMQhelpers.h"
#include "TMessage.h"
#include "TSystem.h"
#include "TApplication.h"
#include "TH1.h"
#include "TH1F.h"
#include <time.h>
#include <string>
#include <map>
#include "TFile.h"
#include "TSystem.h"
#include "TStyle.h"
#include "signal.h"
class MySignalHandler;

//this is meant to become a class, hence the structure with global vars etc.
//Also the code is rather flat - it is a bit of a playground to test ideas.
//TODO structure this at some point, e.g. introduce a SIMPLE unified way of handling
//zmq payloads, maybe a AliZMQmessage class which would by default be multipart and provide
//easy access to payloads based on topic or so (a la HLT GetFirstInputObject() etc...)

//methods
int ProcessOptionString(TString arguments);
int GetData();
int UpdateCanvas();
int ProcessObject(TObject*);
int DumpToFile(TObject* object);
void* run(void* arg);
Int_t countpads(TVirtualPad *pad);

//configuration vars
Bool_t fVerbose = kFALSE;
TString fZMQconfigIN  = "PULL>tcp://localhost:60211";
int fZMQsocketModeIN=-1;
TString fZMQsubscriptionIN = "";
TString fFilter = "";

TString fFileName="";
TFile* fFile=NULL;

int fPollInterval = 0;
int fPollTimeout = 100000; //100s
Bool_t fSort = kTRUE;
stringMap fInfoMap;

//internal state
void* fZMQcontext = NULL;             //ze zmq context
void* fZMQin  = NULL;                 //the in socket - entry point for the data to be merged.

TApplication* gApp;
TCanvas* fCanvas;

std::vector<TObject*> fDrawables;
std::vector<TObject*> fIncomingObjects;

int fNdrawables = 0;
int fNpads = 0;

TString fStatus = "";
Int_t fRunNumber = 0;
TPRegexp* fSelectionRegexp = NULL;
TPRegexp* fUnSelectionRegexp = NULL;
TString fDrawOptions;
Bool_t fScaleLogX = kFALSE;
Bool_t fScaleLogY = kFALSE;
Bool_t fScaleLogZ = kFALSE;
Bool_t fResetOnRequest = kFALSE;
Int_t fHistStats = 0;

Bool_t fAllowResetAtSOR = kTRUE;

ULong64_t iterations=0;

const char* fUSAGE = 
"ZMQhstViewer: Draw() all ROOT drawables in a message\n"
"options: \n"
" -in : data in\n"
" -sleep : how long to sleep in between requests for data in s (if applicable)\n"
" -timeout : how long to wait for the server to reply (s)\n"
" -Verbose : be verbose\n"
" -select : only show selected histograms (by regexp)\n"
" -unselect : as select, only inverted\n"
" -drawoptions : what draw option to use\n"
" -file : dump input to file and exit\n"
" -log[xyz] : use log scale on [xyz] dimension\n"
" -histstats : histogram stat box options (default 0)\n"
" -AllowResetAtSOR : 0/1 to reset at change of run\n"
" -sort : 0/1 sort by title, by default on\n"
;

struct TObjectTitleComparator {
  bool operator()(const TObject* left, const string& right) {
    return right.compare(left->GetTitle())>0;
  }

  bool operator()(const TObject* left, const TObject* right) {
    return strcmp(left->GetTitle(),right->GetTitle())<0;
  }
};

struct TObjectNameComparator {
  bool operator()(const TObject* left, const string& right) {
    return right.compare(left->GetTitle())>0;
  }

  bool operator()(const TObject* left, const TObject* right) {
    return strcmp(left->GetTitle(),right->GetTitle())<0;
  }
};

//_______________________________________________________________________________________
class MySignalHandler : public TSignalHandler
{
  public:
    MySignalHandler(ESignals sig) : TSignalHandler(sig) {}
    Bool_t Notify()
    {
      Printf("signal received, exiting");
      fgTerminationSignaled = true;
      return TSignalHandler::Notify();
    }
    static bool TerminationSignaled() { return fgTerminationSignaled; }
    static bool fgTerminationSignaled;
};
bool MySignalHandler::fgTerminationSignaled = false;

void sig_handler(int signo)
{
  if (signo == SIGINT)
    printf("received SIGINT\n");
  MySignalHandler::fgTerminationSignaled=true;
}

//_______________________________________________________________________________________
void* run(void* arg)
{
  fDrawables.reserve(1000);
  fIncomingObjects.reserve(1000);

  //main loop
  while(!MySignalHandler::TerminationSignaled())
  {
    errno=0;
    //send a request if we are using REQ
    if (fZMQsocketModeIN==ZMQ_REQ)
    {
      string request;

      if (fSelectionRegexp || fUnSelectionRegexp) 
      {
        if (fSelectionRegexp) request += " select="+fSelectionRegexp->GetPattern();
        if (fUnSelectionRegexp) request += " unselect="+fUnSelectionRegexp->GetPattern();
        if (fResetOnRequest) request += " ResetOnRequest";
        alizmq_msg_send("CONFIG", request, fZMQin, ZMQ_SNDMORE);
      }

      if (fVerbose) Printf("sending request CONFIG %s", request.c_str());
      alizmq_msg_send("", "", fZMQin, 0);
    }

    //wait for the data
    zmq_pollitem_t sockets[] = { 
      { fZMQin, 0, ZMQ_POLLIN, 0 }, 
    };
    int rc = zmq_poll(sockets, 1, (fZMQsocketModeIN==ZMQ_REQ)?fPollTimeout:-1);

    if (rc==-1 && errno==ETERM)
    {
      Printf("jumping out of main loop");
      break;
    }

    if (!(sockets[0].revents & ZMQ_POLLIN))
    {
      //server died
      Printf("connection timed out, server %s died?", fZMQconfigIN.Data());
      fZMQsocketModeIN = alizmq_socket_init(fZMQin, fZMQcontext, fZMQconfigIN.Data());
      if (fZMQsocketModeIN < 0) return NULL;
      continue;
    }
    else
    {
      GetData();
      UpdateCanvas();

    }//socket 0
    gSystem->ProcessEvents();
    usleep(fPollInterval);
  }//main loop
  return NULL;
}

//_______________________________________________________________________________________
int GetData()
{
  //get all data (topic+body), possibly many of them
  fIncomingObjects.clear();
  aliZMQmsg message;
  alizmq_msg_recv(&message, fZMQin, 0);

  //process message, deserialize objects, puth them in the container 
  for (aliZMQmsg::iterator i=message.begin(); i!=message.end(); ++i)
  {
    if (alizmq_msg_iter_check(i, "INFO")==0)
    {
      //check if we have a runnumber in the string
      string info;
      alizmq_msg_iter_data(i,info);
      if (fVerbose) Printf("processing INFO %s", info.c_str());

      fCanvas->SetTitle(info.c_str());

      fInfoMap = ParseParamString(info);
      int runnumber = atoi(fInfoMap["run"].c_str());

      if (fVerbose) printf("received run=%i\n",runnumber);

      if (runnumber!=fRunNumber && fAllowResetAtSOR) 
      {
        if (fVerbose) printf("Run changed, resetting!\n");
        vector<TObject*> tmp = fDrawables;
        fDrawables.clear();
        fCanvas->Clear();
        gSystem->ProcessEvents();

        for (vector<TObject*>::iterator i=tmp.begin(); i!=tmp.end(); ++i)
        {
          delete *i;
        }
      }
      fRunNumber = runnumber; 
      continue;
    }

    TObject* object;
    alizmq_msg_iter_data(i, object);
    ProcessObject(object);

    if (!fFileName.IsNull()) 
    {
      if (object) DumpToFile(object);
    }
  } //for iterator i
  alizmq_msg_close(&message);
  return 0;
}

//_______________________________________________________________________________________
int main(int argc, char** argv)
{
  //process args
  int noptions = ProcessOptionString(AliOptionParser::GetFullArgString(argc,argv));
  if (noptions<=0) 
  {
    printf("%s",fUSAGE);
    return 1;
  }

  TH1::AddDirectory(kFALSE);
  TDirectory::AddDirectory(kFALSE);
  gApp = new TApplication("viewer", &argc, argv); 
  gApp->SetReturnFromRun(true);
  //gApp->Run();

  gStyle->SetOptStat(fHistStats);
  fCanvas = new TCanvas();
  gSystem->ProcessEvents();

  int mainReturnCode=0;

  //init stuff
  //globally enable schema evolution for serializing ROOT objects
  TMessage::EnableSchemaEvolutionForAll(kTRUE);
  //ZMQ init
  fZMQcontext = zmq_ctx_new();
  fZMQsocketModeIN = alizmq_socket_init(fZMQin, fZMQcontext, fZMQconfigIN.Data(), -1, 2);
  if (fZMQsocketModeIN < 0) return 1;

  gSystem->ResetSignal(kSigPipe);
  gSystem->ResetSignal(kSigQuit);
  gSystem->ResetSignal(kSigInterrupt);
  gSystem->ResetSignal(kSigTermination);
  //gSystem->AddSignalHandler(new MySignalHandler(kSigPipe));
  //gSystem->AddSignalHandler(new MySignalHandler(kSigQuit));
  //gSystem->AddSignalHandler(new MySignalHandler(kSigInterrupt));
  //gSystem->AddSignalHandler(new MySignalHandler(kSigTermination));

  if (signal(SIGINT, sig_handler) == SIG_ERR)
    printf("\ncan't catch SIGINT\n");

  run(NULL);

  Printf("exiting...");
  if (fFile) fFile->Close();
  delete fFile; fFile=0;

  //destroy ZMQ sockets

  int linger=0;
  zmq_setsockopt(fZMQin, ZMQ_LINGER, &linger, sizeof(linger));
  zmq_close(fZMQin);
  zmq_ctx_term(fZMQcontext);
  return mainReturnCode;
}

//______________________________________________________________________________
int DumpToFile(TObject* object)
{
  Option_t* fileMode="RECREATE";
  if (!fFile) fFile = new TFile(fFileName,fileMode);
  if (fVerbose) Printf("writing object %s to %s",object->GetName(), fFileName.Data());
  int rc = object->Write(object->GetName(),TObject::kOverwrite);
  MySignalHandler::fgTerminationSignaled=true;
  return rc;
}

//______________________________________________________________________________
int ProcessObject(TObject* object)
{
  if (!object) return -1;
  string name = object->GetName();
  string title = object->GetTitle();

  if (fVerbose) Printf("in: %s (%s)", name.c_str(), object->ClassName());

  Bool_t selected = kTRUE;
  Bool_t unselected = kFALSE;
  if (fSelectionRegexp) selected = fSelectionRegexp->Match(name);
  if (fUnSelectionRegexp) unselected = fUnSelectionRegexp->Match(name);
  if (!selected || unselected) 
  {
    delete object;
    return 0;
  }

  fIncomingObjects.push_back(object);

  return 0;
}

//______________________________________________________________________________
Int_t countpads(TVirtualPad *pad) {
  //count the number of pads in pad
  if (!pad) return 0;
  Int_t npads = 0;
  TObject *obj;
  TIter next(pad->GetListOfPrimitives());
  while ((obj = next())) {
    if (obj->InheritsFrom(TVirtualPad::Class())) npads++;
  }
  return npads;
}

//______________________________________________________________________________
int UpdateCanvas()
{
  std::vector<TObject*> fOldDrawables; 
  fOldDrawables.reserve(1000);

  std::vector<int> oldPadNumbers(fIncomingObjects.size());

  bool weHaveNewPlots = kFALSE;

  int ipad = 0;
  for (vector<TObject*>::iterator i=fIncomingObjects.begin();
      i!=fIncomingObjects.end();
      ++i)
  {
    TObject* newObj = *i;
    TObject* oldObj = NULL;

    const char* newObjName = newObj->GetName();

    //check if it will replace something
    bool found=false;
    int jpad = 0;
    for (vector<TObject*>::iterator j=fDrawables.begin();
        j!=fDrawables.end();
        ++j)
    {
      oldObj = *j;
      const char* oldObjName = oldObj->GetName();
      if (strcmp(newObjName,oldObjName)==0)
      {
        if (fVerbose) printf("found! %s\n",oldObjName);
        found = true;
        *j = *i;
        fOldDrawables.push_back(oldObj); //tag for removal
        oldPadNumbers[ipad]=jpad;
        break;
      }
      jpad++;
    }
    if (!found)
    {
      if (fVerbose) printf("new object %s\n",newObjName);
      fDrawables.push_back(newObj); //completely new one
      weHaveNewPlots = true;
    }
    else
    {
      if (fVerbose) printf("found...\n");
    }
    ipad++;
  }

  //re-sort the new list if we have new plots
  if (weHaveNewPlots && fSort) {
    if (fVerbose) { printf("sorting\n");
      printf("before\n");
      for (vector<TObject*>::iterator i=fDrawables.begin(); i!=fDrawables.end(); ++i)
      {printf("  %s",(*i)->GetTitle());}
      printf("\n");
    }

    std::sort(fDrawables.begin(), fDrawables.end(), TObjectTitleComparator());
    
    if (fVerbose) {
      printf("after\n");
      for (vector<TObject*>::iterator i=fDrawables.begin(); i!=fDrawables.end(); ++i)
      {printf("  %s",(*i)->GetTitle());}
      printf("\n");
    }
  }

  //after we clear the canvas, the pads are gone, clear the pad cache as well
  fNdrawables = fDrawables.size();
  fNpads = countpads(fCanvas);

  if (fVerbose) printf("nDrawables: %i, nPads: %i\n", fNdrawables, fNpads);
  if (fNdrawables > fNpads)
  {
    fCanvas->Clear();
    fCanvas->DivideSquare(fNdrawables);
    fNpads = countpads(fCanvas);
    if (fVerbose) printf("reorganizing canvas, now %i pads\n",fNpads);
  }

  //if we have new plots we need to replot the whole thing as the old plots may have been
  //moved
  //if we have no new plots, we just draw the incoming ones at correct locations
  int padIndex = 0;
  if (weHaveNewPlots)
  {
    if (fVerbose) printf("have new plots(%lu), drawing\n", fDrawables.size());
    for (vector<TObject*>::iterator i=fDrawables.begin(); i!=fDrawables.end(); ++i)
    {
      if (fVerbose) printf("plotting %s\n", (*i)->GetName());
      fCanvas->cd(++padIndex);
      gPad->GetListOfPrimitives()->Remove(*i);
      gPad->Clear();
      (*i)->Draw(fDrawOptions);
      gPad->Modified(true);
    }
  }
  else
  {
    if (fVerbose) printf("dont have new plots, redrawing %lu\n",fIncomingObjects.size());
    for (vector<TObject*>::iterator i=fIncomingObjects.begin(); i!=fIncomingObjects.end(); ++i)
    {
      if (fVerbose) printf("plotting %s at pad %i\n", (*i)->GetName(),oldPadNumbers[padIndex]);
      fCanvas->cd(oldPadNumbers[padIndex++]+1);
      gPad->GetListOfPrimitives()->Remove(*i);
      gPad->Clear();
      (*i)->Draw(fDrawOptions);
      gPad->Modified(true);
    }
  }

  //now we can delete old plots
  for (vector<TObject*>::iterator i=fOldDrawables.begin(); i!=fOldDrawables.end(); ++i)
  {
    delete *i;
  }
  fOldDrawables.clear();

  gSystem->ProcessEvents();
  fCanvas->Update();
  return 0;
}

//______________________________________________________________________________
int ProcessOptionString(TString arguments)
{
  //process passed options
  aliStringVec* options = AliOptionParser::TokenizeOptionString(arguments);
  int nOptions = 0;
  for (aliStringVec::iterator i=options->begin(); i!=options->end(); ++i)
  {
    //Printf("  %s : %s", i->first.data(), i->second.data());
    const TString& option = i->first;
    const TString& value = i->second;
    if (option.EqualTo("PollInterval") || option.EqualTo("sleep"))
    {
      fPollInterval = round(value.Atof()*1e6);
    }
    else if (option.EqualTo("PollTimeout") || option.EqualTo("timeout"))
    {
      fPollTimeout = round(value.Atof()*1e3);
    }
    else if (option.EqualTo("ZMQconfigIN") || option.EqualTo("in") )
    {
      fZMQconfigIN = value;
    }
    else if (option.EqualTo("Verbose"))
    {
      fVerbose=kTRUE;
    }
    else if (option.EqualTo("select"))
    {
      delete fSelectionRegexp;
      fSelectionRegexp=new TPRegexp(value);
    }
    else if (option.EqualTo("unselect"))
    {
      delete fUnSelectionRegexp;
      fUnSelectionRegexp=new TPRegexp(value);
    }
    else if (option.EqualTo("ResetOnRequest"))
    {
      fResetOnRequest = kTRUE;
    }
    else if (option.EqualTo("drawoptions"))
    {
      fDrawOptions = value;
    }
    else if (option.EqualTo("logx"))
    {
      fScaleLogX=kTRUE;
    }
    else if (option.EqualTo("logy"))
    {
      fScaleLogY=kTRUE;
    }
    else if (option.EqualTo("logz"))
    {
      fScaleLogZ=kTRUE;
    }
    else if (option.EqualTo("file"))
    {
      fFileName = value;
    }
    else if (option.EqualTo("sort"))
    {
      fSort=value.Contains(0)?kFALSE:kTRUE;
    }
    else if (option.EqualTo("histstats"))
    {
      fHistStats = value.Atoi();
    }
    else if (option.EqualTo("AllowResetAtSOR"))
    {
      fAllowResetAtSOR = (option.Contains("0")||option.Contains("no"))?kFALSE:kTRUE;
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

