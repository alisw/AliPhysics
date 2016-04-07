#include "zmq.h"
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
int UpdatePad(TObject*);
int DumpToFile(TObject* object);
void* run(void* arg);

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

//internal state
void* fZMQcontext = NULL;             //ze zmq context
void* fZMQin  = NULL;                 //the in socket - entry point for the data to be merged.

TApplication* gApp;
TCanvas* fCanvas;
TObjArray fDrawables;

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
    ;
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
      //get all data (topic+body), possibly many of them
      aliZMQmsg message;
      alizmq_msg_recv(&message, fZMQin, 0);
      for (aliZMQmsg::iterator i=message.begin(); i!=message.end(); ++i)
      {
        if (alizmq_msg_iter_check(i, "INFO")==0)
        {
            //check if we have a runnumber in the string
            string info;
            alizmq_msg_iter_data(i,info);
            if (fVerbose) Printf("processing INFO %s", info.c_str());
            
            fCanvas->SetTitle(info.c_str());
            size_t runTagPos = info.find("run");
            if (runTagPos != std::string::npos)
            {
              size_t runStartPos = info.find("=",runTagPos);
              size_t runEndPos = info.find(" ");
              string runString = info.substr(runStartPos+1,runEndPos-runStartPos-1);
              if (fVerbose) printf("received run=%s\n",runString.c_str());
  
              int runnumber = atoi(runString.c_str());
  
              if (runnumber!=fRunNumber && fAllowResetAtSOR) 
              {
                  if (fVerbose) printf("Run changed, resetting!\n");
                  fDrawables.Delete();
                  fCanvas->Clear();
                  gSystem->ProcessEvents();
              }
              fRunNumber = runnumber; 
            }
            continue;
        }

        TObject* object;
        alizmq_msg_iter_data(i, object);
        if (object) UpdatePad(object);

        if (!fFileName.IsNull()) 
        {
          if (object) DumpToFile(object);
        }
      }
      alizmq_msg_close(&message);

    }//socket 0
    gSystem->ProcessEvents();
    usleep(fPollInterval);
  }//main loop
  return NULL;
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
int UpdatePad(TObject* object)
{
  if (!object) return -1;
  const char* name = object->GetName();
  
  TObject* drawable = fDrawables.FindObject(name);
  int padIndex = fDrawables.IndexOf(drawable);

  if (fVerbose) Printf("in: %s (%s)", name, object->ClassName());
  Bool_t selected = kTRUE;
  Bool_t unselected = kFALSE;
  if (fSelectionRegexp) selected = fSelectionRegexp->Match(name);
  if (fUnSelectionRegexp) unselected = fUnSelectionRegexp->Match(name);
  if (!selected || unselected) 
  {
      delete object;
      return 0;
  }
 
  if (drawable)
  {
    //only redraw the one thing
    if (fVerbose) Printf("  redrawing %s in pad %i", name, padIndex);
    fCanvas->cd(padIndex+1);
    gPad->GetListOfPrimitives()->Remove(drawable);
    gPad->Clear();
    fDrawables.RemoveAt(padIndex);
    delete drawable;
    fDrawables.AddAt(object, padIndex);
    object->Draw(fDrawOptions);
    gPad->Modified(kTRUE);
  }
  else
  {
    if (fVerbose) Printf("  new object %s", name);
    //add the new object to the collection, restructure the canvas 
    //and redraw everything
    fDrawables.AddLast(object);

    if (fSort)
    {
      TObjArray sortedTitles(fDrawables.GetEntries());
      for (int i=0; i<fDrawables.GetEntries(); i++)
      { sortedTitles.AddAt(new TNamed(fDrawables[i]->GetTitle(),fDrawables[i]->GetName()),i); }
      sortedTitles.Sort();
      TObjArray sortedDrawables(fDrawables.GetEntries());
      for (int i=0; i<fDrawables.GetEntries(); i++)
      {
        const char* name = sortedTitles[i]->GetTitle();
        TObject* tmp = fDrawables.FindObject(name);
        int index = fDrawables.IndexOf(tmp);
        sortedDrawables.AddAt(fDrawables[index],i); 
      }
      for (int i=0; i<sortedDrawables.GetEntries(); i++)
      { fDrawables.AddAt(sortedDrawables[i],i); }
      sortedTitles.Delete();
    }
    
    //after we clear the canvas, the pads are gone, clear the pad cache as well
    fCanvas->Clear();
    fCanvas->DivideSquare(fDrawables.GetLast()+1);
    
    //redraw all objects at their old places and the new one as last
    for (int i=0; i<fDrawables.GetLast()+1; i++)
    {
      TObject* obj = fDrawables[i];
      fCanvas->cd(i+1);
      if (fScaleLogX) gPad->SetLogx();
      if (fScaleLogY) gPad->SetLogy();
      if (fScaleLogZ) gPad->SetLogz();
      if (fVerbose) Printf("  drawing %s in pad %i", obj->GetName(), i);
      if (obj) obj->Draw(fDrawOptions);
    }
  }
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

