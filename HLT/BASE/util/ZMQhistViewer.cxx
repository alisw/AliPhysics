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
int fPollTimeout = 1000; //1s

//internal state
void* fZMQcontext = NULL;             //ze zmq context
void* fZMQin  = NULL;                 //the in socket - entry point for the data to be merged.

TApplication* gApp;
TCanvas* fCanvas;
TObjArray fDrawables;

TPRegexp* fSelectionRegexp = NULL;
TString fDrawOptions;

ULong64_t iterations=0;

const char* fUSAGE = 
    "ZMQhstViewer: Draw() all ROOT drawables in a message\n"
    "options: \n"
    " -in : data in\n"
    " -sleep : how long to sleep in between requests for data in s (if applicable)\n"
    " -timeout : how long to wait for the server to reply (s)\n"
    " -Verbose : be verbose\n"
    " -select : only show selected histograms (by regexp)\n"
    " -drawoptions : what draw option to use\n"
    " -file : dump input to file and exit\n"
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
      TString request;
      if (fSelectionRegexp) request = "select="+fSelectionRegexp->GetPattern();
      TString requestTopic;
      if (fSelectionRegexp) requestTopic = "CONFIG";

      if (fVerbose) Printf("sending request %s %s",requestTopic.Data(), request.Data());
      zmq_send(fZMQin, requestTopic.Data(), requestTopic.Length(), ZMQ_SNDMORE);
      zmq_send(fZMQin, request.Data(), request.Length(), ZMQ_SNDMORE);
      zmq_send(fZMQin, "", 0, ZMQ_SNDMORE);
      zmq_send(fZMQin, "", 0, 0);
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

    if (!sockets[0].revents & ZMQ_POLLIN)
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
        TObject* object;
        alizmq_msg_iter_data(i, object);
        if (object) UpdatePad(object);

        if (!fFileName.IsNull()) 
        {
          DumpToFile(object);
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
  fFile->Close();
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

  if (fVerbose) Printf("in: %s", name);
  Bool_t selected = kTRUE;
  if (fSelectionRegexp) selected = fSelectionRegexp->Match(name);
  if (!selected) return 0;
 
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
    
    //after we clear the canvas, the pads are gone, clear the pad cache as well
    fCanvas->Clear();
    fCanvas->DivideSquare(fDrawables.GetLast()+1);
    
    //redraw all objects at their old places and the new one as last
    for (int i=0; i<fDrawables.GetLast()+1; i++)
    {
      TObject* obj = fDrawables[i];
      fCanvas->cd(i+1);
      if (fVerbose) Printf("  drawing %s in pad %i", obj->GetName(), i);
      if (obj) obj->Draw();
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
  stringMap* options = AliOptionParser::TokenizeOptionString(arguments);
  int nOptions = 0;
  for (stringMap::iterator i=options->begin(); i!=options->end(); ++i)
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
    else if (option.EqualTo("drawoptions"))
    {
      fDrawOptions = value;
    }
    else if (option.EqualTo("file"))
    {
      fFileName = value;
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

