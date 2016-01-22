#include "zmq.h"
#include "AliZMQhelpers.h"
#include <iostream>
#include "TTimeStamp.h"
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
int DumpToFile(TObject* object);
void* run(void* arg);

//configuration vars
Bool_t fVerbose = kFALSE;
TString fZMQconfigIN  = "PULL>tcp://localhost:60211";
int fZMQsocketModeIN=-1;
TString fZMQsubscriptionIN = "";
TString fFilter = "";

TString fFileNameBase="";
TString fFileName="";
TFile* fFile=NULL;
Int_t fFileNumber = 0;

int fPollInterval = 0;
int fPollTimeout = 1000; //1s
Bool_t fSort = kTRUE;

//internal state
void* fZMQcontext = NULL;             //ze zmq context
void* fZMQin  = NULL;                 //the in socket - entry point for the data to be merged.

TString fStatus = "";
Int_t fRunNumber = 0;
TPRegexp* fSelectionRegexp = NULL;
TPRegexp* fUnSelectionRegexp = NULL;

bool fgTerminationSignaled=false;

ULong64_t iterations=0;

const char* fUSAGE = 
    "ZMQfileSink: dump contents of a multipart message into a file\n"
    "options: \n"
    " -in : data in\n"
    " -sleep : how long to sleep in between requests for data in s (if applicable) (-1 means dump just once)\n"
    " -once : write one file and exit (same as -sleep=-1)"
    " -timeout : how long to wait for the server to reply (s)\n"
    " -Verbose : be verbose\n"
    " -select : select objects (by regexp)\n"
    " -unselect : as select, only inverted\n"
    " -file : file name\n"
    ;

//_______________________________________________________________________________________
void sig_handler(int signo)
{
  if (signo == SIGINT)
    printf("received signal\n");
  fgTerminationSignaled=true;
}

//_______________________________________________________________________________________
void* run(void* arg)
{
  //main loop
  while(!fgTerminationSignaled)
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
      
      //count ROOT objects
      int numberOfTObjects=0;
      for (aliZMQmsg::iterator i=message.begin(); i!=message.end(); ++i)
      {
        if (alizmq_msg_iter_check(i, "ROOT")==0) numberOfTObjects++;
      }

      //process
      for (aliZMQmsg::iterator i=message.begin(); i!=message.end(); ++i)
      {
        if (alizmq_msg_iter_check(i, "INFO")==0)
        {
            //check if we have a runnumber in the string
            string info;
            alizmq_msg_iter_data(i,info);
            if (fVerbose) Printf("processing INFO %s", info.c_str());
            
            size_t runTagPos = info.find("run");
            if (runTagPos != std::string::npos)
            {
              size_t runStartPos = info.find("=",runTagPos);
              size_t runEndPos = info.find(" ");
              string runString = info.substr(runStartPos+1,runEndPos-runStartPos-1);
              if (fVerbose) printf("received run=%s\n",runString.c_str());
  
              int runnumber = atoi(runString.c_str());
  
              fRunNumber = runnumber; 
            }
            continue;
        }

        TObject* object = NULL;
        alizmq_msg_iter_data(i, object);

        if (fVerbose) Printf("got %s %s", object->ClassName(), object->GetName());

        if (object && !fFileNameBase.IsNull()) 
        {
          Option_t* fileMode="RECREATE";
          TTimeStamp time;
          TString timestamp = time.AsString("s");
          timestamp.ReplaceAll(" ","_");

          char runStr[10];
          snprintf(runStr,10,"%.9i",fRunNumber);

          fFileName  = fFileNameBase;
          fFileName += "_";
          fFileName += runStr;
          fFileName += "_"+timestamp;
          if (numberOfTObjects==1) 
          {
            TString objectName = object->GetName();
            objectName.ReplaceAll(" ","_");
            fFileName += "_"+objectName;
          }
          fFileName += ".";
          fFileName += fFileNumber;
          fFileName += ".root";
          if (fVerbose) Printf("opening file: %s", fFileName.Data());
          if (!fFile) fFile = new TFile(fFileName,fileMode);
          DumpToFile(object);
          delete object;
        }
        else
        {
          Printf("no object or no file");
        }
      }
      alizmq_msg_close(&message);
      delete fFile; fFile=NULL;
      fFileNumber++;

    }//socket 0
    if (fPollInterval<0) break;
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

  if (fFileNameBase.IsNull())
  {
    Printf("a filen name must be specified with -file option!");
    return 1;
  }

  int mainReturnCode=0;

  //init stuff
  //globally enable schema evolution for serializing ROOT objects
  TMessage::EnableSchemaEvolutionForAll(kTRUE);
  TDirectory::AddDirectory(kFALSE);
  //ZMQ init
  fZMQcontext = zmq_ctx_new();
  fZMQsocketModeIN = alizmq_socket_init(fZMQin, fZMQcontext, fZMQconfigIN.Data(), -1, 2);
  if (fZMQsocketModeIN < 0) return 1;

  if (signal(SIGHUP, sig_handler) == SIG_ERR)
  printf("\ncan't catch SIGHUP\n");
  if (signal(SIGINT, sig_handler) == SIG_ERR)
  printf("\ncan't catch SIGINT\n");
  if (signal(SIGQUIT, sig_handler) == SIG_ERR)
  printf("\ncan't catch SIGQUIT\n");
  if (signal(SIGTERM, sig_handler) == SIG_ERR)
  printf("\ncan't catch SIGTERM\n");

  run(NULL);

  Printf("closing file(s) and exiting...");
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
  if (fVerbose) Printf("writing object %s to %s",object->GetName(), fFileName.Data());
  int rc = object->Write(object->GetName(),TObject::kOverwrite);
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
    //Printf("  %s : %s", i->first.data(), i->second.data());
    TString option = i->first;
    TString value = i->second;
    if (option.EqualTo("once"))
    {
      fPollInterval = -1;
    }
    else if (option.EqualTo("PollInterval") || option.EqualTo("sleep"))
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
    else if (option.EqualTo("file"))
    {
      fFileNameBase = value.ReplaceAll(".root","");
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

