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
#include "TKey.h"
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
int DumpToFile(TObject* object, std::string name);
int processInfo(aliZMQmsg::iterator& i);
void* run(void* arg);
void* runToFile(void* arg);
void* runFromFile(void* arg);

//configuration vars
TString fConfigIN = "";
TString fConfigOUT = "";
Bool_t fToFile = true;
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
TPRegexp fRunNumberRegex("([-_/]+)(000[0-9][0-9][0-9][0-9][0-9][0-9])([-_/]+)");
Int_t fRunNumber = -1;
std::string fInfo = "";
TPRegexp* fSelectionRegexp = NULL;
TPRegexp* fUnSelectionRegexp = NULL;

bool fgTerminationSignaled=false;
int fNumberOfTObjectsInMessage=0;

ULong64_t iterations=0;

const char* fUSAGE = 
    "ZMQfileProxy: proxy between a (multi-part) message and a ROOT file (sink/source)\n"
    "options: \n"
    " -in : data in (can be file://)\n"
    " -out : data out (can be file://)\n"
    " -sleep : how long to sleep in between requests for data in/out (s) (if applicable) (-1 means send/dump once and exit)\n"
    " -once : write one file (send contents of one file) and exit (same as -sleep=-1)\n"
    " -timeout : how long to wait for the server (when using a REQ socket) to reply before retrying (s)\n"
    " -Verbose : be verbose\n"
    " -select : request selected objects by name with a (perl compatible-) regexp\n"
    " -unselect : as select, only inverted\n"
    ;

//_______________________________________________________________________________________
void sig_handler(int signo)
{
  if (signo == SIGINT)
    printf("received signal\n");
  fgTerminationSignaled=true;
}

//_______________________________________________________________________________________
void* runFromFile(void* arg)
{
  //main loop
  while(!fgTerminationSignaled)
  {
    errno = 0;
    int rc = 0;

    //if we are replying, wait for a request
    if (fZMQsocketModeIN==ZMQ_REP)
    {
      aliZMQmsg request;
      alizmq_msg_recv(&request, fZMQin, 0);
      alizmq_msg_close(&request);
    }

    //init message to be sent
    aliZMQmsg message;

    if (!fFile)
    {
      fFile = new TFile(fFileName.Data(),"read");
      if (fFile->IsZombie())
      {
        Printf("cannot open %s",fFileName.Data());
        return NULL;
      }
    }

    TObject* infoObject = fFile->Get("_ZMQ_internal_fInfo");
    if (infoObject) {
      TObjString* objstr = dynamic_cast<TObjString*>(infoObject);
      if (objstr) {
        fInfo = objstr->String();
        fRunNumber = atoi(GetParamString("run",fInfo).c_str());
      }
      if (fVerbose) printf("read metadata: %s\n", fInfo.c_str());
      delete infoObject;
    }

    //parse the file name for run number
    if (fRunNumber<0) {
      TObjArray* eearr = fRunNumberRegex.MatchS(fFileName);
      TObject* runstr = NULL;
      if (eearr) runstr = eearr->At(2);
      if (runstr) fRunNumber = atoi(runstr->GetName());
      else fRunNumber = 0;
      delete eearr;
      fInfo+=Form("run=%i",fRunNumber);
    }

    alizmq_msg_add(&message, "INFO", fInfo);

    //attach each key as a message part
    TList* listOfKeys = fFile->GetListOfKeys();
    TIter keys(listOfKeys);
    while (TKey* key = (TKey*)keys.Next())
    {
      string objectName = key->GetName();
      
      //skip the info part (handled before)
      if (objectName=="_ZMQ_internal_fInfo") {
        continue;
      }
      
      //do the selection
      Bool_t selected = kTRUE;
      Bool_t unselected = kFALSE;
      if (fSelectionRegexp) selected = fSelectionRegexp->Match(objectName.c_str());
      if (fUnSelectionRegexp) unselected = fUnSelectionRegexp->Match(objectName.c_str());
      if (!selected || unselected)
      {
        if (fVerbose) Printf("skipping %s",objectName.c_str());
        continue;
      }

      //read and attach object
      TObject* object = key->ReadObj();
      if (object)
      {

        if (fVerbose) Printf("attaching %s",objectName.c_str());
        AliHLTDataTopic topic = kAliHLTDataTypeTObject;
        alizmq_msg_add(&message, &topic, object);
      }
    }

    //send the full message
    if (fVerbose) Printf("- sending, runNumber=%i",fRunNumber);
    alizmq_msg_send(&message, fZMQin, 0);
    alizmq_msg_close(&message);
    
    delete fFile; fFile=NULL;
    
    if (fZMQsocketModeIN==ZMQ_REQ)
    {
      //get a request if we are using REQ
      //discard the reply
      aliZMQmsg reply;
      rc = alizmq_msg_recv(&reply, fZMQin, 0);
      alizmq_msg_close(&reply);

      if (rc < 0)
      {
        //server died
        Printf("connection timed out, server %s died?", fZMQconfigIN.Data());
        fZMQsocketModeIN = alizmq_socket_init(fZMQin, fZMQcontext, fZMQconfigIN.Data());
        if (fZMQsocketModeIN < 0) 
        {
          Printf("cannot reinit ZMQ socket %s, %s, exiting...", fZMQconfigIN.Data(), zmq_strerror(errno));
          return NULL;
        }
        continue;
      }


    }
    
    if (fPollInterval<0) break;
    usleep(fPollInterval);
  }//main loop
  return NULL;
}

//_______________________________________________________________________________________
void* runToFile(void* arg)
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

    if (rc == 0 && fZMQsocketModeIN==ZMQ_REQ)
    {
      //server died
      Printf("connection timed out, server %s died?", fZMQconfigIN.Data());
      fZMQsocketModeIN = alizmq_socket_init(fZMQin, fZMQcontext, fZMQconfigIN.Data());
      if (fZMQsocketModeIN < 0) 
      {
        Printf("cannot reinit ZMQ socket %s, %s, exiting...", fZMQconfigIN.Data(), zmq_strerror(errno));
        return NULL;
      }
      continue;
    }
    

    if (sockets[0].revents & ZMQ_POLLIN)
    {
      //get all data (topic+body), possibly many of them
      aliZMQmsg message;
      alizmq_msg_recv(&message, fZMQin, 0);
      
      //count ROOT objects
      for (aliZMQmsg::iterator i=message.begin(); i!=message.end(); ++i)
      {
        if (alizmq_msg_iter_check_id(i, "ROOT")==0) fNumberOfTObjectsInMessage++;
      }

      //process
      for (aliZMQmsg::iterator i=message.begin(); i!=message.end(); ++i)
      {
        TObject* object = NULL;
        string name;
        if (alizmq_msg_iter_check_id(i, "INFO")==0)
        {
          //check if we have a runnumber in the string
          alizmq_msg_iter_data(i,fInfo);
          if (fVerbose) Printf("processing INFO %s", fInfo.c_str());
          fRunNumber = atoi(GetParamString("run",fInfo).c_str());
          object = new TObjString(fInfo.c_str());
          name = "_ZMQ_internal_fInfo";
        }
        else
        {
          alizmq_msg_iter_data(i, object);
          if (object) name = object->GetName();
        }

        if (fVerbose && object) Printf("got %s %s", object->ClassName(), object->GetName());

        if (object && !fFileNameBase.IsNull()) 
        {
          DumpToFile(object,name);
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
void* run(void* arg)
{
  if (fToFile)
    return runToFile(arg);
  else
    return runFromFile(arg);
}

//_______________________________________________________________________________________
int processInfo(aliZMQmsg::iterator& i)
{
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

  //reading FROM file
  if (fConfigIN.BeginsWith("file://") && !fConfigOUT.BeginsWith("file://"))
  {
    fFileName = fConfigIN;
    fFileName.ReplaceAll("file://","");
    fFileNameBase = fFileName;
    fZMQconfigIN = fConfigOUT;
    fToFile=false;
  }
  else
  //dumping TO file
  if (!fConfigIN.BeginsWith("file://") && fConfigOUT.BeginsWith("file://"))
  {
    fFileName = fConfigOUT;
    fFileName.ReplaceAll("file://","");
    fFileNameBase = fFileName;
    fZMQconfigIN = fConfigIN;
    fToFile=true;
  }
  else
  {
    Printf("either -in or -out MUST be \"file://\", not both");
    return 1;
  }
  
  fFileNameBase.ReplaceAll(".root","");
  if (fFileNameBase.IsNull())
  {
    Printf("no file specified!");
    return 1;
  }
  
  int mainReturnCode=0;

  //init stuff
  TDirectory::AddDirectory(kFALSE);
  //ZMQ init
  fZMQcontext = zmq_ctx_new();
  fZMQsocketModeIN = alizmq_socket_init(fZMQin, fZMQcontext, fZMQconfigIN.Data());
  if (fZMQsocketModeIN < 0) return 1;
  printf("in:  (%s) %s\n", alizmq_socket_name(fZMQsocketModeIN), fZMQconfigIN.Data());

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
int DumpToFile(TObject* object, std::string name)
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
  //if we just have one object, append it to file name for clarity)
  if (fNumberOfTObjectsInMessage==1) 
  {
    TString objectName = name.c_str();
    objectName.ReplaceAll(" ","_");
    fFileName += "_"+objectName;
  }
  fFileName += ".";
  fFileName += fFileNumber;
  fFileName += ".root";
  if (fVerbose) Printf("opening file: %s", fFileName.Data());
  if (!fFile) fFile = new TFile(fFileName,fileMode);
  if (fVerbose) Printf("writing object %s to %s",name.c_str(), fFileName.Data());
  int rc = object->Write(name.c_str(),TObject::kOverwrite);
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
    else if ( option.EqualTo("in") )
    {
      fConfigIN = value;
    }
    else if ( option.EqualTo("out") )
    {
      fConfigOUT = value;
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

