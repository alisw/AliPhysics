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

//this is meant to become a class, hence the structure with global vars etc.
//Also the code is rather flat - it is a bit of a playground to test ideas.
//TODO structure this at some point, e.g. introduce a SIMPLE unified way of handling
//zmq payloads, maybe a AliZMQmessage class which would by default be multipart and provide
//easy access to payloads based on topic or so (a la HLT GetFirstInputObject() etc...)

typedef map<std::string,std::string> stringMap;

//methods
TObject* UnpackMessage( zmq_msg_t message);
TString GetFullArgString(int argc, char** argv);
int ProcessOptionString(TString arguments);
stringMap* TokenizeOptionString(const TString str);
int ProcessOption(TString option, TString value);
void deleteObject(void*, void*);
int UpdatePad(TObject*);

//configuration vars
Bool_t fVerbose = kFALSE;
TString fZMQconfigIN  = "PULL>tcp://localhost:60211";
int fZMQsocketModeIN=-1;
TString fZMQsubscriptionIN = "";
TString fFilter = "";

int fPollInterval = 0;
int fPollTimeout = 1000; //1s

//internal state
void* fZMQcontext = NULL;             //ze zmq context
void* fZMQin  = NULL;                 //the in socket - entry point for the data to be merged.

TApplication* gApp;
TCanvas* fCanvas;
TObjArray fDrawables;

ULong64_t iterations=0;

//_______________________________________________________________________________________
int main(int argc, char** argv)
{
  TH1::AddDirectory(kFALSE);
  TDirectory::AddDirectory(kFALSE);
  gApp = new TApplication("viewer", &argc, argv); 
  gApp->SetReturnFromRun(true);
  //gApp->Run();
  fCanvas = new TCanvas();
  gSystem->ProcessEvents();
  
  int mainReturnCode=0;

  //process args
  if (ProcessOptionString(GetFullArgString(argc,argv))!=0) return 1;

  //init stuff
  //globally enable schema evolution for serializing ROOT objects
  TMessage::EnableSchemaEvolutionForAll(kTRUE);
  //ZMQ init
  fZMQcontext = zmq_ctx_new();
  fZMQsocketModeIN = alizmq_socket_init(fZMQin, fZMQcontext, fZMQconfigIN.Data());
  if (fZMQsocketModeIN < 0) return 1;

  //main loop
  while(1)
  {
    errno=0;
    //send a request if we are using REQ
    if (fZMQsocketModeIN==ZMQ_REQ)
    {
      if (fVerbose) Printf("sending request");
      zmq_send(fZMQin, "*", 1, ZMQ_SNDMORE);
      zmq_send(fZMQin, "", 4, 0);
    }
    
    //wait for the data
    zmq_pollitem_t sockets[] = { 
                                 { fZMQin, 0, ZMQ_POLLIN, 0 }, 
                               };
    zmq_poll(sockets, 1, (fZMQsocketModeIN==ZMQ_REQ)?fPollTimeout:-1);

    if (!sockets[0].revents & ZMQ_POLLIN)
    {
      //server died
      Printf("connection timed out");
      fZMQsocketModeIN = alizmq_socket_init(fZMQin, fZMQcontext, fZMQconfigIN.Data());
      if (fZMQsocketModeIN < 0) return 1;
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
      }
      alizmq_msg_close(&message);

    }//socket 0
    gSystem->ProcessEvents();
    usleep(fPollInterval);
  }//main loop

  //destroy ZMQ sockets
  zmq_close(fZMQin);
  zmq_ctx_destroy(fZMQcontext);
  return mainReturnCode;
}

//______________________________________________________________________________
int UpdatePad(TObject* object)
{
  if (!object) return -1;
  const char* name = object->GetName();
  
  TObject* drawable = fDrawables.FindObject(name);
  int padIndex = fDrawables.IndexOf(drawable);
 
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
    object->Draw();
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
void deleteObject(void*, void* hint)
{
  //delete the TMessage, for use in zmq_msg_init_data(...) only.
  TObject* object = static_cast<TObject*>(hint);
  delete object;
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
 
  if (option.EqualTo("PollInterval") || option.EqualTo("sleep"))
  {
    fPollInterval = round(value.Atof()*1e6);
  }
  if (option.EqualTo("PollTimeout"))
  {
    fPollTimeout = round(value.Atof()*1e3);
  }
  else if (option.EqualTo("ZMQconfigIN") || option.EqualTo("in") )
  {
    fZMQconfigIN = value;
  }
  if (option.EqualTo("Verbose"))
  {
    fVerbose=kTRUE;
  }
  return 0; 
}

//_______________________________________________________________________________________
TObject* UnpackMessage( zmq_msg_t message )
{
  size_t size = zmq_msg_size(&message);
  void* data = zmq_msg_data(&message);
  return AliHLTMessage::Extract(data, size);
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
    Printf("  %s : %s", i->first.data(), i->second.data());
    if (ProcessOption(i->first,i->second) != 0) return 1;
  }
  delete options; //tidy up

  return 0; 
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
  //"-option option1=value1 --option2 value2 -option4=\'some string\'"
  
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

