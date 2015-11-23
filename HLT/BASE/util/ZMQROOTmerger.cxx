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
#include "TTimeStamp.h"

//this is meant to become a class, hence the structure with global vars etc.
//Also the code is rather flat - it is a bit of a playground to test ideas.
//TODO structure this at some point, e.g. introduce a SIMPLE unified way of handling
//zmq payloads, maybe a AliZMQmessage class which would by default be multipart and provide
//easy access to payloads based on topic or so (a la HLT GetFirstInputObject() etc...)

typedef std::map<std::string,std::string> stringMap;

//methods
TObject* UnpackMessage(zmq_msg_t* message);
Int_t ProcessOptionString(TString arguments);
Int_t InitZMQ();
void* InitZMQsocket(void* context, Int_t socketMode, const char* configs);
void* work(void* param);
Int_t Run();
long GetMilliSecSince(TTimeStamp* time);

Int_t HandleDataIn(zmq_msg_t* topicMsg, zmq_msg_t* dataMsg, void* /*socket*/=NULL);
Int_t HandleRequest(zmq_msg_t* /*topicMsg*/, zmq_msg_t* /*dataMsg*/, void* /*socket*/=NULL);

Int_t DoReceive(zmq_msg_t* topicMsg, zmq_msg_t* dataMsg, void* socket);
Int_t DoSend(void* socket);
Int_t DoReply(zmq_msg_t* topicMsg, zmq_msg_t* dataMsg, void* socket);
Int_t DoRequest(void* /*socket*/);

//merger private functions
void ResetOutputData();
Int_t Merge(TObject* object, TCollection* list);

//configuration vars
Bool_t  fVerbose = kFALSE;
TString fZMQconfigIN  = "PULL";
TString fZMQsubscriptionIN = "";
TString fZMQconfigOUT  = "PUSH";
TString fZMQconfigMON  = "REP";
Int_t   fZMQmaxQueueSize = 10;
Int_t   fZMQtimeout = 0;

Bool_t  fResetOnSend = kFALSE;      //reset on each send (also on scheduled pushing)
Bool_t  fResetOnRequest = kFALSE;   //reset once after a single request

TPRegexp* fSendSelection = NULL;

//internal state
TMap fMergeObjectMap;        //map of the merged objects, all incoming stuff is merged into these
TMap fMergeListMap;          //map with the lists of objects to be merged in
Int_t fMaxObjects = 1;        //trigger merge after this many messages

long fPushbackPeriod = -1;        //! in milliseconds
TTimeStamp fLastPushBackTime;

//ZMQ stuff
void* fZMQcontext = NULL;    //ze zmq context

void* fZMQmon = NULL;        //the request-reply socket, here we request the merged data
void* fZMQout = NULL;        //the monitoring socket, here we publish a copy of the data
void* fZMQin  = NULL;        //the in socket - entry point for the data to be merged.

const char* fUSAGE = 
    "ZMQROOTmerger options: Merge() all ROOT mergeables in the message.\n"
    "merge based on what GetName() returns, the merged data can be retrieved at any time.\n"
    " -in : data in, zmq config string, e.g. PUSH>tcp://localhost:123123\n"
    " -out : data out\n"
    " -mon : monitoring socket\n"
    " -Verbose : print some info\n"
    " -pushback-period : push the merged data once every n ms\n"
    " -ResetOnSend : always reset after send\n"
    " -ResetOnRequest : reset once after reply\n"
    " -MaxObjects : merge after this many objects are in (default 1)\n"
    " -reset : reset NOW\n"
    " -select : set the selection regex for sending out objects,\n" 
    "           valid for one reply if used in a request,\n"
    ;

void* work(void* /*param*/)
{
  return NULL;
}

//_______________________________________________________________________________________
long GetMilliSecSince(TTimeStamp* last)
{
  TTimeStamp now;
  long s = now.GetSec() - last->GetSec();
  long ms = now.GetNanoSec()/1000000 - last->GetNanoSec()/1000000;
  return s*1000+ms; 
}

//_______________________________________________________________________________________
Int_t Run()
{
  Int_t rc = 0;
  Int_t nSockets=3;
  zmq_pollitem_t sockets[] = { 
    { fZMQin, 0, ZMQ_POLLIN, 0 },
    { fZMQout, 0, ZMQ_POLLIN, 0 },
    { fZMQmon, 0, ZMQ_POLLIN, 0 },
  };
  //main loop
  while(1)
  {
    errno=0;

    Int_t inType=alizmq_socket_type(fZMQin);
    Int_t outType=alizmq_socket_type(fZMQout);
    Int_t monType=alizmq_socket_type(fZMQmon);
    
    //request first
    if (inType==ZMQ_REQ) DoRequest(fZMQin);
    if (outType==ZMQ_REQ) DoRequest(fZMQout);
    if (monType==ZMQ_REQ) DoRequest(fZMQmon);

    //wait for the data
    //poll sockets - we want to take action on one of two conditions:
    //  1 - request comes in - then we merge whatever is not yet merged and send
    //  2 - data comes in - then we add it to the merging list
    rc = zmq_poll(sockets, nSockets, -1); //poll sockets
    if (rc==-1 && errno==ETERM)
    {
      //this can only happen it the context was terminated, one of the sockets are
      //not valid or operation was interrupted
      Printf("zmq_poll was interrupted! rc = %i, %s", rc, zmq_strerror(errno));
      break;
    }

    //data present socket 0 - in
    if (sockets[0].revents & ZMQ_POLLIN)
    {
      aliZMQmsg message;
      alizmq_msg_recv(&message, fZMQin, 0);
      for (aliZMQmsg::iterator i=message.begin(); i!=message.end(); ++i)
      {
        if (alizmq_socket_type(fZMQin)==ZMQ_REP) 
        { HandleRequest(i->first, i->second, fZMQin); }
        else
        { HandleDataIn(i->first, i->second, fZMQout); }
      }
      alizmq_msg_close(&message);
    } //socket 0

    //data present socket 1 - out
    if (sockets[1].revents & ZMQ_POLLIN)
    {
      aliZMQmsg message;
      alizmq_msg_recv(&message, fZMQout, 0);
      for (aliZMQmsg::iterator i=message.begin(); i!=message.end(); ++i)
      {
        if (alizmq_socket_type(fZMQout)==ZMQ_REP) 
        { HandleRequest(i->first, i->second, fZMQout); }
        else
        { HandleDataIn(i->first, i->second, fZMQin); }
      }
      alizmq_msg_close(&message);
    }//socket 1
    
    //data present socket 1 - mon
    if (sockets[2].revents & ZMQ_POLLIN)
    {
      aliZMQmsg message;
      alizmq_msg_recv(&message, fZMQmon, 0);
      for (aliZMQmsg::iterator i=message.begin(); i!=message.end(); ++i)
      {
        if (alizmq_socket_type(fZMQmon)==ZMQ_REP) 
        { HandleRequest(i->first, i->second, fZMQmon); }
        else
        { HandleDataIn(i->first, i->second, fZMQmon); }
      }
      alizmq_msg_close(&message);
    }//socket 1
    
  }//main loop

  return 0;
}

//_____________________________________________________________________
Int_t HandleControlMessage(zmq_msg_t* topicMsg, zmq_msg_t* dataMsg, void* socket)
{
  string tmp;
  tmp.assign((char*)zmq_msg_data(topicMsg),zmq_msg_size(topicMsg));
  if (strncmp((char*)zmq_msg_data(topicMsg),"CONFIG",6)==0)
  {
    //reconfigure (first send a reply to not cause problems on the other end)
    std::string requestBody;
    requestBody.assign(static_cast<char*>(zmq_msg_data(dataMsg)), zmq_msg_size(dataMsg));

    //std::string reply = "Reconfiguring...";
    //zmq_send(socket, "INFO", 4, ZMQ_SNDMORE);
    //zmq_send(socket, reply.c_str(), reply.size(), 0);
    ProcessOptionString(requestBody.c_str());
    return 1;
  }
  else if (strncmp((char*)zmq_msg_data(topicMsg),"INFO",4)==0)
  {
    //do nothing, maybe log, send back an empty info reply
    //zmq_send(socket, "INFO", 4, ZMQ_SNDMORE);
    //zmq_send(socket, 0, 0, 0);
    return 1;
  }
  else
    return 0;
}

//_____________________________________________________________________
Int_t HandleRequest(zmq_msg_t* topicMsg, zmq_msg_t* dataMsg, void* socket)
{
  if (HandleControlMessage(topicMsg, dataMsg, socket)>0) return 0;
  return DoReply(topicMsg, dataMsg, socket);
}

//_____________________________________________________________________
Int_t HandleDataIn(zmq_msg_t* topicMsg, zmq_msg_t* dataMsg, void* socket)
{
  if (HandleControlMessage(topicMsg, dataMsg, socket)>0) return 0;
  return DoReceive(topicMsg, dataMsg, socket);
}

//_____________________________________________________________________
Int_t DoReply(zmq_msg_t* topicMsg, zmq_msg_t* dataMsg, void* socket)
{
  int rc = DoSend(socket);

  //reset the "one shot" options to default values
  fResetOnRequest = kFALSE;
  if (fVerbose && fSendSelection) 
    Printf("unsetting fSendSelection=%s",fSendSelection->GetPattern().Data());
  delete fSendSelection; fSendSelection=NULL;
  return rc;
}

//_____________________________________________________________________
Int_t DoReceive(zmq_msg_t* topicMsg, zmq_msg_t* dataMsg, void* socket)
{
  //handle the message
  //add to the list of objects to merge for each object type (by name)
  //topic
  AliHLTDataTopic dataTopic;
  memcpy(&dataTopic, zmq_msg_data(topicMsg),std::min(zmq_msg_size(topicMsg),sizeof(dataTopic)));
  if (fVerbose) Printf("in: data: %s, size: %zu bytes", dataTopic.Description().c_str(), zmq_msg_size(dataMsg));
  TObject* object = UnpackMessage(dataMsg );
  if (object)
  {
    const char* name = object->GetName();
    TList* mergingList = static_cast<TList*>(fMergeListMap.GetValue(name));
    TObject* mergingObject = fMergeObjectMap.GetValue(name);
    if (!mergingObject)
    {
      if (fVerbose) Printf("adding %s to fMergeObjectMap as first instance", name);
      fMergeObjectMap.Add(new TObjString(name), object);
    }
    else if (!mergingList) 
    {
      if (fVerbose) Printf("adding a new list %s to fMergeObjectMap", name);
      mergingList = new TList();
      mergingList->SetOwner();
      fMergeListMap.Add(new TObjString(name), mergingList);
    }
    else
    {
      //add object and maybe merge
      mergingList->Add(object);

      if (mergingList->GetEntries() >= fMaxObjects)
      {
        if (fVerbose) Printf("%i %s's in, merging",mergingList->GetEntries(),name);
        Merge(mergingObject, mergingList);

      }
    }
  }
  else
  {
    if (fVerbose) Printf("no object!");
  }
  
  if (fPushbackPeriod>=0)
  {
    if (GetMilliSecSince(&fLastPushBackTime)>fPushbackPeriod)
      DoSend(socket);
  }

  return 0;
}

//______________________________________________________________________________
Int_t DoRequest(void*)
{
  return 0;
}

//______________________________________________________________________________
Int_t DoSend(void* socket)
{
  //send back merged data, one object per frame

  aliZMQmsg message;
  Int_t rc = 0;
  TObject* object = NULL;
  TObject* key = NULL;
  
  TIter mapIter(&fMergeObjectMap);
  while ((key = mapIter.Next()))
  {
    //the topic
    AliHLTDataTopic topic = kAliHLTDataTypeTObject;
    //the data
    object = fMergeObjectMap.GetValue(key);

    const char* objectName = object->GetName();
    if (fSendSelection && !fSendSelection->Match(objectName)) 
    {
      if (fVerbose) Printf("     object %s did NOT make the selection %s", 
                           objectName, fSendSelection->GetPattern().Data());
      continue;
    }

    rc = alizmq_msg_add(&message, &topic, object);
    if (fResetOnSend || fResetOnRequest) 
    {
      TPair* pair = fMergeObjectMap.RemoveEntry(key);
      delete pair->Key();
      delete pair->Value();
      delete pair;
    }
  }

  //send
  int sentBytes = alizmq_msg_send(&message, socket, 0);
  if (fVerbose) Printf("merger sent %i bytes", sentBytes);
  alizmq_msg_close(&message);

  //always at least send an empty reply if we are replying
  if (sentBytes==0 && alizmq_socket_type(socket)==ZMQ_REP)
    alizmq_msg_send("INFO","NODATA",socket,0);

  fLastPushBackTime.Set();

  return 0;
}

//______________________________________________________________________________
void ResetOutputData()
{
  fMergeObjectMap.DeleteAll();
}

//_______________________________________________________________________________________
Int_t InitZMQ()
{
  //init or reinit stuff
  Int_t rc = 0;
  rc = alizmq_socket_init(fZMQin,  fZMQcontext, fZMQconfigIN.Data(), fZMQtimeout, fZMQmaxQueueSize);
  printf("in:  (%s) %s\n", alizmq_socket_name(rc), fZMQconfigIN.Data());
  rc = alizmq_socket_init(fZMQout, fZMQcontext, fZMQconfigOUT.Data(), fZMQtimeout, fZMQmaxQueueSize);
  printf("out: (%s) %s\n", alizmq_socket_name(rc), fZMQconfigOUT.Data());
  rc = alizmq_socket_init(fZMQmon, fZMQcontext, fZMQconfigMON.Data(), fZMQtimeout, fZMQmaxQueueSize);
  printf("mon: (%s) %s\n", alizmq_socket_name(rc) , fZMQconfigMON.Data());
  return 0;
}

//_______________________________________________________________________________________
TObject* UnpackMessage(zmq_msg_t* message)
{
  size_t size = zmq_msg_size(message);
  void* data = zmq_msg_data(message);

  TObject* object = NULL;
  object = AliHLTMessage::Extract(data, size);
  return object;
}

//_______________________________________________________________________________________
Int_t Merge(TObject* object, TCollection* mergeList)
{
		if (!object->IsA()->GetMethodWithPrototype("Merge", "TCollection*"))
		{
			Printf("Object does not implement a merge function!");
			return(-1);
		}
		Int_t error = 0;
		TString listHargs;
		listHargs.Form("((TCollection*)0x%lx)", (ULong_t) mergeList);
    //Printf("listHargs: %s", listHargs.Data());
		object->Execute("Merge", listHargs.Data(), &error);
		if (error)
		{
			Printf("Error %i running merge!", error);
			return(-1);
		}
    mergeList->Delete();
  return 0;
}

//______________________________________________________________________________
Int_t ProcessOptionString(TString arguments)
{
  //process passed options
  Int_t nOptions=0;
  stringMap* options = AliOptionParser::TokenizeOptionString(arguments);
  for (stringMap::iterator i=options->begin(); i!=options->end(); ++i)
  {
    //Printf("  %s : %s", i->first.data(), i->second.data());
    const TString& option = i->first; 
    const TString& value = i->second;
    if (option.EqualTo("reset")) 
    {
      ResetOutputData();
    }
    else if (option.EqualTo("ResetOnRequest"))
    {
      fResetOnRequest = value.Contains("0")?kFALSE:kTRUE;
    }
    else if (option.EqualTo("ResetOnSend"))
    {
      fResetOnSend = value.Contains("0")?kFALSE:kTRUE;
    }
    else if (option.EqualTo("MaxObjects"))
    {
      fMaxObjects = value.Atoi();
    }
    else if (option.EqualTo("ZMQconfigIN") || option.EqualTo("in"))
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
    else if (option.EqualTo("Verbose"))
    {
      fVerbose=kTRUE;
    }
    else if (option.EqualTo("pushback-period"))
    {
      fPushbackPeriod=value.Atoi();
    }
    else if (option.EqualTo("ZMQmaxQueueSize"))
    {
      fZMQmaxQueueSize=value.Atoi();
    }
    else if (option.EqualTo("ZMQtimeout"))
    {
      fZMQtimeout=value.Atoi();
    }
    else if (option.EqualTo("select"))
    {
      delete fSendSelection;
      fSendSelection = new TPRegexp(value);
      if (fVerbose) Printf("setting new regex %s",fSendSelection->GetPattern().Data());
    }
    else
    {
      Printf("unrecognized option %s",option.Data());
      nOptions=-1;
      break;
    }
    nOptions++;
  }
  delete options; //tidy up

  return nOptions; 
}

//_______________________________________________________________________________________
int main(Int_t argc, char** argv)
{
  Int_t mainReturnCode=0;

  //process args
  TString argString = AliOptionParser::GetFullArgString(argc,argv);
  if (ProcessOptionString(argString)<=0)
  {
    printf("%s",fUSAGE);
    return 1;
  }

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

