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
TString GetFullArgString(Int_t argc, char** argv);
Int_t ProcessOptionString(TString arguments);
stringMap* TokenizeOptionString(const TString str);
Int_t ProcessOption(TString option, TString value);
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
TString fZMQconfigIN  = "PUSH";
TString fZMQsubscriptionIN = "";
TString fZMQconfigOUT  = "PULL";
TString fZMQconfigMON  = "REP";

Bool_t  fResetOnSend = kFALSE;

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
  if (strncmp((char*)zmq_msg_data(topicMsg),"CONFIG",6))
  {
    //reconfigure (first send a reply to not cause problems on the other end)
    std::string requestBody;
    requestBody.assign(static_cast<char*>(zmq_msg_data(dataMsg)), zmq_msg_size(dataMsg));

    std::string reply = "Reconfiguring...";
    zmq_send(socket, "INFO", 4, ZMQ_SNDMORE);
    zmq_send(socket, reply.c_str(), reply.size(), 0);
    ProcessOptionString(requestBody.c_str());
    return 1;
  }
  else if (strncmp((char*)zmq_msg_data(topicMsg),"INFO",4))
  {
    //do nothing, maybe log, send back an empty info reply
    zmq_send(socket, "INFO", 4, ZMQ_SNDMORE);
    zmq_send(socket, 0, 0, 0);
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
  return DoSend(socket);
}

//_____________________________________________________________________
Int_t DoReceive(zmq_msg_t* topicMsg, zmq_msg_t* dataMsg, void* socket)
{
  //handle the message
  //add to the list of objects to merge for each object type (by name)
  //topic
  AliHLTDataTopic dataTopic;
  memcpy(&dataTopic, zmq_msg_data(topicMsg),std::min(zmq_msg_size(topicMsg),sizeof(dataTopic)));
  if (fVerbose) Printf("in: data: %s", dataTopic.Description().c_str());
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

  Int_t rc = 0;
  TIter mapIter(&fMergeObjectMap);
  TObject* object = NULL;
  Int_t objectNumber=0;
  while ((object = mapIter.Next()))
  {
    //the topic
    AliHLTDataTopic topic = kAliHLTDataTypeTObject;
    //the data
    
    Int_t flags = ( objectNumber < fMergeObjectMap.GetEntries()-1 ) ? ZMQ_SNDMORE : 0;    
    rc = alizmq_msg_send(topic, fMergeObjectMap.GetValue(object), socket, flags, 0);
    if (rc<0)
    {
      Printf("could not send object");
      break;
    }
    objectNumber++;

    if (fResetOnSend)
    {
      ResetOutputData();
      if (fVerbose) Printf("data reset on send");
    }
  }
  //always at least send an empty reply if we are replying
  if (objectNumber==0 && alizmq_socket_type(socket)==ZMQ_REP)
    alizmq_msg_send("INFO","NODATA",socket,0);

  fLastPushBackTime.Set();

  return 0;
}

//______________________________________________________________________________
void ResetOutputData()
{
  fMergeObjectMap.DeleteAll();
}

//______________________________________________________________________________
Int_t ProcessOption(TString option, TString value)
{
  //process option
  //to be implemented by the user
  
  //if (option.EqualTo("ZMQpollIn"))
  //{
  //  fZMQpollIn = (value.EqualTo("0"))?kFALSE:kTRUE;
  //}
 
  if (option.EqualTo("reset")) 
  {
    ResetOutputData();
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

  return 1; 
}

//_______________________________________________________________________________________
Int_t InitZMQ()
{
  //init or reinit stuff
  Int_t rc = 0;
  printf("in:  ");
  rc += alizmq_socket_init(fZMQin,  fZMQcontext, fZMQconfigIN.Data(), 0, 200);
  printf("out: ");
  rc += alizmq_socket_init(fZMQout, fZMQcontext, fZMQconfigOUT.Data(), 0, 200);
  printf("mon: ");
  rc += alizmq_socket_init(fZMQmon, fZMQcontext, fZMQconfigMON.Data(), 0, 200);
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

////////////////////////////////////////////////////////////////////////////////
//_______________________________________________________________________________________
TString GetFullArgString(Int_t argc, char** argv)
{
  TString argString;
  TString argument="";
  if (argc>0) {
    for (Int_t i=1; i<argc; i++) {
      argument=argv[i];
      if (argument.IsNull()) continue;
      if (!argString.IsNull()) argString+=" ";
      argString+=argument;
    }  
  }
  return argString;
}

//______________________________________________________________________________
Int_t ProcessOptionString(TString arguments)
{
  //process passed options
  Int_t nOptions=0;
  stringMap* options = TokenizeOptionString(arguments);
  for (stringMap::iterator i=options->begin(); i!=options->end(); ++i)
  {
    //Printf("  %s : %s", i->first.data(), i->second.data());
    if (ProcessOption(i->first,i->second)>0)
      nOptions++;
    else
      break;
  }
  delete options; //tidy up

  return nOptions; 
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
  // options can be separated by ' ' or ',' arbitrarily combined, e.g:
  //"-option option1=value1 --option2 value2, -option4=\'some string\'"
  
  //optionRE by construction contains a pure option name as 3rd submatch (without --,-, =)
  //valueRE does NOT match options
  TPRegexp optionRE("(?:(-{1,2})|((?='?[^,=]+=?)))"
                    "((?(2)(?:(?(?=')'(?:[^'\\\\]++|\\.)*+'|[^, =]+))(?==?))"
                    "(?(1)[^, =]+(?=[= ,$])))");
  TPRegexp valueRE("(?(?!(-{1,2}|[^, =]+=))"
                   "(?(?=')'(?:[^'\\\\]++|\\.)*+'"
                   "|[^, =]+))");

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

//_______________________________________________________________________________________
int main(Int_t argc, char** argv)
{
  Int_t mainReturnCode=0;

  //process args
  TString argString = GetFullArgString(argc,argv);
  if (ProcessOptionString(argString)<=0)
  {
    printf("options: \n");
    printf("in : data in, zmq config string, e.g. PUSH>tcp://localhost:123123\n");
    printf("out : data out\n");
    printf("mon : monitoring socket\n");
    printf("Verbose : print some info\n");
    printf("pushback-period : push the merged data once every n ms\n");
    printf("ResetOnSend : always reset after send\n");
    printf("reset : reset NOW\n");
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

