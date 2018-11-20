#include "zmq.h"
#include <iostream>
#include "AliHLTDataTypes.h"
#include "AliHLTComponent.h"
#include "AliHLTMessage.h"
#include "TClass.h"
#include "TPRegexp.h"
#include "TObjString.h"
#include "TH1.h"
#include "TList.h"
#include "TMessage.h"
#include "TRint.h"
#include "TApplication.h"
#include <time.h>
#include <sys/time.h>
#include <string>
#include <map>
#include "AliHLTZMQhelpers.h"
#include "AliOptionParser.h"
#include "TCollection.h"
#include "AliLog.h"
#include "AliAnalysisDataContainer.h"
#include "AliAnalysisTask.h"
#include "TFile.h"
#include "TKey.h"
#include "TSystem.h"
#include "signal.h"
#include <unistd.h>
#include <pthread.h>
#include "TMethodCall.h"
#include "TMethod.h"
#include "AliMergeable.h"
#include "TROOT.h"
#include "TInterpreter.h"

//this is meant to become a class, hence the structure with global vars etc.
//Also the code is rather flat - it is a bit of a playground to test ideas.
//TODO structure this at some point, e.g. introduce a SIMPLE unified way of handling
//zmq payloads, maybe a AliZMQmessage class which would by default be multipart and provide
//easy access to payloads based on topic or so (a la HLT GetFirstInputObject() etc...)

using namespace AliZMQhelpers;

//methods
Int_t ProcessOptionString(TString arguments, Bool_t verbose=kFALSE);
Int_t InitZMQ();
Int_t Run();

Int_t HandleDataIn(aliZMQmsg::iterator block, void* /*socket*/=NULL);
Int_t HandleRequest(aliZMQmsg::iterator block, void* /*socket*/=NULL);

Int_t DoReceive(aliZMQmsg::iterator block, void* socket);
Int_t DoSend(void* socket);
Int_t DoReply(aliZMQmsg::iterator block, void* socket);
Int_t DoRequest(void*& socket, TString* config=NULL);
Int_t DoControl(aliZMQmsg::iterator block, void* socket);
Int_t DoProcessAllData(void* outputsocket);

//processor private functions
void SetLastPushBackTime();

//configuration vars
Bool_t  fVerbose = kFALSE;
TString fZMQconfigIN  = "PULL";
TString fZMQsubscriptionIN = "";
TString fZMQconfigOUT  = "PUSH";
TString fZMQconfigMON  = "REP";
TString fZMQconfigSYNC  = "PUB";
TString fOnResetSendTo = "";
Int_t   fZMQmaxQueueSize = 10;
Int_t   fZMQtimeout = -1;
std::string fInitFile = "";

Bool_t  fResetOnSend = kFALSE;      //reset on each send (also on scheduled pushing)
Bool_t  fResetOnRequest = kFALSE;   //reset once after a single request
Bool_t  fRequestResetOnRequest = kFALSE; //when requesting from another processor, request also a reset
Bool_t  fClearOnReset = kFALSE;

Bool_t  fAllowGlobalReset=kTRUE;
Bool_t  fAllowControlSequences=kTRUE;
Bool_t  fAllowResetOnRequest=kTRUE;
Bool_t  fAllowResetAtSOR=kTRUE;
Bool_t  fAllowClearAtSOR=kFALSE;

Bool_t  fUnpackCollections = kTRUE;
Bool_t  fUnpackContainers = kFALSE;
Bool_t  fFullyDestroyAnalysisDataContainer = kFALSE;
std::string fLoadLibs;

TPRegexp* fSendSelection = NULL;
TPRegexp* fUnSendSelection = NULL;
std::string fNameList = "";
TString fTitleAnnotation = "";
TString fNameAnnotation = "";
Bool_t fTitleAnnotationWithContainerName = kFALSE;
Bool_t fIgnoreDefaultNamesWhenUnpacking = kFALSE;

Int_t fRunNumber = 0;
std::string fInfo;           //cache for the info string

std::string fID;  //processor ID/name

//internal state
TObjArray fInputObjects(10);
TObjArray fOutputObjects(10);

std::string fMacroPath;
long fPushbackPeriod = -1;        //in seconds, -1 means never
long fRequestPeriod = -1;        //in seconds, -1 means never
long fRequestTimeout = 10000;    //default 10s
Bool_t fCacheOnly = kFALSE;
aliZMQrootStreamerInfo fSchema;
bool fSchemaOnRequest = false;
bool fSchemaOnSend = false;
int fCompression = 1;
bool fgTerminationSignaled=false;

unsigned long long fLastPushBackTime = 0;
unsigned long long fLastRequestTime = 0;
struct timeval fCurrentTimeStruct;

double fBytesIN = 0;
double fBytesOUT = 0;
unsigned long fSecondsStart = 0;
unsigned long fSecondsStop = 0;
unsigned long long fNumberOfMessagesReceived = 0;
unsigned long long fNumberOfMessagesSent = 0;

//ZMQ stuff
void* fZMQcontext = NULL;    //ze zmq context

void* fZMQmon = NULL;        //the request-reply socket, here we request the merged data
void* fZMQout = NULL;        //the monitoring socket, here we publish a copy of the data
void* fZMQin  = NULL;        //the in socket - entry point for the data to be merged.
void* fZMQsync = NULL;
void* fZMQresetBroadcast = NULL;
bool fReconfigureZMQ = true;

//request trigger
void* fZMQtrig = NULL; //dummy socket to trigger requests at constant time interval
std::string fZMQconfigTRIG = "PAIR@inproc://trigger";
pthread_t fTRIGthread;
void* runTRIGthread(void*);

const char* fUSAGE =
    "ZMQprocessor options:\n"
    "process input, send results on output\n"
    " -id : some string identifier\n"
    " -in : data in, zmq config string, e.g. PULL>tcp://localhost:123123\n"
    " -out : data out\n"
    " -mon : monitoring socket\n"
    " -sync : sync socket, will send the INFO block on run change, has to be PUB or SUB\n"
    " -Verbose : print some info\n"
    " -request-period : request period [ms] - limited by request-timeout\n"
    " -request-timeout : timeout for REQ socket [ms] after which socket is reinitialized\n"
    " -SchemaOnRequest : include streamers ONCE (after a request)\n"
    " -SchemaOnSend : include streamers ALWAYS in each sent message\n"
    " -loadlibs : load ROOT libs, comma separated list\n"
    " -macro : macro to run on input\n"
    "          has to define a function: int process(TCollection* in, TCollection out);\n"
    "          first argument contains the ROOT objects present in the input,\n"
    "          second argument is for registering output objects, takes ownership of the data\n"
    ;
//_______________________________________________________________________________________
void sig_handler(int signo)
{
  if (signo == SIGINT)
    printf(" received SIGINT\n");
  fgTerminationSignaled=true;
}

//_______________________________________________________________________________________
Int_t Run()
{
  //start time (sec)
  gettimeofday(&fCurrentTimeStruct, NULL);
  fSecondsStart = fCurrentTimeStruct.tv_sec;

  //main loop
  while(!fgTerminationSignaled)
  {
    Int_t nSockets=5;
    zmq_pollitem_t sockets[] = {
      { fZMQin, 0, ZMQ_POLLIN, 0 },
      { fZMQout, 0, ZMQ_POLLIN, 0 },
      { fZMQmon, 0, ZMQ_POLLIN, 0 },
      { fZMQsync, 0, ZMQ_POLLIN, 0 },
      { fZMQtrig, 0, ZMQ_POLLIN, 0 },
    };

    Int_t rc = 0;
    errno=0;

    Int_t inType=alizmq_socket_type(fZMQin);
    Int_t outType=alizmq_socket_type(fZMQout);
    Int_t monType=alizmq_socket_type(fZMQmon);
    Int_t syncType=alizmq_socket_type(fZMQsync);

    //wait for the data
    //poll sockets - we want to take action on one of two conditions:
    //  1 - request comes in - then we merge whatever is not yet merged and send
    //  2 - data comes in - then we add it to the merging list
    rc = zmq_poll(sockets, nSockets, fZMQtimeout); //poll sockets
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
      int pushBack = 0;
      aliZMQmsg message;
      fBytesIN += alizmq_msg_recv(&message, fZMQin, 0);
      fNumberOfMessagesReceived++;
      for (aliZMQmsg::iterator i=message.begin(); i!=message.end(); ++i)
      {
        if (alizmq_socket_type(fZMQin)==ZMQ_REP)
        { HandleRequest(i, fZMQin); }
        else
        { pushBack += HandleDataIn(i, fZMQout); }
      }
      alizmq_msg_close(&message);
      DoProcessAllData(fZMQout);

    } //socket 0

    //data present socket 1 - out
    if (sockets[1].revents & ZMQ_POLLIN)
    {
      int pushBack = 0;
      aliZMQmsg message;
      fBytesIN += alizmq_msg_recv(&message, fZMQout, 0);
      fNumberOfMessagesReceived++;
      for (aliZMQmsg::iterator i=message.begin(); i!=message.end(); ++i)
      {
        if (alizmq_socket_type(fZMQout)==ZMQ_REP)
        { HandleRequest(i, fZMQout); }
        else
        { pushBack += HandleDataIn(i, fZMQin); }
      }
      alizmq_msg_close(&message);
      DoProcessAllData(fZMQin);
    }//socket 1

    //data present socket 2 - mon
    if (sockets[2].revents & ZMQ_POLLIN)
    {
      int pushBack = 0;
      aliZMQmsg message;
      fBytesIN += alizmq_msg_recv(&message, fZMQmon, 0);
      fNumberOfMessagesReceived++;
      for (aliZMQmsg::iterator i=message.begin(); i!=message.end(); ++i)
      {
        if (alizmq_socket_type(fZMQmon)==ZMQ_REP)
        { HandleRequest(i, fZMQmon); }
        else
        { pushBack += HandleDataIn(i, fZMQmon); }
      }
      alizmq_msg_close(&message);
      DoProcessAllData(fZMQmon);
    }//socket 2

    //data present socket 3 - sync
    if (sockets[3].revents & ZMQ_POLLIN)
    {
      aliZMQmsg message;
      fBytesIN += alizmq_msg_recv(&message, fZMQsync, 0);
      fNumberOfMessagesReceived++;
      for (aliZMQmsg::iterator i=message.begin(); i!=message.end(); ++i)
      {
        HandleDataIn(i, fZMQsync);
      }
      alizmq_msg_close(&message);
      DoProcessAllData(fZMQsync);
    }//socket 3

    //data present socket 4 - trig
    if (sockets[4].revents & ZMQ_POLLIN)
    {
      zmq_recv(fZMQtrig, 0, 0, 0);
      if (fVerbose) printf("trigger\n");
      // first
      if (inType==ZMQ_REQ) DoRequest(fZMQin, &fZMQconfigIN);
      if (outType==ZMQ_REQ) DoRequest(fZMQout, &fZMQconfigOUT);
      if (monType==ZMQ_REQ) DoRequest(fZMQmon, &fZMQconfigMON);

      //this might trigger internal cleanup
      alizmq_socket_state(fZMQin);
      alizmq_socket_state(fZMQout);
      alizmq_socket_state(fZMQmon);
      alizmq_socket_state(fZMQsync);
    }//socket 3


  }//main loop

  //stop time (sec)
  gettimeofday(&fCurrentTimeStruct, NULL);
  fSecondsStop = fCurrentTimeStruct.tv_sec;

  {
    printf("number of messages received: %llu\n", fNumberOfMessagesReceived);
    printf("number of messages sent    : %llu\n", fNumberOfMessagesSent);
    printf("in  %.2f MB, %.2f MB/s\n", (double)(fBytesIN/1024./1024.),  (float)(fBytesIN/(fSecondsStop-fSecondsStart)/1024./1024.));
    printf("out %.2f MB, %.2f MB/s\n", (double)(fBytesOUT/1024./1024.), (float)(fBytesOUT/(fSecondsStop-fSecondsStart)/1024./1024.));
  }

  return 0;
}

//_____________________________________________________________________
Int_t DoControl(aliZMQmsg::iterator block, void* socket)
{
  AliHLTDataTopic topic;
  alizmq_msg_iter_topic(block, topic);

  if (topic.GetID().compare(0,kAliHLTComponentDataTypefIDsize,kAliHLTDataTypeStreamerInfo.fID,kAliHLTComponentDataTypefIDsize)==0)
  {
    //extract the streamer infos
    if (fVerbose) printf("unpacking ROOT streamer infos... %s\n", topic.GetID().c_str());
    alizmq_msg_iter_init_streamer_infos(block);
    return 1;
  }
  else if (fAllowControlSequences && topic.GetID().compare(0,6,"CONFIG")==0)
  {
    //reconfigure (first send a reply to not cause problems on the other end)
    std::string requestBody;
    alizmq_msg_iter_data(block, requestBody);
    if (fVerbose) printf("received CONFIG %s\n", requestBody.c_str());
    ProcessOptionString(requestBody.c_str());
    return 1;
  }
  else if (topic.GetID().compare(0,4,"INFO")==0)
  {
    //check if we have a runnumber in the string
    alizmq_msg_iter_data(block, fInfo);
    int runnumber = atoi(GetParamString("run",fInfo).c_str());

    if (fVerbose) printf("received run=%i\n",runnumber);

    return 1;
  }

  return 0;
}

//_____________________________________________________________________
Int_t HandleRequest(aliZMQmsg::iterator block, void* socket)
{
  DoControl(block, socket);
  return DoReply(block, socket);
}

//_____________________________________________________________________
Int_t HandleDataIn(aliZMQmsg::iterator block, void* socket)
{
  if (DoControl(block, socket)>0) return 0;
  return DoReceive(block, socket);
}

//_____________________________________________________________________
Int_t DoReply(aliZMQmsg::iterator block, void* socket)
{
  if (fVerbose) printf("replying!\n");
  int rc = DoSend(socket);
  return rc;
}

//_____________________________________________________________________
Int_t DoProcessAllData(void* outputsocket)
{
  if (!fMacroPath.empty() && fInputObjects.GetEntries()>0)
  {
    if (fVerbose) {printf("processing...\n");}
    Int_t error = 0;
    gROOT->ProcessLineFast(Form("process((TCollection*)%p,(TCollection*)%p)",&fInputObjects,&fOutputObjects),&error);
    if (error != TInterpreter::kNoError) {printf("error %i executing process(...)\n");}
    DoSend(outputsocket);
    fOutputObjects.Delete(); //Delete output objects, they are already serialized.
    fInputObjects.Delete(); //Delete input objects after use
  }
  return 0;
}

//_____________________________________________________________________
Int_t DoReceive(aliZMQmsg::iterator block, void* socket)
{
  //handle the message
  //add to the list of objects to merge for each object type (by name)
  //topic
  AliHLTDataTopic dataTopic;
  alizmq_msg_iter_topic(block, dataTopic);

  if (fVerbose) Printf("in: data: %s, size: %zu bytes", dataTopic.Description().c_str(),
                       zmq_msg_size(block->second));
  TObject* object = NULL;
  alizmq_msg_iter_data_hlt(block, object);
  if (!object) {
    if (fVerbose) {
      printf("message has no TObject!\n");
    }
    return 0;
  }

  //if we get a collection, always set ownership to prevent mem leaks
  //if we request unpacking: unpack what was requestd, otherwise just add
  do { //fake loop for easy breaking
    fInputObjects.AddLast(object);
  } while (false);

  if (fPushbackPeriod>=0)
  {
    gettimeofday(&fCurrentTimeStruct, NULL);
    unsigned long long currentTime = 1000*fCurrentTimeStruct.tv_sec+fCurrentTimeStruct.tv_usec/1000;

    unsigned long timeInterval = currentTime - fLastPushBackTime;
    if (timeInterval >= fPushbackPeriod)
    {
      return 1; //signal we will want to send after message is done
    }

  }

  return 0;
}

//______________________________________________________________________________
Int_t DoRequest(void*& socket, TString* config)
{
  //if we cannot send it means we were waiting for a reply.
  //after enough time passed we assume server died and reinit the socket
  gettimeofday(&fCurrentTimeStruct, NULL);
  unsigned long long currentTime = 1000*fCurrentTimeStruct.tv_sec +
                                        fCurrentTimeStruct.tv_usec/1000;
  unsigned long long sinceLastRequest = currentTime - fLastRequestTime;

  if ((alizmq_socket_state(socket)&ZMQ_POLLOUT)==0) {

    //if we are still waiting, do nothing
    if (fLastRequestTime>0 && sinceLastRequest < fRequestTimeout) {
      if (fVerbose) printf("still waiting for a reply....\n");
      return -1;
    }

    //if we are not still waiting and socket is in wrong state, we have to reinit the socket
    if (fVerbose) printf("no reply from %s in %li ms, server died? - reinit socket\n",
        config->Data(), fRequestTimeout);
    alizmq_socket_init(socket, fZMQcontext, config->Data(), fZMQtimeout, fZMQmaxQueueSize);
  }

  //request
  if (fVerbose) Printf("sending an empty request");
  alizmq_msg_send("", "", socket, 0);
  fNumberOfMessagesSent++;

  fLastRequestTime = currentTime;
  return 0;
}

//______________________________________________________________________________
Int_t DoSend(void* socket)
{
  //only send if we actually CAN send
  if ( !(alizmq_socket_state(socket) & ZMQ_POLLOUT) ) {
    if (fVerbose) printf("cannot send to %s, socket in wrong state\n", alizmq_socket_name(alizmq_socket_type(socket)));
    return 0;
  }

  //send back merged data, one object per frame

  DataTopic fInfoTopic = kDataTypeInfo;
  fInfoTopic.SetOrigin(kAliHLTDataOriginOut);

  aliZMQmsg message;
  //forward the (run-)info string
  alizmq_msg_add(&message, &fInfoTopic, fInfo);
  fInfo.clear();
  Bool_t reset = kFALSE;
  Int_t rc = 0;
  int parts = 1;

  TIter next(&fOutputObjects);
  while (TObject* object = next())
  {
    //the topic
    DataTopic topic = kDataTypeTObject;
    topic.SetOrigin(kAliHLTDataOriginOut);
    //the data
    string objectName = object->GetName();

    //add data to be sent
    rc = alizmq_msg_add(&message, &topic, object, fCompression, NULL);
    parts++;
  }

  if ((fSchemaOnRequest || fSchemaOnSend)) {
    alizmq_msg_prepend_streamer_infos(&message, &fSchema);
    parts++;
  }

  int sentBytes = 0;
  //send a copy to the broadcast socket if requested on reset
  if (reset && fZMQresetBroadcast && fZMQresetBroadcast!=socket) {
    aliZMQmsg messageCopy;
    alizmq_msg_copy(&messageCopy, &message);
    sentBytes += alizmq_msg_send(&messageCopy, fZMQresetBroadcast, 0);
    fNumberOfMessagesSent++;
    if (fVerbose) printf("a copy was broadcasted, %i bytes\n", sentBytes);
    alizmq_msg_close(&messageCopy);
  }

  //send
  sentBytes += alizmq_msg_send(&message, socket, 0);
  fNumberOfMessagesSent++;
  fBytesOUT += sentBytes;
  if (fVerbose) Printf("processor sent %i bytes, %i blocks, on socket %p", sentBytes, parts, socket);
  alizmq_msg_close(&message);

  SetLastPushBackTime();

  return sentBytes;
}

//______________________________________________________________________________
void SetLastPushBackTime()
{
  gettimeofday(&fCurrentTimeStruct, NULL);
  fLastPushBackTime = 1000*fCurrentTimeStruct.tv_sec+fCurrentTimeStruct.tv_usec/1000;
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
  rc = alizmq_socket_init(fZMQsync, fZMQcontext, fZMQconfigSYNC.Data(), fZMQtimeout, fZMQmaxQueueSize);
  printf("sync: (%s) %s\n", alizmq_socket_name(rc) , fZMQconfigSYNC.Data());
  rc = alizmq_socket_init(fZMQtrig, fZMQcontext, fZMQconfigTRIG, -1, 1);
  printf("trigger: (%s) %s\n", alizmq_socket_name(rc) , fZMQconfigTRIG.c_str());
  return 0;
}

//______________________________________________________________________________
Int_t ProcessOptionString(TString arguments, Bool_t verbose)
{
  //process passed options
  Int_t nOptions=0;
  aliStringVec* options = AliOptionParser::TokenizeOptionString(arguments);
  for (aliStringVec::iterator i=options->begin(); i!=options->end(); ++i)
  {
    const TString& option = i->first;
    const TString& value = i->second;
    if (option.EqualTo("ResetOnRequest"))
    {
      fResetOnRequest = value.Contains("0")?kFALSE:kTRUE;
    }
    else if (option.EqualTo("RequestResetOnRequest"))
    {
      fRequestResetOnRequest = value.Contains("0")?kFALSE:kTRUE;
    }
    else if (option.EqualTo("ResetOnSend"))
    {
      fResetOnSend = value.Contains("0")?kFALSE:kTRUE;
    }
    else if (option.EqualTo("ClearOnReset"))
    {
      fClearOnReset = value.Contains("0")?kFALSE:kTRUE;
    }
    else if (option.EqualTo("OnResetBroadcastTo"))
    {
      if (value.EqualTo("in")||value.EqualTo("out")||value.EqualTo("mon")) {
        fOnResetSendTo = value;
      } else {
        printf("OnResetBroadcastTo must be [in|out|mon]\n");
        return -1;
      }
    }
    else if (option.EqualTo("ZMQconfigIN") || option.EqualTo("in"))
    {
      fZMQconfigIN = value;
      fReconfigureZMQ = true;
    }
    else if (option.EqualTo("ZMQconfigOUT") || option.EqualTo("out"))
    {
      fZMQconfigOUT = value;
      fReconfigureZMQ = true;
    }
    else if (option.EqualTo("ZMQconfigMON") || option.EqualTo("mon"))
    {
      fZMQconfigMON = value;
      fReconfigureZMQ = true;
    }
    else if (option.EqualTo("ZMQconfigSYNC") || option.EqualTo("sync"))
    {
      int type = alizmq_socket_type(value.Data());
      if (type==ZMQ_PUB || type==ZMQ_SUB) {
        fZMQconfigSYNC = value;
      } else {
        printf("sync socket has to be PUB or SUB!\n");
        return -1;
      }
      fReconfigureZMQ = true;
    }
    else if (option.EqualTo("Verbose"))
    {
      fVerbose=value.EqualTo("0")?kFALSE:kTRUE;
    }
    else if (option.EqualTo("pushback-period"))
    {
      fPushbackPeriod=value.Atoi();
    }
    else if (option.EqualTo("ZMQmaxQueueSize"))
    {
      fZMQmaxQueueSize=value.Atoi();
      fReconfigureZMQ = true;
    }
    else if (option.EqualTo("ZMQtimeout"))
    {
      fZMQtimeout=value.Atoi();
      fReconfigureZMQ = true;
    }
    else if (option.EqualTo("request-period"))
    {
      fRequestPeriod=value.Atoi();
    }
    else if (option.EqualTo("request-timeout"))
    {
      fRequestTimeout=value.Atoi();
    }
    else if (option.EqualTo("select"))
    {
      delete fSendSelection;
      fSendSelection = new TPRegexp(value);
      if (fVerbose) Printf("setting new regex %s",fSendSelection->GetPattern().Data());
    }
    else if (option.EqualTo("unselect"))
    {
      delete fUnSendSelection;
      fUnSendSelection = new TPRegexp(value);
      if (fVerbose) Printf("setting new regex %s",fUnSendSelection->GetPattern().Data());
    }
    else if (option.EqualTo("list"))
    {
      fNameList = value.Data();
      if (fVerbose) Printf("setting a selection list %s", fNameList.c_str());
    }
    else if (option.EqualTo("cache"))
    {
      fCacheOnly = kTRUE;
    }
    else if (option.EqualTo("annotateTitle"))
    {
      fTitleAnnotation = value;
    }
    else if (option.EqualTo("annotateName"))
    {
      fNameAnnotation = value;
    }
    else if (option.EqualTo("annotateTitleWithContainerName"))
    {
      fTitleAnnotationWithContainerName = value.Contains("0")?kFALSE:kTRUE;
    }
    else if (option.EqualTo("AllowGlobalReset"))
    {
      fAllowGlobalReset=(value.Contains("0")||value.Contains("no"))?kFALSE:kTRUE;
    }
    else if (option.EqualTo("AllowControlSequences"))
    {
      fAllowControlSequences = (value.Contains("0")||value.Contains("no"))?kFALSE:kTRUE;
    }
    else if (option.EqualTo("AllowResetOnRequest"))
    {
      fAllowResetOnRequest = (value.Contains("0")||value.Contains("no"))?kFALSE:kTRUE;
    }
    else if (option.EqualTo("AllowResetAtSOR"))
    {
      fAllowResetAtSOR = (value.Contains("0")||value.Contains("no"))?kFALSE:kTRUE;
    }
    else if (option.EqualTo("AllowClearAtSOR"))
    {
      fAllowClearAtSOR = (value.Contains("0")||value.Contains("no"))?kFALSE:kTRUE;
    }
    else if (option.EqualTo("SchemaOnRequest"))
    {
      fSchemaOnRequest = true;
    }
    else if (option.EqualTo("SchemaOnSend"))
    {
      fSchemaOnSend = (value.Contains("0"))?false:true;
    }
    else if (option.EqualTo("UnpackCollections"))
    {
      fUnpackCollections = (value.Contains("0") || value.Contains("no"))?kFALSE:kTRUE;
    }
    else if (option.EqualTo("UnpackContainers"))
    {
      fUnpackContainers = (value.Contains("0") || value.Contains("no"))?kFALSE:kTRUE;
    }
    else if (option.EqualTo("IgnoreDefaultContainerNames"))
    {
      fIgnoreDefaultNamesWhenUnpacking = (value.Contains("0"))?kFALSE:kTRUE;
    }
    else if (option.EqualTo("FullyDestroyAnalysisDataContainer"))
    {
      fFullyDestroyAnalysisDataContainer = (value.Contains("0"))?kFALSE:kTRUE;
    }
    else if (option.EqualTo("statefile"))
    {
      fInitFile = value.Data();
    }
    else if (option.EqualTo("id"))
    {
      fID = value.Data();
    }
    else if (option.EqualTo("loadlibs"))
    {
      fLoadLibs = value.Data();
    }
    else if (option.EqualTo("macro"))
    {
      fMacroPath = value.Data();
      std::string tmp = fMacroPath+"++";
      if (0>gROOT->LoadMacro(tmp.c_str())) {
        printf("ERROR: could not load macro %s\n",fMacroPath.c_str());
      };
    }
    else
    {
      Printf("unrecognized option |%s|",option.Data());
      nOptions=-1;
      break;
    }
    nOptions++;
  }

  if (nOptions<1) fReconfigureZMQ=false;
  if (fReconfigureZMQ && (InitZMQ()<0)) {
    Printf("failed ZMQ init");
    return -1;
  }
  fReconfigureZMQ=false;

  if (!fLoadLibs.empty()) {
    if (LoadROOTlibs(fLoadLibs,verbose)<0) {
      Printf("problem loading libraries %s",fLoadLibs.c_str());
      nOptions=-1;
    }
  }

  if (fOnResetSendTo.IsNull()) fZMQresetBroadcast = NULL;
  else if (fOnResetSendTo.EqualTo("in")) fZMQresetBroadcast = fZMQin;
  else if (fOnResetSendTo.EqualTo("mon")) fZMQresetBroadcast = fZMQmon;
  else if (fOnResetSendTo.EqualTo("out")) fZMQresetBroadcast = fZMQout;

  delete options; //tidy up

  if (verbose)
  {
    if (fRequestTimeout<100) printf("WARNING: setting the socket timeout to %lu ms can be dagerous,\n"
        "         choose something more realistic or leave the default as it is\n", fRequestTimeout);

    if (fZMQresetBroadcast) printf("configured to bradcast resets on %s, socket %p\n", fOnResetSendTo.Data(), fZMQresetBroadcast);
    if (fFullyDestroyAnalysisDataContainer) printf("configured to delete the fProducer/fConsumers of AliAnalysisDataContainer\n");
    if (fIgnoreDefaultNamesWhenUnpacking) printf("ignoring default cont names\n");
    if (fUnpackContainers) printf("unpacking containers\n");
    if (fUnpackCollections) printf("unpacking collections\n");
  }

  return nOptions;
}

//_______________________________________________________________________________________
int main(Int_t argc, char** argv)
{
  Int_t mainReturnCode=0;

  //catch signals
  if (signal(SIGHUP, sig_handler) == SIG_ERR)
  printf("\ncan't catch SIGHUP\n");
  if (signal(SIGINT, sig_handler) == SIG_ERR)
  printf("\ncan't catch SIGINT\n");
  if (signal(SIGQUIT, sig_handler) == SIG_ERR)
  printf("\ncan't catch SIGQUIT\n");
  if (signal(SIGTERM, sig_handler) == SIG_ERR)
  printf("\ncan't catch SIGTERM\n");

  //the context
  fZMQcontext = alizmq_context();
  if (!fZMQcontext) {
    printf("could not init the ZMQ context\n");
    return 1;
  }

  //process args
  TString argString = AliOptionParser::GetFullArgString(argc,argv);
  if (ProcessOptionString(argString, kTRUE)<=0)
  {
    printf("%s",fUSAGE);
    return 1;
  }

  //switch off logging if not verbose
  if (!fVerbose) {
    AliLog::SetGlobalLogLevel(AliLog::kWarning);
  }

  if (fRequestPeriod>=0) {
    pthread_create(&fTRIGthread, NULL, runTRIGthread, NULL);
  }

  Run();

  //destroy ZMQ sockets
  zmq_close(fZMQmon);
  zmq_close(fZMQin);
  zmq_close(fZMQout);
  zmq_close(fZMQsync);
  zmq_close(fZMQtrig);
  zmq_ctx_destroy(fZMQcontext);
  return mainReturnCode;
}

//______________________________________________________________________________
void* runTRIGthread(void*)
{
  if (fVerbose) printf("starting trigger thread\n");
  int delay = fRequestPeriod;
  void* trig = NULL;
  alizmq_socket_init(trig, alizmq_context(), "PAIR+inproc://trigger");
  while (true) {
    if (delay < 0) break;

    zmq_send(trig,0,0,0);

    //poll to handle the delay and exit condition
    zmq_pollitem_t sockets[] = { { trig, 0, ZMQ_POLLIN, 0 } };
    int rc = zmq_poll(sockets, 1, delay);
    if (rc==-1 && errno==ETERM)
    {
      break;
    }
    if (sockets[0].revents & ZMQ_POLLIN) {
      zmq_recv(trig,&delay,sizeof(delay),0);
    }

  }
  alizmq_socket_close(trig, 0);
  if (fVerbose) printf("stopping trigger thread\n");
  return NULL;
}

