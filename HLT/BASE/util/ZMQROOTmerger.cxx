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

//this is meant to become a class, hence the structure with global vars etc.
//Also the code is rather flat - it is a bit of a playground to test ideas.
//TODO structure this at some point, e.g. introduce a SIMPLE unified way of handling
//zmq payloads, maybe a AliZMQmessage class which would by default be multipart and provide
//easy access to payloads based on topic or so (a la HLT GetFirstInputObject() etc...)

class AliZMQROOTMergerEntry;
typedef std::map<std::string, AliZMQROOTMergerEntry> mergeMap_t;

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
Int_t GetObjects(AliAnalysisDataContainer* kont, std::vector<TObject*>* list, std::string prefix="");
Int_t GetObjects(AliMergeable* unpackable, std::vector<TObject*>* list, std::string prefix="");
Int_t GetObjects(TCollection* collection, std::vector<TObject*>* list, std::string prefix="");
Int_t ReadFromFile(std::string file);
Int_t DumpToFile(std::string file);

//merger private functions
int ResetOutputData(Bool_t force=kFALSE);
int ClearOutputData();
int RemoveEntry(string& name, mergeMap_t& map);
Int_t AddObject(TObject* object);
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
Bool_t  fRequestResetOnRequest = kFALSE; //when requesting from another merger, request also a reset
Bool_t  fClearOnReset = kFALSE;

Bool_t  fAllowGlobalReset=kTRUE;
Bool_t  fAllowControlSequences=kTRUE;
Bool_t  fAllowResetOnRequest=kTRUE;
Bool_t  fAllowResetAtSOR=kTRUE;
Bool_t  fAllowClearAtSOR=kFALSE;

Bool_t  fUnpackCollections = kFALSE;
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

AliHLTDataTopic fInfoTopic = kAliHLTDataTypeInfo;

Int_t fRunNumber = 0;
std::string fInfo;           //cache for the info string

std::string fID;  //merger ID/name

//internal state
mergeMap_t fMergeObjectMap;        //map of the merged objects, all incoming stuff is merged into these
Int_t fMaxObjects = 1;        //trigger merge after this many messages
std::vector<TObject*> fListOfObjects;

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
    "ZMQROOTmerger options: Merge() all ROOT mergeables in the message.\n"
    "merge based on what GetName() returns, the merged data can be retrieved at any time.\n"
    " -id : some string identifier\n"
    " -in : data in, zmq config string, e.g. PUSH>tcp://localhost:123123\n"
    " -out : data out\n"
    " -mon : monitoring socket\n"
    " -sync : sync socket, will send the INFO block on run change, has to be PUB or SUB\n"
    " -Verbose : print some info\n"
    " -pushback-period : push the merged data once every n milliseconds (if updated)\n"
    " -request-period : request period [ms] - limited by request-timeout\n"
    " -request-timeout : timeout for REQ socket [ms] after which socket is reinitialized\n"
    " -ResetOnSend : always reset after send\n"
    " -ResetOnRequest : reset once after reply\n"
    " -OnResetBroadcastTo : push (a copy) before resetting to socket [in|out|mon]\n"
    " -RequestResetOnRequest : if requesting form another merger, request a ResetOnRequest\n"
    " -ClearOnReset : instead of destroying data - clear it by calling TH1::Reset() on all histograms\n"
    " -AllowGlobalReset :  allow a global \'reset\' on request\n"
    " -AllowResetOnRequest : allow reset on request\n"
    " -AllowResetAtSOR : allow reset at change of run\n"
    " -AllowClearAtSOR : clear the histograms at change of run, works only if AllowResetAtSOR=0\n"
    " -AllowControlSequences : allow control seqs (CONFIG messages)\n"
    " -MaxObjects : merge after this many objects are in (default 1)\n"
    " -reset : reset NOW\n"
    " -select : set the selection regex for sending out objects,\n"
    "           valid for one reply if used in a request,\n"
    " -unselect : as above, only inverted\n"
    " -list : a list of (full) names to send (arb. delimiter)\n"
    " -cache : don't merge, only cache (i.e. replace)\n"
    " -annotateTitle : prepend string to title (if applicable)\n"
    " -annotateTitleWithContainerName : prepend the container name to title\n"
    " -annotateName : prepend string to the name\n"
    " -ZMQtimeout: when to timeout the sockets (milliseconds)\n"
    " -SchemaOnRequest : include streamers ONCE (after a request)\n"
    " -SchemaOnSend : include streamers ALWAYS in each sent message\n"
    " -UnpackCollections : cache/merge the contents of the collections instead of the collection itself\n"
    " -UnpackContainers : unpack the contents of AliAnalysisDataContainers\n"
    " -IgnoreDefaultContainerNames : don't prefix default container names (TList,TObjArray)\n"
    " -FullyDestroyAnalysisDataContainer : explicitly delete consumers and producer members in the container to work around a memory leak\n"
    " -loadlibs : load ROOT libs, comma separated list\n"
    " -statefile : save/restore state on exit/start\n"
    ;
//_______________________________________________________________________________________
void sig_handler(int signo)
{
  if (signo == SIGINT)
    printf(" received SIGINT\n");
  fgTerminationSignaled=true;
}

//_____________________________________________________________________
class AliZMQROOTMergerEntry {

private:
  TObject* fObject;
  TList fMergeList;
  TMethodCall* fMergeCall;

public:
  AliZMQROOTMergerEntry(): fObject(NULL), fMergeList(), fMergeCall(NULL)
  {
    fMergeList.SetOwner(kTRUE);
  }

  AliZMQROOTMergerEntry(TObject* o)
    :fObject(o)
    ,fMergeList()
    ,fMergeCall(NULL)
  {
    fMergeList.SetOwner(kTRUE);
  }

  AliZMQROOTMergerEntry(const AliZMQROOTMergerEntry& that)
    :fObject(that.fObject)
    ,fMergeList()
    ,fMergeCall(that.fMergeCall)
  {
    fMergeList.SetOwner(kTRUE);
  }

  AliZMQROOTMergerEntry& operator=(const AliZMQROOTMergerEntry& that)
  {
    delete fObject;
    fObject = that.fObject;
    return *this;
  }

  AliZMQROOTMergerEntry& operator=(TObject* o)
  {
    if (!fObject) {
      if (fVerbose) printf("adding first cache instance: %p %s\n",o,o->GetName());
      if (fVerbose) printf("  getting streamers for new object %s\n",o->GetName());
      alizmq_update_streamerlist_from_object(&fSchema, o);
    }
    delete fObject;
    fObject = o;
    return *this;
  }

  AliZMQROOTMergerEntry& operator+=(TObject* o)
  {
    Merge(o);
    return *this;
  }

  TObject* GetObject() {return fObject;}

  Int_t Merge(TObject* o)
  {
    //this method takes ownership of input
    if (!fObject) {
      if (fVerbose) printf("adding first instance: %p %s\n",o,o->GetName());
      if (fVerbose) printf("  getting streamers for new object %s\n",o->GetName());
      alizmq_update_streamerlist_from_object(&fSchema, o);
      fObject = o;
      return 0;
    }

    fMergeList.Add(o);

    int rc=0;
    TH1* hist = dynamic_cast<TH1*>(fObject);

    AliMergeable* mergeable = NULL;
    if (!hist)
    {
      mergeable = dynamic_cast<AliMergeable*>(fObject);
    }
    else if (!fMergeCall && fObject->IsA())
    {
      //use ROOT magic as last resort
      fMergeCall = new TMethodCall;
      fMergeCall->InitWithPrototype(fObject->IsA(), "Merge", "TCollection*");
      if (!fMergeCall->IsValid()) { delete fMergeCall; fMergeCall=NULL; }
    }

    if (hist)
    {
      if (fVerbose) printf("merging TH1\n");
      rc = hist->Merge(&fMergeList);
      if (!(rc<0)) {
        fMergeList.Delete();
        return rc;
      }
    }
    else if (mergeable)
    {
      if (fVerbose) printf("merging via AliMergeable ptr\n");
      rc = mergeable->Merge(&fMergeList);
      if (!(rc<0)) {
        fMergeList.Delete();
        return rc;
      }
    }
    else if (fMergeCall)
    {
      if (fVerbose) printf("merging via TMethodCall\n");
      fMergeCall->SetParam((Long_t) &fMergeList);
      fMergeCall->Execute(fObject);
      fMergeList.Delete();
      return 0;
    }
    else
    {
      if (fVerbose) printf("not mergeable, caching instead: %p %s\n",o,o->GetName());
      delete fObject;
      fObject = fMergeList.Remove(o);
      return 0;
    }
    return -1;
  }

  ~AliZMQROOTMergerEntry() {
    delete fObject;
    fMergeList.Delete();
    delete fMergeCall;
  }
};

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

      if (pushBack>0)
      {
        if (fVerbose) printf("pushback!\n");
        DoSend(fZMQout);
        SetLastPushBackTime();
      }
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
      if (pushBack>0)
      {
        if (fVerbose) printf("pushback!\n");
        DoSend(fZMQin);
        SetLastPushBackTime();
      }
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
      if (pushBack>0)
      {
        if (fVerbose) printf("pushback!\n");
        DoSend(fZMQmon);
        SetLastPushBackTime();
      }
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

  if (topic.GetID().compare(0,kAliHLTComponentDataTypefIDsize,kAliHLTDataTypeCDBEntry.fID,kAliHLTComponentDataTypefIDsize)==0)
  {
    //dont merge CDB entries, just cache them
  }
  else if (topic.GetID().compare(0,kAliHLTComponentDataTypefIDsize,kAliHLTDataTypeStreamerInfo.fID,kAliHLTComponentDataTypefIDsize)==0)
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

    //on run change
    if (runnumber != fRunNumber)
    {
      fRunNumber = runnumber;
      DoSend(socket);
      if (fAllowResetAtSOR)
      {
        if (ResetOutputData(fAllowResetAtSOR)>0)
        {
          if (fVerbose) printf("Run changed, merger reset!\n");
        }
      }
      else if (fAllowClearAtSOR)
      {
        if (ClearOutputData()>0)
        {
          if (fVerbose) printf("Run changed, objects cleared!\n");
        }
      }
      if (fZMQsync)
      {
        fBytesOUT += alizmq_msg_send(fInfoTopic, fInfo, fZMQsync, ZMQ_DONTWAIT);
        fNumberOfMessagesSent++;
      }
    }//on run change

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

  //reset the "one shot" options to default values
  fResetOnRequest = kFALSE;
  fSchemaOnRequest = false;
  fNameList.clear();
  if (fVerbose && (fSendSelection || fUnSendSelection))
  {
    Printf("unsetting include=%s, exclude=%s",
        (fSendSelection)?fSendSelection->GetPattern().Data():"",
        (fUnSendSelection)?fUnSendSelection->GetPattern().Data():"");
  }
  delete fSendSelection; fSendSelection=NULL;
  delete fUnSendSelection; fUnSendSelection=NULL;
  return rc;
}

//_____________________________________________________________________
int RemoveEntry(string& entry, mergeMap_t& map)
{
  map.erase(entry);
  return 0;
}

//_____________________________________________________________________
Int_t AddObject(TObject* object)
{
  if (!object)
  {
    if (fVerbose) Printf("no object!");
    return -1;
  }

  //add annotations
  string name = object->GetName();
  if (!fTitleAnnotation.IsNull())
  {
    TNamed* named = dynamic_cast<TNamed*>(object);
    if (!named) return 0;

    string title = named->GetTitle();
    title = fTitleAnnotation + " " + title;
    named->SetTitle(title.c_str());
  }
  if (!fNameAnnotation.IsNull())
  {
    TNamed* named = dynamic_cast<TNamed*>(object);
    if (!named) return 0;

    string name = named->GetName();
    name = fNameAnnotation + " " + name;
    named->SetName(name.c_str());
  }

  //do the merge or cache
  AliZMQROOTMergerEntry& mergerEntry = fMergeObjectMap[name];
  if (fCacheOnly)
  {
    mergerEntry = object;
  }
  else
  {
    mergerEntry += object;
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
  alizmq_msg_iter_data(block, object);

  //if we get a collection, always set ownership to prevent mem leaks
  //if we request unpacking: unpack what was requestd, otherwise just add
  do { //fake loop for easy breaking
    if (fUnpackContainers) {
      AliMergeable* unpackable = dynamic_cast<AliMergeable*>(object);
      if (unpackable) {
        if (fVerbose) printf("unpacking custom AliMergeable container %p\n", unpackable);
        GetObjects(unpackable, &fListOfObjects);
        if (fVerbose) printf("deleting custom AliMergeable container %p\n", unpackable);
        delete unpackable;
        break;
      }

      AliAnalysisDataContainer* container = dynamic_cast<AliAnalysisDataContainer*>(object);
      if (container) {
        //unpack an analysis data container
        if (fVerbose) printf("unpacking analysis container with name %s %p\n", container->GetName(), container);
        GetObjects(container, &fListOfObjects);
        if (fVerbose) printf("deleting analysis container with name %s %p\n", container->GetName(), container);
        delete container;
        break;
      }
      AliMergeable* unpackable = dynamic_cast<AliMergeable*>(object);
      if (unpackable) {
        if (fVerbose) printf("unpacking custom AliMergeable container %p\n", unpackable);
        GetObjects(unpackable, &fListOfObjects);
        if (fVerbose) printf("deleting custom AliMergeable container %p\n", unpackable);
        delete unpackable;
        break;
      }
    }

    if (TCollection* collection = dynamic_cast<TCollection*>(object)) {
      //unpack a collection
      if (fUnpackCollections) {
        if (fVerbose) printf("unpacking collection with name %s %p\n", collection->GetName(), collection);
        GetObjects(collection, &fListOfObjects);
        if (fVerbose) printf("deleting collection with name %s %p\n", collection->GetName(), collection);
        delete collection;
        break;
      } else {
        collection->SetOwner(kTRUE);
        fListOfObjects.push_back(collection);
        break;
      }
    } else {
      fListOfObjects.push_back(object);
      break;
    }
  } while (false);

  //add all extracted objects to the list of mergeable objects
  for (std::vector<TObject*>::iterator i=fListOfObjects.begin();
       i!=fListOfObjects.end(); ++i) {
    AddObject(*i);
  }
  fListOfObjects.clear();

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
  if (fRequestResetOnRequest) {
    if (fVerbose) Printf("sending an ResetOnRequest request");
    alizmq_msg_send(kAliHLTDataTypeConfig, "ResetOnRequest",socket,0);
    fNumberOfMessagesSent++;
  } else {
    if (fVerbose) Printf("sending an empty request");
    alizmq_msg_send("", "", socket, 0);
    fNumberOfMessagesSent++;
  }

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

  aliZMQmsg message;
  //forward the (run-)info string
  alizmq_msg_add(&message, &fInfoTopic, fInfo);
  Bool_t reset = kFALSE;
  Int_t rc = 0;
  int parts = 0;

  for (mergeMap_t::iterator i = fMergeObjectMap.begin(); i!=fMergeObjectMap.end(); ++i)
  {
    //the topic
    AliHLTDataTopic topic = kAliHLTDataTypeTObject|kAliHLTDataOriginOut;
    //the data
    string objectName = i->first;
    AliZMQROOTMergerEntry& entry = i->second;
    TObject* object = entry.GetObject();

    Bool_t selected = kTRUE;
    Bool_t unselected = kFALSE;
    if (fSendSelection) selected = fSendSelection->Match(objectName);
    if (fUnSendSelection) unselected = fUnSendSelection->Match(objectName);
    if (!fNameList.empty()) unselected = unselected ||
                                         (fNameList.find(objectName)==std::string::npos);
    if (!selected || unselected)
    {
      if (fVerbose) Printf("     object %s did NOT make the selection [%s] && ![%s]",
                           objectName.c_str(), (fSendSelection)?fSendSelection->GetPattern().Data():"",
                           (fUnSendSelection)?fUnSendSelection->GetPattern().Data():"");
      continue;
    }

    //add data to be sent
    rc = alizmq_msg_add(&message, &topic, object, fCompression, NULL);
    parts++;

    if (fResetOnSend || ( fResetOnRequest && fAllowResetOnRequest ))
    {
      reset = kTRUE;
      if (fClearOnReset) {
        if (fVerbose) {printf("resetting %s\n",objectName.c_str());}
        TH1* hist = dynamic_cast<TH1*>(object);
        if (hist) hist->Reset();
      } else {
        if (fVerbose) {printf("destroying %s\n",objectName.c_str());}
        fMergeObjectMap.erase(objectName);
      }
    }
  }

  if ((fSchemaOnRequest || fSchemaOnSend)) {
    alizmq_msg_prepend_streamer_infos(&message, &fSchema);
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
  if (fVerbose) Printf("merger sent %i bytes, %i parts, on socket %p", sentBytes, parts, socket);
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

//______________________________________________________________________________
int ResetOutputData(Bool_t force)
{
  if (fAllowGlobalReset || force)
  {
    if (fZMQresetBroadcast) {
      if (fVerbose) Printf("broadcasting data before reset");
      DoSend(fZMQresetBroadcast);
    }
    if (fClearOnReset) {
      if (fVerbose) Printf("Clearing the merger by calling Reset() on all data");
      return ClearOutputData();
    } else {
      if (fVerbose) Printf("Resetting the merger");
      fMergeObjectMap.clear();
      return 1;
    }
  }
  return 0;
}

//______________________________________________________________________________
int ClearOutputData()
{
  if (fZMQresetBroadcast) {
    if (fVerbose) Printf("broadcasting data before clear");
    DoSend(fZMQresetBroadcast);
  }

  for (mergeMap_t::iterator i = fMergeObjectMap.begin(); i!=fMergeObjectMap.end(); ++i)
  {
    //the data
    AliZMQROOTMergerEntry& entry = i->second;
    TObject* object = entry.GetObject();

    TH1* hist = dynamic_cast<TH1*>(object);
    if (!hist) continue;
    if (fVerbose) printf("clearing %s\n",hist->GetName());
    hist->Reset();
  }
  return 1;
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
    if (option.EqualTo("reset"))
    {
      ResetOutputData();
    }
    else if (option.EqualTo("ResetOnRequest"))
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
    else if (option.EqualTo("MaxObjects"))
    {
      fMaxObjects = value.Atoi();
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

  //init other stuff
  fListOfObjects.reserve(1000);

  ReadFromFile(fInitFile);

  Run();

  DumpToFile(fInitFile);

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
Int_t GetObjects(AliAnalysisDataContainer* kont, std::vector<TObject*>* list, std::string prefix)
{
  std::string kontName = kont->GetName();
  if (fIgnoreDefaultNamesWhenUnpacking && ( kontName=="" )) {
    kontName = "";
  } else {
    kontName+="/";
  }
  std::string kontPrefix = prefix + kontName;

  TObject* analData = kont->GetData();

  if (fFullyDestroyAnalysisDataContainer) {
    AliAnalysisTask* producer = kont->GetProducer();
    TObjArray* consumers = kont->GetConsumers();
    if (fVerbose) {
      printf("container %s has producer: %s(%p), consumers: %p, data: %p\n,  "
             "deleting producer/consumers\n",kontName.c_str(),
             (producer)?producer->GetName():"", producer, consumers, analData);
    }
    if (consumers) {
      consumers->SetOwner(kTRUE);
      delete consumers;
    }
    if (producer) {
      delete producer;
    }
  }

  if (TCollection* collection = dynamic_cast<TCollection*>(analData)) {
    //a collection
    if (fVerbose) Printf("  have a collection %s %p",collection->GetName(),collection);
    GetObjects(collection, list, kontPrefix);
    if (fVerbose) printf("  destroying collection %s %p\n",collection->GetName(),collection);
    delete collection;
    kont->SetDataOwned(kFALSE);

  } else {
    //just an object
    TNamed* named = dynamic_cast<TNamed*>(analData);
    if (named) {
      std::string name = kontPrefix + analData->GetName();
      named->SetName(name.c_str());
      if (fTitleAnnotationWithContainerName) {
        std::string title = kontPrefix + analData->GetTitle();
        named->SetTitle(title.c_str());
      }
      if (fVerbose) Printf("--in (from analysis container): %s (%s), %p",
          named->GetName(),
          named->ClassName(),
          named );
    }
    kont->SetDataOwned(kFALSE);
    list->push_back(analData);
  }
  return 0;
}

//______________________________________________________________________________
Int_t GetObjects(AliMergeable* unpackable, std::vector<TObject*>* list, std::string prefix)
{
  TCollection* collection = unpackable->GetListOfDrawableObjects();
  if (!collection) {
    if (fVerbose) {
      TObject* tmp = dynamic_cast<TObject*>(unpackable);
      printf("requested unpacking, but %s unpacks to NULL, skipping\n",(tmp)?tmp->GetName():"UNKNOWN");
    }
    return -1;
  }
  std::string collName = collection->GetName();
  if (fIgnoreDefaultNamesWhenUnpacking && ( collName=="TObjArray" || collName=="TList" || 
                                            collName=="AliHLTObjArray" || collName=="AliHLTList")) {
    collName = "";
  } else {
    collName+="/";
  }
  std::string collPrefix = prefix + collName;

  TIter next(collection);
  while (TObject* tmp = next()) {
    collection->Remove(tmp);
    TCollection* subcollection = dynamic_cast<TCollection*>(tmp);
    if (subcollection) {
      if (fVerbose) Printf("  have a subcollection %s %p",subcollection->GetName(), subcollection);
      GetObjects(subcollection, list, collPrefix);
      if (fVerbose) Printf("  destroying a subcollection %s %p",subcollection->GetName(), subcollection);
      delete subcollection;
    } else {
      //..or just an object
      TNamed* named = dynamic_cast<TNamed*>(tmp);
      if (named) {
        std::string name = collPrefix + named->GetName();
        named->SetName(name.c_str());
        if (fTitleAnnotationWithContainerName) {
          std::string title = collPrefix + named->GetTitle();
          named->SetTitle(title.c_str());
        }
        if (fVerbose) Printf("--in: %s (%s), %p",
                             named->GetName(),
                             named->ClassName(),
                             named );
      }
      list->push_back(tmp);
    }
  }
  collection->SetOwner(kTRUE);
  delete collection;
  return 0;
}

//______________________________________________________________________________
Int_t GetObjects(TCollection* collection, std::vector<TObject*>* list, std::string prefix)
{
  std::string collName = collection->GetName();
  if (fIgnoreDefaultNamesWhenUnpacking && ( collName=="TObjArray" || collName=="TList" || 
                                            collName=="AliHLTObjArray" || collName=="AliHLTList")) {
    collName = "";
  } else {
    collName+="/";
  }
  std::string collPrefix = prefix +collName;
  TIter next(collection);
  while (TObject* tmp = next()) {
    collection->Remove(tmp);

    AliAnalysisDataContainer* analKont = NULL;
    TCollection* subcollection = NULL;
    AliMergeable* unpackable = NULL;

    if (fUnpackContainers) {
      analKont = dynamic_cast<AliAnalysisDataContainer*>(tmp);
    }

    if (fUnpackContainers && !analKont) {
      unpackable = dynamic_cast<AliMergeable*>(tmp);
    }

    if (!analKont && !unpackable) {
      subcollection = dynamic_cast<TCollection*>(tmp);
    }

    if (analKont) {
      //analysis container
      if (fVerbose) Printf("  have an analysis container %s %p",analKont->GetName(), analKont);
      GetObjects(analKont,list,collPrefix);
      if (fVerbose) printf("  destroying anal container %s %p\n",analKont->GetName(), analKont);
      delete analKont;

    } else if (subcollection) {
      //embedded collection
      if (fVerbose) Printf("  have a subcollection %s %p",subcollection->GetName(), subcollection);
      GetObjects(subcollection, list, collPrefix);
      if (fVerbose) Printf("  destroying a subcollection %s %p",subcollection->GetName(), subcollection);
      delete subcollection;

    } else if (unpackable) {
        //something implementing a custom method to unpack into a list
        if (fVerbose) Printf("  have a AliMergeable object %p",unpackable);
        GetObjects(unpackable, list, collPrefix);
        if (fVerbose) Printf("  destroying the AliMergeable object %p",unpackable);
        delete unpackable;

    } else {
      //..or just an object
      TNamed* named = dynamic_cast<TNamed*>(tmp);
      if (named) {
        std::string name = collPrefix + named->GetName();
        named->SetName(name.c_str());
        if (fTitleAnnotationWithContainerName) {
          std::string title = collPrefix + named->GetTitle();
          named->SetTitle(title.c_str());
        }
        if (fVerbose) Printf("--in: %s (%s), %p",
                             named->GetName(),
                             named->ClassName(),
                             named );
      }
      list->push_back(tmp);
    }
    collection->SetOwner(kTRUE);
  } //while
  return 0;
}

//______________________________________________________________________________
Int_t ReadFromFile(std::string file)
{
  TDirectory::AddDirectory(kFALSE);
  TH1::AddDirectory(kFALSE);

  if (gSystem->AccessPathName(file.c_str())) { return -1; }
  TFile f (file.c_str(),"read");
  if (f.IsZombie()) { return -1; }

  TList* listOfKeys = f.GetListOfKeys();
  TIter keys(listOfKeys);
  while (TKey* key = static_cast<TKey*>(keys.Next()))
  {
    string objectName = key->GetName();

    //read and attach object
    TObject* object = key->ReadObj();

    //restore the metadata (runnumberr)
    if (objectName=="_ZMQ_internal_fInfo") {
      TObjString* objstr = dynamic_cast<TObjString*>(object);
      if (objstr) {
        fInfo = objstr->String();
        fRunNumber = atoi(GetParamString("run",fInfo).c_str());
      }
      if (fVerbose) printf("restoring metadata: %s\n", fInfo.c_str());
      delete objstr;
      continue;
    }

    if (object)
    {
      if (fVerbose) Printf("file (%s): attaching %s",file.c_str(),objectName.c_str());
      //prevent annotations to be added on top of the old ones
      TString oldTitleAnn = fTitleAnnotation;
      TString oldNameAnn = fNameAnnotation;
      fTitleAnnotation="";
      fNameAnnotation="";
      AddObject(object);
      fTitleAnnotation=oldTitleAnn;
      fNameAnnotation=oldNameAnn;
    }
  }
  f.Close();
  return 0;
}

//______________________________________________________________________________
Int_t DumpToFile(std::string file)
{
  int rc = 0;
  if (file.empty()) { return 0; }
  TFile f(file.c_str(),"recreate");
  if (f.IsZombie()) { return -1; }

  for (mergeMap_t::iterator i = fMergeObjectMap.begin(); i!=fMergeObjectMap.end(); ++i)
  {
    //the data
    AliZMQROOTMergerEntry& entry = i->second;
    TObject* object = entry.GetObject();
    if (!object) continue;
    if (fVerbose) printf("dumping %s to file %s\n", object->GetName(), file.c_str());
    rc = object->Write(object->GetName(),TObject::kOverwrite);
    if (rc<0) { rc=-1; }
  }

  TObjString info(fInfo.c_str());
  info.Write("_ZMQ_internal_fInfo", TObject::kOverwrite);

  f.Close();
  return 0;
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

