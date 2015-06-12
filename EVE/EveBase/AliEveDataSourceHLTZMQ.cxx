//
//  AliEveDataSourceHLTZMQ
//
//  blame: Mikolaj Krzewicki, mikolaj.krzewicki@cern.ch
//
//

#include "AliEveDataSourceHLTZMQ.h"

#include "AliEveConfigManager.h"
#include "AliCDBManager.h"
#include "AliCDBStorage.h"
#include "AliGRPPreprocessor.h"
#include <TEnv.h>
#include <TInterpreter.h>
#include <iostream>

#include "AliHLTComponent.h"
#include "AliHLTMessage.h"

#ifdef ZMQ
#include "zmq.h"
#endif


using namespace std;

AliEveDataSourceHLTZMQ::AliEveDataSourceHLTZMQ(bool storageManager) :
    AliEveDataSource("HLT"),
    fEventListenerThreadHLT(0),
    fCurrentRun(-1),
    fZMQContext(NULL),
    fZMQeventQueue(NULL),
    fHLTPublisherAddress("tcp://localhost:60201")
{
  //ctor
  Init();
}

void AliEveDataSourceHLTZMQ::Init()
{
#ifdef ZMQ
  //get the address of the HLT proxy from the environment
  if (gSystem->Getenv("HLT_ZMQ_proxy")) 
    fHLTPublisherAddress=gSystem->Getenv("HLT_ZMQ_proxy");
  //single ZMQ context for inter thread comm. etc.
  if (!fZMQContext) fZMQContext = zmq_ctx_new();
  //single ZMQ socket for gathering the events form various listening threads
  //must be bound before threads can connect
  if (! fZMQeventQueue) 
  {
    fZMQeventQueue = zmq_socket(fZMQContext, ZMQ_PULL);
    zmq_bind(fZMQeventQueue, "inproc://fCurrentEvent");
  }

  AliInfo("Starting HLT subscriber thread.");
  if(fEventListenerThreadHLT)
  {
    fEventListenerThreadHLT->Join();
    fEventListenerThreadHLT->Kill();
    delete fEventListenerThreadHLT;
    AliInfo("HLT listener thread killed and deleted");
  }
  fEventListenerThreadHLT = new TThread("fEventListenerThreadHLT",DispatchEventListenerHLT,(void*)this);
  fEventListenerThreadHLT->Run();
#endif
}

AliEveDataSourceHLTZMQ::~AliEveDataSourceHLTZMQ()
{
#ifdef ZMQ
  if (fZMQeventQueue)
  {
    int lingerValue = 0;
    int rc = zmq_setsockopt(fZMQeventQueue, ZMQ_LINGER, &lingerValue, sizeof(int));
    if (rc<0) printf("error setting linger on fZMQeventQueue\n");
    rc = zmq_close(fZMQeventQueue);
    if (rc<0) printf("error closing socket fZMQeventQueue\n");
  }

  printf("trying to destroy fZMQContext\n");
  if (fZMQContext) zmq_ctx_destroy(fZMQContext);
  printf("destroyed fZMQContext\n");
#endif

  if(fEventListenerThreadHLT)
  {
    fEventListenerThreadHLT->Join();
    fEventListenerThreadHLT->Kill();
    delete fEventListenerThreadHLT;
    cout<<"HLT listener thread killed and deleted"<<endl;
  }


}

void AliEveDataSourceHLTZMQ::PullEventFromHLT()
{
#ifdef ZMQ
  int rc = 0;
  
  if (!fZMQContext) return;
  
  //connection to the HLT data
  void* listenerSocket = zmq_socket(fZMQContext, ZMQ_SUB);
  if (!listenerSocket) {return;}
  printf("connecting to ZMQ socket: %s\n", fHLTPublisherAddress.Data());
  rc = zmq_connect(listenerSocket, fHLTPublisherAddress.Data());
  rc = zmq_setsockopt (listenerSocket, ZMQ_SUBSCRIBE, NULL, 0);
  if (rc < 0) {return;}

  //internal publisher
  void* internalPublisher = zmq_socket(fZMQContext, ZMQ_PUSH);
  if (!internalPublisher) {return;}
  rc = zmq_connect(internalPublisher, "inproc://fCurrentEvent");
  if (rc < 0) {return;}
  
  //create a default selection of any data:
  char requestedTopic[kAliHLTComponentDataTypeTopicSize+1]; memset(&requestedTopic, 0, sizeof(requestedTopic));
  memset(&requestedTopic, '*', kAliHLTComponentDataTypeTopicSize);
  strncpy(requestedTopic, "ALIESDV0HLT**************", kAliHLTComponentDataTypeTopicSize);
  char requestedTopic1[kAliHLTComponentDataTypeTopicSize+1]; memset(&requestedTopic1, 0, sizeof(requestedTopic1));
  memset(&requestedTopic1, '*', kAliHLTComponentDataTypeTopicSize);
  strncpy(requestedTopic1, "FLATESD**************", kAliHLTComponentDataTypeTopicSize);

  Bool_t done=false;
  while (!done)
  {
    int more = 0;
    size_t more_size = sizeof more;
    do
    {
      //here the HLT data comes in
      //receive the topic
      char receivedTopic[kAliHLTComponentDataTypeTopicSize+1]; memset(&receivedTopic, 0, sizeof(receivedTopic));
      int receivedTopicSize = zmq_recv(listenerSocket, &receivedTopic, sizeof(receivedTopic), 0);
      if (receivedTopicSize < 0 && errno==ETERM) {AliInfo("listenerSocket received ETERM, exiting main loop"); done=true; break;}

      //reive the message if there is any
      zmq_msg_t message;
      rc = zmq_msg_init(&message);
      rc = zmq_getsockopt(listenerSocket, ZMQ_RCVMORE, &more, &more_size);
      if (more) 
      {
        rc = zmq_msg_recv(&message, listenerSocket, 0);
        if (rc < 0 && errno==ETERM) { zmq_msg_close(&message); done=true; break; }
      }
      rc = zmq_getsockopt(listenerSocket, ZMQ_RCVMORE, &more, &more_size);
      
      //after receiving a single part (topic+message) publish it internally
      //in general this will be an HLT ESD.
      //do check first if this is what we wanted
      bool isSelected = Topicncmp(requestedTopic, receivedTopic, kAliHLTComponentDataTypeTopicSize, receivedTopicSize) ||
                        Topicncmp(requestedTopic1, receivedTopic, kAliHLTComponentDataTypeTopicSize, receivedTopicSize);
                      
      if (!isSelected)
      {
        //clean up the memory and continue
        //printf("  dropped %s\n",receivedTopic);
        zmq_msg_close(&message);
        continue;
      }

      //publish the data for further use
      rc = -1;
      //printf("pushing %s\n", receivedTopic);
      rc = zmq_send(internalPublisher, &receivedTopic, receivedTopicSize, ZMQ_SNDMORE);
      rc = zmq_msg_send(&message, internalPublisher, 0);
      if (rc < 0) zmq_msg_close(&message);

    } while (more != 0);
  }
  //cleanly exit the thread
  int lingerValue = 10;
  rc = zmq_setsockopt(listenerSocket, ZMQ_LINGER, &lingerValue, sizeof(lingerValue));
  if (rc<0) printf("error setting linger on listenerSocket, errno: %i\n", errno);
  rc = zmq_setsockopt(internalPublisher, ZMQ_LINGER, &lingerValue, sizeof(lingerValue));
  if (rc<0) printf("error setting linger on internalPublisher, errno: %i\n",errno);
  rc = zmq_close(listenerSocket);
  if (rc<0) printf("error closing listenerSocket, errno: %i\n", errno);
  rc = zmq_close(internalPublisher);
  if (rc<0) printf("error closing internalPublisher, errno: %i\n", errno);
#endif
}

void AliEveDataSourceHLTZMQ::InitOCDB(int runNo)
{  
  AliInfo(Form("Initializing OCDB for run: %i",runNo));
  AliCDBManager* cdb = AliCDBManager::Instance();
  if (!cdb) return;

  if (cdb->IsDefaultStorageSet() == kTRUE)
  {
    AliInfo(Form("CDB already set - using the old storage:  '%s'", cdb->GetDefaultStorage()->GetURI().Data()));
    return;
  }
  else {
    if (gSystem->Getenv("ocdbStorage")) {
    cdb->SetDefaultStorage(gSystem->Getenv("ocdbStorage"));
    if (runNo!=AliEveEventManager::GetMaster()->GetCurrentRun()){
        AliEveEventManager::GetMaster()->SetCurrentRun(runNo);
      cdb->SetRun(runNo);
    }
}

void AliEveDataSourceHLTZMQ::GotoEvent(Int_t /*event*/)
{
    NextEvent();
    return;
}

void AliEveDataSourceHLTZMQ::NextEvent()
{
    static const TEveException kEH("AliEveDataSourceHLTZMQ::NextEvent ");
  //read event from queue

#ifdef ZMQ
  //init some stuff
  int rc = 0;
  AliESDEvent* esdObject = NULL;
  char topic[kAliHLTComponentDataTypeTopicSize+1]; memset(topic, 0, sizeof(topic));
  memset(topic, 0, kAliHLTComponentDataTypeTopicSize);
  zmq_msg_t message;
  rc = zmq_msg_init(&message);
  int64_t more = 0;
  size_t more_size = sizeof(more);

  //try to receive a new topic+message from the queue
  int topicSize = zmq_recv(fZMQeventQueue, &topic, kAliHLTComponentDataTypeTopicSize, ZMQ_DONTWAIT);
  int errnoTopic = errno;
  rc = zmq_getsockopt(fZMQeventQueue, ZMQ_RCVMORE, &more, &more_size);
  if (more) rc = zmq_msg_recv(&message, fZMQeventQueue, ZMQ_DONTWAIT);

  //decode the message into an AliESDEvent,
  //for HLT messages AliHLTMessage needs to be used
  if( topicSize > 0 /*&& errnoTopic != EAGAIN*/ && zmq_msg_size(&message) > 0 )
  {
    if (Topicncmp(topic, "ALIESDV0HLT", topicSize, 11))
    {
      printf("trying to decode %s with AliHLTMessage\n", topic);
      AliHLTMessage msg(zmq_msg_data(&message), zmq_msg_size(&message));
      TClass* objclass = msg.GetClass();
      esdObject = static_cast<AliESDEvent*>(msg.ReadObject(AliESDEvent::Class()));
    }
    else if (Topicncmp(topic, "ALIESDV0POR", topicSize, 11))
    {
      printf("trying to decode %s with TBufferFile\n", topic);
      TBufferFile *mess = new TBufferFile(TBuffer::kRead,
          zmq_msg_size(&message)+sizeof(UInt_t),
          zmq_msg_data(&message));
      mess->InitMap();
      mess->ReadClass();// get first the class stored in message
      mess->SetBufferOffset(sizeof(UInt_t) + sizeof(kMESS_OBJECT));
      mess->ResetMap();
      esdObject = (AliESDEvent*)(mess->ReadObjectAny(AliESDEvent::Class()));
    } 
    else
    {
      printf("Don't know what to do with %s\n",topic);
    }
  }
  if (esdObject) 
  {
    printf("setting new event\n");
    esdObject->GetStdContent();
    int runNumber = esdObject->GetRunNumber();
    InitOCDB(runNumber);

    printf("deleting current data\n");
    //replace the ESD
    fCurrentData.Clear();
    fCurrentData.fESD = esdObject;
    
    AliEveEventManager::GetMaster()->SetHasEvent(true);
    AliEveEventManager::GetMaster()->AfterNewEventLoaded();
    //if (AliEveEventManager::GetMaster()->GetAutoLoad()) {
    //  AliEveEventManager::GetMaster()->StartAutoLoadTimer();
    //}
  }
  else
  {
    cout<<"No new event is avaliable."<<endl;
    //AliEveEventManager::GetMaster()->NoEventLoaded();
  }
  zmq_msg_close(&message);
#endif
}

void AliEveDataSourceHLTZMQ::SetCdbUri(const TString& cdb)
{
    // Set path to CDB, there is no default.
    if ( ! cdb.IsNull()) fgCdbUri = cdb;
    AliCDBManager* cdbManager = AliCDBManager::Instance();
    cdbManager->SetDefaultStorage(cdb);
}


