//
//  AliEveDataSourceHLTZMQ
//
//  blame: Mikolaj Krzewicki, mikolaj.krzewicki@cern.ch
//
//

#include "AliEveDataSourceHLTZMQ.h"

#include "AliEveConfigManager.h"
#include "AliEveInit.h"
#include "AliGRPPreprocessor.h"
#include <TEnv.h>
#include <TInterpreter.h>
#include <iostream>

#include "AliHLTDataTypes.h"
#include "TSystem.h"
#include "AliESDEvent.h"
#include "AliCDBStorage.h"
#include "AliCDBEntry.h"

#ifdef ZMQ
#include "AliZMQhelpers.h"
#include "zmq.h"
#endif


using namespace std;

AliEveDataSourceHLTZMQ::AliEveDataSourceHLTZMQ(bool storageManager) :
    AliEveDataSource("HLT"),
    fZMQContext(NULL),
    fZMQin(NULL),
    fZMQtimeout(5000)
{
  //ctor
  TEnv settings;
  AliEveInit::GetConfig(&settings);
  fSourceURL=settings.GetValue("HLT.ZMQ.proxy","SUB>tcp://localhost:60201");
  Init();
}

void AliEveDataSourceHLTZMQ::Init()
{
  int rc=0;
#ifdef ZMQ
  //get the address of the HLT proxy from the environment if not set
  if (gSystem->Getenv("HLT_ZMQ_proxy"))
    fSourceURL=gSystem->Getenv("HLT_ZMQ_proxy");
  //single ZMQ context for inter thread comm. etc.
  if (!fZMQContext) fZMQContext = alizmq_context();
  rc = alizmq_socket_init(fZMQin, fZMQContext, fSourceURL.Data(), 10000, 10);
  printf("Initialized ZMQ socket %s (%s), rc=%i %s\n",
      alizmq_socket_name(alizmq_socket_type(fZMQin)),fSourceURL.Data(),rc,(rc<0)?zmq_strerror(errno):"");
#endif
}

AliEveDataSourceHLTZMQ::~AliEveDataSourceHLTZMQ()
{
#ifdef ZMQ
  int rc = alizmq_socket_close(fZMQin);

  printf("trying to destroy fZMQContext\n");
  if (fZMQContext) zmq_ctx_destroy(fZMQContext);
  printf("destroyed fZMQContext\n");
#endif
}

void AliEveDataSourceHLTZMQ::GotoEvent(Int_t /*event*/)
{
    NextEvent();
    return;
}

void AliEveDataSourceHLTZMQ::RequestData()
{
  aliZMQmsg request;
  if ((alizmq_socket_state(fZMQin) & ZMQ_POLLOUT)==ZMQ_POLLOUT)
  {
    alizmq_msg_add(&request, "", "");
  }
  alizmq_msg_send(&request, fZMQin, 0);
  printf("sent data request on %s\n", fSourceURL.Data());
}

void AliEveDataSourceHLTZMQ::NextEvent()
{
    static const TEveException kEH("AliEveDataSourceHLTZMQ::NextEvent ");
  //read event from queue

#ifdef ZMQ
  //init some stuff
  int rc = 0;
  TObject* object = NULL;

  //if we are requesting and socket is not already waiting for reply
  if (alizmq_socket_type(fZMQin)==ZMQ_REQ)
  {
    RequestData();
  }

  //wait for the data
  Int_t nSockets=1;
  zmq_pollitem_t sockets[] = { 
    { fZMQin, 0, ZMQ_POLLIN, 0 },
  };
  rc = zmq_poll(sockets, nSockets, fZMQtimeout); //poll sockets
  if (rc==-1 && errno==ETERM)
  {
    //this can only happen if the context was terminated, one of the sockets are
    //not valid or operation was interrupted
    Printf("zmq_poll was interrupted! rc = %i, %s", rc, zmq_strerror(errno));
    return;
  }

  //if we time out (waiting for a response) reinit the REQ socket(s)
  if (rc==0 && alizmq_socket_type(fZMQin)==ZMQ_REQ)
  {
    printf("no reply from %s in %i ms, server died?\n",
        fSourceURL.Data(), fZMQtimeout);
    rc = alizmq_socket_init(fZMQin, fZMQContext, fSourceURL.Data(), 5000, 10);
    return;
  }

  aliZMQmsg message;
  rc = alizmq_msg_recv(&message,fZMQin,ZMQ_DONTWAIT);
  printf("received rc=%i data on %s\n",rc,fSourceURL.Data());
  
  for (aliZMQmsg::iterator i=message.begin(); i!=message.end(); ++i)
  {
    AliHLTDataTopic topic;
    alizmq_msg_iter_topic(i,topic);
    AliHLTDataTopic esdTopic = kAliHLTDataTypeESDObject;
    if (topic.GetID() == esdTopic.GetID())
    {
      printf("  unpacking %s\n", topic.Description().c_str());
      alizmq_msg_iter_data(i, object);
    }
  }  
  alizmq_msg_close(&message);

  AliESDEvent* esdObject = dynamic_cast<AliESDEvent*>(object);
  if (esdObject) 
  {
    printf("setting new event\n");
    esdObject->GetStdContent();
    int runNumber = esdObject->GetRunNumber();
    printf("  run: %i, number of tracks: %i\n", runNumber, esdObject->GetNumberOfTracks());

    //replace the ESD
    printf("deleting current data\n");
    fCurrentData.Clear();
    fCurrentData.fESD = esdObject;

    AliEveEventManager::GetMaster()->DestroyElements();

    if(esdObject->GetRunNumber() != AliEveEventManager::GetMaster()->GetCurrentRun()){
      AliEveEventManager::GetMaster()->ResetMagneticField();
      AliEveEventManager::GetMaster()->SetCurrentRun(esdObject->GetRunNumber());
    }

    AliEveEventManager::GetMaster()->SetHasEvent(true);
    AliEveEventManager::GetMaster()->AfterNewEventLoaded();
  }
  else
  {
    printf("No new event is avaliable on %s\n", fSourceURL.Data());
  }
#endif
}

Bool_t AliEveDataSourceHLTZMQ::ReceivePromptRecoParameters(Int_t runNo)
{
#ifdef ZMQ
  //makes no sense to do anything here if we cannot request anything
  if (alizmq_socket_type(fZMQin)!=ZMQ_REQ) return kFALSE;

  aliZMQmsg request;
  if ((alizmq_socket_state(fZMQin) & ZMQ_POLLOUT)==ZMQ_POLLOUT)
  {
    alizmq_msg_add(&request, "CDBENTRY", "GRP/GRP/Data");
    alizmq_msg_add(&request, "ROOTSTRI", "");
  }
  alizmq_msg_send(&request, fZMQin, 0);
  printf("requested GRP and streamers\n");

  aliZMQmsg message;
  int rc = alizmq_msg_recv(&message,fZMQin,0);

  TObject* cdbGRPEntryObj = NULL;
  AliHLTDataTopic schemaTopic = kAliHLTDataTypeStreamerInfo;
  AliHLTDataTopic cdbTopic = kAliHLTDataTypeCDBEntry;
  for (aliZMQmsg::iterator i=message.begin(); i!=message.end(); ++i)
  {
    AliHLTDataTopic topic;
    alizmq_msg_iter_topic(i,topic);
    if (topic.GetID() == schemaTopic.GetID())
    {
      printf("received schema, initializing...\n");
      alizmq_msg_iter_init_streamer_infos(i);
    }
    //GRP
    else if (topic.GetID() == cdbTopic.GetID())
    {
      printf("received the GRP entry, unpacking...\n");
      alizmq_msg_iter_data(i, cdbGRPEntryObj);
    }
  }  
  alizmq_msg_close(&message);  

  AliCDBEntry* cdbGRPEntry = dynamic_cast<AliCDBEntry*>(cdbGRPEntryObj);
  if (cdbGRPEntry) {
    printf("storing the GRP in local://OCDB\n");
    AliCDBStorage* localStorage = AliCDBManager::Instance()->GetStorage("local://OCDB");
    localStorage->Put(cdbGRPEntry);
    AliCDBManager* man = AliCDBManager::Instance();
    man->SetSpecificStorage("GRP/GRP/Data","local://OCDB");
    delete cdbGRPEntry;
    return kTRUE;
  }
#endif
  return kFALSE;
}
