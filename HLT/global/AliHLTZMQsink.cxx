/**************************************************************************
 * This file is property of and copyright by the ALICE HLT Project        * 
 * ALICE Experiment at CERN, All rights reserved.                         *
 *                                                                        *
 * Primary Author: Mikolaj Krzewicki, mikolaj.krzewicki@cern.ch           *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

#include "AliHLTZMQsink.h"
#include "AliHLTErrorGuard.h"
#include "TDatime.h"
#include "TRandom3.h"
#include <TObject.h>
#include "zmq.h"
#include "AliZMQhelpers.h"
#include "AliCDBManager.h"
#include "AliCDBStorage.h"
#include "AliCDBEntry.h"
#include "AliHLTMessage.h"
#include "TStreamerInfo.h"
#include "TCollection.h"
#include "TList.h"
#include "AliZMQhelpers.h"

using namespace std;

ClassImp(AliHLTZMQsink)

//______________________________________________________________________________
AliHLTZMQsink::AliHLTZMQsink() :
  AliHLTComponent()
  , fZMQcontext(NULL)
  , fZMQout(NULL)
  , fZMQsocketType(-1)
  , fZMQoutConfig("PUB")
  , fZMQpollIn(kFALSE)
  , fPushbackDelayPeriod(-1)
  , fIncludePrivateBlocks(kFALSE)
  , fZMQneverBlock(kTRUE)
  , fSendRunNumber(kTRUE)
  , fNskippedErrorMessages(0)
  , fZMQerrorMsgSkip(100)
  , fSendECSparamString(kFALSE)
  , fECSparamString()
  , fSendStreamerInfos(kFALSE)
  , fCDBpattern("^/*[a-zA-Z0-9_.-]+/[a-zA-Z0-9_.-]+/[a-zA-Z0-9_.-]+/*$")
{
  //ctor
}

//______________________________________________________________________________
AliHLTZMQsink::~AliHLTZMQsink()
{
  //dtor
  zmq_close(fZMQout);
  zmq_ctx_destroy(fZMQcontext);
}

//______________________________________________________________________________
const Char_t* AliHLTZMQsink::GetComponentID()
{
  //id
  return "ZMQsink";
}

//______________________________________________________________________________
AliHLTComponentDataType AliHLTZMQsink::GetOutputDataType()
{
  // default method as sink components do not produce output
  AliHLTComponentDataType dt =
    {sizeof(AliHLTComponentDataType),
     kAliHLTVoidDataTypeID,
     kAliHLTVoidDataOrigin};
  return dt;
}

//______________________________________________________________________________
void AliHLTZMQsink::GetInputDataTypes( vector<AliHLTComponentDataType>& list)
{
  //what data types do we accept
  list.clear();
  list.push_back(kAliHLTAnyDataType);
}

//______________________________________________________________________________
void AliHLTZMQsink::GetOutputDataSize( unsigned long& constBase, double& inputMultiplier )
{
  // default method as sink components do not produce output
  constBase=0;
  inputMultiplier=0;
}

//______________________________________________________________________________
AliHLTComponent* AliHLTZMQsink::Spawn()
{
  //Spawn a new instance
  return new AliHLTZMQsink();
}

//______________________________________________________________________________
Int_t AliHLTZMQsink::DoInit( Int_t /*argc*/, const Char_t** /*argv*/ )
{
  // see header file for class documentation
  Int_t retCode=0;

  //process arguments
  if (ProcessOptionString(GetComponentArgs())<0) 
  {
    HLTFatal("wrong config string! %s", GetComponentArgs().c_str());
    return -1;
  }

  int rc = 0;
  //init ZMQ context
  fZMQcontext = alizmq_context();
  HLTMessage(Form("ctx create ptr %p %s",fZMQcontext,(rc<0)?zmq_strerror(errno):""));
  if (!fZMQcontext) return -1;

  //init ZMQ socket
  rc = alizmq_socket_init(fZMQout, fZMQcontext, fZMQoutConfig.Data(), 0, 10 ); 
  if (!fZMQout || rc<0) 
  {
    HLTError("cannot initialize ZMQ socket %s, %s",fZMQoutConfig.Data(),zmq_strerror(errno));
    return -1;
  }
  
  HLTMessage(Form("socket create ptr %p %s",fZMQout,(rc<0)?zmq_strerror(errno):""));
  HLTImportant(Form("ZMQ connected to: %s (%s(id %i)) rc %i %s",
               fZMQoutConfig.Data(),alizmq_socket_name(fZMQsocketType),
               fZMQsocketType,rc,(rc<0)?zmq_strerror(errno):""));
  
  return retCode;
}

//______________________________________________________________________________
Int_t AliHLTZMQsink::DoDeinit()
{
  // see header file for class documentation
  return 0;
}

//______________________________________________________________________________
int AliHLTZMQsink::DoProcessing( const AliHLTComponentEventData& evtData,
                                const AliHLTComponentBlockData* blocks, 
                                AliHLTComponentTriggerData& /*trigData*/,
                                AliHLTUInt8_t* /*outputPtr*/, 
                                AliHLTUInt32_t& /*size*/,
                                AliHLTComponentBlockDataList& outputBlocks,
                                AliHLTComponentEventDoneData*& /*edd*/ )
{ 
  // see header file for class documentation
  Int_t retCode=0;
  
  //create a default selection of any data:
  int requestTopicSize=-1;
  char requestTopic[kAliHLTComponentDataTypeTopicSize];
  memset(requestTopic, '*', kAliHLTComponentDataTypeTopicSize);

  int rc = 0;
  Bool_t doSend = kTRUE;
  Bool_t doSendECSparamString = kFALSE;
  Bool_t doSendStreamerInfos = kFALSE;
  Bool_t doSendCDB = kFALSE;
  AliCDBEntry* cdbEntry = NULL;
  
  //cache an ECS param topic
  char ecsParamTopic[kAliHLTComponentDataTypeTopicSize];
  DataType2Topic(kAliHLTDataTypeECSParam, ecsParamTopic);
  TString requestedCDBpath;

  //in case we reply to requests instead of just pushing/publishing
  //we poll for requests
  if (fZMQpollIn)
  {
    zmq_pollitem_t items[] = { { fZMQout, 0, ZMQ_POLLIN, 0 } };
    zmq_poll(items, 1, 0);

    if (items[0].revents & ZMQ_POLLIN)
    {
      int64_t more=0;
      size_t moreSize=sizeof(more);
      do //request could be multipart, get all parts
      {
        requestTopicSize = zmq_recv (fZMQout, requestTopic, kAliHLTComponentDataTypeTopicSize, 0);
        zmq_getsockopt(fZMQout, ZMQ_RCVMORE, &more, &moreSize);
        zmq_msg_t requestMsg;
        int requestSize=-1;
        if (more) {
          zmq_msg_init(&requestMsg);
          requestSize = zmq_msg_recv(&requestMsg, fZMQout, 0);
          zmq_getsockopt(fZMQout, ZMQ_RCVMORE, &more, &moreSize);
        }
        //if request is for ECS params, set the flag
        if (*reinterpret_cast<const AliHLTUInt64_t*>(requestTopic) ==
            *reinterpret_cast<const AliHLTUInt64_t*>(kAliHLTDataTypeECSParam.fID))
        {
          doSendECSparamString = kTRUE;
        }
        //if request is for streamer infos, set the flag
        else if (*reinterpret_cast<const AliHLTUInt64_t*>(requestTopic) ==
                 *reinterpret_cast<const AliHLTUInt64_t*>(kAliHLTDataTypeStreamerInfo.fID))
        {
          doSendStreamerInfos = kTRUE;
        }
        //if request is for an OCDB object, set the flag
        else if (requestSize>0 &&
            *reinterpret_cast<const AliHLTUInt64_t*>(requestTopic) ==
            *reinterpret_cast<const AliHLTUInt64_t*>(kAliHLTDataTypeCDBEntry.fID))
        {
          requestedCDBpath.Append(static_cast<char*>(zmq_msg_data(&requestMsg)),requestSize);
          
          //get the CDB entry in a safe way
          do {
            if (!requestedCDBpath.Contains(fCDBpattern)) {
              HLTWarning("malformed CDB path: %s", requestedCDBpath.Data());
              break;
            }
            AliCDBManager* cdbMan = AliCDBManager::Instance();
            if (!cdbMan) break;
            AliCDBStorage* cdbStorage = cdbMan->GetDefaultStorage();
            if (!cdbStorage) break;
            AliCDBId* cdbId = cdbStorage->GetId(requestedCDBpath.Data(),GetRunNo());
            if (!cdbId) {
              HLTWarning("cannot get CDB entry: %s", requestedCDBpath.Data());
              break;
            }
            cdbEntry = cdbMan->Get(requestedCDBpath.Data());
            if (!cdbEntry) break;
            doSendCDB = kTRUE;
          } while (false);
        }
        zmq_msg_close(&requestMsg);
      } while (more==1);
    }
    else { doSend = kFALSE; }
  }
 
  //if enabled (option -pushback-period), send at most so often
  if (fPushbackDelayPeriod>0)
  {
    TDatime time;
    if ((Int_t)time.Get()-fLastPushbackDelayTime<fPushbackDelayPeriod) 
    {
      doSend=kFALSE;
    }
  }  

  //caching the ECS param string has to happen for non data event
  if (!IsDataEvent())
  {
    const AliHLTComponentBlockData* inputBlock = NULL;
    for (int iBlock = 0;
         iBlock < evtData.fBlockCnt;
         iBlock++) 
    {
      inputBlock = &blocks[iBlock];
      //cache the ECS param string
      if (*reinterpret_cast<const AliHLTUInt64_t*>(inputBlock->fDataType.fID) ==
          *reinterpret_cast<const AliHLTUInt64_t*>(kAliHLTDataTypeECSParam.fID))
      {
        const char* ecsparamstr = reinterpret_cast<const char*>(inputBlock->fPtr);
        int ecsparamsize = inputBlock->fSize;
        if (ecsparamstr[ecsparamsize-1]!=0)
        {
          fECSparamString.Insert(0, ecsparamstr, ecsparamsize);
          fECSparamString += "";
        }
        else
        {
          fECSparamString = ecsparamstr;
        }
        break;
      }
    }
  }

  //always cache streamer infos
  {
    const AliHLTComponentBlockData* inputBlock = NULL;
    for (int iBlock = 0;
        iBlock < evtData.fBlockCnt;
        iBlock++) 
    {
      inputBlock = &blocks[iBlock];
      //cache the streamer info
      if (*reinterpret_cast<const AliHLTUInt64_t*>(inputBlock->fDataType.fID) ==
          *reinterpret_cast<const AliHLTUInt64_t*>(kAliHLTDataTypeStreamerInfo.fID))
      {
        TObject* obj = NULL;
        AliHLTUInt32_t firstWord=*((AliHLTUInt32_t*)inputBlock->fPtr);
        if (firstWord==inputBlock->fSize-sizeof(AliHLTUInt32_t))
        {
          HLTDebug("create object from block %d size %d", iBlock, inputBlock->fSize);
          AliHLTMessage msg(inputBlock->fPtr, inputBlock->fSize);
          TClass* objclass=msg.GetClass();
          obj=msg.ReadObject(objclass);
        }
        TCollection* coll = dynamic_cast<TCollection*>(obj);
        if (coll)
        {
          HLTMessage("updating streamer infos");
          UpdateSchema(coll);
        }
        //delete the remaining infos and destroy the collection
        coll->SetOwner(kTRUE);
        delete coll;
      } //if kAliHLTDataTypeStreamerInfo
    } //for iBlock
  } //dummy scope

  if (doSend)
  {
    //set the time of current push
    if (fPushbackDelayPeriod>0)
    {
      TDatime time;
      fLastPushbackDelayTime=time.Get();
    }

    //first make a map of selected blocks, and identify the last one
    //so we can properly mark the last block for multipart ZMQ sending later
    const AliHLTComponentBlockData* inputBlock = NULL;
    std::vector<int> selectedBlockIdx;
    for (int iBlock = 0;
         iBlock < evtData.fBlockCnt;
         iBlock++) 
    {
      inputBlock = &blocks[iBlock];

      //don't include provate data unless explicitly asked to
      if (!fIncludePrivateBlocks && 
          *reinterpret_cast<const AliHLTUInt32_t*>(inputBlock->fDataType.fOrigin) ==
          *reinterpret_cast<const AliHLTUInt32_t*>(kAliHLTDataOriginPrivate))
      {
        continue;
      }

      //check if the data type matches the request
      char blockTopic[kAliHLTComponentDataTypeTopicSize];
      DataType2Topic(inputBlock->fDataType, blockTopic);
      if (Topicncmp(requestTopic, blockTopic, requestTopicSize))
      {
        selectedBlockIdx.push_back(iBlock);
      }
    }
    int nSelectedBlocks = selectedBlockIdx.size();
    int nSentBlocks = 0;

    aliZMQmsg message;

    //only send the INFO block if there is some data to send
    if (fSendRunNumber && nSelectedBlocks>0)
    {
      string runNumberString = "run=";
      char tmp[34];
      snprintf(tmp,34,"%i",GetRunNo()); 
      runNumberString+=tmp;
      rc = alizmq_msg_add(&message, "INFO", runNumberString);
      if (rc<0) {
        HLTWarning("ZMQ error adding INFO");
      }
    }

    //maybe send the ECS param string
    //once if requested or always if so configured
    if ((fSendECSparamString && nSelectedBlocks>0) || doSendECSparamString)
    {
      AliHLTDataTopic topic = kAliHLTDataTypeECSParam;
      rc = alizmq_msg_add(&message, &topic, fECSparamString.Data());
      if (rc<0) {
        HLTWarning("ZMQ error adding ECS param string");
      }
      doSendECSparamString = kFALSE;
    }

    //send the streamer infos if requested
    if ((fSendStreamerInfos && nSelectedBlocks>0)|| doSendStreamerInfos)
    {
      AliHLTDataTopic topic = kAliHLTDataTypeStreamerInfo;
      rc = alizmq_msg_add(&message, &topic, GetSchema(), GetCompressionLevel());
      if (rc<0) {
        HLTWarning("ZMQ error adding schema infos");
      }
      doSendStreamerInfos = kFALSE;
    }

    //send the CDB entry if requested (on request only)
    if (doSendCDB && cdbEntry)
    {
      AliHLTDataTopic topic = kAliHLTDataTypeCDBEntry;
      rc = alizmq_msg_add(&message, &topic, cdbEntry, GetCompressionLevel());
      if (rc<0) {
        HLTWarning("ZMQ error adding CDB entry %s", requestedCDBpath.Data());
      }
      doSendCDB = kFALSE;
    }

    //send the selected blocks
    for (int iSelectedBlock = 0;
         iSelectedBlock < selectedBlockIdx.size();
         iSelectedBlock++) 
    {
      inputBlock = &blocks[selectedBlockIdx[iSelectedBlock]];
      AliHLTDataTopic blockTopic = *inputBlock;

      rc = alizmq_msg_add(&message, &blockTopic, inputBlock->fPtr, inputBlock->fSize);
      if (rc<0) {
        HLTWarning("ZMQ error adding block %s", blockTopic.Description().c_str());
      }
      HLTMessage(Form("send data rc %i %s",rc,(rc<0)?zmq_strerror(errno):""));
    }

    
    //send an empty message if we really need a reply (ZMQ_REP mode)
    //only in case no blocks were sent
    if (message.size()==0 && fZMQsocketType==ZMQ_REP)
    { 
      rc = alizmq_msg_add(&message, "", "");
      if (rc<0) {
        HLTWarning("ZMQ error adding dummy rep data");
      }
    }
    rc = alizmq_msg_send(&message, fZMQout, 0);
    if (rc<0){ 
      HLTWarning("ZMQ error sending message: %s", zmq_strerror(errno));
    }
    alizmq_msg_close(&message);
  }

  outputBlocks.clear();
  return retCode;
}

//______________________________________________________________________________
int AliHLTZMQsink::ProcessOption(TString option, TString value)
{
  //process option
  //to be implemented by the user
  
  if (option.EqualTo("out"))
  {
    fZMQoutConfig = value;
    fZMQsocketType = alizmq_socket_type(value.Data());
    switch (fZMQsocketType)
    {
      case ZMQ_REP:
        fZMQpollIn=kTRUE;
        break;
      case ZMQ_PUSH:
        fZMQpollIn=kFALSE;
        break;
      case ZMQ_PUB:
        fZMQpollIn=kFALSE;
        break;
      default:
        HLTFatal("use of socket type %s for a sink is currently unsupported! (config: %s)", alizmq_socket_name(fZMQsocketType), fZMQoutConfig.Data());
        return -EINVAL;
    }
  }
  else if (option.EqualTo("SendRunNumber"))
  {
    fSendRunNumber=(value.EqualTo("0") || value.EqualTo("no") || value.EqualTo("false"))?kFALSE:kTRUE;
  }

  else if (option.EqualTo("SendECSparamString"))
  {
    fSendECSparamString=(value.EqualTo("0") || value.EqualTo("no") || value.EqualTo("false"))?kFALSE:kTRUE;
  }

  else if (option.EqualTo("pushback-period"))
  {
    HLTMessage(Form("Setting pushback delay to %i", atoi(value.Data())));
    fPushbackDelayPeriod = atoi(value.Data());
  }

  else if (option.EqualTo("IncludePrivateBlocks"))
  {
    fIncludePrivateBlocks=kTRUE;
  }

  else if (option.EqualTo("ZMQneverBlock"))
  {
    if (value.EqualTo("0") || value.EqualTo("no") || value.Contains("false",TString::kIgnoreCase))
      fZMQneverBlock = kFALSE;
    else if (value.EqualTo("1") || value.EqualTo("yes") || value.Contains("true",TString::kIgnoreCase) )
      fZMQneverBlock = kTRUE;
  }
  
  else if (option.EqualTo("ZMQerrorMsgSkip"))
  {
    fZMQerrorMsgSkip = value.Atoi();
  }

  else if (option.EqualTo("schema"))
  {
    fSendStreamerInfos = kTRUE;
  }

  else
  {
    HLTError("unrecognized option %s", option.Data());
    return -1;
  }

  return 1; 
}

