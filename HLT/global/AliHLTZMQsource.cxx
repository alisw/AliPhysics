// $Id$

///**************************************************************************
///* This file is property of and copyright by the                          *
///* ALICE Experiment at CERN, All rights reserved.                         *
///*                                                                        *
///* Primary Authors: Mikolaj Krzewicki, mikolaj.krzewicki@cern.ch          *
///*                                                                        *
///* Permission to use, copy, modify and distribute this software and its   *
///* documentation strictly for non-commercial purposes is hereby granted   *
///* without fee, provided that the above copyright notice appears in all   *
///* copies and that both the copyright notice and this permission notice   *
///* appear in the supporting documentation. The authors make no claims     *
///* about the suitability of this software for any purpose. It is          *
///* provided "as is" without express or implied warranty.                  *
///**************************************************************************

/// @file   AliHLTZMQsource.cxx
/// @author Mikolaj Krzewicki
/// @date
/// @brief  HLT ZMQ component implementation. */
///

#include "AliHLTZMQsource.h"
#include "AliHLTErrorGuard.h"
#include "AliLog.h"
#include <TPRegexp.h>
#include "zmq.h"
#include "AliHLTZMQhelpers.h"

using namespace std;

/** ROOT macro for the implementation of ROOT specific class methods */
ClassImp(AliHLTZMQsource)

//______________________________________________________________________________
AliHLTZMQsource::AliHLTZMQsource()
  : AliHLTComponent()
  , fOutputDataTypes()
  , fZMQcontext(NULL)
  , fZMQin(NULL)
  , fZMQsocketType(-1)
  , fZMQinConfig("SUB")
  , fMessageFilter("")
  , fZMQrequestTimeout(1000000)
  , fZMQneverBlock(kTRUE)
  , fForwardHLTinput(false)
  , fOutputBufferSize(10000000)
  , fZMQinit(NULL)
  , fZMQinitConfig("REQ")
  , fIncomingData()
  , fSkipSOR(kFALSE)
  , fOnlyOnDataEvents(kFALSE)
{
}

//______________________________________________________________________________
AliHLTZMQsource::~AliHLTZMQsource()
{
  //dtor
  zmq_close(fZMQin);
}

//______________________________________________________________________________
const char* AliHLTZMQsource::GetComponentID()
{
  // overloaded from AliHLTComponent
  return "ZMQsource";
}

//______________________________________________________________________________
void AliHLTZMQsource::GetInputDataTypes(AliHLTComponentDataTypeList& list)
{
  // overloaded from AliHLTComponent
  list.clear();
  list.push_back(kAliHLTAllDataTypes);
}

//______________________________________________________________________________
AliHLTComponentDataType AliHLTZMQsource::GetOutputDataType()
{
  // overloaded from AliHLTComponent
  return kAliHLTAllDataTypes;
}

////______________________________________________________________________________
//int AliHLTZMQsource::GetOutputDataTypes(AliHLTComponentDataTypeList& tgtList)
//{
//  // overloaded from AliHLTComponent
//  tgtList.assign(fOutputDataTypes.begin(), fOutputDataTypes.end());
//  HLTMessage("%s %p provides %d output data types", GetComponentID(), this, fOutputDataTypes.size());
//  return fOutputDataTypes.size();
//}

//______________________________________________________________________________
void AliHLTZMQsource::GetOutputDataSize( unsigned long& constBase, double& inputMultiplier )
{
  // overloaded from AliHLTComponent
  constBase=fOutputBufferSize;
  inputMultiplier=1.0;
}

//______________________________________________________________________________
AliHLTComponent* AliHLTZMQsource::Spawn()
{
  // overloaded from AliHLTComponent
  return new AliHLTZMQsource;
}

//______________________________________________________________________________
int AliHLTZMQsource::DoInit( int argc, const char** argv )
{
  // overloaded from AliHLTComponent: initialization
  int retCode=0;
  //process arguments
  if (ProcessOptionString(GetComponentArgs())<0)
  {
    HLTFatal("wrong options %s", GetComponentArgs().c_str());
    return -1;
  }

  int rc = 0;
  //init ZMQ stuff
  fZMQcontext = alizmq_context();
  HLTMessage(Form("ctx create rc %i %p errno %i",rc,fZMQcontext,errno));

  //in ZMQ socket
  rc = alizmq_socket_init(fZMQin, fZMQcontext, fZMQinConfig.Data(), 0, 10 );
  if (!fZMQin || rc<0)
  {
    HLTError("cannot initialize the in ZMQ socket %s, %s",fZMQinConfig.Data(),zmq_strerror(errno));
    return -1;
  }

  //subscribe
  rc = zmq_setsockopt(fZMQin, ZMQ_SUBSCRIBE, fMessageFilter.Data(), fMessageFilter.Length());

  HLTMessage(Form("in socket create ptr %p %s",fZMQin,(rc<0)?zmq_strerror(errno):""));
  HLTMessage(Form("  ZMQ connected to: %s rc %i %s",fZMQinConfig.Data(),rc,(rc<0)?zmq_strerror(errno):""));

  if (fZMQinitConfig.Length()>3)
  {
    //init ZMQ socket
    rc = alizmq_socket_init(fZMQinit, fZMQcontext, fZMQinitConfig.Data(), 0, 10 );
    if (!fZMQinit || rc<0)
    {
      HLTError("cannot initialize the init ZMQ socket %s, %s",fZMQinConfig.Data(),zmq_strerror(errno));
      return -1;
    }

    HLTMessage(Form("init socket create ptr %p %s",fZMQinit,(rc<0)?zmq_strerror(errno):""));
    HLTMessage(Form("  ZMQ connected to: %s rc %i %s",fZMQinitConfig.Data(),rc,(rc<0)?zmq_strerror(errno):""));

    //request data
    HLTMessage("requesting init data");
    aliZMQmsg request;
    alizmq_msg_add(&request, "", "");
    alizmq_msg_send(&request, fZMQinit, 0);
    alizmq_msg_close(&request);

    //get reply
    //need to use the poller here to control the interruptions and timeout behaviour
    zmq_pollitem_t sockets[] = { {fZMQinit, 0, ZMQ_POLLIN, 0} };
    rc = zmq_poll( sockets, 1, fZMQrequestTimeout );
    if (rc==-1)
    {
      HLTError("init request interrupted");
      //interrupted, stop processing
      return -1;
    }
    if (! (sockets[0].revents & ZMQ_POLLIN))
    {
      //if we got no reply reset the connection, probably source died
      HLTError("init request timed out after %i us",fZMQrequestTimeout);
      return -1;
    }
    HLTMessage("receiving init data");
    alizmq_msg_recv(&fIncomingData, fZMQinit, 0);
  }

  return retCode;
}

//______________________________________________________________________________
int AliHLTZMQsource::DoDeinit()
{
  // overloaded from AliHLTComponent: cleanup
  int retCode=0;
  retCode = alizmq_socket_close(fZMQin);
  retCode+= alizmq_socket_close(fZMQinit);
  return retCode;
}

//______________________________________________________________________________
int AliHLTZMQsource::DoProcessing( const AliHLTComponentEventData& evtData,
                  const AliHLTComponentBlockData* blocks,
                  AliHLTComponentTriggerData& /*trigData*/,
                  AliHLTUInt8_t* outputBuffer,
                  AliHLTUInt32_t& outputBufferSize,
                  AliHLTComponentBlockDataList& outputBlocks,
                  AliHLTComponentEventDoneData*& edd )
{
  // overloaded from AliHLTComponent: event processing
  int retCode=0;

  // process data events only
  //if (!IsDataEvent()) return 0;

  //handle non data events
  AliHLTUInt32_t eventType = 0;
  if (!IsDataEvent(&eventType)) {
    if (eventType==gkAliEventTypeStartOfRun && fSkipSOR)
    {
      return 0;
    }
    if (fOnlyOnDataEvents) { return 0; }
  }

  //if we want to forward HLT input:
  if (fForwardHLTinput)
  {
    for (int iBlock = 0;
        iBlock < evtData.fBlockCnt;
        iBlock++)
    {
      outputBlocks.push_back(blocks[iBlock]);
    }
  }

  int rc = -1;

  //in case we do requests: request first and poll for replies
  //if no reply arrives after a timeout period, reset the connection
  if (fZMQsocketType==ZMQ_REQ)
  {
    //send request (header + an empty body for good measure)
    HLTMessage("sending request");
    zmq_send(fZMQin, fMessageFilter.Data(), fMessageFilter.Length(), ZMQ_SNDMORE);
    zmq_send(fZMQin, 0, 0, 0);
    //wait for reply
    zmq_pollitem_t sockets[] = { {fZMQin, 0, ZMQ_POLLIN, 0} };
    rc = zmq_poll( sockets, 1, fZMQrequestTimeout );
    if (rc==-1)
    {
      HLTImportant("request interrupted");
      //interrupted, stop processing
      return 0;
    }
    if (! (sockets[0].revents & ZMQ_POLLIN))
    {
      //if we got no reply reset the connection, probably source died
      HLTImportant("request timed out after %i us",fZMQrequestTimeout);
      rc = alizmq_socket_init(fZMQin, fZMQcontext, fZMQinConfig.Data(), 0, 10 );
      if (rc<0)
      {
        HLTError("cannot reinitialize ZMQ socket %s, %s",fZMQinConfig.Data(),zmq_strerror(errno));
        return -1;
      }

      //just return normally
      return 0;
    }
  }

  //get data via ZMQ
  alizmq_msg_recv(&fIncomingData, fZMQin, (fZMQneverBlock)?ZMQ_DONTWAIT:0);

  //forward all parts to the HLT chain
  //init internal
  AliHLTUInt32_t initialOutputBufferCapacity = outputBufferSize;
  AliHLTUInt32_t outputBufferCapacity = outputBufferSize;
  outputBufferSize=0;
  int blockSize = 0;
  void* block = NULL;
  for (aliZMQmsg::iterator i=fIncomingData.begin(); i!=fIncomingData.end(); ++i)
  {
    outputBufferCapacity -= blockSize;
    outputBufferSize += blockSize;
    block = outputBuffer + outputBufferSize;

    //get the topic
    AliHLTDataTopic blockTopic = kAliHLTAnyDataType | kAliHLTDataOriginAny;
    if (alizmq_msg_iter_topic(i, blockTopic)!=0) continue; //no valid topic - skip also data

    //get the data
    void* dataIN = NULL;
    size_t dataINsize = 0;
    alizmq_msg_iter_data(i, dataIN, dataINsize);
    if (dataINsize > outputBufferCapacity)
    {
      HLTWarning("output buffer too small: %i, doubling size", initialOutputBufferCapacity);
      fOutputBufferSize = 2*fOutputBufferSize;
      retCode = -ENOSPC; break;
    }

    //copy data into the HLT block and set size
    std::memcpy(block, dataIN, dataINsize);
    blockSize = dataINsize;

    //prepare the component block data descriptor
    AliHLTComponentBlockData blockHeader; FillBlockData(blockHeader);
    blockHeader.fPtr      = outputBuffer;
    blockHeader.fOffset   = outputBufferSize;
    blockHeader.fSize     = blockSize;
    blockTopic.Fill(blockHeader.fDataType);
    blockHeader.fSpecification = blockTopic.fSpecification;

    //register the block in the output buffer list
    std::string datatypestr; datatypestr.assign(blockHeader.fDataType.fID,kAliHLTComponentDataTypefIDsize);
    HLTMessage(Form("pushing back %s, %i bytes", datatypestr.c_str(), blockSize));
    outputBlocks.push_back(blockHeader);

  }

  // after we sent the data, discard the ZMQ data
  alizmq_msg_close(&fIncomingData);

  edd=NULL;
  return retCode;
}

//______________________________________________________________________________
int AliHLTZMQsource::ProcessOption(TString option, TString value)
{
  //process option
  //to be implemented by the user

  if (option.EqualTo("in"))
  {
    fZMQinConfig = value;
    fZMQsocketType = alizmq_socket_type(value.Data());
    switch (fZMQsocketType)
    {
      case ZMQ_REQ:
        break;
      case ZMQ_PULL:
        break;
      case ZMQ_SUB:
        break;
      default:
        HLTWarning("use of socket type %s for a source is currently unsupported! (config: %s)", alizmq_socket_name(fZMQsocketType), fZMQinConfig.Data());
        return -EINVAL;
    }
  }
  else if (option.EqualTo("init"))
  {
    fZMQinitConfig = value;
    int socketType = alizmq_socket_type(value.Data());
    switch (socketType)
    {
      case ZMQ_REQ:
        break;
      default:
        HLTError("use of socket type %s for an init source is currently unsupported! (config: %s), must be REQ", alizmq_socket_name(fZMQsocketType), fZMQinConfig.Data());
        return -EINVAL;
    }
  }
  else if (option.EqualTo("MessageFilter") || option.EqualTo("subscription"))
  {
    fMessageFilter = value;
  }
  else if (option.EqualTo("ZMQrequestTimeout")) //this is in ms
  {
    fZMQrequestTimeout = round(value.Atoi()*1e3);
  }
  else if (option.EqualTo("ZMQneverBlock"))
  {
    if (value.EqualTo("0") || value.EqualTo("no") || value.Contains("false",TString::kIgnoreCase))
      fZMQneverBlock = kFALSE;
    else if (value.EqualTo("1") || value.EqualTo("yes") || value.Contains("true",TString::kIgnoreCase) )
      fZMQneverBlock = kTRUE;
  }
  else if (option.EqualTo("OutputBufferSize"))
  {
    fOutputBufferSize = value.Atoi();
  }
  else if (option.EqualTo("forwardHLTinput"))
  {
    fForwardHLTinput = true;
  }
  else if (option.EqualTo("SkipSOR"))
  {
    fSkipSOR = kTRUE;
  }
  else if (option.EqualTo("OnlyOnDataEvents"))
  {
    fOnlyOnDataEvents = !value.EqualTo("0");
  }
  else
  {
    HLTError("unrecognized option %s", option.Data());
    return -1;
  }

  return 1;
}

