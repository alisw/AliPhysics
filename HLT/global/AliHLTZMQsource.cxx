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
#include "AliZMQhelpers.h"

using namespace std;

/** ROOT macro for the implementation of ROOT specific class methods */
ClassImp(AliHLTZMQsource)

//______________________________________________________________________________
AliHLTZMQsource::AliHLTZMQsource()
  : AliHLTComponent()
  , fOutputDataTypes()
  , fZMQcontext(NULL)
  , fZMQin(NULL)
  , fZMQsocketType(ZMQ_SUB)
  , fZMQendpoint(">tcp://localhost:60201")
  , fMessageFilter("")
  , fZMQrequestTimeout(1000)
  , fZMQneverBlock(kTRUE)
{
}

//______________________________________________________________________________
AliHLTZMQsource::~AliHLTZMQsource()
{
  //dtor
  zmq_close(fZMQin);
  zmq_ctx_destroy(fZMQcontext);
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
  constBase=10000000;
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
  ProcessOptionString(GetComponentArgs());
  HLTMessage("options processed\n");

  int rc = 0;
  //init ZMQ stuff
  fZMQcontext = zmq_ctx_new();
  HLTMessage(Form("ctx create rc %i errno %i",rc,errno));
  fZMQin = zmq_socket(fZMQcontext, fZMQsocketType); 
  HLTMessage(Form("socket create rc %i errno %i",rc,errno));

  //set socket options
  //default subscription to all
  rc = zmq_setsockopt(fZMQin, ZMQ_SUBSCRIBE, fMessageFilter.Data(), fMessageFilter.Length());
  HLTMessage(Form("setopt ZMQ_SUBSCRIBE=\'%s\' arglength=%i rc=%i errno=%i", fMessageFilter.Data(), fMessageFilter.Length(), rc, errno));
  
  int lingerValue = 10;
  rc = zmq_setsockopt(fZMQin, ZMQ_LINGER, &lingerValue, sizeof(lingerValue));
  HLTMessage(Form("setopt ZMQ_LINGER=%i rc=%i errno=%i", lingerValue, rc, errno));
  int highWaterMarkSend = 20;
  rc = zmq_setsockopt(fZMQin, ZMQ_SNDHWM, &highWaterMarkSend, sizeof(highWaterMarkSend));
  HLTMessage(Form("setopt ZMQ_SNDHWM=%i rc=%i errno=%i",highWaterMarkSend, rc, errno));
  int highWaterMarkRecv = 20;
  rc = zmq_setsockopt(fZMQin, ZMQ_RCVHWM, &highWaterMarkRecv, sizeof(highWaterMarkRecv));
  HLTMessage(Form("setopt ZMQ_RCVHWM=%i rc=%i errno=%i",highWaterMarkRecv, rc, errno));
  int rcvtimeo = 0;
  rc = zmq_setsockopt(fZMQin, ZMQ_RCVTIMEO, &rcvtimeo, sizeof(rcvtimeo));
  HLTMessage(Form("setopt ZMQ_RCVTIMEO=%i rc=%i errno=%i",rcvtimeo, rc, errno));
  int sndtimeo = 0;
  rc = zmq_setsockopt(fZMQin, ZMQ_SNDTIMEO, &sndtimeo, sizeof(sndtimeo));
  HLTMessage(Form("setopt ZMQ_SNDTIMEO=%i rc=%i errno=%i",sndtimeo, rc, errno));

  //connect or bind, after setting socket options
  HLTMessage(Form("ZMQ connect to %s",fZMQendpoint.Data()));
  rc = alizmq_attach(fZMQin,fZMQendpoint.Data());
  if (rc==-1) retCode=-1;
  HLTMessage(Form("  connect rc %i errno %i",rc,errno));

  return retCode;
}

//______________________________________________________________________________
int AliHLTZMQsource::DoDeinit()
{
  // overloaded from AliHLTComponent: cleanup
  int retCode=0;
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
  if (!IsDataEvent()) return 0;

  //init internal
  AliHLTUInt32_t outputBufferCapacity = outputBufferSize;
  outputBufferSize=0;
  int blockSize = 0;
  void* block = NULL;

  int blockTopicSize=-1;
  AliHLTDataTopic blockTopic;
  int rc = -1;
 
  //in case we do requests: request first and poll for replies
  //if no reply arrives after a timeout period, reset the connection
  if (fZMQsocketType==ZMQ_REQ)
  {
    //send request (header + an empty body for good measure)
    zmq_send(fZMQin, fMessageFilter.Data(), fMessageFilter.Length(), ZMQ_SNDMORE);
    zmq_send(fZMQin, 0, 0, 0);
    //wait for reply
    zmq_pollitem_t sockets[] = { {fZMQin, 0, ZMQ_POLLIN, 0} };
    rc = zmq_poll( sockets, 1, fZMQrequestTimeout );
    if (rc==-1) 
    {
      //interrupted, stop processing
      return 0;
    }
    if (! (sockets[0].revents & ZMQ_POLLIN))
    {
      //if we got no reply reset the connection, probably source died
      int lingerValue = 0;
      rc = zmq_setsockopt(fZMQin, ZMQ_LINGER, &lingerValue, sizeof(lingerValue));
      rc = zmq_close(fZMQin);
      fZMQin = zmq_socket(fZMQcontext, fZMQsocketType); 
      rc = alizmq_attach(fZMQin,fZMQendpoint.Data());
      
      //just return normally
      return 0;
    }
  }

  int64_t more=0;
  size_t moreSize=sizeof(more);
  do //multipart, get all parts
  {
    outputBufferCapacity -= blockSize;
    outputBufferSize += blockSize;
    block = outputBuffer + outputBufferSize;
    
    //get (fill) the block topic
    blockTopicSize = zmq_recv (fZMQin, &blockTopic, sizeof(blockTopic), (fZMQneverBlock)?ZMQ_DONTWAIT:0);
    if (blockTopicSize<0 && errno==EAGAIN) break; //nothing on the socket
    zmq_getsockopt(fZMQin, ZMQ_RCVMORE, &more, &moreSize);
    if (more) {
      //get (fill) the block data
      blockSize = zmq_recv(fZMQin, block, outputBufferCapacity, (fZMQneverBlock)?ZMQ_DONTWAIT:0);
      if (blockSize < 0 && errno == EAGAIN) break; //nothing on the socket
      if (blockSize > outputBufferCapacity) {retCode = ENOSPC; break;}//no space for message
      zmq_getsockopt(fZMQin, ZMQ_RCVMORE, &more, &moreSize);
    }

    if (blockTopicSize <= 0) continue; //empty header, dont push back

    HLTMessage(Form("pushing back %s, %i bytes", blockTopic.Description().c_str(), blockSize));
    
    //prepare the component block data descriptor
    AliHLTComponentBlockData blockHeader; FillBlockData(blockHeader);
    blockHeader.fPtr      = outputBuffer;
    blockHeader.fOffset   = outputBufferSize;
    blockHeader.fSize     = blockSize;
    blockHeader.fDataType = blockTopic.fTopic;
    blockHeader.fSpecification = blockTopic.fSpecification;
    
    //register the block in the output buffer list
    outputBlocks.push_back(blockHeader);
  } while (more==1);
  
  edd=NULL;
  return retCode;
}

//______________________________________________________________________________
int AliHLTZMQsource::ProcessOption(TString option, TString value)
{
  //process option
  //to be implemented by the user
  
  if (option.EqualTo("ZMQsocketMode")) 
  {
    if (value.EqualTo("SUB"))  fZMQsocketType=ZMQ_SUB;
    if (value.EqualTo("PULL")) fZMQsocketType=ZMQ_PULL;
    if (value.EqualTo("REQ"))  fZMQsocketType=ZMQ_REQ;
  }
 
  if (option.EqualTo("ZMQendpoint"))
  {
    fZMQendpoint = value;
  }

  if (option.EqualTo("MessageFilter"))
  {
    fMessageFilter = value;
  }

  if (option.EqualTo("ZMQrequestTimeout"))
  {
    fZMQrequestTimeout = value.Atoi();
  }

  if (option.EqualTo("ZMQneverBlock"))
  {
    if (value.EqualTo("0") || value.EqualTo("no") || value.Contains("false",TString::kIgnoreCase))
      fZMQneverBlock = kFALSE;
    else if (value.EqualTo("1") || value.EqualTo("yes") || value.Contains("true",TString::kIgnoreCase) )
      fZMQneverBlock = kTRUE;
  }

  return 1; 
}

////////////////////////////////////////////////////////////////////////////////

//______________________________________________________________________________
int AliHLTZMQsource::ProcessOptionString(TString arguments)
{
  //process passed options
  HLTMessage("Argument string: %s", arguments.Data());
  stringMap* options = TokenizeOptionString(arguments);
  for (stringMap::iterator i=options->begin(); i!=options->end(); ++i)
  {
    HLTMessage("  %s : %s", i->first.data(), i->second.data());
    ProcessOption(i->first,i->second);
  }
  delete options; //tidy up

  return 1; 
}

//______________________________________________________________________________
AliHLTZMQsource::stringMap* AliHLTZMQsource::TokenizeOptionString(const TString str)
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

  for (stringMap::iterator i=options->begin(); i!=options->end(); ++i)
  {
    Printf("%s : %s", i->first.data(), i->second.data());
  }
  return options;
}
