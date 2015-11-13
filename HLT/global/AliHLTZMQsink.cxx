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
#include <TPRegexp.h>
#include "zmq.h"
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
  list.push_back(kAliHLTAllDataTypes);
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
  ProcessOptionString(GetComponentArgs());

  int rc = 0;
  //init ZMQ context
  fZMQcontext = zmq_ctx_new();
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
  HLTImportant(Form("ZMQ connected to: %s (%s(id %i)) rc %i %s",fZMQoutConfig.Data(),alizmq_socket_name(fZMQsocketType),fZMQsocketType,rc,(rc<0)?zmq_strerror(errno):""));
  
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
  int requestSize=-1;
  char request[kAliHLTComponentDataTypeTopicSize];
  memset(request, '*', kAliHLTComponentDataTypeTopicSize);

  int rc = 0;
  Bool_t doSend = kTRUE;
  
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
        if (more) {
          requestSize = zmq_recv(fZMQout, request, kAliHLTComponentDataTypeTopicSize, 0);
          zmq_getsockopt(fZMQout, ZMQ_RCVMORE, &more, &moreSize);
        }
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
      if (!fIncludePrivateBlocks)
      {
        if (!memcmp(inputBlock->fDataType.fOrigin, &kAliHLTDataOriginPrivate, kAliHLTComponentDataTypefOriginSize))
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

    //send the selected blocks
    for (int iSelectedBlock = 0;
         iSelectedBlock < selectedBlockIdx.size();
         iSelectedBlock++) 
    {
      inputBlock = &blocks[selectedBlockIdx[iSelectedBlock]];
      AliHLTDataTopic blockTopic = *inputBlock;

      //send:
      //  first part : AliHLTComponentDataType in string format
      //  second part: Payload
      rc = zmq_send(fZMQout, &blockTopic, sizeof(blockTopic), ZMQ_SNDMORE);
      HLTMessage(Form("send topic rc %i %s",rc,(rc<0)?zmq_strerror(errno):""));
      if (rc<0) HLTWarning("error sending topic frame %s, %s", blockTopic.Description().c_str(),zmq_strerror(errno));
      int flags = 0;
      if (fZMQneverBlock) flags = ZMQ_DONTWAIT;
      if (iSelectedBlock < (selectedBlockIdx.size()-1)) flags = ZMQ_SNDMORE;
      rc = zmq_send(fZMQout, inputBlock->fPtr, inputBlock->fSize, flags);
      if (rc<0) HLTWarning("error sending data frame %s, %s", blockTopic.Description().c_str(),zmq_strerror(errno));
      HLTMessage(Form("send data rc %i %s",rc,(rc<0)?zmq_strerror(errno):""));
    }
    
    //send an empty message if we really need a reply (ZMQ_REP mode)
    //only in case no blocks were selected
    if (selectedBlockIdx.size() == 0 && fZMQsocketType==ZMQ_REP)
    { 
      rc = zmq_send(fZMQout, 0, 0, ZMQ_SNDMORE);
      HLTMessage(Form("send endframe rc %i %s",rc,(rc<0)?zmq_strerror(errno):""));
      if (rc<0) HLTWarning("error sending dummy REP topic");
      rc = zmq_send(fZMQout, 0, 0, 0);
      HLTMessage(Form("send endframe rc %i %s",rc,(rc<0)?zmq_strerror(errno):""));
      if (rc<0) HLTWarning("error sending dummy REP data");
    }
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
        HLTWarning("use of socket type %s for a sink is currently unsupported! (config: %s)", alizmq_socket_name(fZMQsocketType), fZMQoutConfig.Data());
        return -EINVAL;
    }
  }

  if (option.EqualTo("pushback-period"))
  {
    HLTMessage(Form("Setting pushback delay to %i", atoi(value.Data())));
    fPushbackDelayPeriod = atoi(value.Data());
  }

  if (option.EqualTo("IncludePrivateBlocks"))
  {
    fIncludePrivateBlocks=kTRUE;
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
int AliHLTZMQsink::ProcessOptionString(TString arguments)
{
  //process passed options
  HLTMessage("Argument string: %s\n", arguments.Data());
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
AliHLTZMQsink::stringMap* AliHLTZMQsink::TokenizeOptionString(const TString str)
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
  //"-option option1=value1 --option2 value2, -option4=\'some string\'"
  
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
