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

using namespace std;

ClassImp(AliHLTZMQsink)

//______________________________________________________________________________
AliHLTZMQsink::AliHLTZMQsink() :
  AliHLTDataSink()
  , fZMQcontext(NULL)
  , fZMQout(NULL)
  , fZMQsocketType(ZMQ_PUB)
  , fZMQconnectMode("bind")
  , fZMQendpoint("tcp://*:60201")
  , fZMQpollIn(kFALSE)
  , fPushbackDelayPeriod(-1)
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
void AliHLTZMQsink::GetInputDataTypes( vector<AliHLTComponentDataType>& list)
{
  //what data types do we accept
  list.clear();
  list.push_back(kAliHLTAllDataTypes);
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

  //process arguments
  ProcessOptionString(GetComponentArgs());

  int rc = 0;
  //init ZMQ stuff
  fZMQcontext = zmq_ctx_new();
  HLTMessage(Form("ctx create ptr %p errno %i",fZMQcontext,errno));
  if (!fZMQcontext) return -1;
  fZMQout = zmq_socket(fZMQcontext, fZMQsocketType); 
  HLTMessage(Form("socket create ptr %p errno %i",fZMQout,errno));
  if (!fZMQout) return -1;

  //set socket options
  int lingerValue = 10;
  rc = zmq_setsockopt(fZMQout, ZMQ_LINGER, &lingerValue, sizeof(int));
  HLTMessage(Form("setopt ZMQ_LINGER=%i   rc=%i errno=%i", lingerValue, rc, errno));
  int highWaterMarkSend = 100;
  rc = zmq_setsockopt(fZMQout, ZMQ_SNDHWM, &highWaterMarkSend, sizeof(int));
  HLTMessage(Form("setopt ZMQ_SNDHWM=%i   rc=%i errno=%i",highWaterMarkSend, rc, errno));
  int highWaterMarkRecv = 100;
  rc = zmq_setsockopt(fZMQout, ZMQ_RCVHWM, &highWaterMarkRecv, sizeof(int));
  HLTMessage(Form("setopt ZMQ_RCVHWM=%i   rc=%i errno=%i",highWaterMarkRecv, rc, errno));
  int rcvtimeo = 1000;
  rc = zmq_setsockopt(fZMQout, ZMQ_RCVTIMEO, &rcvtimeo, sizeof(int));
  HLTMessage(Form("setopt ZMQ_RCVTIMEO=%i rc=%i errno=%i",rcvtimeo, rc, errno));
  int sndtimeo = 1000;
  rc = zmq_setsockopt(fZMQout, ZMQ_SNDTIMEO, &sndtimeo, sizeof(int));
  HLTMessage(Form("setopt ZMQ_SNDTIMEO=%i rc=%i errno=%i",sndtimeo, rc, errno));

  //connect or bind, after setting socket options
  if (fZMQconnectMode.Contains("connect")) 
  {
    HLTMessage(Form("ZMQ connect to %s",fZMQendpoint.Data()));
    rc = zmq_connect(fZMQout,fZMQendpoint.Data());
    HLTMessage(Form("connect rc %i errno %i",rc,errno));
  }
  else 
  {
    HLTMessage(Form("ZMQ bind to %s",fZMQendpoint.Data()));
    rc = zmq_bind(fZMQout,fZMQendpoint.Data());
    HLTMessage(Form("bind rc=%i errno=%i",rc,errno));
  }
  return 0;
}

//______________________________________________________________________________
Int_t AliHLTZMQsink::DoDeinit()
{
  // see header file for class documentation
  return 0;
}

//______________________________________________________________________________
int AliHLTZMQsink::DumpEvent( const AliHLTComponentEventData& evtData,
                              const AliHLTComponentBlockData* blocks, 
                              AliHLTComponentTriggerData& trigData )
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

    //Bool_t nothingSelected=kTRUE; //if we select nothing, send a void reply
    int iBlock=0;
    const AliHLTComponentBlockData* inputBlock=NULL;
    for (inputBlock=GetFirstInputBlock();
        inputBlock!=NULL;
        inputBlock=GetNextInputBlock(), iBlock++) 
    {
      //check if the block is selected
      char blockTopic[kAliHLTComponentDataTypeTopicSize];
      DataType2Topic(inputBlock->fDataType, blockTopic);
      if (!Topicncmp(requestTopic, blockTopic, requestTopicSize)) continue;
      //nothingSelected=kFALSE;

      //send:
      //  first part : AliHLTComponentDataType in string format
      //  second part: Payload
      rc = zmq_send(fZMQout, &blockTopic, kAliHLTComponentDataTypeTopicSize, ZMQ_SNDMORE);
      HLTMessage(Form("send topic rc %i errno %i",rc,errno));
      rc = zmq_send(fZMQout, inputBlock->fPtr, inputBlock->fSize, ZMQ_SNDMORE);
      HLTMessage(Form("send data rc %i errno %i",rc,errno));
    }
    
    //send an empty message if we really need a reply (ZMQ_REP mode)
    //if (nothingSelected && fZMQsocketType==ZMQ_REP)
    //{ 
      rc = zmq_send(fZMQout, 0, 0, ZMQ_SNDMORE);
      HLTMessage(Form("send endframe rc %i errno %i",rc,errno));
      rc = zmq_send(fZMQout, 0, 0, 0);
      HLTMessage(Form("send endframe rc %i errno %i",rc,errno));
    //}
  }

  return retCode;
}

//______________________________________________________________________________
int AliHLTZMQsink::ProcessOption(TString option, TString value)
{
  //process option
  //to be implemented by the user
  
  //if (option.Contains("ZMQpollIn"))
  //{
  //  fZMQpollIn = (value.Contains("0"))?kFALSE:kTRUE;
  //}
 
  if (option.Contains("ZMQsocketMode")) 
  {
    if (value.Contains("PUB"))  fZMQsocketType=ZMQ_PUB;
    if (value.Contains("REP"))  fZMQsocketType=ZMQ_REP;
    if (value.Contains("PUSH")) fZMQsocketType=ZMQ_PUSH;
    
    //always poll when REPlying
    fZMQpollIn=(fZMQsocketType==ZMQ_REP)?kTRUE:kFALSE;
  }
 
  if (option.Contains("ZMQconnectMode"))
  {
    if (! (
          value.Contains("connect") ||
          value.Contains("bind")
          )
       ) {return 1;}
    fZMQconnectMode = value;
  }
 
  if (option.Contains("ZMQendpoint"))
  {
    fZMQendpoint = value;
  }

  if (option.Contains("pushback-period"))
  {
    HLTMessage(Form("Setting pushback delay to %i", atoi(value.Data())));
    fPushbackDelayPeriod = atoi(value.Data());
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
    HLTMessage("  %s : %s\n", i->first.data(), i->second.data());
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
