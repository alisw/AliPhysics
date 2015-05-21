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
#include <TObject.h>
#include <TPRegexp.h>
#include "zmq.h"
#include "AliLog.h"

using namespace std;

ClassImp(AliHLTZMQsink)

//______________________________________________________________________________
AliHLTZMQsink::AliHLTZMQsink() :
  AliHLTDataSink()
  , fZMQcontext(zmq_ctx_new())
  , fZMQout(NULL)
  , fZMQsocketType(ZMQ_PUB)
  , fZMQconnectMode("bind")
  , fZMQendpoint("tcp://*:60201")
  , fZMQpollIn(kTRUE)
  , fZMQsendAllInOne(kTRUE)
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
  list.push_back(kAliHLTAnyDataType);
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

  //init sockets
  fZMQout = zmq_socket(fZMQcontext, fZMQsocketType); 

  int rc = 0;
  
  //set socket options
  int lingerValue = 10;
  rc = zmq_setsockopt(fZMQout, ZMQ_LINGER, &lingerValue, sizeof(lingerValue));
  int highWaterMarkSend = 2;
  rc = zmq_setsockopt(fZMQout, ZMQ_SNDHWM, &highWaterMarkSend, sizeof(highWaterMarkSend));
  int highWaterMarkRecv = 2;
  rc = zmq_setsockopt(fZMQout, ZMQ_RCVHWM, &highWaterMarkRecv, sizeof(highWaterMarkRecv));

  //connect or bind, after setting socket options
  if (fZMQconnectMode.Contains("connect")) 
  {
    rc = zmq_connect(fZMQout,fZMQendpoint);
  }
  else if(fZMQconnectMode.Contains("bind"))
  {
    rc = zmq_bind(fZMQout,fZMQendpoint);
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
  Int_t iResult=0;
  
  Bool_t doSend = kTRUE;

  //if enabled (option -pushback-period), send at most so often
  //if (fPushbackPeriod>0)
  //{
  //  TDatime time;
  //  if (fLastPushBackTime<0 || (int)time.Get()-fLastPushBackTime<fPushbackPeriod) 
  //    doSend=kFALSE;
  //}

  //create a default selection of any data:
  char requestedTopic[kAliHLTComponentDataTypeTopicSize];
  memset(requestedTopic, '*', kAliHLTComponentDataTypeTopicSize);

  //in case we reply to requests instead of just pushing/publishing
  //we poll for requests
  int requestedTopicSize=0;
  if (fZMQpollIn)
  {
    zmq_pollitem_t items[] = { { fZMQout, 0, ZMQ_POLLIN, 0 } };
    zmq_poll(items, 1, 0);

    if (items[0].revents & ZMQ_POLLIN)
    {
      requestedTopicSize = zmq_recv (fZMQout, requestedTopic, kAliHLTComponentDataTypeTopicSize, 0);
    }
    else { doSend = kFALSE; }
  } 
 
  if (doSend)
  {
    Bool_t nothingSelected=kTRUE; //if we select nothing, send a void reply
    int iBlock=0;
    const AliHLTComponentBlockData* inputBlock=NULL;
    for (inputBlock=GetFirstInputBlock();
        inputBlock!=NULL;
        inputBlock=GetNextInputBlock(), iBlock++) 
    {
      //check if the block is selected
      char blockTopic[kAliHLTComponentDataTypeTopicSize];
      DataType2Topic(inputBlock->fDataType, blockTopic);
      if (!Topicncmp(requestedTopic, blockTopic, requestedTopicSize)) continue;
      nothingSelected=kFALSE;

      char topic[kAliHLTComponentDataTypeTopicSize];
      DataType2Topic(inputBlock->fDataType, topic);

      //send:
      //  first part : AliHLTComponentDataType in string format
      //  second part: Payload
      zmq_send(fZMQout, &topic, kAliHLTComponentDataTypeTopicSize, ZMQ_SNDMORE);
      zmq_send(fZMQout, inputBlock->fPtr, inputBlock->fSize, (fZMQsendAllInOne)?ZMQ_SNDMORE:0);
    }
    if (nothingSelected && fZMQsocketType==ZMQ_REP)
    { 
      //empty frame of type void if we really need a reply (ZMQ_REP mode)
      char delimiter[kAliHLTComponentDataTypeTopicSize];
      memset(delimiter, 0, kAliHLTComponentDataTypeTopicSize);
      zmq_send(fZMQout, &delimiter, kAliHLTComponentDataTypeTopicSize, ZMQ_SNDMORE);
      zmq_send(fZMQout, 0, 0, 0);
    }
  }

  return iResult;
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
    if (value.Contains("REP")) fZMQsocketType=ZMQ_PUB;
    if (value.Contains("PUB")) fZMQsocketType=ZMQ_REP;
    
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

  if (option.Contains("ZMQsendAllInOne"))
  {
    fZMQsendAllInOne = (value.Contains("0"))?kFALSE:kTRUE;
  }

  return 1; 
}

////////////////////////////////////////////////////////////////////////////////

//______________________________________________________________________________
int AliHLTZMQsink::ProcessOptionString(TString arguments)
{
  //process passed options
  HLTInfo("Argument string: %s\n", arguments.Data());
  stringMap* options = TokenizeOptionString(arguments);
  for (stringMap::iterator i=options->begin(); i!=options->end(); ++i)
  {
    HLTInfo("  %s : %s\n", i->first.data(), i->second.data());
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

  for (stringMap::iterator i=options->begin(); i!=options->end(); ++i)
  {
    Printf("%s : %s", i->first.data(), i->second.data());
  }
  return options;
}
