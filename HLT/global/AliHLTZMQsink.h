#ifndef ALIHLTZMQSINK_H
#define ALIHLTZMQSINK_H

/* This file is property of and copyright by the ALICE HLT Project        *
 * ALICE Experiment at CERN, All rights reserved.                         *
 * See cxx source for full Copyright notice                               */

/** @file    AliHLTZMQsink.h
    @author  Mikolaj Krzewicki
    @brief   ZeroMQ sink
*/

#include "AliHLTComponent.h"
#include "AliZMQhelpers.h"
#include <map>
#include <string>
#include "TRegexp.h"
#include "AliOptionParser.h"

using namespace AliZMQhelpers;

class AliHLTZMQsink : public AliHLTComponent, public AliOptionParser {
public:

  AliHLTZMQsink();
  virtual ~AliHLTZMQsink();

  //interface functions
  const char* GetComponentID();
  void GetInputDataTypes( vector<AliHLTComponentDataType>& list);
  AliHLTComponentDataType GetOutputDataType();
  TComponentType GetComponentType() { return AliHLTComponent::kSink;}
  void GetOutputDataSize( unsigned long& constBase, double& inputMultiplier );
  AliHLTComponent* Spawn();

  //overload from AliOptionParser
  int ProcessOption(TString option, TString value);

protected:
  //interface functions
  Int_t DoInit( Int_t /*argc*/, const Char_t** /*argv*/ );
  Int_t DoDeinit();
  int DoProcessing( const AliHLTComponentEventData& evtData,
                    const AliHLTComponentBlockData* blocks,
                    AliHLTComponentTriggerData& trigData,
                    AliHLTUInt8_t* outputPtr,
                    AliHLTUInt32_t& size,
                    AliHLTComponentBlockDataList& outputBlocks,
                    AliHLTComponentEventDoneData*& edd );

private:

  AliHLTZMQsink(const AliHLTZMQsink&);
  AliHLTZMQsink& operator=(const AliHLTZMQsink&);

  void* fZMQcontext;       //!ZMQ context pointer
  void* fZMQout;           //!the output socket
  int fZMQsocketType;      //ZMQ_REP,ZMQ_PUB,ZMQ_PUSH
  std::string fZMQoutConfig;   //config the ZMQ socket: e.g. SUB+tcp://localhost:123123
  Bool_t fZMQpollIn;       //do we poll fo incoming requests?
  Int_t fPushbackDelayPeriod;   //how often do we send?
  Int_t fLastPushbackDelayTime; //last push back time
  Bool_t fIncludePrivateBlocks; //include private blocks?
  Bool_t fZMQneverBlock;        //never block, even with a PUSH sock.
  Bool_t fSendRunNumber; //send the run number in each message
  unsigned long fNskippedErrorMessages; //keep the number of pushed messages
  Int_t fZMQerrorMsgSkip; //downscale ZMQ send warnings
  Bool_t fSendECSparamString; //send the full ECS param string
  std::string fECSparamString; //the ECS param string from SOR event
  stringMap fECSparamMap; //parsed map of ECS param strings
  Bool_t fSendStreamerInfos;   //send the cached streamer infos
  TRegexp fCDBpattern;          //keep the pattern for cdb string checks
  std::string fInfoString;  //cache the string with the run information

  ClassDef(AliHLTZMQsink, 1)
};
#endif
