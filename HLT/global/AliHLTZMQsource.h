#ifndef __AliHLTZMQsource__
#define __AliHLTZMQsource__
//* This file is property of and copyright by the                          * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//* See cxx source for full Copyright notice                               *

//  @file   AliHLTZMQsource.h
//  @author Mikolaj Krzewicki, mikolaj.krzewicki@cern.ch
//  @date   
//  @brief  An HLT ZMQ data source component.
// 

#include "AliHLTComponent.h"
#include "AliZMQhelpers.h"
#include <TList.h>

class TFile;

class AliHLTZMQsource : public AliHLTComponent, public AliOptionParser  {
  public:
    AliHLTZMQsource();
    virtual ~AliHLTZMQsource();

    const char* GetComponentID();
    AliHLTComponentDataType GetOutputDataType();
    void GetInputDataTypes( AliHLTComponentDataTypeList& list);
    void GetOutputDataSize( unsigned long& constBase, double& inputMultiplier );
    TComponentType GetComponentType() { return AliHLTComponent::kSource;}
    AliHLTComponent* Spawn();

    //overload AliOptionParser
    int ProcessOption(TString option, TString value);

  protected:
    virtual int DoInit( int argc, const char** argv );
    int DoDeinit();

    int DoProcessing( const AliHLTComponentEventData& evtData,
                      const AliHLTComponentBlockData* blocks, 
                      AliHLTComponentTriggerData& trigData,
                      AliHLTUInt8_t* outputPtr, 
                      AliHLTUInt32_t& size,
                      AliHLTComponentBlockDataList& outputBlocks,
                      AliHLTComponentEventDoneData*& edd );

  private:
    /** prohibit copy constructor */
    AliHLTZMQsource(const AliHLTZMQsource&);
    /** prohibit assignment operator */
    AliHLTZMQsource& operator=(const AliHLTZMQsource&);

  protected:
    /** output data types  */
    AliHLTComponentDataTypeList fOutputDataTypes;                    //! transient
    void* fZMQcontext;       //!ZMQ context pointer
    void* fZMQin;           //!the output socket
    int fZMQsocketType;     //!cache the value of the socket type
    TString fZMQinConfig;    //e.g. "PUSH@tcp://*:60100" "PULL>tcp://ecs0:60100"
    TString fMessageFilter;   //ZMQ subscription
    ULong_t fZMQrequestTimeout;  //timeout in ms
    Bool_t fZMQneverBlock;    //dont block on receive
    Bool_t fForwardHLTinput;  //forward everything on HLT input as well as ZMQ input

    ClassDef(AliHLTZMQsource, 0)
};
#endif
