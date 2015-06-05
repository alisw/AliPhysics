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

#include "AliHLTDataSource.h"
#include <TList.h>
#include <map>
#include <string>

class TFile;

class AliHLTZMQsource : public AliHLTDataSource  {
  public:
    AliHLTZMQsource();
    virtual ~AliHLTZMQsource();

    typedef map<std::string,std::string> stringMap;

    const char* GetComponentID();
    AliHLTComponentDataType GetOutputDataType();
    int GetOutputDataTypes(AliHLTComponentDataTypeList& tgtList);
    void GetOutputDataSize( unsigned long& constBase, double& inputMultiplier );
    AliHLTComponent* Spawn();

    //new option parser
    static stringMap* TokenizeOptionString(const TString str);
    int ProcessOptionString(TString arguments);
    int ProcessOption(TString option, TString value);

  protected:
    virtual int DoInit( int argc, const char** argv );
    int DoDeinit();

    int GetEvent( const AliHLTComponentEventData& evtData,
        AliHLTComponentTriggerData& trigData,
        AliHLTUInt8_t* outputPtr, 
        AliHLTUInt32_t& size,
        AliHLTComponentBlockDataList& outputBlocks );

    using AliHLTDataSource::GetEvent;

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
    int fZMQsocketType;      //ZMQ_REP,ZMQ_PUB,ZMQ_PUSH
    TString fZMQconnectMode; //"connect" or "bind"
    TString fZMQendpoint;    //e.g. "tcp://*:60100" "tcp://ecs0:60100"
    TString fMessageFilter;   //ZMQ subscription

    ClassDef(AliHLTZMQsource, 0)
};
#endif
