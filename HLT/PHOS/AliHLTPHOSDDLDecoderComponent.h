// 1
// 2
// 3
// 4
// 5
#ifndef ALIHLTPHOSDDLDECODERCOMPONENT_H
#define ALIHLTPHOSDDLDECODERCOMPONENT_H

/* Copyright(c) 2006, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                          */

#include "AliHLTPHOSRcuProcessor.h"



class AliHLTPHOSRcuCellEnergyDataStruct;
class AliHLTPHOSRcuChannelDataStruct;
class AliHLTPHOSPulseGenerator;
class AliHLTPHOSDataCorruptor;
class AliHLTDDLDecoder;
class AliHLTAltroData;
class AliHLTPHOSMapper;



class AliHLTPHOSDDLDecoderComponent:public AliHLTPHOSRcuProcessor
{
 public:
  AliHLTPHOSDDLDecoderComponent();
  virtual ~AliHLTPHOSDDLDecoderComponent();
  virtual int DoInit( int argc, const char** argv );
  virtual int Deinit();
  virtual const char* GetComponentID();
  virtual void GetInputDataTypes( vector <AliHLTComponentDataType>& list);
  virtual AliHLTComponentDataType GetOutputDataType();
  virtual void GetOutputDataSize(unsigned long& constBase, double& inputMultiplier);
  virtual AliHLTComponent* Spawn();
  
  virtual int DoEvent( const AliHLTComponentEventData& evtData, const AliHLTComponentBlockData* blocks, 
		     AliHLTComponentTriggerData& trigData, AliHLTUInt8_t* outputPtr, 
		     AliHLTUInt32_t& size, vector<AliHLTComponentBlockData>& outputBlocks ); 


 private:
  AliHLTPHOSDataCorruptor *fDataCorruptorPtr;                  /**<Pointer to data corruptor*/
  Double_t fTmpChannelData[ALTRO_MAX_SAMPLES];                 /**<Temporary storage for altro dat from a single channel*/
  AliHLTPHOSRcuChannelDataStruct*  fOutPtr;                    /**<Pointer to outputbuffer to write results from the component into shared memory*/
  AliHLTDDLDecoder *fDecoderPtr; //comment
  AliHLTAltroData *fAltroDataPtr; //comment
  AliHLTPHOSMapper   *fMapperPtr; //comment
};
#endif

