#ifndef ALIHLTPHOSDDLDECODERCOMPONENT_H
#define ALIHLTPHOSDDLDECODERCOMPONENT_H

/* Copyright(c) 2006, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                          */

#include "AliHLTPHOSProcessor.h"
#include "AliHLTPHOSDefinitions.h"
#include "AliHLTPHOSCommonDefs.h"


class AliHLTPHOSRcuCellEnergyDataStruct;
class AliHLTPHOSRcuChannelDataStruct;
class AliHLTPHOSPulseGenerator;
class AliHLTPHOSDataCorruptor;
class AliRawReaderMemory;
class AliCaloRawStream;

class AliHLTPHOSDDLDecoderComponent:public AliHLTPHOSProcessor
{
 public:
  AliHLTPHOSDDLDecoderComponent();
  virtual ~AliHLTPHOSDDLDecoderComponent();
  AliHLTPHOSDDLDecoderComponent(const AliHLTPHOSDDLDecoderComponent & );
  AliHLTPHOSDDLDecoderComponent & operator = (const AliHLTPHOSDDLDecoderComponent &)
    {
      return *this;
    };
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
  AliCaloRawStream *fPHOSRawStream;                            /**<Streamer for PHOS raw data, used by fPHOSRawMemory reader*/ 
  AliRawReaderMemory *fRawMemoryReader;                        /**<Decoder to read PHOS raw data on the altro format*/ 
  AliHLTPHOSRcuChannelDataStruct*  fOutPtr;                    /**<Pointer to outputbuffer to write results from the component into shared memory*/
};
#endif

