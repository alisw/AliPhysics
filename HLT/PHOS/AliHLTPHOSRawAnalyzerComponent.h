#ifndef ALIHLTPHOSRAWANALYZERCOMPONENT_H
#define ALIHLTPHOSRAWANALYZERCOMPONENT_H

/* Copyright(c) 2006, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                          */

//
//Base class for PHOS HLT raw data analysis components
// see cxx file for more details

#include "AliHLTPHOSDefinitions.h"
#include "AliHLTPHOSCommonDefs.h"
#include "AliHLTPHOSProcessor.h"


class AliRawReaderMemory;
class AliCaloRawStream;
class AliHLTPHOSRawAnalyzer;
class AliHLTPHOSRcuCellEnergyDataStruct;
class AliHLTPHOSRcuChannelDataStruct;
class AliHLTPHOSMapper;
class AliHLTDDLDecoder;
class AliHLTAltroData;

class AliHLTPHOSRawAnalyzerComponent: public AliHLTPHOSProcessor
{
 public:
  AliHLTPHOSRawAnalyzerComponent();
  virtual ~AliHLTPHOSRawAnalyzerComponent();
  AliHLTPHOSRawAnalyzerComponent(const AliHLTPHOSRawAnalyzerComponent & );
  AliHLTPHOSRawAnalyzerComponent & operator = (const AliHLTPHOSRawAnalyzerComponent &)
   {
      return *this;
   };

  virtual int DoInit(int argc =0, const char** argv  = 0);
  virtual int Deinit();

  //  void DumpData(int gain =0) const;
  //  void DumpChannelData(Double_t *data =0) const; 

  virtual const char* GetComponentID() = 0;
  virtual void GetInputDataTypes( vector <AliHLTComponentDataType>& list);
  virtual AliHLTComponentDataType GetOutputDataType();
  virtual void GetOutputDataSize(unsigned long& constBase, double& inputMultiplier);
  virtual AliHLTComponent* Spawn() = 0; 
 protected:
  AliHLTPHOSRawAnalyzer *fAnalyzerPtr;  /**<Pointer to an analyzer object used for raw data anlysis*/ 

 private:
  virtual int DoEvent( const AliHLTComponentEventData& evtData, const AliHLTComponentBlockData* blocks, 
		     AliHLTComponentTriggerData& trigData, AliHLTUInt8_t* outputPtr, 
		     AliHLTUInt32_t& size, vector<AliHLTComponentBlockData>& outputBlocks ); 
  void Reset();
  void ResetDataPtr(int startindex = 0, int sampleCnt = 0);
  Bool_t fSendChannelData;       /**<wether or not to send raw data from the component into shared memory*/
  Double_t fTmpChannelData[ALTRO_MAX_SAMPLES];                        /**<temporary variable to store raw samples from a single altro channel*/
  Double_t fMaxValues[N_MODULES][N_ROWS_MOD][N_COLUMNS_MOD][N_GAINS]; /**<array to store cell energies*/
  AliCaloRawStream *fPHOSRawStream;                   /**<Streamer for PHOS raw data, used by fPHOSRawMemory reader*/ 
  AliRawReaderMemory *fRawMemoryReader;               /**<Decoder to read PHOS raw data on the altro format*/  
  AliHLTPHOSRcuCellEnergyDataStruct* fOutPtr;         /**<Pointer to outputbuffer to write results from the component into shared memory*/
 
  AliHLTPHOSMapper *fMapperPtr;

  AliHLTDDLDecoder *fDecoderPtr;
  AliHLTAltroData  *fAltroDataPtr;
};
#endif

