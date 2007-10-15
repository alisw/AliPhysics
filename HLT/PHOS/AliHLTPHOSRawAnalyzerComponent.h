#ifndef ALIHLTPHOSRAWANALYZERCOMPONENT_H
#define ALIHLTPHOSRAWANALYZERCOMPONENT_H

/* Copyright(c) 2006, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                          */
#include "AliHLTPHOSRcuProcessor.h"

class AliHLTPHOSRawAnalyzer;
class AliHLTPHOSRcuCellEnergyDataStruct;
class AliHLTPHOSRcuChannelDataStruct;
class AliHLTPHOSMapper;
class AliHLTDDLDecoder;
class AliHLTAltroData;
class AliHLTAltroBunch;
class AliHLTPHOSSanityInspector;

class AliHLTPHOSRawAnalyzerComponent: public AliHLTPHOSRcuProcessor
{
 public:
  AliHLTPHOSRawAnalyzerComponent();
  virtual ~AliHLTPHOSRawAnalyzerComponent();
  virtual int DoInit(int argc =0, const char** argv  = 0);
  virtual int Deinit();
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
  void SetBaselines(const char* baselineFile);
  Bool_t fSendChannelData;       /**<wether or not to send raw data from the component into shared memory*/
  Double_t fTmpChannelData[ALTRO_MAX_SAMPLES];                        /**<temporary variable to store raw samples from a single altro channel*/
  Double_t fMaxValues[N_MODULES][N_ZROWS_MOD][N_XCOLUMNS_MOD][N_GAINS]; /**<array to store cell energies*/
  // AliHLTPHOSRcuCellEnergyDataStruct* fOutPtr;         /**<Pointer to outputbuffer to write results from the component into shared memory*/

  AliHLTPHOSRcuCellEnergyDataStruct* fOutPtr;
  AliHLTPHOSMapper *fMapperPtr;
  AliHLTDDLDecoder *fDecoderPtr;
  AliHLTAltroData  *fAltroDataPtr;
  AliHLTAltroBunch *fAltroBunchPtr;
  AliHLTPHOSSanityInspector *fSanityInspectorPtr;
  Bool_t fUseBaselineSubtraction;
  Float_t fBaselines[N_XCOLUMNS_RCU][N_ZROWS_RCU][N_GAINS];
  int fDebugCnt;
  
};
#endif

