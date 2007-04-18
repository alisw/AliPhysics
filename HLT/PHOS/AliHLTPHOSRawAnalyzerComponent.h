#ifndef ALIHLTPHOSRAWANALYZERCOMPONENT_H
#define ALIHLTPHOSRAWANALYZERCOMPONENT_H

/* Copyright(c) 2006, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                          */

#include "AliHLTProcessor.h"
#include "AliHLTPHOSRawAnalyzer.h"
#include "AliRawReaderMemory.h"
#include "AliCaloRawStream.h"
#include "AliHLTPHOSDefinitions.h"
#include "AliHLTPHOSCommonDefs.h"


class AliHLTPHOSRcuCellEnergyDataStruct;

class AliHLTPHOSRawAnalyzerComponent: public AliHLTProcessor
{
 public:


  AliHLTPHOSRawAnalyzerComponent();
  virtual ~AliHLTPHOSRawAnalyzerComponent();
  AliHLTPHOSRawAnalyzerComponent(const AliHLTPHOSRawAnalyzerComponent & );
  AliHLTPHOSRawAnalyzerComponent & operator = (const AliHLTPHOSRawAnalyzerComponent &)
   {
      return *this;
   };

  virtual int DoInit( int argc, const char** argv );
  virtual int Deinit();
  virtual int DoDeinit();
  void DumpData(int gain);
  void DumpData();
  void DumpChannelData(Double_t *data); 
  void SetEquippmentID(AliHLTUInt16_t id);
  AliHLTUInt16_t  GetEquippmentID();
  void SetCoordinates(AliHLTUInt16_t equippmentID);
  virtual const char* GetComponentID() = 0;
  virtual void GetInputDataTypes(std::vector<AliHLTComponentDataType, std::allocator<AliHLTComponentDataType> >&);
  virtual AliHLTComponentDataType GetOutputDataType();
  virtual void GetOutputDataSize(unsigned long& constBase, double& inputMultiplier);
  virtual AliHLTComponent* Spawn() = 0;
  virtual int DoEvent(const AliHLTComponentEventData&, const AliHLTComponentBlockData*, AliHLTComponentTriggerData&, AliHLTUInt8_t*, AliHLTUInt32_t&, std::vector<AliHLTComponentBlockData, std::allocator<AliHLTComponentBlockData> >&);

 protected:
  AliHLTPHOSRawAnalyzer *fAnalyzerPtr; 

 private:
  void Reset();
  void ResetDataPtr();
  void ResetDataPtr(int sampleCnt);
  void ResetDataPtr(int startindex, int sampleCnt);
  static int fgEventCount;
  AliHLTUInt16_t fEquippmentID;
  AliHLTUInt8_t  fRcuX;
  AliHLTUInt8_t  fRcuZ;
  AliHLTUInt8_t  fRcuZOffset;
  AliHLTUInt8_t  fRcuXOffset;
  AliHLTUInt8_t  fModuleID;
  Double_t fTmpChannelData[ALTRO_MAX_SAMPLES];
  Double_t fMaxValues[N_MODULES][N_ROWS_MOD][N_COLUMNS_MOD][N_GAINS];
  AliCaloRawStream *fPHOSRawStream;
  AliRawReaderMemory *fRawMemoryReader;
  AliHLTPHOSRcuCellEnergyDataStruct* fOutPtr;
  static const AliHLTComponentDataType fgkInputDataTypes[];
};
#endif

