#ifndef ALIHLTPHOSRAWANALYZERCOMPONENT_H
#define ALIHLTPHOSRAWANALYZERCOMPONENT_H

/* Copyright(c) 2006, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                          */

#include "AliHLTProcessor.h"
#include "AliHLTPHOSRawAnalyzer.h"
#include "AliRawReaderMemory.h"
#include "AliCaloRawStream.h"
#include "AliHLTPHOSDefinitions.h"

class AliHLTPHOSRawAnalyzerComponent: public AliHLTProcessor
{
 public:
  AliHLTPHOSRawAnalyzerComponent();
  ~AliHLTPHOSRawAnalyzerComponent();
  AliHLTPHOSRawAnalyzerComponent(const AliHLTPHOSRawAnalyzerComponent & );
  AliHLTPHOSRawAnalyzerComponent & operator = (const AliHLTPHOSRawAnalyzerComponent &)
   {
      return *this;
   };

  virtual int DoInit( int argc, const char** argv );
  virtual int Deinit();
  virtual int DoDeinit();
  void DumpData();
  void SetEquippmentID(AliHLTUInt32_t id);
  int GetEquippmentID();
  void SetCoordinates( AliHLTUInt32_t equippmentID);
  virtual const char* GetComponentID() = 0;
  virtual void GetInputDataTypes(std::vector<AliHLTComponentDataType, std::allocator<AliHLTComponentDataType> >&);
  virtual AliHLTComponentDataType GetOutputDataType();
  virtual void GetOutputDataSize(unsigned long& constBase, double& inputMultiplier);
  virtual AliHLTComponent* Spawn() = 0;
  virtual int DoEvent(const AliHLTComponentEventData&, const AliHLTComponentBlockData*, AliHLTComponentTriggerData&, AliHLTUInt8_t*, AliHLTUInt32_t&, std::vector<AliHLTComponentBlockData, std::allocator<AliHLTComponentBlockData> >&);

 protected:
  AliHLTPHOSRawAnalyzer *analyzerPtr; 
  void Reset();
  void ResetDataPtr();

 private:
  int fEventCount;
  //  int fEquippmentId;
  //  AliHLTUInt16_t fEquippmentId;
  AliHLTUInt32_t fEquippmentID;
  AliHLTUInt8_t fRcuX;
  AliHLTUInt8_t fRcuY;
  AliHLTUInt8_t fModuleID;

  Double_t fTmpChannelData[1008];
  Double_t fMaxValues[5][64][56][2];
  AliCaloRawStream *fPHOSRawStream;
  AliRawReaderMemory *fRawMemoryReader;

  static const AliHLTComponentDataType inputDataTypes[];
  //  static const AliHLTComponentDataType outputDataType;
  
};
#endif
