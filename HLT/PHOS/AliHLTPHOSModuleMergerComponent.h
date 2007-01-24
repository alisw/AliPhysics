#ifndef ALIHLTPHOSMODULEMERGERCOMPONENT_H
#define ALIHLTPHOSMODULEMERGERCOMPONENT_H

/* Copyright(c) 2006, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice  */ 

#include "AliHLTProcessor.h"
#include "AliHLTPHOSDefinitions.h"


class AliHLTPHOSModuleMergerComponent: public AliHLTProcessor
{
 public:
  AliHLTPHOSModuleMergerComponent();
  ~AliHLTPHOSModuleMergerComponent();
  AliHLTPHOSModuleMergerComponent(const AliHLTPHOSModuleMergerComponent & );
  AliHLTPHOSModuleMergerComponent & operator = (const AliHLTPHOSModuleMergerComponent &)
   {
      return *this;
   };



  virtual int DoInit( int argc, const char** argv );
  virtual int Deinit();
  virtual int DoDeinit();
  void DumpData();
  void SetEquippmentId(int id);
  int GetEquippmentId();
  //  virtual const char* GetComponentID() = 0;
  virtual const char* GetComponentID();

  virtual void GetInputDataTypes(std::vector<AliHLTComponentDataType, std::allocator<AliHLTComponentDataType> >&);
  virtual AliHLTComponentDataType GetOutputDataType();
  virtual void GetOutputDataSize(unsigned long& constBase, double& inputMultiplier);
  //  virtual AliHLTComponent* Spawn() = 0;

  virtual int DoEvent(const AliHLTComponentEventData&, const AliHLTComponentBlockData*, AliHLTComponentTriggerData&, AliHLTUInt8_t*, AliHLTUInt32_t&, std::vector<AliHLTComponentBlockData, std::allocator<AliHLTComponentBlockData> >&);
  virtual AliHLTComponent* Spawn();
 protected:
  //  AliHLTPHOSRawAnalyzer *analyzerPtr; 
  void Reset();
  void ResetDataPtr();

 private:
  int fEventCount;
  AliHLTUInt32_t fEquippmentID;
  Double_t fTmpChannelData[1008];
  Double_t fMaxValues[5][64][56][2];
  //  TH2S *legoPlotPtr;

  //  AliCaloRawStream *fPHOSRawStream;
  //  AliRawReaderMemory *fRawMemoryReader;
  static const AliHLTComponentDataType inputDataTypes[];
  static const AliHLTComponentDataType outputDataType;
  
};

#endif
