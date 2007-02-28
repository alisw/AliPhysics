#ifndef ALIHLTPHOSMODULEMERGERCOMPONENT_H
#define ALIHLTPHOSMODULEMERGERCOMPONENT_H

/* Copyright(c) 2006, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice  */ 

#include "AliHLTProcessor.h"
#include "AliHLTPHOSDefinitions.h"
#include "AliHLTPHOSCommonDefs.h"

class AliHLTPHOSModuleMergerComponent:public AliHLTProcessor
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
  virtual int DoEvent(const AliHLTComponentEventData&, const AliHLTComponentBlockData*, AliHLTComponentTriggerData&, AliHLTUInt8_t*, AliHLTUInt32_t&, std::vector<AliHLTComponentBlockData, std::allocator<AliHLTComponentBlockData> >&);
  
  void DumpData(int gain);
 
  int GetEquippmentId();
  virtual const char* GetComponentID();
  virtual void GetInputDataTypes(std::vector<AliHLTComponentDataType, std::allocator<AliHLTComponentDataType> >&);
  virtual AliHLTComponentDataType GetOutputDataType();
  virtual void GetOutputDataSize(unsigned long& constBase, double& inputMultiplier);
  void SetEquippmentId(int id);
  virtual AliHLTComponent* Spawn();

 protected:
  void Reset();
  void ResetDataPtr();

 private:
  int fEventCount;
  AliHLTUInt32_t fEquippmentID;
  Double_t fTmpChannelData[ALTRO_MAX_SAMPLES];
  Double_t fMaxValues[N_MODULES][N_ROWS_MOD][N_COLUMNS_MOD][N_GAINS];
  static const AliHLTComponentDataType inputDataTypes[];
  static const AliHLTComponentDataType outputDataType;
};

#endif
