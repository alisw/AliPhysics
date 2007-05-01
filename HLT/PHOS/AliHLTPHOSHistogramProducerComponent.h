#ifndef ALIHLTPHOSHISTOGRAMPRODUCERCOMPONENT_H
#define ALIHLTPHOSHISTOGRAMPRODUCERCOMPONENT_H 


/* Copyright(c) 2006, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice  */ 

#include "AliHLTProcessor.h"
#include "AliHLTPHOSDefinitions.h"
#include "AliHLTPHOSCommonDefs.h"
//#include "TH2.h"
//#include "AliHLTPHOSCommonDefs.h"
//#include "AliHLTPHOSModuleCellAccumulatedEnergyDataStruct.h"

//class TH2;
class AliHLTPHOSModuleCellAccumulatedEnergyDataStruct;

class AliHLTPHOSHistogramProducerComponent:public AliHLTProcessor
{
 public:
  AliHLTPHOSHistogramProducerComponent();
  
  //  ~AliHLTPHOSHistogramProducerComponent();
  virtual ~AliHLTPHOSHistogramProducerComponent();
  AliHLTPHOSHistogramProducerComponent(const AliHLTPHOSHistogramProducerComponent & );
  AliHLTPHOSHistogramProducerComponent & operator = (const AliHLTPHOSHistogramProducerComponent &)
   {
      return *this;
   };
  virtual int DoInit( int argc = 0, const char** argv = 0);
  virtual int Deinit();
  virtual int DoDeinit();
  virtual int DoEvent(const AliHLTComponentEventData&, const AliHLTComponentBlockData*, AliHLTComponentTriggerData&, AliHLTUInt8_t*, AliHLTUInt32_t&, std::vector<AliHLTComponentBlockData, std::allocator<AliHLTComponentBlockData> >&);
  void DumpData(int gain = 0);
  int GetEquippmentId();
  virtual const char* GetComponentID();
  virtual void GetInputDataTypes(std::vector<AliHLTComponentDataType, std::allocator<AliHLTComponentDataType> >&);
  virtual AliHLTComponentDataType GetOutputDataType();
  virtual void GetOutputDataSize(unsigned long& constBase, double& inputMultiplier);
  void SetEquippmentId(int id = 0);
  virtual AliHLTComponent* Spawn();
  
 protected:
  void Reset();
  void ResetDataPtr();

 private:
  Double_t fEnergyAverageValues[N_ZROWS_MOD][N_XCOLUMNS_MOD][N_GAINS];  
  Double_t fAccumulatedValues[N_ZROWS_MOD][N_XCOLUMNS_MOD][N_GAINS];
  Double_t fTimingAverageValues[N_ZROWS_MOD][N_XCOLUMNS_MOD][N_GAINS]; 
  AliHLTUInt32_t fHits[N_ZROWS_MOD][N_XCOLUMNS_MOD][N_GAINS];
  int fEventCount;
  AliHLTUInt32_t fEquippmentID;
  Double_t fTmpChannelData[ALTRO_MAX_SAMPLES];
  AliHLTPHOSModuleCellAccumulatedEnergyDataStruct*  fOutPtr;
  static const AliHLTComponentDataType fgkInputDataTypes[];
  static const AliHLTComponentDataType fgkOutputDataType;
};

#endif
