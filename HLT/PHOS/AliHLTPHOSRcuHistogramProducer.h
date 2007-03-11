#ifndef ALIHLTPHOSRCUHISTOGRAMPRODUCER_H
#define ALIHLTPHOSRCUHISTOGRAMPRODUCER_H 

/* Copyright(c) 2006, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice  */ 

#include "AliHLTPHOSDefinitions.h"
#include "AliHLTPHOSCommonDefs.h"
#include "TH1.h"
//#include "AliHLTPHOSModuleCellAccumulatedEnergyDataStruct.h"
#include "AliHLTPHOSRcuCellAccumulatedEnergyDataStruct.h"

#define XBIN_LOW  0
#define XBIN_UP   1023
#define N_BINS    1023

class AliHLTPHOSRcuHistogramProducer
{
 public:
  AliHLTPHOSRcuHistogramProducer();
  AliHLTPHOSRcuHistogramProducer(AliHLTUInt8_t moduleID, AliHLTUInt8_t rcuX, AliHLTUInt8_t rcuZ);

  virtual ~AliHLTPHOSRcuHistogramProducer();
  int GetEquippmentId();
  const AliHLTPHOSRcuCellAccumulatedEnergyDataStruct& GetCellAccumulatedEnergies(); 
  int IncrementEventCounter();
  void Init();
  void SetEquippmentId(int id);
  void SetRcuX(AliHLTUInt8_t X);
  void SetRcuZ(AliHLTUInt8_t Z);
  void SetModuleID(AliHLTUInt8_t moduleID); 
  void FillEnergy(AliHLTUInt8_t x, AliHLTUInt8_t z,  AliHLTUInt8_t gain, float energy);
  void FillTime(AliHLTUInt8_t x,   AliHLTUInt8_t z,  AliHLTUInt8_t gain, float time); 
  void Reset();
  void WriteEnergyHistograms();


 protected:
 
 private:
  AliHLTPHOSRcuHistogramProducer(const AliHLTPHOSRcuHistogramProducer & );
  AliHLTPHOSRcuHistogramProducer & operator = (const AliHLTPHOSRcuHistogramProducer &)
   {
      return *this;
   };

  TH1F *fEnergyHistogramPtrs[N_XCOLUMNS_RCU][N_ZROWS_RCU][N_GAINS];
  TH1F *fTimingHistogramPtrs[N_XCOLUMNS_RCU][N_ZROWS_RCU][N_GAINS];

  Float_t fEnergyAverageValues[N_XCOLUMNS_RCU][N_ZROWS_RCU][N_GAINS];  
  Double_t fAccumulatedValues[N_XCOLUMNS_RCU][N_ZROWS_RCU][N_GAINS];
  Float_t fTimingAverageValues[N_XCOLUMNS_RCU][N_ZROWS_RCU][N_GAINS]; 
  AliHLTUInt32_t fHits[N_XCOLUMNS_RCU][N_ZROWS_RCU][N_GAINS];
  int fEventCount;
  AliHLTUInt32_t fEquippmentID;
  Double_t fTmpChannelData[ALTRO_MAX_SAMPLES];
  //  AliHLTPHOSModuleCellAccumulatedEnergyDataStruct*  fOutPtr;
  AliHLTPHOSRcuCellAccumulatedEnergyDataStruct fCellAccEnergy;
  //  AliHLTPHOSRcuCellAccumulatedEnergyDataStruct*  fOutPtr;
  AliHLTUInt8_t fModuleID;
  AliHLTUInt8_t fRcuX;
  AliHLTUInt8_t fRcuZ;
};

#endif
