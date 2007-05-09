#ifndef ALIHLTPHOSRCUHISTOGRAMPRODUCER_H
#define ALIHLTPHOSRCUHISTOGRAMPRODUCER_H 

/* Copyright(c) 2006, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice  */ 

#include "AliHLTPHOSDefinitions.h"
#include "AliHLTPHOSCommonDefs.h"
#include "TH1.h"
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
  const AliHLTPHOSRcuCellAccumulatedEnergyDataStruct& GetCellAccumulatedEnergies(); 
  void Init();
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

  TH1F *fEnergyHistogramPtrs[N_XCOLUMNS_RCU][N_ZROWS_RCU][N_GAINS];    /**<Array to store energy distribution per channel for one rcu*/
  TH1F *fTimingHistogramPtrs[N_XCOLUMNS_RCU][N_ZROWS_RCU][N_GAINS];    /**<Array to store timing distribution per channel for one rcu*/

  Float_t fEnergyAverageValues[N_XCOLUMNS_RCU][N_ZROWS_RCU][N_GAINS];  /**<Accumulated energy divided by  hits*/
  Double_t fAccumulatedValues[N_XCOLUMNS_RCU][N_ZROWS_RCU][N_GAINS];   /**<Array to store accumulated energy per channel for one rcu during run*/
  Float_t fTimingAverageValues[N_XCOLUMNS_RCU][N_ZROWS_RCU][N_GAINS];  /**<Avereage TOF*/
  AliHLTUInt32_t fHits[N_XCOLUMNS_RCU][N_ZROWS_RCU][N_GAINS];
  Double_t fTmpChannelData[ALTRO_MAX_SAMPLES];        /**<temporary variable to store raw samples from a single altro channel*/
  AliHLTPHOSRcuCellAccumulatedEnergyDataStruct fCellAccEnergy;
  AliHLTUInt8_t fModuleID; /**<ID of the module this component read data from (0-4)*/
  AliHLTUInt8_t fRcuX;     /**<X position of RCU the data from this Equippment comes from (0 or 1)*/
  AliHLTUInt8_t fRcuZ;     /**<Z position of RCU the data from this Equippment comes from (0 or 1)*/
};

#endif
