// 1
// 2
// 3
// 4
// 5

#ifndef ALIHLTPHOSRCUHISTOGRAMPRODUCER_H
#define ALIHLTPHOSRCUHISTOGRAMPRODUCER_H 

/* Copyright(c) 2006, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice  */ 

#include "AliHLTPHOSDefinitions.h"
#include "TH1.h"
#include "TH2D.h"
#include "AliHLTPHOSRcuCellAccumulatedEnergyDataStruct.h"
//#include "AliHLTPHOSBase.h"
#include "AliHLTPHOSRcuProcessor.h"

class TH1;
class TH2D;
class AliHLTPHOSRcuCellAccumulatedEnergyDataStruct;

#define XBIN_LOW  0
#define XBIN_UP   1023
#define N_BINS    1023

class AliHLTPHOSRcuHistogramProducer : public AliHLTPHOSBase
//class AliHLTPHOSRcuHistogramProducer : public AliHLTPHOSRcuProcessor
{
 public:
  //  AliHLTPHOSRcuHistogramProducer();
  AliHLTPHOSRcuHistogramProducer(AliHLTUInt8_t moduleID, AliHLTUInt8_t rcuX, AliHLTUInt8_t rcuZ);
  virtual ~AliHLTPHOSRcuHistogramProducer();
  const AliHLTPHOSRcuCellAccumulatedEnergyDataStruct& GetCellAccumulatedEnergies(); 
  void Init();
  void SetRcuX(AliHLTUInt8_t X);
  void SetRcuZ(AliHLTUInt8_t Z);
  void SetModuleID(AliHLTUInt8_t moduleID); 
  void SetHistoOutDir(char *outDir);
  void FillEnergy(AliHLTUInt8_t x, AliHLTUInt8_t z,  AliHLTUInt8_t gain, float energy);
  void FillTime(AliHLTUInt8_t x,   AliHLTUInt8_t z,  AliHLTUInt8_t gain, float time); 

  void FillLiveChannels(Int_t data[], int size, Int_t x, Int_t z, Int_t gain);
  void FillLiveChannelHistograms();

  void Reset();
  void WriteAllHistograms(char opt[] = "update");

 protected:
 
 private:
  void SetDefaultHistoOutDir(); 
  void ScanTimeString(char *timeString);

  AliHLTPHOSRcuHistogramProducer();
  char fHistoOutDir[512];


  TH1F *fEnergyHistogramPtrs[N_XCOLUMNS_RCU][N_ZROWS_RCU][N_GAINS];    /**<Array to store energy distribution per channel for one rcu*/
  TH1F *fTimingHistogramPtrs[N_XCOLUMNS_RCU][N_ZROWS_RCU][N_GAINS];    /**<Array to store timing distribution per channel for one rcu*/
  //  TH1D *fDeadChannelMapHistogramPtrs[N_XCOLUMNS_RCU][N_ZROWS_RCU][N_GAINS];
  TH2D *fDeadChannelMapHistogramPtrs[N_GAINS];

  Float_t fEnergyAverageValues[N_XCOLUMNS_RCU][N_ZROWS_RCU][N_GAINS];  /**<Accumulated energy divided by  hits*/
  Double_t fAccumulatedValues[N_XCOLUMNS_RCU][N_ZROWS_RCU][N_GAINS];   /**<Array to store accumulated energy per channel for one rcu during run*/
  Float_t fTimingAverageValues[N_XCOLUMNS_RCU][N_ZROWS_RCU][N_GAINS];  /**<Avereage TOF*/
  AliHLTUInt32_t fHits[N_XCOLUMNS_RCU][N_ZROWS_RCU][N_GAINS]; //comment
  AliHLTUInt32_t fDeadChannelMap[N_XCOLUMNS_RCU][N_ZROWS_RCU][N_GAINS]; //comment
  Double_t fTmpChannelData[ALTRO_MAX_SAMPLES];        /**<temporary variable to store raw samples from a single altro channel*/
  AliHLTPHOSRcuCellAccumulatedEnergyDataStruct fCellAccEnergy; //comment
  AliHLTUInt8_t fModuleID; /**<ID of the module this component read data from (0-4)*/
  AliHLTUInt8_t fRcuX;     /**<X position of RCU the data from this Equippment comes from (0 or 1)*/
  AliHLTUInt8_t fRcuZ;     /**<Z position of RCU the data from this Equippment comes from (0 or 1)*/
};

#endif
