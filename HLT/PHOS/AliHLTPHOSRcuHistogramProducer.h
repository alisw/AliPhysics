//-*- Mode: C++ -*-
// $Id$

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


class TH1;
class TH2D;
class AliHLTPHOSRcuCellAccumulatedEnergyDataStruct;
class AliHLTPHOSUtilities; 

#define XBINLOW  0
#define XBINUP   1023
#define NBINS    1023


class AliHLTPHOSRcuHistogramProducer 
{
 public:
  //  AliHLTPHOSRcuHistogramProducer();
  AliHLTPHOSRcuHistogramProducer(AliHLTUInt8_t moduleID, AliHLTUInt8_t rcuX, AliHLTUInt8_t rcuZ);
  virtual ~AliHLTPHOSRcuHistogramProducer();
  const AliHLTPHOSRcuCellAccumulatedEnergyDataStruct& GetCellAccumulatedEnergies(); 
  void Init();
  void SetHistoOutDir(char *outDir);
  void FillEnergy(AliHLTUInt8_t x, AliHLTUInt8_t z,  AliHLTUInt8_t gain, float energy);
  void FillTime(AliHLTUInt8_t x,   AliHLTUInt8_t z,  AliHLTUInt8_t gain, float time); 
  void FillLiveChannels(Int_t data[], int size, Int_t x, Int_t z, Int_t gain);
  void FillLiveChannelHistograms();
  void Reset();
  void WriteAllHistograms(char opt[] = "update");

 
 private:
  AliHLTPHOSRcuHistogramProducer(const AliHLTPHOSRcuHistogramProducer & );
  AliHLTPHOSRcuHistogramProducer & operator = (const AliHLTPHOSRcuHistogramProducer &);
  void SetDefaultHistoOutDir(); 
  void ScanTimeString(char *timeString);
  AliHLTPHOSRcuHistogramProducer();
  char fHistoOutDir[512];
  TH1F *fEnergyHistogramPtrs[NXCOLUMNSRCU][NZROWSRCU][NGAINS];    /**<Array to store energy distribution per channel for one rcu*/
  TH1F *fTimingHistogramPtrs[NXCOLUMNSRCU][NZROWSRCU][NGAINS];    /**<Array to store timing distribution per channel for one rcu*/
  //  TH1D *fDeadChannelMapHistogramPtrs[NXCOLUMNSRCU][NZROWSRCU][NGAINS];
  TH2D *fDeadChannelMapHistogramPtrs[NGAINS];
  Float_t fEnergyAverageValues[NXCOLUMNSRCU][NZROWSRCU][NGAINS];  /**<Accumulated energy divided by  hits*/
  Double_t fAccumulatedValues[NXCOLUMNSRCU][NZROWSRCU][NGAINS];   /**<Array to store accumulated energy per channel for one rcu during run*/
 

  Float_t fTimingAverageValues[NXCOLUMNSRCU][NZROWSRCU][NGAINS];  /**<Avereage TOF*/
  AliHLTUInt32_t fHits[NXCOLUMNSRCU][NZROWSRCU][NGAINS]; //comment
  AliHLTUInt32_t fDeadChannelMap[NXCOLUMNSRCU][NZROWSRCU][NGAINS]; //comment
  Double_t fTmpChannelData[ALTROMAXSAMPLES];        /**<temporary variable to store raw samples from a single altro channel*/
  AliHLTPHOSRcuCellAccumulatedEnergyDataStruct fCellAccEnergy; //comment
  AliHLTUInt8_t fModuleID; /**<ID of the module this component read data from (0-4)*/
  AliHLTUInt8_t fRcuX;     /**<X position of RCU the data from this Equippment comes from (0 or 1)*/
  AliHLTUInt8_t fRcuZ;     /**<Z position of RCU the data from this Equippment comes from (0 or 1)*/

  AliHLTPHOSUtilities *fUtilitiesPtr;

};

#endif
