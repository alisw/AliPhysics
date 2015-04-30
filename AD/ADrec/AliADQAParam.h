/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

#ifndef ALIADQAPARAM_H
#define ALIADQAPARAM_H

#include "TNamed.h"

class AliADQAParam : public TNamed
{
 public: 
  AliADQAParam();
  virtual ~AliADQAParam();
  
  //HPTDC time
  void SetNTdcTimeBins (Int_t val) { fNTdcTimeBins = val; }
  void SetTdcTimeMin (Float_t val) { fTdcTimeMin = val; }
  void SetTdcTimeMax (Float_t val) { fTdcTimeMax = val; }
  //Width
  void SetNTdcWidthBins (Int_t val) { fNTdcWidthBins = val; }
  void SetTdcWidthMin (Float_t val) { fTdcWidthMin = val; }
  void SetTdcWidthMax (Float_t val) { fTdcWidthMax = val; }
  //Charge per channel
  void SetNChargeChannelBins (Int_t val) { fNChargeChannelBins = val; }
  void SetChargeChannelMin (Int_t val) { fChargeChannelMin = val; }
  void SetChargeChannelMax (Int_t val) { fChargeChannelMax = val; }
  void SetChargeChannelZoomMin (Int_t val) { fChargeChannelZoomMin = val; }
  void SetChargeChannelZoomMax (Int_t val) { fChargeChannelZoomMax = val; }
  //Charge per side
  void SetNChargeSideBins (Int_t val) { fNChargeSideBins = val; }
  void SetChargeSideMin (Int_t val) { fChargeSideMin = val; }
  void SetChargeSideMax (Int_t val) { fChargeSideMax = val; }
  //Chharge correlation
  void SetNChargeCorrBins (Int_t val) { fNChargeCorrBins = val; }
  void SetChargeCorrMin (Int_t val) { fChargeCorrMin = val; }
  void SetChargeCorrMax (Int_t val) { fChargeCorrMax = val; }
  //Pair time correlation
  void SetNPairTimeCorrBins(Int_t val) { fNPairTimeCorrBins = val;}
  void SetPairTimeCorrMin(Float_t val) { fPairTimeCorrMin = val;}
  void SetPairTimeCorrMax(Float_t val) { fPairTimeCorrMax = val;} 
  //Pair time diference
  void SetNPairTimeDiffBins(Int_t val) { fNPairTimeDiffBins = val;}
  void SetPairTimeDiffMin(Float_t val) { fPairTimeDiffMin = val;}
  void SetPairTimeDiffMax(Float_t val) { fPairTimeDiffMax = val;} 
  //Mean time correlation
  void SetNMeanTimeCorrBins(Int_t val) { fNMeanTimeCorrBins = val;}
  void SetMeanTimeCorrMin(Float_t val) { fMeanTimeCorrMin = val;}
  void SetMeanTimeCorrMax(Float_t val) { fMeanTimeCorrMax = val;} 
  
  //QA checker limits
  void SetSatMed(Float_t val) {fSatMed = val;}
  void SetSatHigh(Float_t val) {fSatHigh = val;}
  void SetSatHuge(Float_t val) {fSatHuge = val;} 
  
  void SetMaxPedDiff(Int_t val) {fMaxPedDiff = val;} 
  void SetMaxPedWidth(Float_t val) {fMaxPedWidth = val;} 
  
  void SetMaxNoTimeRate(Float_t val) {fMaxNoTimeRate = val;}
  void SetMaxNoFlagRate(Float_t val) {fMaxNoFlagRate = val;}
  void SetMaxBBVariation(Float_t val) {fMaxBBVariation = val;}
  void SetMaxBGVariation(Float_t val) {fMaxBGVariation = val;}

  //HPTDC time
  Int_t GetNTdcTimeBins() const { return fNTdcTimeBins; }
  Float_t GetTdcTimeMin() const { return fTdcTimeMin; }
  Float_t GetTdcTimeMax() const { return fTdcTimeMax; }
  //Width
  Int_t GetNTdcWidthBins() const { return fNTdcWidthBins; }
  Float_t GetTdcWidthMin() const { return fTdcWidthMin; }
  Float_t GetTdcWidthMax() const { return fTdcWidthMax; }
  //Charge per channel
  Int_t GetNChargeChannelBins() const { return fNChargeChannelBins; }
  Int_t GetChargeChannelMin() const { return fChargeChannelMin; }
  Int_t GetChargeChannelMax() const { return fChargeChannelMax; }
  Int_t GetChargeChannelZoomMin() const { return fChargeChannelZoomMin; }
  Int_t GetChargeChannelZoomMax() const { return fChargeChannelZoomMax; }
  //Charge per side
  Int_t GetNChargeSideBins() const { return fNChargeSideBins; }
  Int_t GetChargeSideMin() const { return fChargeSideMin; }
  Int_t GetChargeSideMax() const { return fChargeSideMax; }
  //Charge correlation - be careful with nBins
  Int_t GetNChargeCorrBins() const { return fNChargeCorrBins; }
  Int_t GetChargeCorrMin() const { return fChargeCorrMin; }
  Int_t GetChargeCorrMax() const { return fChargeCorrMax; }
  //Pair time correlation 
  Int_t GetNPairTimeCorrBins() const { return fNPairTimeCorrBins;}
  Float_t GetPairTimeCorrMin() const { return fPairTimeCorrMin;}
  Float_t GetPairTimeCorrMax() const { return fPairTimeCorrMax;} 
  //Pair time difference
  Int_t GetNPairTimeDiffBins() const { return fNPairTimeDiffBins;}
  Float_t GetPairTimeDiffMin() const { return fPairTimeDiffMin;}
  Float_t GetPairTimeDiffMax() const { return fPairTimeDiffMax;} 
  //Mean time correlation
  Int_t  GetNMeanTimeCorrBins() const { return fNMeanTimeCorrBins;}
  Float_t GetMeanTimeCorrMin() const { return fMeanTimeCorrMin;}
  Float_t GetMeanTimeCorrMax() const { return fMeanTimeCorrMax;}
  
  //QA checker limits
  Float_t GetSatMed() const {return fSatMed;}
  Float_t GetSatHigh() const {return fSatHigh;}
  Float_t GetSatHuge() const {return fSatHuge;} 
  
  Int_t GetMaxPedDiff() const {return fMaxPedDiff;} 
  Float_t GetMaxPedWidth() const {return fMaxPedWidth;}
  
  Float_t GetMaxNoTimeRate() const {return fMaxNoTimeRate;}
  Float_t GetMaxNoFlagRate() const {return fMaxNoFlagRate;}
  Float_t GetMaxBBVariation() const {return fMaxBBVariation;}
  Float_t GetMaxBGVariation() const {return fMaxBGVariation;}

 private:
 
  //QA histogram bins and limits
  Int_t	fNTdcTimeBins;	//Time bining
  Float_t fTdcTimeMin;	
  Float_t fTdcTimeMax;
  Int_t fNTdcWidthBins; //Width binning
  Float_t fTdcWidthMin;
  Float_t fTdcWidthMax;
  Int_t fNChargeChannelBins;	//Charge binnings
  Int_t fChargeChannelMin;
  Int_t fChargeChannelMax;
  Int_t fChargeChannelZoomMin;
  Int_t fChargeChannelZoomMax;
  Int_t fNChargeSideBins;
  Int_t fChargeSideMin;
  Int_t fChargeSideMax;
  Int_t fNChargeCorrBins;
  Int_t fChargeCorrMin;
  Int_t fChargeCorrMax;
  
  Int_t   fNPairTimeCorrBins; 
  Float_t fPairTimeCorrMin;
  Float_t fPairTimeCorrMax; 
   
  Int_t   fNPairTimeDiffBins; 
  Float_t fPairTimeDiffMin;
  Float_t fPairTimeDiffMax;
  
  Int_t   fNMeanTimeCorrBins; 
  Float_t fMeanTimeCorrMin;
  Float_t fMeanTimeCorrMax;
  
  //QA checker limits
  Float_t fSatMed;
  Float_t fSatHigh;
  Float_t fSatHuge;
  
  Int_t fMaxPedDiff;
  Float_t fMaxPedWidth;
  
  Float_t fMaxNoTimeRate;
  Float_t fMaxNoFlagRate;
  Float_t fMaxBBVariation;
  Float_t fMaxBGVariation;
  

  ClassDef(AliADQAParam,3)
};
#endif
