/**************************************************************************
 * This file is property of and copyright by the ALICE HLT Project        * 
 * All rights reserved.                                                   *
 *                                                                        *
 * Primary Authors: Oystein Djuvsland                                     *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          * 
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

#ifndef ALIHLTPHOSBASELINEANALYZER_H
#define ALIHLTPHOSBASELINEANALYZER_H

#include "AliHLTPHOSConstants.h"
#include "AliHLTPHOSBase.h"

class AliHLTPHOSRcuCellEnergyDataStruct;
class AliHLTPHOSValidCellDataStruct;
class AliHLTPHOSSanityInspector;
class TTree;
class TClonesArray;
class TH1F;
class TH2F;

using namespace PhosHLTConst;

class AliHLTPHOSBaselineAnalyzer : public AliHLTPHOSBase
{
  
public:
  AliHLTPHOSBaselineAnalyzer();
  
  virtual ~AliHLTPHOSBaselineAnalyzer();

  void CalculateRcuBaselines(AliHLTPHOSRcuCellEnergyDataStruct* rcuData);
  Float_t CalculateChannelBaseline(AliHLTPHOSValidCellDataStruct *cellData, Int_t xOff, Int_t zOff);
  void CalculateAccumulatedChannelBaseline(Int_t x, Int_t z, Int_t gain, Float_t baseline);
  void CalculateChannelsBaselineRMS();
  // void CalculateAccumulatedBaselines();
  void SetRootObjects(TTree *tree, TClonesArray *array);
  void FillTree();

  void WriteAccumulatedBaselines(const Char_t* filename);
  void WriteChannelHistograms(const Char_t* filename);
  void WriteRMSHistogram(const Char_t* filename);
  void ResetBaselines();
  void ResetChannelCount();
  void ResetAccumulatedBaselines();
  void SetNumberOfSamples(Int_t nSamples) { fNSamples = nSamples; }
  void SetMaxCrazyDifference(Int_t diff);
  void SetMaxSignal(Int_t max) { fMaxSignal = max; }
  
 // void SetChannelsHistograms(TH1F *channelLowGainHistArray[N_XCOLUMNS_MOD][N_ZROWS_MOD], TH1F *channelLowGainHistArray[N_XCOLUMNS_MOD][N_ZROWS_MOD]);
 
     
private:
  Int_t fSampleStart; //comment
  Int_t fMaxCrazyDifference; //comment
  Int_t fMaxSignal; //comment
  Int_t fChannelCount; //comment
  Float_t fBaselines[N_XCOLUMNS_MOD][N_ZROWS_MOD][N_GAINS]; //comment
  Float_t fAccumulatedBaselines[N_XCOLUMNS_MOD][N_ZROWS_MOD][N_GAINS][2]; //comment
  TTree *fTreePtr; //comment
  TClonesArray *fBaselineArrayPtr; //comment
 // TH1F *fChannelLowGainHistogramsPtr[N_XCOLUMNS_MOD][N_ZROWS_MOD]; //comment
 // TH1F *fChannelHighGainHistogramsPtr[N_XCOLUMNS_MOD][N_ZROWS_MOD]; //comment
  TH1F *fChannelHistogramsPtr[N_XCOLUMNS_MOD][N_ZROWS_MOD][N_GAINS]; //comment
  TH1F *fFixedChannelHistogramsPtr[N_XCOLUMNS_MOD][N_ZROWS_MOD][N_GAINS]; //comment
  TH1F *fRMSHistogramPtr; //comment
  TH2F *fRMSMapHighGainHistogramPtr; //comment
  TH2F *fRMSMapLowGainHistogramPtr; //comment
  TH1F *fFixedRMSHistogramPtr; //comment
  TH2F *fFixedRMSMapHighGainHistogramPtr; //comment
  TH2F *fFixedRMSMapLowGainHistogramPtr; //comment
  AliHLTPHOSSanityInspector *fSanityInspector;   //comment
 
  ClassDef(AliHLTPHOSBaselineAnalyzer, 1);
};

#endif
 
