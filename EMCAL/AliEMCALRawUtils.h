// -*- mode: c++ -*-
#ifndef ALIEMCALRAWUTILS_H
#define ALIEMCALRAWUTILS_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice     */

/* $Id$ */


//_________________________________________________________________________
//  Utility Class for handling Raw data
//  Does all transitions from Digits to Raw and vice versa, 
//  for simu and reconstruction
//
//  Note: the current version is still simplified. Only 
//    one raw signal per digit is generated; either high-gain or low-gain
//    Need to add concurrent high and low-gain info in the future
//    No pedestal is added to the raw signal.
//
//*-- Author: Marco van Leeuwen (LBL)
//
#include "TObject.h" // for ROOT types
#include <TString.h>


class AliAnalysisManager;
class AliCaloRawStreamV3;
class AliAltroMapping;
class TGraph;
class AliRawReader;
class AliEMCALGeometry;
class AliCaloCalibPedestal;
class AliCaloRawAnalyzer;
class AliCaloRawAnalyzerLMSOffline;
class AliEMCALTriggerRawDigitMaker;
class AliEMCALTriggerData;

#include "AliCaloConstants.h"


class AliEMCALRawUtils : public TObject {
 public:
  AliEMCALRawUtils(Algo::fitAlgorithm fitAlgo = Algo::kStandard);
  virtual ~AliEMCALRawUtils();
  void Digits2Raw();
  void Raw2Digits(AliRawReader *reader, TClonesArray *digitsArr, const AliCaloCalibPedestal* pedbadmap,
				  TClonesArray *digitsTRG=0x0, AliEMCALTriggerData* trgData = 0x0);
  void AddDigit(TClonesArray *digitsArr, Int_t id, Int_t lowGain, Float_t amp, Float_t time, Float_t chi2, Int_t ndf);
  void TrimDigits(TClonesArray *digitsArr);
  Int_t    GetNoiseThreshold()             const { return fNoiseThreshold; }
  Int_t    GetNPedSamples()                const { return fNPedSamples; }
  Bool_t   GetRemoveBadChannels() const {return fRemoveBadChannels;}
  Int_t    GetFittingAlgorithm()  const {return fFittingAlgorithm; }
  Float_t  GetTimeMax()           const {return fTimeMax ;}
  Float_t  GetTimeMin()           const {return fTimeMin ;}
  Bool_t   UseFALTRO()            const {return fUseFALTRO; }
  void SetNoiseThreshold(Int_t val)                {fNoiseThreshold=val; }
  void SetNPedSamples(Int_t val)                   {fNPedSamples=val; }
  void SetRemoveBadChannels(Bool_t val)            {fRemoveBadChannels=val; }
  void SetFittingAlgorithm(Int_t val) ;             
  void SetTimeMin(Float_t t)                       {fTimeMin   = t          ;}
  void SetTimeMax(Float_t t)                       {fTimeMax   = t          ;}
  void SetFALTROUsage(Bool_t val)                  {fUseFALTRO=val; }
  AliCaloRawAnalyzer *GetRawAnalyzer()  const { return fRawAnalyzer;}
  virtual Option_t* GetOption() const { return fOption.Data(); }
  void SetOption(const Option_t* opt) { fOption = opt; }
  
private:
  AliEMCALRawUtils(const AliEMCALRawUtils& rawUtils);  //copy ctor
  AliEMCALRawUtils& operator =(const AliEMCALRawUtils& rawUtils);
  Int_t fNoiseThreshold;                // threshold to consider signal or noise
  Int_t fNPedSamples;                   // number of samples to use in pedestal calculation
  AliEMCALGeometry* fGeom;              // geometry
  AliAltroMapping*  fMapping[4];        // only two for now
  TString fOption;                      // option passed from Reconstructor
  Bool_t  fRemoveBadChannels;           // select if bad channels are removed before fitting
  Int_t   fFittingAlgorithm;            // select the fitting algorithm
  Float_t fTimeMin;                     // minimum threshold for the time of the signal
  Float_t fTimeMax;                     // maximum threshold for the time of the signal
  Bool_t  fUseFALTRO;			// use FALTRO and pass it to the digits
  AliCaloRawAnalyzer *fRawAnalyzer;     // e.g. for sample selection for fits
  AliEMCALTriggerRawDigitMaker* fTriggerRawDigitMaker;	
 
  ClassDef(AliEMCALRawUtils,7)          // utilities for raw signal fitting

};

#endif
