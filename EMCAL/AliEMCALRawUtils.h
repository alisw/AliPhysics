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
//#include "AliCaloRawStreamV3.h"
class AliCaloRawStreamV3;
class AliAltroMapping;
class TGraph;
class TF1;
class AliRawReader;
class AliEMCALGeometry;
class AliCaloCalibPedestal;

class AliEMCALRawUtils : public TObject {
 public:
  AliEMCALRawUtils();
  AliEMCALRawUtils(AliEMCALGeometry *pGeometry);
  virtual ~AliEMCALRawUtils();

  AliEMCALRawUtils(const AliEMCALRawUtils& rawUtils);  //copy ctor
  AliEMCALRawUtils& operator =(const AliEMCALRawUtils& rawUtils);

  void Digits2Raw();
  void Raw2Digits(AliRawReader *reader,TClonesArray *digitsArr, AliCaloCalibPedestal* pedbadmap);

  void AddDigit(TClonesArray *digitsArr, Int_t id, Int_t lowGain, Int_t amp, Float_t time);

  // Signal shape parameters
  Double_t GetRawFormatHighLowGainFactor() const { return fHighLowGainFactor ;}
  Int_t GetRawFormatOrder()                const { return fOrder ; }   
  Double_t GetRawFormatTau()               const { return fTau ; }    
  Int_t GetNoiseThreshold()                const { return fNoiseThreshold; }
  Int_t GetNPedSamples()                   const { return fNPedSamples; }
  // get methods for fast fit simulation
  Double_t GetPedestalValue()  const {return fgPedestalValue;}
  Double_t GetFEENoise()       const {return fgFEENoise;}

  void SetRawFormatHighLowGainFactor(Double_t val) {fHighLowGainFactor=val;}
  void SetRawFormatOrder(Int_t val)                {fOrder=val; }   
  void SetRawFormatTau(Double_t val)               {fTau=val; }    
  void SetNoiseThreshold(Int_t val)                {fNoiseThreshold=val; }
  void SetNPedSamples(Int_t val)                   {fNPedSamples=val; }

  // set methods for fast fit simulation
  void SetFEENoise(Double_t val)                   {fgFEENoise = val;}
  void SetRawFormatTimeBins(Int_t val)             {fgTimeBins = val;}
  
  static Int_t GetRawFormatTimeBins() { return fgTimeBins ; }    
  static Double_t GetRawFormatTimeMax() { return fgTimeBins*fgTimeBinWidth; }   
  static Double_t GetRawFormatTimeBinWidth() { return fgTimeBinWidth; }   
  static Double_t GetRawFormatTimeBin() 
  { return GetRawFormatTimeMax() / GetRawFormatTimeBins(); }   
  Double_t GetRawFormatTimeTrigger() const { return fgTimeTrigger ; }
  Int_t GetRawFormatThreshold() const { return fgThreshold ; }       
  Int_t GetRawFormatDDLPerSuperModule() const { return fgDDLPerSuperModule ; } 

  virtual Option_t* GetOption() const { return fOption.Data(); }
  void SetOption(Option_t* opt) { fOption = opt; }

  // Signal shape functions
  void FitRaw(TGraph * gSig, TF1* signalF, const Int_t lastTimeBin, Float_t & amp, Float_t & time, Float_t & ped, Float_t & ampEstimate, Float_t & timeEstimate, Float_t & pedEstimate, const Float_t cut = 0) const ;
  static Double_t RawResponseFunction(Double_t *x, Double_t *par); 
  Bool_t   RawSampledResponse(Double_t dtime, Double_t damp, Int_t * adcH, Int_t * adcL) const;  


 private:
  Double_t fHighLowGainFactor ;         // high to low gain factor for the raw RO signal
  Int_t fOrder ;                        // order of the gamma function for the RO signal
  Double_t fTau ;                       // tau parameter of gamma function for the RO signal
  Int_t fNoiseThreshold;                // threshold to consider signal or noise
  Int_t fNPedSamples;                   // number of samples to use in pedestal calculation

  static const Int_t fgkOverflowCut = 950;  // cut to discriminate overflowed channels
  static const Int_t fgkRawSignalOverflow = 0x3FF; // maximum signal (10 bits)
  static Int_t fgTimeBins; // number of sampling bins of the raw RO signal

  static Double_t fgTimeTrigger ;       // time of the trigger for the RO signal 
  static Double_t fgTimeBinWidth;       // maximum sampled time of the raw RO signal                             
  static Int_t fgThreshold;             // threshold
  static Int_t fgDDLPerSuperModule;     // number of DDL per SuperModule
  static Int_t fgPedestalValue;         // pedestal value for Digits2Raw
  static Double_t fgFEENoise;           // electronics noise in ADC units

  AliEMCALGeometry* fGeom;         //geometry
  AliAltroMapping*  fMapping[4];   //only two for now

  TString fOption;                      //! option passed from Reconstructor

  ClassDef(AliEMCALRawUtils,3)          // utilities for raw signal fitting
};

#endif
