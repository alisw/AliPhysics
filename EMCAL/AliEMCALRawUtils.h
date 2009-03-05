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
//#include "AliCaloRawStream.h"
class AliCaloRawStream;
class AliAltroMapping;
class TGraph;
class TF1;
class AliRawReader;
class AliEMCALGeometry;

class AliEMCALRawUtils : public TObject {
 public:
  AliEMCALRawUtils();
  AliEMCALRawUtils(AliEMCALGeometry *pGeometry);
  virtual ~AliEMCALRawUtils();

  AliEMCALRawUtils(const AliEMCALRawUtils& rawUtils);  //copy ctor
  AliEMCALRawUtils& operator =(const AliEMCALRawUtils& rawUtils);

  void Digits2Raw();
  void Raw2Digits(AliRawReader *reader,TClonesArray *digitsArr);

  void AddDigit(TClonesArray *digitsArr, Int_t id, Int_t lowGain, Int_t amp, Float_t time);

  // Signal shape parameters
  Double_t GetRawFormatHighLowGainFactor() const { return fHighLowGainFactor ;}
  Int_t GetRawFormatOrder()                const { return fOrder ; }   
  Double_t GetRawFormatTau()               const { return fTau ; }    
  Int_t GetNoiseThreshold()                const { return fNoiseThreshold; }
  Int_t GetNPedSamples()                   const { return fNPedSamples; }

  void SetRawFormatHighLowGainFactor(Double_t val) {fHighLowGainFactor=val;}
  void SetRawFormatOrder(Int_t val)                {fOrder=val; }   
  void SetRawFormatTau(Double_t val)               {fTau=val; }    
  void SetNoiseThreshold(Int_t val)                {fNoiseThreshold=val; }
  void SetNPedSamples(Int_t val)                   {fNPedSamples=val; }
  
  static Int_t GetRawFormatTimeBins() { return fgkTimeBins ; }    
  static Double_t GetRawFormatTimeMax() { return fgkTimeBins*fgTimeBinWidth; }   
  static Double_t GetRawFormatTimeBinWidth() { return fgTimeBinWidth; }   
  static Double_t GetRawFormatTimeBin() 
  { return GetRawFormatTimeMax() / GetRawFormatTimeBins(); }   
  Double_t GetRawFormatTimeTrigger() const { return fgTimeTrigger ; }
  Int_t GetRawFormatThreshold() const { return fgThreshold ; }       
  Int_t GetRawFormatDDLPerSuperModule() const { return fgDDLPerSuperModule ; } 

  virtual Option_t* GetOption() const { return fOption.Data(); }
  void SetOption(Option_t* opt) { fOption = opt; }

  // Signal shape functions
  void FitRaw(TGraph * gSig, TF1* signalF, Float_t & amp, Float_t & time) const ;
  static Double_t RawResponseFunction(Double_t *x, Double_t *par); 
  Bool_t   RawSampledResponse(Double_t dtime, Double_t damp, Int_t * adcH, Int_t * adcL) const;  


 private:
  Double_t fHighLowGainFactor ;         // high to low gain factor for the raw RO signal
  Int_t fOrder ;                        // order of the gamma function for the RO signal
  Double_t fTau ;                       // tau parameter of gamma function for the RO signal
  Int_t fNoiseThreshold;                // threshold to consider signal or noise
  Int_t fNPedSamples;                   // number of samples to use in pedestal calculation

  static const Int_t fgkOverflowCut = 950;  // cut to discriminate overflowed channels
  static const Int_t fgkTimeBins = 256 ; // number of sampling bins of the raw RO signal  
  static const Int_t fgkRawSignalOverflow = 0x3FF; // maximum signal (10 bits)

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
