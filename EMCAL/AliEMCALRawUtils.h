#ifndef ALIEMCALRAWUTILS_H
#define ALIEMCALRAWUTILS_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice     */

/* $Id$ */
/* History of cvs commits:
 *
 * $Log$
 * Revision 1.2  2007/09/03 20:55:35  jklay
 * EMCAL e-by-e reconstruction methods from Cvetan
 *
 * Revision 1.1  2007/03/17 19:56:38  mvl
 * Moved signal shape routines from AliEMCAL to separate class AliEMCALRawUtils to streamline raw data reconstruction code.
 *
 *
 */
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

class TGraph;
class TF1;
class AliRawReader;

class AliEMCALRawUtils : public TObject {
 public:
  AliEMCALRawUtils();
  virtual ~AliEMCALRawUtils();

  void Digits2Raw();
  void Raw2Digits(AliRawReader *reader,TClonesArray *digitsArr);
  void AddDigit(TClonesArray *digitsArr, Int_t id, Int_t lowGain, Int_t amp, Float_t time);

  // Signal shape parameters
  Double_t GetRawFormatHighLowGainFactor() const {return fHighLowGainFactor ;}
  static Int_t GetRawFormatOrder() { return fgOrder ; }   
  static Int_t GetRawFormatTimeBins() { return fgkTimeBins ; }    
  static Double_t GetRawFormatTimeMax() { return fgkTimeBins*fgTimeBinWidth; }   
  static Double_t GetRawFormatTimeBinWidth() { return fgTimeBinWidth; }   
  Double_t GetRawFormatTau() const { return fgTau ; }    
  Double_t GetRawFormatTimeTrigger() const { return fgTimeTrigger ; }
  Int_t GetRawFormatThreshold() const { return fgThreshold ; }       
  Int_t GetRawFormatDDLPerSuperModule() const { return fgDDLPerSuperModule ; } 

  // Signal shape functions
  void FitRaw(TGraph * gSig, TF1* signalF, Float_t & amp, Float_t & time);
  static Double_t RawResponseFunction(Double_t *x, Double_t *par); 
  Bool_t   RawSampledResponse(Double_t dtime, Double_t damp, Int_t * adcH, Int_t * adcL) const;  

  ClassDef(AliEMCALRawUtils,0)

 private:
  Double_t fHighLowGainFactor ;         // high to low gain factor for the raw RO signal
  static Int_t fgOrder ;                // order of the gamma function for the RO signal
  static Double_t fgTau ;                // tau parameter of gamma function for the RO signal
  static Double_t fgTimeTrigger ;       // time of the trigger for the RO signal 

  static const Int_t fgkTimeBins = 256 ; // number of sampling bins of the raw RO signal  
  static Double_t fgTimeBinWidth;       // maximum sampled time of the raw RO signal                             
  static Int_t fgThreshold;             // threshold
  static Int_t fgDDLPerSuperModule;     // number of DDL per SuperModule
};

#endif
