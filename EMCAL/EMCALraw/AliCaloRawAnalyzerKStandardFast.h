#ifndef ALICALORAWANALYZERKSTANDARDFAST_H
#define ALICALORAWANALYZERKSTANDARDFAST_H

/* Copyright(c) 1998-2010, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice     */

//_________________________________________________________________________
/// \class AliCaloRawAnalyzerKStandardFast
/// \ingroup EMCALraw
/// \brief  Raw data fitting: standard TMinuit fit
///
/// Extraction of amplitude and peak position
/// from CALO raw data using
/// least square fit for the
/// Moment assuming identical and 
/// independent errors (equivalent with chi square)
///
/// Extracted from AliEMCALRawUtils
///
/// \author Per Thomas Hille <p.t.hille@fys.uio.no>, Yale. 
/// \author David Silvermyr <David.Silvermyr@cern.ch>, ORNL
/// \author Rudiger Haake <ruediger.haake@cern.ch>, Yale

//_________________________________________________________________________

#include "AliCaloRawAnalyzerFitter.h"

class TH1;

class  AliCaloRawAnalyzerKStandardFast : public AliCaloRawAnalyzerFitter
{
  friend class AliCaloRawAnalyzerFactory; // Factory for creation of raw analyzer (rule checker request)
  
 public:
  
  virtual ~AliCaloRawAnalyzerKStandardFast();
  
  virtual AliCaloFitResults  Evaluate( const std::vector<AliCaloBunchInfo> &bunchvector,
                                       UInt_t altrocfg1,
                                       UInt_t altrocfg2 );
  
  void FitRaw( Int_t firstTimeBin, Int_t lastTimeBin,
               Float_t & amp,  Float_t & time,
               Float_t & chi2, Bool_t & fitDone) const ;
 
 private:
 
  TGraph*                     fSignal; ///< histogram for signal that will be fit
  static Double_t             RawResponseFunction(Double_t *x, Double_t *par); 

  AliCaloRawAnalyzerKStandardFast();
  AliCaloRawAnalyzerKStandardFast(               const AliCaloRawAnalyzerKStandardFast & );
  AliCaloRawAnalyzerKStandardFast  & operator = (const AliCaloRawAnalyzerKStandardFast & );
  
  /// \cond CLASSIMP
  ClassDef(AliCaloRawAnalyzerKStandardFast, 1);
  /// \endcond

};

#endif //ALICALORAWANALYZERKSTANDARDFAST_H
