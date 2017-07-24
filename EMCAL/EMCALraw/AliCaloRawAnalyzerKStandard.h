#ifndef ALICALORAWANALYZERKSTANDARD_H
#define ALICALORAWANALYZERKSTANDARD_H

/* Copyright(c) 1998-2010, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice     */

//_________________________________________________________________________
/// \class AliCaloRawAnalyzerKStandard
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
//_________________________________________________________________________

#include "AliCaloRawAnalyzerFitter.h"
class  TGraph;

class  AliCaloRawAnalyzerKStandard : public AliCaloRawAnalyzerFitter
{
  friend class AliCaloRawAnalyzerFactory; // Factory for creation of raw analyzer (rule checker request)
  
 public:
  
  virtual ~AliCaloRawAnalyzerKStandard();
  
  virtual AliCaloFitResults  Evaluate( const std::vector<AliCaloBunchInfo> &bunchvector,
                                       UInt_t altrocfg1,
                                       UInt_t altrocfg2 );
  
  void FitRaw( Int_t firstTimeBin, Int_t lastTimeBin,
               Float_t & amp,  Float_t & time,
               Float_t & chi2, Bool_t & fitDone) const ;
 
 private:
  
  AliCaloRawAnalyzerKStandard();
  AliCaloRawAnalyzerKStandard(               const AliCaloRawAnalyzerKStandard & );
  AliCaloRawAnalyzerKStandard  & operator = (const AliCaloRawAnalyzerKStandard & );
  
  /// \cond CLASSIMP
  ClassDef(AliCaloRawAnalyzerKStandard, 2);
  /// \endcond

};

#endif //ALICALORAWANALYZERKSTANDARD_H
