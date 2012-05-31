#ifndef ALICALORAWANALYZERKSTANDARD_H
#define ALICALORAWANALYZERKSTANDARD_H
/**************************************************************************
 * This file is property of and copyright by                              *
 * the Relativistic Heavy Ion Group (RHIG), Yale University, US, 2009     *
 *                                                                        *
 * Primary Author: Per Thomas Hille <p.t.hille@fys.uio.no>                *
 *                                                                        *
 * Contributors are mentioned in the code where appropriate.              *
 * Please report bugs to p.t.hille@fys.uio.no                             *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

// Extraction of amplitude and peak position
// FRom CALO raw data using
// Chi square fit

#include "AliCaloRawAnalyzerFitter.h"
class  TGraph;

class  AliCaloRawAnalyzerKStandard : public AliCaloRawAnalyzerFitter
{
  friend class AliCaloRawAnalyzerFactory; // Factory for creation of raw analyzer (= shutting up the rule checker )
  
 public:
  virtual ~AliCaloRawAnalyzerKStandard();
  virtual AliCaloFitResults  Evaluate( const std::vector<AliCaloBunchInfo> &bunchvector, const UInt_t altrocfg1,  const UInt_t altrocfg2 );
  void FitRaw(const Int_t firstTimeBin, const Int_t lastTimeBin, Float_t & amp, Float_t & time, 
	      Float_t & chi2, Bool_t & fitDone) const ;
 
 private:
  AliCaloRawAnalyzerKStandard();
  AliCaloRawAnalyzerKStandard(const AliCaloRawAnalyzerKStandard & );
  AliCaloRawAnalyzerKStandard  & operator = (const AliCaloRawAnalyzerKStandard  &);
  ClassDef(AliCaloRawAnalyzerKStandard, 2)
};

#endif
