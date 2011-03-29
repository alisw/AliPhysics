#ifndef ALICALORAWANALYZERFAKEALTRO_H
#define ALICALORAWANALYZERFAKEALTRO_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */
/*

 
Author: R. GUERNANE LPSC Grenoble CNRS/IN2P3
*/

#include "AliCaloRawAnalyzerFitter.h"

class  AliCaloRawAnalyzerFakeALTRO : public AliCaloRawAnalyzerFitter
{
  friend class AliCaloRawAnalyzerFactory;

 public:
  virtual ~AliCaloRawAnalyzerFakeALTRO();
  virtual AliCaloFitResults  Evaluate( const std::vector<AliCaloBunchInfo> &bunchvector, const UInt_t altrocfg1,  const UInt_t altrocfg2 );
  
 private:
  AliCaloRawAnalyzerFakeALTRO();
  AliCaloRawAnalyzerFakeALTRO(const AliCaloRawAnalyzerFakeALTRO & );
  AliCaloRawAnalyzerFakeALTRO  & operator = (const AliCaloRawAnalyzerFakeALTRO  &);
  ClassDef(AliCaloRawAnalyzerFakeALTRO,1)

};

#endif
