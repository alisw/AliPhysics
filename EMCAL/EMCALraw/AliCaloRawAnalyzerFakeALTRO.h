#ifndef ALICALORAWANALYZERFAKEALTRO_H
#define ALICALORAWANALYZERFAKEALTRO_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//_________________________________________________________________________
/// \class AliCaloRawAnalyzerFakeALTRO
/// \ingroup EMCALraw
/// \brief  Raw data fitting: FALTRO
///
/// Extraction of Amplitude and peak
/// position of FastALTRO 
///
/// \author Rachid Guernane, < guernane@lpsc.in2p3.fr>, LPSC-IN2P3-CNRS. 
//_________________________________________________________________________

#include "AliCaloRawAnalyzerFitter.h"

class  AliCaloRawAnalyzerFakeALTRO : public AliCaloRawAnalyzerFitter
{
  friend class AliCaloRawAnalyzerFactory;

 public:
  
  virtual ~AliCaloRawAnalyzerFakeALTRO();
  
  virtual AliCaloFitResults  Evaluate( const std::vector<AliCaloBunchInfo> &bunchvector,
                                       UInt_t altrocfg1,
                                       UInt_t altrocfg2 );
  
 private:
  
  AliCaloRawAnalyzerFakeALTRO();
  AliCaloRawAnalyzerFakeALTRO(               const AliCaloRawAnalyzerFakeALTRO & );
  AliCaloRawAnalyzerFakeALTRO  & operator = (const AliCaloRawAnalyzerFakeALTRO & );
  
  /// \cond CLASSIMP
  ClassDef(AliCaloRawAnalyzerFakeALTRO,1) ;
  /// \endcond

};

#endif //ALICALORAWANALYZERFAKEALTRO_H
