#ifndef ALICALORAWANALYZERFASTFIT_H
#define ALICALORAWANALYZERFASTFIT_H

/* Copyright(c) 1998-2010, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice     */

#include "AliCaloRawAnalyzerFitter.h"
#include "Rtypes.h"

//_________________________________________________________________________
/// \class AliCaloRawAnalyzerFastFit
/// \ingroup EMCALraw
/// \brief  Raw data fitting: special fast fit
///
/// Extraction of Amplitude and peak
/// position using special algorithm
///
/// \author Alexei Palinov 
//_________________________________________________________________________

class  AliCaloRawAnalyzerFastFit : public AliCaloRawAnalyzerFitter
{
  friend class  AliCaloRawAnalyzerFactory; // RuleChecker request
  
public:
  
  virtual ~AliCaloRawAnalyzerFastFit() { ; }
  
  virtual AliCaloFitResults Evaluate( const std::vector<AliCaloBunchInfo> &bunchvector,
                                     UInt_t altrocfg1, UInt_t altrocfg2 );
  
private:
  
  AliCaloRawAnalyzerFastFit();
  
  /// \cond CLASSIMP
  ClassDef( AliCaloRawAnalyzerFastFit, 1 ) ;
  /// \endcond
  
};

#endif //ALICALORAWANALYZERFASTFIT_H
