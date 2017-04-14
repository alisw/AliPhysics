#ifndef ALICALORAWANALYZERCRUDE_H
#define ALICALORAWANALYZERCRUDE_H

/* Copyright(c) 1998-2010, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice     */

//_________________________________________________________________________
/// \class AliCaloRawAnalyzerCrude
/// \ingroup EMCALraw
/// \brief  Raw data fitting: crude fit
///
/// Evaluation of amplitude
/// as max sample value - pedestal
/// Not very accurate, but very robust
///
/// \author Per Thomas Hille <p.t.hille@fys.uio.no>, Yale. 
//_________________________________________________________________________

#include "AliCaloRawAnalyzer.h"

class AliCaloFitResults;
class AliCaloBunchInfo;

class  AliCaloRawAnalyzerCrude : public  AliCaloRawAnalyzer
{
  friend class AliCaloRawAnalyzerFactory;

 public:
   AliCaloRawAnalyzerCrude(); 

  virtual ~AliCaloRawAnalyzerCrude() { ; }
  
  virtual AliCaloFitResults Evaluate( const std::vector<AliCaloBunchInfo> &bunchvector,
				       const UInt_t altrocfg1,  const UInt_t altrocfg2 );
  
 private:
  
  /// \cond CLASSIMP
  ClassDef(AliCaloRawAnalyzerCrude, 1) ;
  /// \endcond

};

#endif //ALICALORAWANALYZERCRUDE_H
