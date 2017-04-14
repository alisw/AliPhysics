// -*- mode: c++ -*-

#ifndef ALICALORAWANALYZERFITTER_H
#define ALICALORAWANALYZERFITTER_H

/* Copyright(c) 1998-2010, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice     */

//_________________________________________________________________________
/// \class AliCaloRawAnalyzerFitter
/// \ingroup EMCALraw
/// \brief  Raw data fitters base class
///
/// Raw data fitters base class
///
/// \author Per Thomas Hille <p.t.hille@fys.uio.no>, Yale. 
//_________________________________________________________________________

#include "AliCaloRawAnalyzer.h"
#include "AliCaloConstants.h"

using namespace ALTRO;
using namespace EMCAL;

class  TF1;
class  TGraph;

class  AliCaloRawAnalyzerFitter : public AliCaloRawAnalyzer
{

public:
  
  AliCaloRawAnalyzerFitter( const char *name, const char *nameshort );

  virtual ~AliCaloRawAnalyzerFitter();

  TF1 * GetFit() const { return fTf1; };
  
  void PrintFitResult(const TF1 *f) const;

protected:
  
  const double fkEulerSquared;          ///< e^2 = 7.389056098930650227
  TF1        * fTf1;                    ///< Analytical formula of the Semi Gaussian to be fitted
  double       fXaxis[ALTROMAXSAMPLES]; ///< Axis if time bins, ( used by TGraph )

private:
  
  AliCaloRawAnalyzerFitter(               const AliCaloRawAnalyzerFitter & );
  AliCaloRawAnalyzerFitter  & operator = (const AliCaloRawAnalyzerFitter & );
  AliCaloRawAnalyzerFitter();

};

#endif //ALICALORAWANALYZERFITTER_H
