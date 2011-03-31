// -*- mode: c++ -*-

#ifndef ALICALORAWANALYZERFITTER_H
#define ALICALORAWANALYZERFITTER_H

/**************************************************************************
 * This file is property of and copyright by the Experimental Nuclear     *
 * Physics Group, Yale University, US 2011                                *
 *                                                                        *
 * Author: Per Thomas Hille <perthomas.hille@yale.edu> for the ALICE      *
 * experiment. Contributors are mentioned in the code where appropriate.  *
 * Please report bugs to  perthomas.hille@yale.edu                        *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

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
  
  //Bool_t GetFixTau() const { return fFixTau; }; 
  Bool_t GetFixTau() const; 
  void SetFixTau(Bool_t b) { fFixTau = b; };
  TF1 * GetFit() const { return fTf1; };
  
  void PrintFitResult(const TF1 *f) const;

protected: 
  const double fkEulerSquared; //e^2 = 7.389056098930650227
  TF1 *fTf1;  // Analytical formula of the Semi Gaussian to be fitted 
  double fXaxis[ALTROMAXSAMPLES]; //Axis if time bins, ( used by TGraph )
  Bool_t fFixTau; // flag if tau should be fix

private:
  AliCaloRawAnalyzerFitter(const AliCaloRawAnalyzerFitter & );
  AliCaloRawAnalyzerFitter  & operator = (const AliCaloRawAnalyzerFitter  &);
  AliCaloRawAnalyzerFitter();
};

#endif
