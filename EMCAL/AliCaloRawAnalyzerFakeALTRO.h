#ifndef ALICALORAWANALYZERFAKEALTRO_H
#define ALICALORAWANALYZERFAKEALTRO_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/*

 
Author: R. GUERNANE LPSC Grenoble CNRS/IN2P3
*/

#include "AliCaloRawAnalyzer.h"


class  TF1;
class  TGraph;

class  AliCaloRawAnalyzerFakeALTRO : public AliCaloRawAnalyzer
{
  friend class AliCaloRawAnalyzerFactory;

 public:
  //AliCaloRawAnalyzerFakeALTRO();
  virtual ~AliCaloRawAnalyzerFakeALTRO();
  
  virtual AliCaloFitResults  Evaluate( const std::vector<AliCaloBunchInfo> &bunchvector, const UInt_t altrocfg1,  const UInt_t altrocfg2 );
  void PrintFitResult(const TF1 *f) const;
  
  // shaper tau value, in time-bins, and flag for keeping tau fixed
  Float_t GetTau() const { return fTau;};
  void SetTau(Float_t f) { fTau = f; }; 
  Bool_t GetFixTau() const { return fFixTau; }; 
  void SetFixTau(Bool_t b) { fFixTau = b; }; 

  // extra interfaces
  TF1 * GetFit() const { return fTf1; };

 private:
  AliCaloRawAnalyzerFakeALTRO();
  AliCaloRawAnalyzerFakeALTRO(const AliCaloRawAnalyzerFakeALTRO & );
  AliCaloRawAnalyzerFakeALTRO  & operator = (const AliCaloRawAnalyzerFakeALTRO  &);
 
  double fXaxis[ALTROMAXSAMPLES]; //Axis if time bins, ( used by TGraph )
  const double fkEulerSquared; //e^2 = 7.389056098930650227
  TF1 *fTf1;     // Analytical formula of the Semi Gaussian to be fitted

  Float_t fTau; // shaper tau, in time bins
  Bool_t fFixTau; // flag if tau should be fix

  ClassDef(AliCaloRawAnalyzerFakeALTRO,1)

};

#endif
