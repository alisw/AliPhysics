#ifndef ALIBLASTWAVEFITSPECTRA2_H
#define ALIBLASTWAVEFITSPECTRA2_H


/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id: AliBlastwaveFitSpectra2.h 49869 2012-05-17 04:49:51Z fnoferin $ */

/////////////////////////////////////////////////
//                                             //
//      Blastwave Fit Class in Spectra2        //
//           noferini@bo.infn.it               //
/////////////////////////////////////////////////

#include "TF1.h"
#include "TF2.h"

#include "AliBlastwaveFit.h"

class AliBlastwaveFitSpectra2 : public AliBlastwaveFit
{
public:
  enum{
      kParT,
      kParBetaAv,
      kParGamma
  };
  
  static const char *fgParName[5];
  static Float_t fgStartValues[3];
  static const Float_t fgStepValues[3];
  static const Float_t fgMinValues[3];
  static const Float_t fgMaxValues[3];
  
  AliBlastwaveFitSpectra2(const char *name,Double_t mass=0);
  AliBlastwaveFitSpectra2();
  ~AliBlastwaveFitSpectra2();

  // required methods (pure virtual in the base class)
  void Initialize();
  Int_t SetParameter(Int_t ipar,Double_t val);
  Int_t SetNormalization();
  Int_t GetNpar() {return 3;};
  
  const Float_t GetParStart(Int_t ipar) {if(ipar>=0 && ipar <3) return fgStartValues[ipar]; else return 0.0;};
  const Float_t GetParStep(Int_t ipar) {if(ipar>=0 && ipar <3) return fgStepValues[ipar]; else return 0.0;};
  const Float_t GetParMin(Int_t ipar) {if(ipar>=0 && ipar <3) return fgMinValues[ipar]; else return 0.0;};
  const Float_t GetParMax(Int_t ipar) {if(ipar>=0 && ipar <3) return fgMaxValues[ipar]; else return 0.0;};
  
  const Float_t GetMeanBeta();
  const Float_t GetMeanBeta(Double_t par[]);

  // method suggested to be replace (virtual in the base class)
  void SetMass(Double_t mass);
  const char *GetParName(Int_t i) {if(i >= 0 && i < 5) return fgParName[i]; else return "";};
  
  // to allow double integration on R and phi
  static Double_t FunctionIntYield(Double_t x[],Double_t par[]); // integrand (phi,r)
  
  static Double_t Pt(Double_t x[],Double_t par[]); // integrated on phi and r
  
  
private:
  static TF1 *fgFuncIntYield; // function used to integrate FunctionIntYield
  
  ClassDef(AliBlastwaveFitSpectra2,1)  // blast wave fit Spectra beta parameterization
};

#endif


