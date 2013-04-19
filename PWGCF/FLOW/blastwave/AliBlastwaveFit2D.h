#ifndef ALIBLASTWAVEFIT2D_H
#define ALIBLASTWAVEFIT2D_H


/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id: AliBlastwaveFit2D.h 49869 2012-05-17 04:49:51Z fnoferin $ */

/////////////////////////////////////////////////
//                                             //
//        Blastwave Fit Class in 2D            //
//           noferini@bo.infn.it               //
/////////////////////////////////////////////////

#include "TF1.h"
#include "TF2.h"

#include "AliBlastwaveFit.h"

class AliBlastwaveFit2D : public AliBlastwaveFit
{
 public:
  enum{
    kParT,
    kParS2,
    kParRho0,
    kParRho2,
    kParGamma
  };

  static const char *fgParName[7];
  static Float_t fgStartValues[5];
  static const Float_t fgStepValues[5];
  static const Float_t fgMinValues[5];
  static const Float_t fgMaxValues[5];

  AliBlastwaveFit2D(const char *name,Double_t mass=0);
  AliBlastwaveFit2D();
  ~AliBlastwaveFit2D();

  // required methods (pure virtual in the base class)
  void Initialize();
  Int_t SetParameter(Int_t ipar,Double_t val);
  Int_t SetNormalization();
  Int_t GetNpar() {return 5;};

  Float_t GetParStart(Int_t ipar) {if(ipar>=0 && ipar <5) return fgStartValues[ipar]; else return 0.0;};
  Float_t GetParStep(Int_t ipar) {if(ipar>=0 && ipar <5) return fgStepValues[ipar]; else return 0.0;};
  Float_t GetParMin(Int_t ipar) {if(ipar>=0 && ipar <5) return fgMinValues[ipar]; else return 0.0;};
  Float_t GetParMax(Int_t ipar) {if(ipar>=0 && ipar <5) return fgMaxValues[ipar]; else return 0.0;};

  Float_t GetMeanBeta();
  Float_t GetMeanBeta(Double_t par[]);

  // method suggested to be replace (virtual in the base class)
  void SetMass(Double_t mass);
  const char *GetParName(Int_t i) {if(i >= 0 && i < 7) return fgParName[i]; else return "";};
  void SwitchOffFlow(TMinuit *m) const;

  // to allow double integration on R and phi
  static Double_t FunctionIntYield(Double_t x[],Double_t par[]); // integrand (phi,r)
  static Double_t FunctionIntV2(Double_t x[],Double_t par[]); // integrand (phi,r)
 
  static Double_t V2(Double_t x[],Double_t par[]); // integrated on phi and r
  static Double_t Pt(Double_t x[],Double_t par[]); // integrated on phi and r


 private:
  AliBlastwaveFit2D(const AliBlastwaveFit2D & old);
  AliBlastwaveFit2D& operator=(const AliBlastwaveFit2D &/*source*/); // ass. op.


  static TF2 *fgFuncIntYield; // function used to integrate FunctionIntYield
  static TF2 *fgFuncIntV2; // function used to integrate FunctionIntV2

  ClassDef(AliBlastwaveFit2D,1)  // blast wave fit 2D
};

#endif


