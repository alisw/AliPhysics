#ifndef ALIGENFUNCTION_H
#define ALIGENFUNCTION_H
/* Copyright(c) 1998-2007, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */
// Generator for particles according generic functions
//  
//  TF1 *   fFMomentum;           // momentum distribution function inGeV
//  TF1 *   fFPhi;                // phi distribution function in rad
//  TF1 *   fFTheta;              // theta distribution function in rad
//  TF3 *   fFPosition;           // position distribution function in cm
//  TF1 *   fFPdg;                // pdg distribution function  
//  We assume that the moment, postion and PDG code of particles are independent  
//  Only tracks/particle crossing the reference radius at given z range
//
// Origin: marian.ivanov@cern.ch


#include "AliGenerator.h"
class TF3;

class AliGenFunction : public AliGenerator
{
public:

  AliGenFunction();
  AliGenFunction(const AliGenFunction& func);
  AliGenFunction &operator=(const AliGenFunction& func);
  virtual ~AliGenFunction() {}
  virtual void Generate();
  virtual void Init();
  void SetFunctions(TF1 * momentum, TF1 *fphi, TF1 *ftheta, TF3 * position, TF1* pdg);
  void SetCylinder(Double_t refR, Double_t zmin, Double_t zmax);
  void SetBkG(Float_t b) { fBkG=b; return; }
private:
  Bool_t IntersectCylinder(Float_t r,Float_t zmin, Float_t zmax,Int_t pdg,
			   Float_t o[3],Float_t p[3]) const;  

  Float_t fBkG;                 // field in kGauss
  TF1 *   fFMomentum;           // momentum distribution function
  TF1 *   fFPhi;                // phi distribution function
  TF1 *   fFTheta;              // theta distribution function
  TF3 *   fFPosition;           // position distribution function 
  TF1 *   fFPdg;                // pdg distribution function  
  //
  Double_t fRefRadius;          // reference radius to be crossed
  Double_t fZmin;               // minimal z at reference radius
  Double_t fZmax;               // maximal z at reference radius
  Int_t    fMaxTrial;           // maximal number of attempts
  //
  ClassDef(AliGenFunction,1) // parametrized cosmics generator
};

#endif
