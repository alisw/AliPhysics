#ifndef ALIGENPROMPTPHOTONS_H
#define ALIGENPROMPTPHOTONS_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//-------------------------------------------------------------------------
// author: Sergey Kiselev, ITEP, Moscow
// e-mail: Sergey.Kiselev@cern.ch
// tel.: 007 495 129 95 45
//-------------------------------------------------------------------------
// Generator of prompt photons for the reaction A+B, sqrt(S)
//
// main assumptions:
// 1. flat rapidity distribution
// 2. all existing p+p(pbar) data at y_{c.m.} can be described by the function
//           F(x_T) = (sqrt(s))^5 Ed^3sigma/d^3p, x_T = 2p_t/sqrt(s)
//           all data points cover the region x_T: 0.01 - 0.6
//    see Nucl.Phys.A783:577-582,2007, hep-ex/0609037
// 3. binary scaling: for A+B at the impact parameter b
//    Ed^3N^{AB}(b)/d^3p = Ed^3sigma^{pp}/d^3p A B T_{AB}(b),
//    T_{AB}(b) - nuclear overlapping fuction, calculated in the Glauber approach,
//                nuclear density is parametrized by a Woods-Saxon with nuclear radius
//                R_A = 1.19 A^{1/3} - 1.61 A^{-1/3} fm and surface thickness a=0.54 fm
// 4. nuclear effects (Cronin, shadowing, ...) are ignored
//
// input parameters:
//       fAProjectile, fATarget - number of nucleons in a nucleus A and B
//       fMinImpactParam - minimal impct parameter, fm
//       fMaxImpactParam - maximal impct parameter, fm
//       fEnergyCMS - sqrt(S) per nucleon pair, AGeV
//
//       fYMin - minimal rapidity of photons 
//       fYMax - maximal rapidity of photons
//       fPtMin - minimal p_t value of gamma, GeV/c
//       fPtMax - maximal p_t value of gamma, GeV/c
//-------------------------------------------------------------------------
// comparison with SPS and RHIC data, prediction for LHC can be found in
// arXiv:0811.2634 [nucl-th]
//-------------------------------------------------------------------------

class TF1;

#include "AliGenerator.h"

class AliGenPromptPhotons : public AliGenerator
{
 public:

  AliGenPromptPhotons();
  AliGenPromptPhotons(Int_t npart);
  virtual ~AliGenPromptPhotons();
  virtual void Generate();
  virtual void Init();
  virtual void SetPtRange(Float_t ptmin = 0.1, Float_t ptmax=10.);
  virtual void SetYRange(Float_t ymin = -1., Float_t ymax=1.);

// Setters
    virtual void SetAProjectile(Float_t a = 208) {fAProjectile = a;}
    virtual void SetATarget(Float_t a = 208)     {fATarget     = a;}
    virtual void SetEnergyCMS(Float_t energy = 5500.) {fEnergyCMS = energy;}
    virtual void SetImpactParameterRange(Float_t bmin = 0., Float_t bmax = 0.)
	{fMinImpactParam=bmin; fMaxImpactParam=bmax;}

 protected:
  Float_t fAProjectile;     // Projectile nucleus mass number
  Float_t fATarget;         // Target nucleus mass number
  Float_t fEnergyCMS;       // Center of mass energy
  Float_t fMinImpactParam;  // minimum impact parameter
  Float_t fMaxImpactParam;  // maximum impact parameter	
  
  static Double_t FitData      (const Double_t *xx, const Double_t *par);
  static Double_t WSforNorm    (const Double_t *xx, const Double_t *par);
  static Double_t WSz          (const Double_t *xx, const Double_t *par);
  static Double_t TA           (const Double_t *xx, const Double_t *par);
  static Double_t TB           (const Double_t *xx, const Double_t *par);
  static Double_t TAxTB        (const Double_t *xx, const Double_t *par);
  static Double_t TAB          (const Double_t *xx, const Double_t *par);

  static TF1 *fgDataPt;             // d^{2}#sigma^{pp}/(dp_t dy) from data fit 
  static TF1 *fgWSzA;               // Wood Saxon parameterisation for nucleus A 
  static TF1 *fgWSzB;               // Wood Saxon parameterisation for nucleus B 
  static TF1 *fgTA;                 // nuclear thickness function T_A(b) (1/fm**2) 
  static TF1 *fgTB;                 // nuclear thickness function T_B(phi)=T_B(sqtr(s**2+b**2-2*s*b*cos(phi))) 
  static TF1 *fgTAxTB;              // s * TA(s) * 2 * Integral(0,phiMax) TB(phi(s,b)) 
  static TF1 *fgTAB;                // overlap function T_AB(b) (1/fm**2) 
  
 private:

  AliGenPromptPhotons(const AliGenPromptPhotons & PromptPhotons);
  AliGenPromptPhotons& operator = (const AliGenPromptPhotons & PromptPhotons) ;


  ClassDef(AliGenPromptPhotons, 1) // prompt photon generator
};
#endif
