#ifndef ALIGENTHERMALPHOTONS_H
#define ALIGENTHERMALPHOTONS_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//-------------------------------------------------------------------------
// author: Sergey Kiselev, ITEP, Moscow
// e-mail: Sergey.Kiselev@cern.ch
// tel.: 007 495 129 95 45
//-------------------------------------------------------------------------
// Generator of direct thermal photons for the reaction A+B, sqrt(S)
// main assumptions:
// 1+1 Bjorken scaling hydrodinamics.
// 1st order phase transition
// QGP + Mixed (QGP+HHG) + HHG (Hot Hadron Gas) phases, 
// an ideal massless parton gas and ideal massless HHG 
// see 
// F.D.Steffen, nucl-th/9909035
// F.D.Steffen and M.H.Thoma, Phys.Lett. B510, 98 (2001)
// T.Peitzmann and M.H.Thoma, Phys.Rep., 364, 175 (2002) 
//
// photon rates for QGP: Phys.Rep., 364, 175 (2002), section 2.1.1
//
// photon rates for HHG
// prates for i rho --> pi gamma, pi pi --> rho gamma and rho --> pi pi gamma:
// Song and Fai, Phys.Rev. C58, 1689 (1998)
// rates for omega --> pi gamma: Phys.Rev. D44, 2774 (1991)
//
// input parameters:
//       fAProjectile, fATarget - number of nucleons in a nucleus A and B
//       fMinImpactParam - minimal impct parameter, fm
//       fMaxImpactParam - maximal impct parameter, fm
//       fEnergyCMS - sqrt(S) per nucleon pair, AGeV
//       fTau0 - initial proper time, fm/c
//       fT0 - initial temperature, GeV
//       fTc - critical temperature, GeV
//       fTf - freeze-out temperature, GeV
//       fGhhg - effective number of degrees of freedom in HHG
//
//       fYMin - minimal rapidity of photons 
//       fYMax - maximal rapidity of photons
//              in [fYMin,fYMax] uniform distribution of gamma is assumed
//       fPtMin - minimal p_t value of gamma, GeV/c
//       fPtMax - maximal p_t value of gamma, GeV/c
//-------------------------------------------------------------------------
// comparison with SPS and RHIC data, prediction for LHC can be found in
// arXiv:0811.2634 [nucl-th]
//-------------------------------------------------------------------------

class TH1F;

#include "AliGenerator.h"

class AliGenThermalPhotons : public AliGenerator
{
 public:

  AliGenThermalPhotons();
  AliGenThermalPhotons(Int_t npart);
  virtual ~AliGenThermalPhotons();
  virtual void Generate();
  virtual void Init();
  virtual void SetPtRange(Float_t ptmin = 0.1, Float_t ptmax=10.);
  virtual void SetYRange(Float_t ymin = -1., Float_t ymax=1.);

// Setters
    virtual void SetAProjectile(Float_t a = 208) {fAProjectile = a;}
    virtual void SetATarget(Float_t a = 208)     {fATarget     = a;}
    virtual void SetEnergyCMS(Float_t energy = 5500.) {fEnergyCMS = energy;}
    virtual void    SetImpactParameterRange(Float_t bmin = 0., Float_t bmax = 0.)
	{fMinImpactParam=bmin; fMaxImpactParam=bmax;}
    virtual void    SetTau0(Float_t tau0 = 0.1)             {fTau0   = tau0;}
    virtual void    SetT0(Float_t   T0   = 0.650)           {fT0     = T0;}
    virtual void    SetTc(Float_t   Tc   = 0.170)           {fTc     = Tc;}
    virtual void    SetTf(Float_t   Tf   = 0.100)           {fTf     = Tf;}
    virtual void    SetGhhg(Int_t   Ghhg = 8)               {fGhhg   = Ghhg;}

 protected:
  Float_t fAProjectile;     // Projectile nucleus mass number
  Float_t fATarget;         // Target nucleus mass number
  Float_t fEnergyCMS;       // Center of mass energy
  Float_t fMinImpactParam;  // minimum impact parameter
  Float_t fMaxImpactParam;  // maximum impact parameter	
  Float_t fTau0;            // initial proper time, fm	
  Float_t fT0;              // initial temperature, GeV	
  Float_t fTc;              // critical temperature, GeV	
  Float_t fTf;              // freeze-out temperature, GeV	
  Int_t   fGhhg;            // number of degrees of freedom in HHG	

  TH1F *fSumPt;             // histo with pt from all origins

 private:

  AliGenThermalPhotons(const AliGenThermalPhotons & ThermalPhotons);
  AliGenThermalPhotons& operator = (const AliGenThermalPhotons & ThermalPhotons) ;


  ClassDef(AliGenThermalPhotons, 1) // Direct thermal photon generator
};
#endif
