#ifndef ALIGENTUNEDONPBPB_H
#define ALIGENTUNEDONPBPB_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id: AliGenTunedOnPbPb.h 51126 2013-08-19 13:37:49Z fnoferin $ */

// Parameterisation based on 5.5 ATeV PbPb data
// pi,K, p , K0, lambda, phi, Xi, Omega spectra, v2, v3 (no jets!)
// used for the ALICE TDRs.
// Author: fnoferin@cern.ch

class TF1;

#include "AliGenerator.h"
#include "TMath.h"
#include "TH1F.h"
#include "TH2F.h"

class AliGenTunedOnPbPb : public AliGenerator
{
 public:

  AliGenTunedOnPbPb();
  virtual ~AliGenTunedOnPbPb();
  virtual void Generate();
  virtual void Init();
  virtual void SetPtRange(Float_t ptmin = 0., Float_t ptmax=15.);
  virtual void SetCentralityRange(Float_t cmin,Float_t cmax){fCmin=TMath::Max(cmin,Float_t(0));fCmax=TMath::Min(cmax,Float_t(100));};  

  void SetSpectrum(Int_t species,TH1 *histo){fgHSpectrum[species] = histo;};
  TH1 *GetSpectrum(Int_t species){return fgHSpectrum[species];};
  void SetV2(Int_t species,TH1 *histo){fgHv2[species] = histo;};
  TH1 *GetV2(Int_t species){return fgHv2[species];};

  static Float_t GetEventPlane(){return fgEventplane;};

  static TH1F *GetMultVsCentrality(Int_t species);

  void SetCentralityDependence(Bool_t flag=kTRUE){fChangeWithCentrality=flag;};
  void SetYmax(Float_t value){fYMaxValue=value;};
  void SetYmaxFlatness(Float_t value=2.0){fYlimitForFlatness=value;};
  void SetDecreaseSp(Float_t value=0.2){fYdecreaseSp=value;};
  void SetDecreaseV2(Float_t value=0.2){fYdecreaseV2=value;};

  enum Particles {
    kPiPlus, kPiMinus, kPi0   , kKaonPlus, kKaonMinus, kProton , kAntiProton , kKaon0, kLambda, kAntiLambda,
    kPhi   , kXi     , kAntiXi, kOmega   , kAntiOmega, kNeutron, kAntiNeutron
  };
 private:
  AliGenTunedOnPbPb(const AliGenTunedOnPbPb &para);
  AliGenTunedOnPbPb& operator = (const AliGenTunedOnPbPb &para) ;

  static void DefineSpectra();

  static void SetParameters(Float_t centrality);

  static const Int_t fgNspecies = 17; // number of species available
  static Int_t fgPdgInput[fgNspecies]; // pdgs available
  static Float_t fgMult[fgNspecies]; // current multiplicity  (fixed as a function of centrality)
  static Float_t fgV3Overv2; // v3 / v2 (fixed as a function of centrality)
  static Float_t fgEventplane; // event plane (Psi2)

  static TF1 *fgV2; // function to model the anisotropy

  TH1 *fgHSpectrum[fgNspecies]; // pt distributions (should be passed in Config.C)
  TH1 *fgHv2[fgNspecies]; // v2 vs pt (should be passed in Config.C)

  Float_t fCmin; // min centrality
  Float_t fCmax; // max centrality

  Bool_t fChangeWithCentrality;  // flag to apply a centrality dependence to pt-distr and v2
  Float_t fYMaxValue;            // max value for rapidity (abs)
  Float_t fYlimitForFlatness;    // starting from this value y-profile starts to decrease both for spectra and v2
  Float_t fYdecreaseSp;          // angular coefficient for the decrease above fYlimitForFlatness (spectra)
  Float_t fYdecreaseV2;          // angular coefficient for the decrease above fYlimitForFlatness (v2)

  ClassDef(AliGenTunedOnPbPb,3) // Hijing parametrisation generator
};
#endif
