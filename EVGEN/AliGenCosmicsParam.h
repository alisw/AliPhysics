#ifndef ALIGENCOSMICSPARAM_H
#define ALIGENCOSMICSPARAM_H
/* Copyright(c) 1998-2007, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// Generator for muons according to kinematic parametrizations at ALICE
// (not at the surface).
// Origin: andrea.dainese@lnl.infn.it


#include "AliLog.h"
#include "AliGenerator.h"

class AliGenCosmicsParam : public AliGenerator
{
public:

  AliGenCosmicsParam();
  virtual ~AliGenCosmicsParam() {}
  virtual void Generate();
  virtual void Init();
  void SetParamMI() { fParamMI=kTRUE; fParamACORDE=kFALSE; return; }
  void SetParamACORDE() { fParamMI=kFALSE; fParamACORDE=kTRUE; return; }
  void SetYOrigin(Float_t y=600.) { fYOrigin=y; return; }
  void SetMaxAngleWRTVertical(Float_t max=45.) { 
      if(max<0. || max>90.) AliFatal("angle must be in [0,pi/2]");
      fMaxAngleWRTVertical=max; return; }
  void SetBkG(Float_t b) { fBkG=b; return; }
  void SetInTPC() { fTPC=kTRUE; return; }
  void SetInITS() { fITS=kTRUE; return; }
  void SetInSPDouter() { fSPDouter=kTRUE; return; }
  void SetInSPDinner() { fSPDinner=kTRUE; return; }

private:

  Bool_t IntersectCylinder(Float_t r,Float_t z,Int_t pdg,
			   Float_t o[3],Float_t p[3]) const;  

  Bool_t fParamMI;              // parametrization from M.Ivanov
  Bool_t fParamACORDE;          // parametrization from AliGenACORDE 
                                // (parametrized at ALICE y)
  Float_t fYOrigin;             // y of muon origin
  Float_t fMaxAngleWRTVertical; // maximum angle between momentum and y axis
  Float_t fBkG;                 // field in kGauss
  Bool_t fTPC;                  // acceptance cuts
  Bool_t fITS;                  // acceptance cuts
  Bool_t fSPDouter;             // acceptance cuts
  Bool_t fSPDinner;             // acceptance cuts

  ClassDef(AliGenCosmicsParam,1) // parametrized cosmics generator
};

#endif
