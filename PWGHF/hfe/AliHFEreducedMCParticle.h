/**************************************************************************
* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
*                                                                        *
* Author: The ALICE Off-line Project.                                    *
* Contributors are mentioned in the code where appropriate.              *
*                                                                        *
* Permission to use, copy, modify and distribute this software and its   *
* documentation strictly for non-commercial purposes is hereby granted   *
* without fee, provided that the above copyright notice appears in all   *
* copies and that both the copyright notice and this permission notice   *
* appear in the supporting documentation. The authors make no claims     *
* about the suitability of this software for any purpose. It is          *
* provided "as is" without express or implied warranty.                  *
**************************************************************************/
//
// Debug tree to look at the distribution of the variable we are cutting on
//
//
#ifndef ALIHFEREDUCEDMCPARTICLE_H
#define ALIHFEREDUCEDMCPARTICLE_H

#include <TObject.h>
#include <TMath.h>

class AliHFEreducedMCParticle : public TObject{
 public:
  AliHFEreducedMCParticle();
  AliHFEreducedMCParticle(const AliHFEreducedMCParticle &ref);
  AliHFEreducedMCParticle &operator=(const AliHFEreducedMCParticle &ref);
  ~AliHFEreducedMCParticle() {}
  
  Double_t Pt() const { return TMath::Abs(fSignedPt); }
  Double_t P() const { return fP; }
  Double_t Eta() const { return fEta; }
  Double_t Phi() const { return fPhi; }
  Int_t Charge() const { 
    if(fSignedPt > 0) return 1; 
    else return -1;
  }
  Int_t Pdg() const { return fPdg; }
  Int_t MotherPdg() const { return fMotherPdg; }
  Int_t Source() const { return static_cast<Int_t>(fSource); }
  Bool_t IsSignal() const { return fSignal; }
  Double_t RadialProductionVertex() const { return TMath::Abs(fProductionVertex[0]*fProductionVertex[0]+fProductionVertex[1]*fProductionVertex[1]); }
  Double_t VX() const { return fProductionVertex[0]; }
  Double_t VY() const { return fProductionVertex[1]; }
  Double_t VZ() const { return fProductionVertex[2]; }
  
  void SetSignedPt(Double_t pt, Bool_t positiveCharge){
    double chargesign = positiveCharge ? 1. : -1.;
    fSignedPt = pt * chargesign;
  }
  void SetP(Double_t p) { fP = p; }
  void SetEta(Double_t eta) { fEta = eta; }
  void SetPhi(Double_t phi) { fPhi = phi; }
  void SetPdg(Int_t pdg) { fPdg = pdg; }
  void SetMotherPdg(Int_t pdg) { fPdg = pdg; }
  void SetSource(Int_t source) { fSource = static_cast<Char_t>(source); }
  void SetSignal() { fSignal = kTRUE; }
  void SetProductionVertex(Double_t vx, Double_t vy, Double_t vz) {
    fProductionVertex[0] = vx;
    fProductionVertex[1] = vy;
    fProductionVertex[2] = vz;
  }
  
 private:
  Double_t  fSignedPt;              // signed pt
  Double_t  fP;                     // p
  Double_t  fEta;                   // eta
  Double_t  fPhi;                   // phi
  Int_t     fPdg;                   // pdg
  Int_t     fMotherPdg;             // mother pdg
  Char_t    fSource;                // source
  Bool_t    fSignal;                // signal
  Double_t  fProductionVertex[3];   // production vertex
  
  ClassDef(AliHFEreducedMCParticle, 1)
  
};
#endif
