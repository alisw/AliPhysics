#ifndef ALIGENFIXED_H
#define ALIGENFIXED_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

#include "AliGenerator.h"

class AliGenFixed : public AliGenerator
{
 public:
  AliGenFixed();
  AliGenFixed(Int_t npart);
  virtual ~AliGenFixed() {}
  virtual void Generate();
  virtual void Init() {}
  virtual void SetSigma(Float_t sx, Float_t sy, Float_t sz);
  virtual void SetMomentum(Float_t pmom) {fPMin=pmom; fPMax=pmom;}
  virtual void SetPhi(Float_t phi) {fPhiMin=phi*TMath::Pi()/180; fPhiMax=phi*TMath::Pi()/180;}
  virtual void SetTheta(Float_t theta) {fThetaMin=theta*TMath::Pi()/180; fThetaMax=theta*TMath::Pi()/180;}
  virtual void SetPart(Int_t part) {fIpart=part;}
 
protected:

  Int_t fIpart; // Particle type

  ClassDef(AliGenFixed,1) // Single particle generator
};
#endif
