#ifndef CORRELPARTICLE_H
#define CORRELPARTICLE_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice     */
/* $Id:  $ */

//________________________________________________
// Main container class - stores generic particle.
//-- Author: Paul Constantin

#include "CorrelDefs.h"

class CorrelParticle_t {        
 public:
  CorrelParticle_t();
  CorrelParticle_t(Float_t pt, Float_t p, Float_t t, Float_t m, PartType_t i);
  CorrelParticle_t(const CorrelParticle_t& p);
  virtual ~CorrelParticle_t() {;}
  virtual CorrelParticle_t* Copy();
  
  // data setters
  void SetPt(Float_t v)    {fPt=v;}
  void SetPhi(Float_t v)   {fPhi=v;}
  void SetEta(Float_t v)   {fEta=v;}
  void SetMass(Float_t v)  {fMass=v;}
  void SetID(PartType_t v) {fID=v;}
  
  // data getters
  virtual Float_t Pt()    const {return TMath::Abs(fPt);}
  virtual Float_t Phi()   const {return fPhi;}
  virtual Float_t Eta()   const {return fEta;}    
  virtual Float_t M()     const {return fMass;}
  virtual PartType_t ID() const {return fID;}
  
  // derived data getters
  virtual Short_t Q()     const {return (Pt()>0)?(Short_t(fPt/Pt())):(-99);}
  virtual Float_t Theta() const {return 2.*TMath::ATan(TMath::Exp(-fEta));}
  virtual Float_t SinT()  const {return TMath::Sin(Theta());}
  virtual Float_t TanT()  const {return TMath::Tan(Theta());}
  virtual Float_t Px()    const {return Pt()*TMath::Cos(fPhi);}
  virtual Float_t Py()    const {return Pt()*TMath::Sin(fPhi);}
  virtual Float_t Pz()    const {return (TMath::Abs(TanT())>0)?(Pt()/TanT()):(-99.);}
  virtual Float_t P()     const {return (TMath::Abs(SinT())>0)?(Pt()/SinT()):(-99.);}
  virtual Float_t E()     const {return TMath::Sqrt(P()*P()+fMass*fMass);}
  virtual Float_t Et()    const {return TMath::Sqrt(Pt()*Pt()+fMass*fMass);}
  virtual Float_t Y()     const {return 0.5*TMath::Log((E()+Pz())/(E()-Pz()));}
  
  virtual void Show() const;
  
 protected:
  Float_t fPt;    // pt; store charge as its sign
  Float_t fPhi;   // phi
  Float_t fEta;   // eta
  Float_t fMass;  // mass
  PartType_t fID; // ID
};

#endif
