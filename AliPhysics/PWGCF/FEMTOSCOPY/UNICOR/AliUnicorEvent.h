#ifndef ALIUNICOREVENT_H
#define ALIUNICOREVENT_H

/* Copyright(c) 1998-2048, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice */
/* $Id$ */

// Author: Dariusz Miskowiec <mailto:d.miskowiec@gsi.de> 2007

//=============================================================================
// parent class of all events; analyzers access data via this class
//=============================================================================

#include <cmath>
#include <TObject.h>

class TTree;

//=============================================================================
class AliUnicorEvent : public TObject {
 public:
  AliUnicorEvent() : TObject(){}                      // constructor
  virtual ~AliUnicorEvent()   {}                      // destructor

  // interface part

  virtual void        AttachTree(TTree *tr) = 0;
  virtual Double_t    Etamin() const = 0;         // experiment's acceptance
  virtual Double_t    Etamax() const = 0;
  virtual Bool_t      Good() const = 0;  
  virtual Double_t    Centrality() const = 0;     // centrality (0,1); 0 is most central
  virtual void        RP(Double_t &qx, Double_t &qy) const = 0;
  virtual Double_t    RPphi() const = 0;
  virtual Double_t    Zver() const = 0;           // z-vertex (-1,1)
  virtual Int_t       NParticles() const = 0;     // number of tracks

  virtual Bool_t      ParticleGood(Int_t i, Int_t pidi=0) const = 0;
  virtual Double_t    ParticleP(Int_t i) const = 0;
  virtual Double_t    ParticleTheta(Int_t i) const = 0;
  virtual Double_t    ParticlePhi(Int_t i) const = 0;
  virtual Double_t    ParticleDedx(Int_t i) const = 0;
  virtual Bool_t      PairGood(double p0, double the0, double phi0, double z0, 
			       double p1, double the1, double phi1, double z1) const = 0;

  // toolkit part

  Int_t    NGoodParticles() const {int n=0; for (int i=0; i<NParticles(); i++) if (ParticleGood(i)) n++; return n;}
  void     RP(Double_t &qx, Double_t &qy, Int_t harmonic) const;
  Double_t ParticlePt(Int_t i) const {return ParticleP(i)*sin(ParticleTheta(i));}
  Double_t ParticlePz(Int_t i) const {return ParticleP(i)*cos(ParticleTheta(i));}
  Double_t ParticleEta(Int_t i) const;
  Double_t ParticleY(Int_t i, Double_t mass) const; 

  ClassDef(AliUnicorEvent,0)
};
#endif 
//=============================================================================
