#ifndef AliAODJet_H
#define AliAODJet_H
/* Copyright(c) 1998-2007, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//-------------------------------------------------------------------------
//     AOD track base class
//     Author: Markus Oldenburg, CERN
//-------------------------------------------------------------------------

#include <TRef.h>

#include "AliVirtualParticle.h"
#include "AliAODVertex.h"

class AliAODJet : public AliVirtualParticle {

 public:
  
  AliAODJet();

  virtual ~AliAODJet();
  AliAODJet(const AliAODJet& trk); 
  AliAODJet& operator=(const AliAODJet& trk);

  virtual Double_t Px() const { return 0.;}
  virtual Double_t Py() const { return 0.;}
  virtual Double_t Pz() const { return 0.;}
  virtual Double_t Pt() const { return 0.;}
  virtual Double_t P() const { return 0.;}
  virtual Double_t OneOverPt() const { return 0.;}
  virtual Double_t Phi() const { return 0.;}
  virtual Double_t Theta() const { return 0.;}
  virtual Double_t E() const { return 0.;}
  virtual Double_t M() const { return 0.;}
  virtual Double_t Eta() const { return 0.;}
  virtual Double_t Y() const { return 0.;}
  virtual Short_t Charge() const { return 0;}
  virtual const Double_t* PID() const { return NULL;}

 private :

  // Momentum & position

  ClassDef(AliAODJet,1);
};

#endif
