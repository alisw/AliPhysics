#ifndef AliVirtualParticle_H
#define AliVirtualParticle_H
/* Copyright(c) 1998-2007, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//-------------------------------------------------------------------------
//     base class for ESD and AOD particles
//     Author: Markus Oldenburg, CERN
//-------------------------------------------------------------------------

#include <Rtypes.h>
#include <TObject.h>

class AliVirtualParticle: public TObject {

public:
  AliVirtualParticle() { }
  virtual ~AliVirtualParticle() { }

  // kinematics
  virtual Double_t Px() const = 0;
  virtual Double_t Py() const = 0;
  virtual Double_t Pz() const = 0;
  virtual Double_t Pt() const = 0;
  virtual Double_t P() const = 0;

  virtual Double_t OneOverPt() const = 0;
  virtual Double_t Phi() const = 0;
  virtual Double_t Theta() const = 0;


  virtual Double_t E() const = 0;
  virtual Double_t M() const = 0;
  
  virtual Double_t Eta() const = 0;
  virtual Double_t Y() const = 0;
  
  virtual Short_t Charge() const = 0;

  // PID
  virtual const Double_t *PID() const = 0; // return PID object (to be defined, still)


  ClassDef(AliVirtualParticle,0)  // base class for particles
};

#endif
