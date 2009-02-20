#ifndef AliVParticle_H
#define AliVParticle_H
/* Copyright(c) 1998-2007, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//-------------------------------------------------------------------------
//     base class for ESD and AOD particles
//     Author: Markus Oldenburg, CERN
//-------------------------------------------------------------------------

#include <Rtypes.h>
#include <TObject.h>

#include <float.h>

const Double_t kAlmost1=1. - Double_t(FLT_EPSILON);
const Double_t kAlmost0=Double_t(FLT_MIN);

const Double_t kB2C=0.299792458e-3;
const Double_t kAlmost0Field=1.e-13;

class AliVParticle: public TObject {

public:
  AliVParticle() { }
  virtual ~AliVParticle() { }
  AliVParticle(const AliVParticle& vPart); 
  AliVParticle& operator=(const AliVParticle& vPart);

  // kinematics
  virtual Double_t Px() const = 0;
  virtual Double_t Py() const = 0;
  virtual Double_t Pz() const = 0;
  virtual Double_t Pt() const = 0;
  virtual Double_t P() const = 0;
  virtual Bool_t   PxPyPz(Double_t p[3]) const = 0;

  virtual Double_t Xv() const = 0;
  virtual Double_t Yv() const = 0;
  virtual Double_t Zv() const = 0;
  virtual Bool_t   XvYvZv(Double_t x[3]) const = 0;  

  virtual Double_t OneOverPt() const = 0;
  virtual Double_t Phi() const = 0;
  virtual Double_t Theta() const = 0;


  virtual Double_t E() const = 0;
  virtual Double_t M() const = 0;
  
  virtual Double_t Eta() const = 0;
  virtual Double_t Y() const = 0;
  
  virtual Short_t Charge() const = 0;
  virtual Int_t   GetLabel() const = 0;
  // PID
  virtual const Double_t *PID() const = 0; // return PID object (to be defined, still)

  // coordinate system conversions
  Bool_t   Local2GlobalMomentum(Double_t p[3], Double_t alpha) const;
  Bool_t   Local2GlobalPosition(Double_t r[3], Double_t alpha) const;
  Bool_t   Global2LocalMomentum(Double_t p[3], Short_t charge, Double_t &alpha) const;
  Bool_t   Global2LocalPosition(Double_t r[3], Double_t alpha) const;

  ClassDef(AliVParticle,0)  // base class for particles
};

#endif
