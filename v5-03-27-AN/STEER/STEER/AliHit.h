#ifndef ALIHIT_H
#define ALIHIT_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

// Base class for hits
// This class is used as a base class for 
// hits in the different detectors

#include "TObject.h"

class AliHit : public TObject {
public:
  AliHit();
  AliHit(Int_t shunt, Int_t track);
  virtual ~AliHit() {}
  virtual Int_t GetTrack() const {return fTrack;}
  virtual void SetTrack(Int_t track) {fTrack=track;}
  virtual Float_t X() const {return fX;}
  virtual Float_t Y() const {return fY;}
  virtual Float_t Z() const {return fZ;}
  virtual Int_t Track() const {return fTrack;}
  
protected:
  Int_t     fTrack;  // Track number
  Float_t   fX;      // X position of the hit
  Float_t   fY;      // Y position of the hit
  Float_t   fZ;      // Z position of the hit

  ClassDef(AliHit,1)  //Base class for all Alice hits
};
#endif
