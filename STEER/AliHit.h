#ifndef AliHit_H
#define AliHit_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

#include "TObject.h"

class AliHit : public TObject {
public:
  Int_t     fTrack;      //track number
  //Position of the hit
  Float_t fX;
  Float_t fY;
  Float_t fZ;

public:
  AliHit();
  AliHit(Int_t shunt, Int_t track);
  virtual ~AliHit() {}
  virtual int GetTrack() const {return fTrack;}
  inline virtual void SetTrack(int track) {fTrack=track;}
  
  ClassDef(AliHit,1)  //Base class for all Alice hits
};
#endif
