#ifndef TRDgeometryFull_h
#define TRDgeometryFull_h
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

#include "AliTRDgeometry.h"

class AliTRDgeometryFull : public AliTRDgeometry {

 public:

  AliTRDgeometryFull();
  ~AliTRDgeometryFull();

  virtual void    CreateGeometry(Int_t *);
  virtual Int_t   IsVersion() const { return 1; };
  virtual void    Init();

 protected:

  Float_t         fClengthI[kNplan];               // Length of the inner chambers
  Float_t         fClengthM[kNplan];               // Length of the middle chambers
  Float_t         fClengthO[kNplan];               // Length of the outer chambers
  Float_t         fCwidth[kNplan];                 // Width of the chambers

  ClassDef(AliTRDgeometryFull,1)                   // TRD geometry without hole

};

#endif
