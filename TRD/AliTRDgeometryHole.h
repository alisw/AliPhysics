#ifndef ALITRDGEOMETRYHOLE_H
#define ALITRDGEOMETRYHOLE_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

#include "AliTRDgeometry.h"

class AliTRDgeometryHole : public AliTRDgeometry {

 public:

  AliTRDgeometryHole();
  virtual ~AliTRDgeometryHole();

          void    CreateGeometry(Int_t *idtmed); 
          Int_t   IsVersion() const   { return 0; };
          void    Init();

          void    SetPHOShole()       { };
          void    SetRICHhole()       { };

          void    SetNRowPad();
  virtual void    SetNRowPad(Int_t p, Int_t c, Int_t npad);

          Bool_t  GetPHOShole() const { return kTRUE; };
          Bool_t  GetRICHhole() const { return kTRUE; };

 protected:

  Float_t         fClengthI[kNplan];               // Length of the inner chambers
  Float_t         fClengthM1[kNplan];              // Length of the middle chambers
  Float_t         fClengthM2[kNplan];              // Length of the middle chambers
  Float_t         fClengthO1[kNplan];              // Length of the outer chambers
  Float_t         fClengthO2[kNplan];              // Length of the outer chambers
  Float_t         fClengthO3[kNplan];              // Length of the outer chambers

  ClassDef(AliTRDgeometryHole,1)                   // TRD geometry with hole

};

#endif
