#ifndef TRDgeometryHole_h
#define TRDgeometryHole_h
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

#include "AliTRDgeometry.h"

class AliTRDgeometryHole : public AliTRDgeometry {

 public:

  AliTRDgeometryHole();
  ~AliTRDgeometryHole();

          void    CreateGeometry(Int_t *); 
          Int_t   IsVersion() const { return 0; };
          void    Init();

          void    SetPHOShole()     { };
          void    SetRICHhole()     { };

          Bool_t  GetPHOShole()     { return kTRUE; };
          Bool_t  GetRICHhole()     { return kTRUE; };

 protected:

  Float_t         fClengthI[kNplan];               // Length of the inner chambers
  Float_t         fClengthM1[kNplan];              // Length of the middle chambers
  Float_t         fClengthM2[kNplan];              // 
  Float_t         fClengthO1[kNplan];              // Length of the outer chambers
  Float_t         fClengthO2[kNplan];              // 
  Float_t         fClengthO3[kNplan];              // 

  ClassDef(AliTRDgeometryHole,1)                   // TRD geometry with hole

};

#endif
