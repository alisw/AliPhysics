#ifndef ALITRDGEOMETRYFULL_H
#define ALITRDGEOMETRYFULL_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

#include "AliTRDgeometry.h"

class AliTRDgeometryFull : public AliTRDgeometry {

 public:

  AliTRDgeometryFull();
  virtual ~AliTRDgeometryFull();

          void    CreateGeometry(Int_t *idtmed);
          Int_t   IsVersion() const   { return 1; };
          void    Init();

          void    SetPHOShole()       { fPHOShole = kTRUE; };
          void    SetRICHhole()       { fRICHhole = kTRUE; };

  virtual void    SetRowPadSize(Float_t size);

          Bool_t  GetPHOShole() const { return fPHOShole;  };
          Bool_t  GetRICHhole() const { return fRICHhole;  };

 protected:

  Bool_t          fPHOShole;                       // Switch for the hole in front of the PHOS
  Bool_t          fRICHhole;                       // Switch for the hole in front of the RICH

  Float_t         fClengthI[kNplan];               // Length of the inner chambers
  Float_t         fClengthM1[kNplan];              // Length of the middle chambers
  Float_t         fClengthM2[kNplan];              // Length of the middle chambers
  Float_t         fClengthO1[kNplan];              // Length of the outer chambers
  Float_t         fClengthO2[kNplan];              // Length of the outer chambers
  Float_t         fClengthO3[kNplan];              // Length of the outer chambers

  ClassDef(AliTRDgeometryFull,1)                   // TRD geometry without hole

};

#endif
