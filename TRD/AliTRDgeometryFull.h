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

          void    CreateGeometry(Int_t *);
          Int_t   IsVersion() const { return 1; };
          void    Init();

          void    SetPHOShole()     { fPHOShole = kTRUE; };
          void    SetRICHhole()     { fRICHhole = kTRUE; };

          Bool_t  GetPHOShole()     { return fPHOShole;  };
          Bool_t  GetRICHhole()     { return fRICHhole;  };

 protected:

  Bool_t          fPHOShole;                       // Switch for the hole in front of the PHOS
  Bool_t          fRICHhole;                       // Switch for the hole in front of the RICH

  Float_t         fClengthI[kNplan];               // Length of the inner chambers
  Float_t         fClengthM1[kNplan];              // Length of the middle chambers
  Float_t         fClengthM2[kNplan];              // 
  Float_t         fClengthO1[kNplan];              // Length of the outer chambers
  Float_t         fClengthO2[kNplan];              // 
  Float_t         fClengthO3[kNplan];              // 

  ClassDef(AliTRDgeometryFull,1)                   // TRD geometry without hole

};

#endif
