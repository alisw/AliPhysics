#ifndef ALITRDHIT_H
#define ALITRDHIT_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id: AliTRDhit.h,v */

////////////////////////////////////////////////
//  Hit class for the TRD                     //
////////////////////////////////////////////////

#include "AliHit.h"

//_____________________________________________________________________________
class AliTRDhit : public AliHit {
 
 public:

  AliTRDhit();
  AliTRDhit(Int_t shunt, Int_t track, Int_t *det, Float_t *hits);
  virtual ~AliTRDhit();

          Int_t   GetDetector() { return fDetector; };
          Float_t GetCharge()   { return fQ;        };

 protected:

  Int_t        fDetector;   // TRD detector number
  Float_t      fQ;          // Charge created by a hit (slow simulator only)
 
  ClassDef(AliTRDhit,2)     // Hit for the Transition Radiation Detector

};

#endif
