#ifndef TRDhit_H
#define TRDhit_H
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

  Int_t        fDetector;   // TRD detector number
  Float_t      fQ;          // Charge created by a hit (slow simulator only)
 
 public:

  AliTRDhit() {}
  AliTRDhit(Int_t shunt, Int_t track, Int_t *det, Float_t *hits);
  virtual ~AliTRDhit() {};
 
  ClassDef(AliTRDhit,2)     // Hit for the Transition Radiation Detector

};

#endif
