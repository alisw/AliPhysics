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
  AliTRDhit(Int_t shunt, Int_t track, Int_t det, Float_t *hits, Int_t q);
  virtual ~AliTRDhit();

          Int_t GetDetector() const { return fDetector; };
          Int_t GetCharge() const   { return fQ;        };

 protected:

  UShort_t     fDetector;   // TRD detector number
  Short_t      fQ;          // Charge created by a hit. TR signals are negative.

  ClassDef(AliTRDhit,3)     // Hit for the Transition Radiation Detector

};

#endif
