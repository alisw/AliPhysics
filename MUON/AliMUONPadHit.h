#ifndef ALIMUONPADHIT_H
#define ALIMUONPADHIT_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */
// Revision of includes 07/05/2004

#include <TObject.h>

class AliMUONPadHit : public TObject {
 
public:
   AliMUONPadHit();
   AliMUONPadHit(Int_t *clhits);
   virtual ~AliMUONPadHit() {;}

   Int_t   HitNumber() const {return fHitNumber;}
   Int_t   Cathode()   const {return fCathode;}
   Int_t   Q()         const {return fQ;}
   Int_t   PadX()      const {return fPadX;}   
   Int_t   PadY()      const {return fPadY;}
   Int_t   QPad()      const {return fQpad;}
   Int_t   RSec()      const {return fRSec;}
   
 private:
   Int_t     fHitNumber;    // Hit number
   Int_t     fCathode;      // Cathode number
   Int_t     fQ  ;          // Total charge      
   Int_t     fPadX  ;       // Pad number along X
   Int_t     fPadY  ;       // Pad number along Y
   Int_t     fQpad  ;       // Charge per pad
   Int_t     fRSec  ;       // R -sector of pad
   ClassDef(AliMUONPadHit,1)  // MUON Pad Hit
};
#endif
