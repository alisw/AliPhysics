#ifndef ALIMUONPADHIT_H
#define ALIMUONPADHIT_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

#include <TObject.h>


class AliMUONPadHit : public TObject {
 
public:
   AliMUONPadHit() {
      fHitNumber=fQ=fPadX=fPadY=fQpad=fRSec=0;   
}
   AliMUONPadHit(Int_t *clhits);
   virtual ~AliMUONPadHit() {;}
   Int_t   HitNumber()  {return fHitNumber;}
   Int_t   Cathode()    {return fCathode;}
   Int_t   Q()          {return fQ;}
   Int_t   PadX()       {return fPadX;}   
   Int_t   PadY()       {return fPadY;}
   Int_t   QPad()       {return fQpad;}
   Int_t   RSec()       {return fRSec;}
   
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
