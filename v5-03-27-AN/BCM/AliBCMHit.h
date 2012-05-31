#ifndef ALIBCMHIT_H
#define ALIBCMHIT_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */
//
// Beam Condition Monitor (BCM) Hit Class
//
// andreas.morsch@cern.ch
//
#include "AliHit.h"

class AliBCMHit : public AliHit 
{
 public:
    AliBCMHit();
    AliBCMHit(Int_t shunt, Int_t track, Double_t x[4], Int_t isens, Float_t edep);
    Int_t    Id()   const { return fId;   }
    Float_t  Edep() const { return fEdep; }
    Float_t  Time() const { return fTime; }
 private:
    Int_t fId;      // ID of the sensor (11, 12, 13, 14) for z > 0 and (21, 22, 23, 24) for z < 0
    Float_t fEdep;  // Deposited energy [GeV]
    Float_t fTime;  // Time of hit
    
    ClassDef(AliBCMHit, 1) // BCM hit class 
};
      
#endif
