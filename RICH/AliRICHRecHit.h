#ifndef ALIRICHRECHIT_H
#define ALIRICHRECHIT_H


/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */
#include <TObject.h>
class AliRICHRecHit : public TObject {
public:
  Float_t     fTheta  ;            //Incidence Angle theta
  Float_t     fPhi  ;              //Incidence Angle phi
  Float_t     fOmega;              //Cherenkov angle omega
  Float_t     fX;                  //Impact coordinate x
  Float_t     fY;                  //Impact coordinate y
 public:
    AliRICHRecHit() {
      fTheta=fPhi=fOmega=0;
    }
    AliRICHRecHit(Int_t id, Float_t* rechit);
    virtual ~AliRICHRecHit() {}
    ClassDef(AliRICHRecHit,1)  //Reconstructed hit object for set:RICH
};

#endif






