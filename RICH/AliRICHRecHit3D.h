#ifndef ALIRICHRECHIT3D_H
#define ALIRICHRECHIT3D_H


/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */
#include <TObject.h>
class AliRICHRecHit3D : public TObject {
public:
  Float_t     fTheta  ;            //Incidence Angle theta
  Float_t     fPhi  ;              //Incidence Angle phi
  Float_t     fOmega;              //Cherenkov angle omega
  Float_t     fX;                  //Impact coordinate x
  Float_t     fY;                  //Impact coordinate y
 
 public:
    AliRICHRecHit3D() {
      fTheta=fPhi=fOmega=0;
    }
    AliRICHRecHit3D(Int_t id, Float_t* rechit);
    virtual ~AliRICHRecHit3D() {}
    ClassDef(AliRICHRecHit3D,1)  //Reconstructed hit object for set:RICH
};

#endif






