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
  Float_t     fMeanRadius;         //Mean radius from input digits
  Float_t     fOriginalOmega;      //Particle real Cerenkov angle
  Float_t     fOriginalTheta;      //Particle real incidence angle
  Float_t     fOriginalPhi;        //Particle real azimuthal angle


 public:
    AliRICHRecHit3D() {
      fTheta=fPhi=fOmega=0;
    }
    AliRICHRecHit3D(Int_t id, Float_t* rechit, Float_t omega, Float_t theta, Float_t phi);
    virtual ~AliRICHRecHit3D() {}
    ClassDef(AliRICHRecHit3D,1)  //Reconstructed hit object for set:RICH
};

#endif






