#ifndef ALIRICHRECHIT1D_H
#define ALIRICHRECHIT1D_H


/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */
#include <TObject.h>
class AliRICHRecHit1D : public TObject {
public:
  Float_t     fTheta  ;            //Incidence Angle theta
  Float_t     fPhi  ;              //Incidence Angle phi
  Float_t     fOmega;              //Cherenkov angle omega
  Float_t     fX;                  //Impact coordinate x
  Float_t     fY;                  //Impact coordinate y
  Int_t       fGoodPhotons;        //Number of photons used for reconstruction
  Float_t     fEmissPoint;         //Emission point of the cherenkov photons
  Float_t     fCerPerPhoton[100];  //Recontructed cerenkov angle per photon
  Int_t     fPadsUsedX[100];      //List of pads used for reconstruction (x)
  Int_t     fPadsUsedY[100];      //List of pads used for reconstruction (y)

 public:
    AliRICHRecHit1D() {
      fTheta=fPhi=fOmega=0;
    }
    AliRICHRecHit1D(Int_t id, Float_t* rechit, Float_t* photons, Int_t* padsX, Int_t* padsY);
    virtual ~AliRICHRecHit1D() {}
    ClassDef(AliRICHRecHit1D,1)  //Reconstructed hit object for set:RICH
};

#endif






