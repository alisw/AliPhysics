#ifndef ALIZDCHIT_H
#define ALIZDCHIT_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

////////////////////////////////////////////////
//  Hits classes for set ZDC                  //
////////////////////////////////////////////////

#include "AliHit.h"

class AliZDCHit : public AliHit {

public:
  AliZDCHit() {}
  AliZDCHit(Int_t shunt, Int_t track, Int_t *vol, Float_t *hits);
  AliZDCHit(AliZDCHit* oldhit) {*this=*oldhit;}
  virtual ~AliZDCHit() {}

  // Getters 
  virtual Int_t   GetVolume(Int_t i) {return fVolume[i];}
  virtual Float_t GetLightPMQ() {return fLightPMQ;}
  virtual Float_t GetLightPMC() {return fLightPMC;}
  virtual Float_t GetEnergy() {return fEnergy;}


  // Data members
  Int_t      fVolume[2];    //Array of volumes
  Float_t    fX;	    //X-coord. in the hall RS
  Float_t    fY;	    //Y-coord. in the hall RS
  Float_t    fZ;	    //Z-coord. in the hall RS
  Float_t    fPrimKinEn;    //Primary particle energy
  Float_t    fXImpact;      //x-coord. of the impact point over the ZDC
  Float_t    fYImpact;      //y-coord. of the impact point over the ZDC
  Float_t    fSFlag;        //Secondary flag
  Float_t    fLightPMQ;     //Cerenkov light produced in each quadrant
  Float_t    fLightPMC;     //Cerenkov light seen by the common PM
  Float_t    fEnergy;       //Total energy deposited in eV
 

  // Operators
  Int_t operator == (AliZDCHit &quad) {
     Int_t i;
//      Superfluo finche' c'e' shunt = 1 !?!?
//     if(fTrack!=quad.GetTrack()) return 0;
     for(i=0; i<2; i++) if(fVolume[i]!=quad.GetVolume(i)) return 0;
     return 1;
  }
  
  virtual AliZDCHit& operator + (AliZDCHit &quad) {
     fLightPMQ+=quad.GetLightPMQ();
     fLightPMC+=quad.GetLightPMC();
     fEnergy+=quad.GetEnergy();
     return *this;
  }

  // Print method
  virtual void Print(Option_t *) {
     printf(" -> HIT: vol[0] =  %d vol[1] =  %d Track: %d \n" 
            "  Primary E = %f, Ximpact = %f, Yimpact = %f, SFlag = %f\n"
	    "  PMQLight = %f, PMCLight = %f,  Deposited E = %f\n ", 
	    fVolume[0],fVolume[1],fTrack,fPrimKinEn,fXImpact,fYImpact,
	    fSFlag,fLightPMQ,fLightPMC,fEnergy);
  }

  ClassDef(AliZDCHit,1)  // Hits for the Zero Degree Calorimeters
};
 
#endif
