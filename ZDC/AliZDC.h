#ifndef ALIZDC_H
#define ALIZDC_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

////////////////////////////////////////////////
//  Manager and hits classes for set ZDC      //
////////////////////////////////////////////////
 
#include "AliDetector.h"
#include "AliHit.h"

 
class AliZDC : public AliDetector {

public:
  AliZDC();
  AliZDC(const char *name, const char *title);
  virtual      ~AliZDC() {}
  virtual void  AddHit(Int_t track, Int_t *vol, Float_t *hits);
  virtual void  BuildGeometry();
  virtual void  CreateGeometry() {}
  virtual void  CreateMaterials() {}
  Int_t          DistancetoPrimitive(Int_t px, Int_t py);
  virtual Int_t IsVersion() const =0;
  virtual void  StepManager();
 
protected:
  // Parameters for ZDCs geometry
  Float_t fDimZN[3];  // Dimensions of neutron detector
  Float_t fDimZP[3];  // Dimensions of proton detector
  Float_t fPosZN[3];  // Position of neutron detector
  Float_t fPosZP[3];  // Position of proton detector
  Float_t fFibZN[3];  // Fibers for neutron detector
  Float_t fFibZP[3];  // Fibers for proton detector
  Float_t fGrvZN[3];  // Grooves for neutron detector
  Float_t fGrvZP[3];  // Grooves for proton detector
  Int_t   fDivZN[3];  // Division for neutron detector
  Int_t   fDivZP[3];  // Division for proton detector
  Int_t   fTowZN[2];  // Tower for neutron detector
  Int_t   fTowZP[2];  // Tower for proton detector

   ClassDef(AliZDC,1)  // Zero Degree Calorimeter base class
};

//_____________________________________________________________________________
class AliZDChit : public AliHit {

public:
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
 
public:
  AliZDChit() {}
  AliZDChit(Int_t shunt, Int_t track, Int_t *vol, Float_t *hits);
  AliZDChit(AliZDChit* oldhit) {*this=*oldhit;}
  virtual ~AliZDChit() {}
  virtual Int_t GetVolume(Int_t i) {return fVolume[i];}
  virtual Float_t GetLightPMQ() {return fLightPMQ;}
  virtual Float_t GetLightPMC() {return fLightPMC;}
  virtual Float_t GetEnergy() {return fEnergy;}
  Int_t operator == (AliZDChit &quad) {
     Int_t i;
     if(fTrack!=quad.GetTrack()) return 0;
     for(i=0; i<2; i++) if(fVolume[i]!=quad.GetVolume(i)) return 0;
     return 1;
  }
  virtual AliZDChit& operator + (AliZDChit &quad) {
     fLightPMQ+=quad.GetLightPMQ();
     fLightPMC+=quad.GetLightPMC();
     fEnergy+=quad.GetEnergy();
     return *this;
  }
  virtual void Print(Option_t *) {
     printf(" -> HIT: vol[0] =  %d vol[1] =  %d Track: %d \n" 
            "  hit[3] = %f, hit[4] = %f, hit[5] = %f, SFlag = %f\n"
	    "  PMQLight = %f, PMCLight = %f, Energy %f\n ", 
	    fVolume[0],fVolume[1],fTrack,fPrimKinEn,fXImpact,fYImpact,
	    fSFlag,fLightPMQ,fLightPMC,fEnergy);
  }

  ClassDef(AliZDChit,1)  // Hits for the Zero Degree Calorimeters
};
 
#endif
